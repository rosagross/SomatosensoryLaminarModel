import jax
import jax.numpy as jnp
from jax import jit, lax
from functools import partial

@jit
def sigmoid(v, r, v_thr, m_max):
    return m_max / (1 + jnp.exp(r * (v_thr - v)))

def compute_rate(v_current, sigm_params):
    v_sum = jnp.sum(v_current, axis=-1) * 1e3  # Convert to mV
    return jax.vmap(sigmoid)(v_sum, sigm_params[:, 0], sigm_params[:, 1], sigm_params[:, 2])

@partial(jit, static_argnums=(0))
def run_jax_simulation(nPop, sigm, W, H, tau, Iext, Ib, step_size, steps):

    def step_fn(carry, t_idx):
        v_current, u_t, rate_hist, pot_hist = carry

        rate_current = compute_rate(v_current, sigm)

        def update_i(i, carry_ij):
            v_i, u_i = carry_ij

            def update_j(j, acc):
                v_row, u_row = acc
                tau_ij = tau[i, j]
                h_ij = H[i, j]
                w_ij = W[i, j]
                u_dot = (h_ij / tau_ij) * (w_ij * rate_current[j]) - 2 * u_row[j] / tau_ij - v_row[j] / (tau_ij ** 2)
                u_row = u_row.at[j].add(u_dot * step_size)
                v_row = v_row.at[j].add(u_row[j] * step_size)
                return v_row, u_row

            v_i, u_i = jax.lax.fori_loop(0, nPop, update_j, (v_i, u_i))
            return v_i, u_i

        v_current, u_t = jax.lax.fori_loop(0, nPop, update_i, (v_current, u_t))

        for i in range(nPop):
            # External input
            v_current = v_current.at[i, -1].add(u_t[i, -1] * step_size)
            u_dot_ext = (H[i, -1] / tau[i, -1]) * (W[i, -1] * Iext[i, t_idx]) - 2 * u_t[i, -1] / tau[i, -1] - v_current[i, -1] / (tau[i, -1] ** 2)
            u_t = u_t.at[i, -1].add(u_dot_ext * step_size)

            # Background input
            v_current = v_current.at[i, -2].add(u_t[i, -2] * step_size)
            u_dot_bg = (H[i, -2] / tau[i, -2]) * (W[i, -2] * Ib[i, t_idx]) - 2 * u_t[i, -2] / tau[i, -2] - v_current[i, -2] / (tau[i, -2] ** 2)
            u_t = u_t.at[i, -2].add(u_dot_bg * step_size)

        rate_hist = rate_hist.at[:, t_idx].set(rate_current)
        pot_hist = pot_hist.at[:, :, t_idx].set(v_current)

        return (v_current, u_t, rate_hist, pot_hist), None

    v_init = jnp.zeros((nPop, nPop + 2))
    u_init = jnp.zeros((nPop, nPop + 2))
    rate_hist = jnp.zeros((nPop, len(steps)))
    pot_hist = jnp.zeros((nPop, nPop + 2, len(steps)))

    (v_final, u_final, rate_out, pot_out), _ = lax.scan(
        step_fn,
        (v_init, u_init, rate_hist, pot_hist),
        jnp.arange(len(steps))
    )

    return rate_out, pot_out
