"""
Jax simulation function to run simulations faster!
"""

import jax
import jax.numpy as jnp
from jax import jit, lax
from functools import partial

@jit
def sigmoid(v, r, v_thr, m_max):
    return m_max / (1 + jnp.exp(r * (v_thr - v)))

def compute_rate(v_current, sigm_params):
    v_sum = jnp.sum(v_current, axis=-1) * 1e3  # mV
    return jax.vmap(sigmoid)(v_sum, sigm_params[:,0], sigm_params[:,1], sigm_params[:,2])


# the partial function is part of the jit package, allowing for compilation
@partial(jit, static_argnums=(0,))
def run_simulation_jax(params, W, H, tau, Iext, Ib, step_size, steps):

    nPop = params['nPop']
    sigm = params['sigm']

    def step_fn(carry, t):
        v_current, u_t, rate_hist, pot_hist = carry

        # Compute rate from current potential
        rate_current = compute_rate(v_current, sigm)

        # Update v and u (second derivative)
        def update(i, args):
            v_i, u_i = args
            
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

            # external input 

            # background input

        v_current, u_t = jax.lax.fori_loop(0, nPop, update, (v_current, u_t))


        # Save history
        rate_hist = rate_hist.at[:, t].set(rate_current)
        pot_hist = pot_hist.at[:, :, t].set(v_current)

        return (v_current, u_t, rate_hist, pot_hist), None

    # Init
    v_init = jnp.zeros((nPop, nPop+2))
    u_init = jnp.zeros((nPop, nPop+2))
    rate_hist = jnp.zeros((nPop, len(steps)))
    pot_hist = jnp.zeros((nPop, nPop+2, len(steps)))

    (v_final, u_final, rate_out, pot_out), _ = lax.scan(
        step_fn,
        (v_init, u_init, rate_hist, pot_hist),
        jnp.arange(len(steps))
    )

    return rate_out, pot_out