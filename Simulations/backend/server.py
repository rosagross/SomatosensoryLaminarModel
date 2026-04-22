import os
import sys

from fastapi import FastAPI, WebSocket
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import json
import time
import logging

_SIM_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_MODEL_DIR = os.path.join(_SIM_DIR, "model")
for _p in (_MODEL_DIR, _SIM_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from somato_model import SomatoModel

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

model = None
running = False
timestep = 0
_last_log_t = 0.0
speed_multiplier = 1.0
_last_stat_t = 0.0
_last_stat_step = 0
_log = None
_reset_armed = False
_reset_log_steps = 0
_ws_clients = 0
_ws_set = set()
_pending_params = {}


def _get_logger():
    global _log
    if _log is not None:
        return _log
    log_path = os.path.join(os.path.dirname(__file__), "server.log")
    logger = logging.getLogger("sim_server")
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    _log = logger
    return _log

PARAM_KEYS = [
    "simulation_dur",
    "step_size",
    "input_onset",
    "thal_connect",
    "extI_cellcounts",
    "strength_I",
    "bI_cellcounts",
    "thal_cellcounts",
    "sI_thal",
    "g_thal",
    "input_type",
    "area",
    "g_intercortical",
    "coupling_strength",
    "Ib_strength",
    "Iext_strength",
    "Iext_duration",
]


@app.post("/start")
async def start_simulation():

    global model, running, timestep, _last_stat_t, _last_stat_step, _reset_armed, _reset_log_steps, _pending_params

    if model is None:
        model = SomatoModel()
        if _pending_params:
            model.apply_params(_pending_params)
            _get_logger().info("[start] applied pending params: %s", _pending_params)
            _pending_params = {}
    # Always reset dynamic state on start, but keep current parameters.
    model.initialize_state()

    running = True
    timestep = 0
    _last_stat_t = 0.0
    _last_stat_step = 0
    _reset_armed = False
    _reset_log_steps = 0
    print(f"[start] simulation started at {time.strftime('%H:%M:%S')}")
    _get_logger().info("[start] simulation started")

    return {"status": "started"}

@app.get("/params")
async def get_params():
    global model
    if model is None:
        temp = SomatoModel()
        base = {k: getattr(temp, k) for k in PARAM_KEYS}
        base.update(_pending_params)
        return base
    return {k: getattr(model, k) for k in PARAM_KEYS}


@app.post("/stop")
async def stop_simulation():

    global running
    running = False
    print(f"[stop] simulation stopped at {time.strftime('%H:%M:%S')}")
    _get_logger().info("[stop] simulation stopped")

    return {"status": "stopped"}


@app.post("/reset")
async def reset_simulation():

    global model, timestep, _last_stat_t, _last_stat_step, running, _last_log_t, _reset_armed, _reset_log_steps, _ws_set

    if model is None:
        return {"error": "simulation not started"}

    running = False
    # Close any active WebSocket clients to avoid stepping during reset.
    if _ws_set:
        for client in list(_ws_set):
            try:
                await client.close()
            except Exception:
                pass
        _ws_set.clear()
    current_params = {k: getattr(model, k) for k in PARAM_KEYS}
    model = SomatoModel(params=current_params)
    model.initialize_state()
    timestep = 0
    _last_stat_t = 0.0
    _last_stat_step = 0
    _last_log_t = 0.0
    _reset_armed = True
    _reset_log_steps = 5
    print(f"[reset] state reset at {time.strftime('%H:%M:%S')}")
    print("v_current sum", float(model.v_current.sum()), "u_t sum", float(model.u_t.sum()), "t", model.t)
    _get_logger().info(
        "[reset] t=0 v_max=%.6f u_max=%.6f rate_max=%.6f",
        float(model.v_current.max()),
        float(model.u_t.max()),
        float(model.rate_current.max()),
    )
    _get_logger().info(
        "[reset] step_size=%.6f Iext0=%.6f Ib0=%.6f W_max=%.6f W_min=%.6f strength_I=%.6f g=%.6f Ib=%.6f Iext=%.6f Iext_dur=%.6f",
        float(model.step_size),
        float(model.Iext[0, 0]) if model.Iext.size else 0.0,
        float(model.Ib[0, 0]) if model.Ib.size else 0.0,
        float(model.W.max()),
        float(model.W.min()),
        float(model.strength_I),
        float(model.coupling_strength),
        float(model.Ib_strength),
        float(model.Iext_strength),
        float(model.Iext_duration),
    )

    return {"status": "reset"}

@app.post("/speed")
async def update_speed(payload: dict):
    global speed_multiplier
    try:
        speed_multiplier = float(payload.get("speed_multiplier", 1.0))
    except (TypeError, ValueError):
        return {"error": "invalid speed_multiplier"}
    print(f"[speed] speed_multiplier set to {speed_multiplier}")
    _get_logger().info("[speed] speed_multiplier set to %s", speed_multiplier)
    return {"status": "updated", "speed_multiplier": speed_multiplier}


@app.post("/params")
async def update_params(params: dict):

    global model, _pending_params

    if model is None:
        _pending_params.update(params)
        _get_logger().info("[params] queued (no model yet): %s", params)
        return {"status": "queued"}

    model.apply_params(params)
    print(f"[params] updated: {params}")
    _get_logger().info("[params] updated: %s", params)

    return {"status": "updated"}


@app.websocket("/ws")
async def websocket_endpoint(ws: WebSocket):

    global timestep, _last_log_t, running, speed_multiplier, _last_stat_t, _last_stat_step, _reset_armed, _reset_log_steps, _ws_clients

    await ws.accept()
    print("[ws] client connected")
    _ws_clients += 1
    _ws_set.add(ws)
    _get_logger().info("[ws] client connected (clients=%s)", _ws_clients)

    try:
        while True:

            if running and model is not None:
                step_idx = timestep % len(model.steps)
                rates = model.simulate_step(step_idx)
                potential = model.v_current.sum(axis=1)

                if _reset_armed:
                    _get_logger().info(
                        "[post-reset-step] t=%s v_max=%.6f u_max=%.6f rate_max=%.6f pot_max=%.6f",
                        timestep,
                        float(model.v_current.max()),
                        float(model.u_t.max()),
                        float(model.rate_current.max()),
                        float(potential.max()),
                    )
                    _reset_armed = False
                elif _reset_log_steps > 0:
                    _get_logger().info(
                        "[post-reset-step+%s] t=%s v_max=%.6f u_max=%.6f rate_max=%.6f pot_max=%.6f",
                        5 - _reset_log_steps + 1,
                        timestep,
                        float(model.v_current.max()),
                        float(model.u_t.max()),
                        float(model.rate_current.max()),
                        float(potential.max()),
                    )
                    _reset_log_steps -= 1

                await ws.send_text(json.dumps({
                    "rates": rates.tolist(),
                    "potential": potential.tolist(),
                    "t": timestep
                }))

                timestep += 1

                now = time.time()
                if now - _last_log_t > 2.0:
                    _last_log_t = now
                    print(f"[ws] sent timestep={timestep} rate0={rates[0]:.4f}")
                    _get_logger().info(
                        "[ws] sent timestep=%s rate0=%.6f pot_max=%.6f",
                        timestep,
                        float(rates[0]),
                        float(potential.max()),
                    )
                if now - _last_stat_t > 2.0:
                    dt = now - _last_stat_t if _last_stat_t > 0 else None
                    if dt:
                        steps = timestep - _last_stat_step
                        steps_per_s = steps / dt if dt > 0 else 0.0
                        sim_s_per_s = steps_per_s * model.step_size
                        print(f"[perf] {steps_per_s:.1f} steps/s, {sim_s_per_s:.3f} sim s/s (speed x{speed_multiplier})")
                    _last_stat_t = now
                    _last_stat_step = timestep
                sleep_s = model.step_size
                await asyncio.sleep(sleep_s)
            else:
                await asyncio.sleep(0.01)
    except Exception as e:
        print(f"[ws] error: {e}")
    finally:
        print("[ws] client disconnected")
        _ws_clients = max(0, _ws_clients - 1)
        _ws_set.discard(ws)
        _get_logger().info("[ws] client disconnected (clients=%s)", _ws_clients)
