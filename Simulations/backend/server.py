import os
import sys

from fastapi import FastAPI, WebSocket
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import json
import time

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

PARAM_KEYS = [
    "simulation_dur",
    "step_size",
    "input_onset",
    "thal_connect",
    "extI_cellcounts",
    "balance_EI",
    "bI_cellcounts",
    "thal_cellcounts",
    "bEI_thal",
    "g_thal",
    "input_type",
    "area",
    "coupling_strength",
    "Ib_strength",
    "Iext_strength",
    "Iext_duration",
]


@app.post("/start")
async def start_simulation():

    global model, running, timestep

    model = SomatoModel()
    model.initialize_state()

    running = True
    timestep = 0
    print(f"[start] simulation started at {time.strftime('%H:%M:%S')}")

    return {"status": "started"}

@app.get("/params")
async def get_params():
    global model
    if model is None:
        temp = SomatoModel()
        return {k: getattr(temp, k) for k in PARAM_KEYS}
    return {k: getattr(model, k) for k in PARAM_KEYS}


@app.post("/stop")
async def stop_simulation():

    global running
    running = False
    print(f"[stop] simulation stopped at {time.strftime('%H:%M:%S')}")

    return {"status": "stopped"}


@app.post("/reset")
async def reset_simulation():

    global model, timestep

    if model is None:
        return {"error": "simulation not started"}

    model.initialize_state()
    timestep = 0
    print(f"[reset] state reset at {time.strftime('%H:%M:%S')}")

    return {"status": "reset"}

@app.post("/speed")
async def update_speed(payload: dict):
    global speed_multiplier
    try:
        speed_multiplier = float(payload.get("speed_multiplier", 1.0))
    except (TypeError, ValueError):
        return {"error": "invalid speed_multiplier"}
    print(f"[speed] speed_multiplier set to {speed_multiplier}")
    return {"status": "updated", "speed_multiplier": speed_multiplier}


@app.post("/params")
async def update_params(params: dict):

    global model

    if model is None:
        return {"error": "simulation not started"}

    model.apply_params(params)
    print(f"[params] updated: {params}")

    return {"status": "updated"}


@app.websocket("/ws")
async def websocket_endpoint(ws: WebSocket):

    global timestep, _last_log_t, running, speed_multiplier

    await ws.accept()
    print("[ws] client connected")

    try:
        while True:

            if running and model is not None:
                step_idx = timestep % len(model.steps)
                rates = model.simulate_step(step_idx)
                potential = model.v_current.sum(axis=1)

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
                sleep_s = model.step_size / max(speed_multiplier, 1e-6)
                await asyncio.sleep(sleep_s)
            else:
                await asyncio.sleep(0.01)
    except Exception as e:
        print(f"[ws] error: {e}")
    finally:
        print("[ws] client disconnected")
