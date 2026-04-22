# Interactive plotting page

Use two terminals: one for the backend API/websocket server, one for the frontend static page.

## 1) Start backend (Terminal 1)

```bash
cd /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel
source bat/start.sh
pip install fastapi uvicorn numpy matplotlib seaborn pandas h5py pyyaml
cd Simulations/backend
python -m uvicorn server:app --host 127.0.0.1 --port 8000 --reload
```

## 2) Start frontend (Terminal 2)

```bash
cd /data/p_02989/Modelling/grossmannr_wd/SomatosensoryLaminarModel/Simulations/frontend
python -m http.server 8080
```

## 3) Open the page

Open:

`http://127.0.0.1:8080/index.html`

The frontend is configured to connect to `127.0.0.1:8000`, so keep both processes running.
