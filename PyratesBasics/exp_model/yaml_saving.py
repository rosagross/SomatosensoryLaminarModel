# %%
import numpy as np
from ruamel.yaml import YAML
from pyrates.frontend.template import CircuitTemplate
from pyrates.frontend.fileio.yaml import dump_to_yaml

def circuit_to_yaml(circuit: CircuitTemplate, path: str):
    """
    Save a Pyrates CircuitTemplate to YAML safely, converting NumPy types to Python types.
    """
    # Patch ruamel.yaml to handle numpy scalars
    yaml = YAML()
    def represent_numpy_scalar(dumper, data):
        return dumper.represent_float(float(data))
    yaml.representer.add_representer(np.float64, represent_numpy_scalar)
    yaml.representer.add_representer(np.float32, represent_numpy_scalar)
    yaml.representer.add_representer(np.int32, lambda d, x: d.represent_int(int(x)))
    yaml.representer.add_representer(np.int64, lambda d, x: d.represent_int(int(x)))

    dump_to_yaml(circuit, path=path)
    print(f"CircuitTemplate successfully saved to {path}")
# %%
