# %% Looking at preimplemented PyRates models
path = "model_templates.neural_mass_models.jansenrit.JRC" # PyRates path to the model
from pyrates.frontend.file import parse_path
from pprint import pprint
template_name, filename, directory = parse_path(path)

    # test if file can be found (and potentially add extension)
import os

if "." in filename:
    filepath = os.path.join(directory, filename)
else:
        # this is actually the default case for the internal interface
    for ext in ["yaml", "yml"]:
        filepath = os.path.join(directory, ".".join((filename, ext)))
        if os.path.exists(filepath):
            break
        else:
            raise FileNotFoundError(f"Could not identify file with name {filename} in directory {directory}.")

    # load as yaml file
from ruamel.yaml import YAML

yaml = YAML(typ="safe", pure=True)
with open(filepath, "r") as file:
        file_dict = yaml.load(file)

pprint(file_dict)

# %% intermediate level
from pyrates.frontend.fileio.yaml import dict_from_yaml
template_dict = dict_from_yaml(path)
pprint(template_dict)
# %%
