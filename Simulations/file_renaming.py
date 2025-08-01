import os
import re
from pathlib import Path

# Directory containing the files (change this to your folder path)
directory = Path.cwd().joinpath('output')

# Regex pattern to match the old filename format
old_pattern = re.compile(
    r"potentials_g(\d+)_bEI([\d\.]+)_(.+)_IbStrength([\d\.]+)_Iduration([\d\.]+)_(.+)IextStrength([\d\.]+)_Ionset([\d\.]+)_tauVisual_(.+)_thalEI0_S1S2\.csv"
)  # rates_g200_bEI0.5_somato_IbStrength0_Iduration1_stepIextStrength140_Ionset1.001_tauVisual_thalJiang_thalEI0_S1S2
    #params_g100_bEI0.2_somato_IbStrength0_Iduration0.5_stepIextStrength100_Ionset1.001_tauVisual_thalJiang_thalEI0_S1S2
# New filename format
new_format = "potentials_g{g}_bEI{bEI}_Ib{bI}_Iextd{d}_{input_type}Iexts{s}_Ionset{input_onset}_tauVisual_{thalamus_source}_thalEI0_S1S2.csv"

# Loop through all files in the directory
for filename in os.listdir(directory):
    match = old_pattern.fullmatch(filename)
    if match:
        print('MATCH')
        # Extract matched groups
        g, bEI, cortex_type, bI, d, input_type, s, input_onset, thalamus_source = match.groups()
        
        # Construct new filename
        new_filename = new_format.format(
            g=g,
            bEI=bEI,
            bI=bI,
            d=d,
            input_type=input_type,
            s=s,
            input_onset=input_onset,
            thalamus_source=thalamus_source
        )
        
        # Full paths for old and new names
        old_path = os.path.join(directory, filename)
        new_path = os.path.join(directory, new_filename)
        
        # Rename the file
        os.rename(old_path, new_path)
        print(f"Renamed: {filename} → {new_filename}")