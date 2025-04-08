"""
Perform forward modelling with simulated EEG data. 

"""

# %% 
import mne
from mne.datasets import sample

# %%
data_path = sample.data_path()

# the raw file containing the channel location + types
sample_dir = data_path / "MEG" / "sample"
raw_fname = sample_dir / "sample_audvis_raw.fif"
# The paths to Freesurfer reconstructions
subjects_dir = data_path / "subjects"
subject = "sample"

# %%
################ Computing the forward operator #####################
# We need:
# coregistration info (-trans.fif file)
# source space 
# BEM surfaces (tissue surfaces)


######################## Coregistration #######################
# set up coregistration
# use mne.gui.coregistration() or mne coreg 



# visualize coreg 
# The transformation file obtained by coregistration
trans = sample_dir / "sample_audvis_raw-trans.fif"

info = mne.io.read_info(raw_fname)
# Here we look at the dense head, which isn't used for BEM computations but
# is useful for coregistration.
mne.viz.plot_alignment(
    info,
    trans,
    subject=subject,
    dig=True,
    meg=["helmet", "sensors"],
    subjects_dir=subjects_dir,
    surfaces="head-dense",
)

# %%
########################  BEM surfaces ##############################


# visualize BEM
plot_bem_kwargs = dict(
    subject=subject,
    subjects_dir=subjects_dir,
    brain_surfaces="white",
    orientation="coronal",
    slices=[50, 100, 150, 200], 
)

mne.viz.plot_bem(**plot_bem_kwargs)

# %%
###################### Compute source space ########################
# two different options here
# surface based
# volumetric / discrete

src = mne.setup_source_space(
    subject, spacing="oct4", add_dist="patch", subjects_dir=subjects_dir
)
print(src)

# %% visualize the sources on top of the BEM images
mne.viz.plot_bem(src=src, **plot_bem_kwargs)

# %% 
# VOLUME based 

sphere = (0.0, 0.0, 0.04, 0.09)
vol_src = mne.setup_volume_source_space(
    subject,
    subjects_dir=subjects_dir,
    sphere=sphere,
    sphere_units="m",
    add_interpolator=False,
)  # just for speed!
print(vol_src)

mne.viz.plot_bem(src=vol_src, **plot_bem_kwargs)

# %%
# visualization volume based 

fig = mne.viz.plot_alignment(
    subject=subject,
    subjects_dir=subjects_dir,
    surfaces="white",
    coord_frame="mri",
    src=src,
)
mne.viz.set_3d_view(
    fig,
    azimuth=173.78,
    elevation=101.75,
    distance=0.30,
    focalpoint=(-0.03, -0.01, 0.03),
) 

# %%

conductivity = (0.3,)  # for single layer
# conductivity = (0.3, 0.006, 0.3)  # for three layers
model = mne.make_bem_model(
    subject="sample", ico=4, conductivity=conductivity, subjects_dir=subjects_dir
)
bem = mne.make_bem_solution(model)

# %%
########################## Compute Forward Solution ###################

# specify how many layers to compute with the conductivity parameter
conductivity = (0.3,)  # for single layer
# conductivity = (0.3, 0.006, 0.3)  # for three layers

# compute BEM model of head geometry and tissue conductivities 
model = mne.make_bem_model(
    subject="sample", ico=4, conductivity=conductivity, subjects_dir=subjects_dir
)
bem = mne.make_bem_solution(model)

# %%

# computation of the lead field matrix 
fwd = mne.make_forward_solution(
    raw_fname,
    trans=trans,
    src=src,
    bem=bem,
    meg=True,
    eeg=False,
    mindist=5.0,
    n_jobs=None,
    verbose=True,
)
print(fwd)

# %%
leadfield = fwd["sol"]["data"]
print(f"Leadfield size : {leadfield.shape[0]} sensors x {leadfield.shape[1]} dipoles")


# %%
