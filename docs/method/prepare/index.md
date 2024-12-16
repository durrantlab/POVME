# Trajectory preparation

Before measuring pocket volumes, it is essential to ensure that all protein structures in the simulation are aligned.
Proteins are flexible, and their binding pockets can move, rotate, or change shape during a simulation.
To make sure POVME measures the same pocket across all frames, the structures need to be aligned to a reference.
Tools like [Visual Molecular Dynamics (VMD)](https://www.ks.uiuc.edu/Research/vmd/) or [MDAnalysis](https://www.mdanalysis.org/) are often used to align the protein structures.

!!! note
    POVME currently only supports loading data from PDB files.

## Aligning with MDAnalysis

TODO:

```python
import os

import MDAnalysis as mda
from MDAnalysis import transformations

# Specifying where our simulations are
dir_simulations = "~/simulations/example/trajectories"
path_topology = os.path.join(dir_simulations, "mol.prmtop")
path_trajectories = ["prod-1.nc", "prod-2.nc", "prod-3.nc"]
path_trajectories = [os.path.join(dir_simulations, p) for p in path_trajectories]

# Specifying where to store our PDB file
dir_save = "~/simulations/example/analyses"
path_pdb = os.path.join(dir_save, "traj_aligned.pdb")
os.makedirs(dir_save, exist_ok=True)

# Loading simulation
u = mda.Universe(path_topology, path_trajectories)

# Selecting atoms for alignment
atoms_str = "protein"
atoms = u.select_atoms(atoms_str)

ref_u = u.copy()
atoms_ref = ref_u.select_atoms(atoms_str)

# Preparing alignment transformations
# These are done on-the-fly for each trajectory step.
ag = u.atoms
workflow = (
    transformations.unwrap(ag),
    transformations.center_in_box(atoms, center="mass"),
    transformations.fit_rot_trans(atoms, atoms_ref),
)
u.trajectory.add_transformations(*workflow)

# Write PDB file for every `stride` steps
stride = 5
with mda.Writer(path_pdb, atoms.n_atoms) as W:
    for ts in u.trajectory[None:None:stride]:
        W.write(u.select_atoms(atoms_str))
```
