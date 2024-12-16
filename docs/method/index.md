# Methodology

The POVME (POcket Volume MEasurer) algorithm is an advanced computational framework designed to quantify and characterize the volumes and shapes of ligand-binding pockets within macromolecular structures.

## General workflow

POVME works step by step to provide accurate measurements of binding pockets.
The process begins by aligning the protein structures, so the binding pocket stays in the same position across all frames of a molecular dynamics simulation.
Then, regions around the binding pocket are defined.
These regions are filled with points, and points that don’t belong in the pocket are removed.
Finally, the volume of the pocket is calculated based on the remaining points.
