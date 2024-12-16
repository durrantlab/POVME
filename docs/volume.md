# Computing volume

The `povme volume` command in POVME is designed to compute and analyze the volume of molecular pockets based on user-specified inclusion and exclusion regions.
This command can process both static protein structures and molecular dynamics (MD) trajectories to identify binding pocket volumes and characteristics.

## YAML configuration

The command is configured through a YAML input file that will be used to update [`PocketVolumeConfig`][config.volume.PocketVolumeConfig] that defines the parameters for pocket volume calculations.

```yaml
grid_spacing: 1.0
load_points_path: null

points_inclusion_sphere:
  - center: [65.0, 98.0, 50.0]
    radius: 16.0
  - center: [-100.0, -100.0, -100.0]
    radius: 10.0

save_points: true
distance_cutoff: 1.09
convex_hull_exclusion: true

contiguous_pocket_seed_sphere:
  - center: [67.0, 102.0, 57.0]
    radius: 4.0
contiguous_points_criteria: 3

use_ray: true
n_cores: 8

save_individual_pocket_volumes: true
save_pocket_volumes_trajectory: true
output_equal_num_points_per_frame: true
save_volumetric_density_map: true
compress_output: false
```

::: config.volume.PocketVolumeConfig
    options:
      show_bases: false
      heading_level: 3
      show_root_heading: false
      show_root_members_full_path: false
      parameter_headings: true
      show_symbol_type_heading: false
      show_symbol_type_toc: false
      show_root_toc_entry: false
      filters:
      - "!^log$"

## Running the command

Once the configuration file is ready, you can use the `povme detect` command to identify pockets.
The command syntax is:

```bash
povme volume -c config.yaml -i input.pdb -o output/
```

-   `-c`: Path to the YAML configuration file.
-   `-i`: Input file, typically a PDB structure or trajectory.
-   `-o`: Output directory.

## Outputs

By default, POVME writes several files to the disk.
POVME will save the volume measurements to a CSV file that can be easily pasted into popular spreadsheet programs.

POVME also optionally saves the pocket-occupying points of each frame to PDB file(s) on this disk.
The user can instruct the program to save these points to separate files and/or to a single PDB trajectory.
Some visualization programs (e.g., VMD) are only compatible with trajectories with the same
number of atoms in each frame.
To satisfy this requirement, POVME can optionally write extra points to the origin (0.0, 0.0, 0.0) on a frame-by-frame basis.

Finally, POVME also optionally saves a volumetric density map in the Data Explorer (DX) format.
A volumetric density value is associated with each pocket-occupying point by calculating the fraction of all trajectory pocket volumes that include the given point.
If the density map is displayed as an isosurface, the value of the isosurface expresses the fraction of time (e.g., throughout the simulation) that the pocket included the displayed volume.
