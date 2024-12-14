# Pocket detection

Pocket detection is an essential process in computational biology and drug design, enabling researchers to identify and characterize binding sites on protein surfaces.
These pockets or cavities are critical for understanding molecular interactions, particularly in the context of ligand binding, drug discovery, and protein functionality.
By accurately detecting pockets, researchers can explore their shapes, volumes, and properties, facilitating the design of effective inhibitors or activators.

This guide focuses on using `POVME` for pocket detection, detailing its purpose, the methodology behind the detection process, and how to use the `povme detect` command-line interface (CLI) for real-world scenarios.

## Purpose of Pocket Detection

Proteins interact with small molecules, ions, and other macromolecules through specific regions called pockets.
These pockets are often formed by the three-dimensional arrangement of amino acid residues and are highly significant for biological function.
Identifying pockets is the first step toward analyzing potential druggable sites and understanding how ligands interact with proteins.

Pocket detection involves identifying regions on the protein surface that:

-   Have sufficient depth and volume to accommodate a ligand.
-   Are spatially enclosed by surrounding residues.
-   Could serve as potential binding sites for drugs or other molecules.

## Pocket Detection Methodology

Pocket detection in `POVME` is a multi-step process that systematically identifies biologically significant pockets on protein surfaces.
This methodology is grounded in computational geometry and configurable algorithms that balance precision with computational efficiency.

### Initial Point Selection

The first step involves selecting initial points around the protein structure.
These points are spaced according to the `pocket_detection_resolution` parameter in the configuration file.
This coarse resolution allows for efficient initial exploration of potential pocket regions.

Example Configuration:

```yaml
pocket_detection_resolution: 4.0
```

Code Example:

```python
from povme.points import GridMesh

# Generate a grid of points surrounding the protein
protein_box = molecule.information.get_bounding_box()
initial_points = GridMesh(protein_box, config.pocket_detection_resolution * 4)
```

### Convex Hull Calculation

The convex hull is computed to identify the smallest convex shape that encloses the alpha carbons of the protein.
This step excludes shallow surface indentations and ensures that only biologically relevant regions are considered.

Key algorithm: **Gift Wrapping Algorithm**.
This algorithm iteratively selects points to form triangular facets of the convex hull.

Code Example:

```python
from povme.points.hull import ConvexHull

# Calculate convex hull of alpha carbon coordinates
alpha_carbon_coords = molecule.selections.get_coordinates("alpha_carbons")
convex_hull = ConvexHull(alpha_carbon_coords)
```

### Removing Redundant Points

Points that fall outside the convex hull are removed, reducing computational overhead.
Additionally, points within a specified distance (`clashing_cutoff`) of protein atoms are excluded.

Example Configuration:

```yaml
clashing_cutoff: 3.0
```

Code Example:

```python
# Remove points outside the convex hull
points_inside_hull = convex_hull.filter_points(initial_points)

# Exclude points too close to protein atoms
filtered_points = points_inside_hull.remove_close_points(
    molecule.get_coordinates(), config.clashing_cutoff
)
```

### Refining the Point Field

The remaining points are refined to improve pocket detection accuracy.
This includes:

1.  **Increasing Point Density**: Points are added around existing ones for higher resolution, controlled by `pocket_measuring_resolution`.
2.  **Filtering Isolated Points**: Points with fewer than `n_neighbors` are removed to eliminate noise.

Example Configuration:

```yaml
pocket_measuring_resolution: 1.0
n_neighbors: 4
```

### Identifying and Saving Pockets

Finally, the refined points are grouped into discrete pockets using clustering algorithms. Inclusion spheres are generated for each pocket, and the results are saved in PDB format.

Example Configuration:

```yaml
n_spheres: 5
sphere_padding: 5.0
```

Code Example:

```python
from povme.io import write_pdbs

# Separate points into pockets
pockets = final_points.separate_into_pockets()

# Save each pocket as a PDB file
write_pdbs(pockets, output_prefix="output/pocket")
```

## Using the `povme detect` CLI

`POVME` simplifies pocket detection through its `povme detect` CLI, which provides an intuitive and efficient way to run the detection process. Below, we outline the steps for using this command.

### Preparing the Configuration File

To begin, create a YAML configuration file that defines the parameters for pocket detection. Here is an example configuration:

```yaml
pocket_detection_resolution: 4.0
pocket_measuring_resolution: 1.0
clashing_cutoff: 3.0
n_neighbors: 4
n_spheres: 5
sphere_padding: 5.0
use_ray: false
n_cores: 2
```

The parameters include:

-   `pocket_detection_resolution`: The spacing between probe points for initial detection.
-   `pocket_measuring_resolution`: A finer spacing for detailed measurement.
-   `clashing_cutoff`: Minimum distance between probe points and protein atoms.
-   `n_neighbors`: Minimum number of neighboring points for a point to be retained.
-   `n_spheres`: Number of inclusion spheres per pocket.
-   `sphere_padding`: Additional padding for inclusion spheres.
-   `use_ray` and `n_cores`: Enable parallel computation for improved performance.

### Running the Command

Once the configuration file is ready, you can use the `povme detect` command to identify pockets. The command syntax is:

```bash
povme detect -c config.yaml -i protein.pdb -o output/
```

Here:

-   `-c` specifies the path to the YAML configuration file.
-   `-i` specifies the input PDB file.
-   `-o` specifies the output directory or file prefix.

### Example Workflow

1.  **Input PDB File**: Ensure that your protein structure file (`protein.pdb`) is prepared. The PDB file should contain all relevant atoms and be free from major structural errors.
2.  **Configuration File**: Save the YAML configuration as `config.yaml`.
3.  **Run the Command**: Execute the following command:

   ```bash
   povme detect -c config.yaml -i protein.pdb -o ./output/
   ```

4.  **Review Output**: Navigate to the specified output directory to find the results.

### Analyzing the Results

The output of `povme detect` includes several key files:

-   **Pocket PDB Files**: Files such as `output/pocket1.pdb` contain the three-dimensional points defining each detected pocket.
-   **Log File**: Detailed logs summarize the detection process, including any warnings or errors encountered.
-   **Inclusion Sphere Parameters**: Each pocket PDB file includes header comments listing the inclusion sphere parameters for that pocket.

For example, the header of `pocket1.pdb` might look like this:

```text
REMARK Pocket #1
REMARK CHAIN A: PointsInclusionSphere 10.5 20.7 15.3 5.0
```

This indicates the center and radius of the inclusion sphere representing the pocket.

### Visualization

The output files can be visualized using molecular visualization tools such as PyMOL or VMD.
These tools allow you to inspect the detected pockets and verify their biological relevance.

## A Complete Example

Here is a step-by-step example demonstrating pocket detection on a sample protein.

### Example Configuration File (`config.yaml`)

```yaml
pocket_detection_resolution: 4.0
pocket_measuring_resolution: 1.0
clashing_cutoff: 3.0
n_neighbors: 4
n_spheres: 5
sphere_padding: 5.0
use_ray: false
n_cores: 2
```

### Command to Run

```bash
povme detect -c config.yaml -i sample_protein.pdb -o results/
```

### Output Files

After running the command, the `results/` directory will contain:

1.  `pocket1.pdb`: A PDB file defining the first detected pocket.
2.  `pocket2.pdb`: Additional PDB files for other pockets.
3.  `povme.log`: A log file summarizing the pocket detection process.

### Visualization

Open `pocket1.pdb` in PyMOL:

```bash
pymol results/pocket1.pdb
```

This will display the first pocket, along with its inclusion spheres, allowing you to evaluate the detection results.
