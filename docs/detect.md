# Pocket detection

`POVME` simplifies pocket detection through its `povme detect` CLI, which provides an intuitive and efficient way to run the detection process.

## Preparing the Configuration File

To begin, create a YAML configuration file that will be used to update [`PocketDetectConfig`][config.detect.PocketDetectConfig] that defines the parameters for pocket detection.
Here is an example configuration:

```yaml
pocket_detection_resolution: 4.0
pocket_measuring_resolution: 1.0
clashing_cutoff: 3.0
n_neighbors: 4
n_spheres: 5
sphere_padding: 5.0
use_ray: false
n_cores: 4
```

::: config.detect.PocketDetectConfig
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

## Running the Command

Once the configuration file is ready, you can use the `povme detect` command to identify pockets.
The command syntax is:

```bash
povme detect -c config.yaml -i protein.pdb -o output/
```

Here:

-   `-c` specifies the path to the YAML configuration file.
-   `-i` specifies the input PDB file.
-   `-o` specifies the output directory or file prefix.

## Analyzing the results

The `povme detect` command produces a set of PDB files that comprehensively describe the detected pockets and their associated parameters.
Each PDB file contains a set of points that represent a POVME-detected pocket.
For example,

```text
ATOM      1    A AAA A   1      35.000  19.000  30.000                       A
ATOM      1    B AAA B   1      34.000  19.000  29.000                       B
ATOM      2    B AAB B   2      34.000  20.000  29.000                       B
ATOM      1    C AAA C   1      35.000  21.000  28.000                       C
ATOM      2    C AAB C   2      36.000  20.000  27.000                       C
ATOM      3    C AAC C   3      36.000  21.000  27.000                       C
ATOM      4    C AAD C   4      36.000  21.000  28.000                       C
ATOM      1    D AAA D   1      36.000  20.000  28.000                       D
ATOM      1    E AAA E   1      35.000  20.000  29.000                       E
ATOM      2    E AAB E   2      35.000  20.000  30.000                       E
ATOM      3    E AAC E   3      35.000  21.000  29.000                       E
ATOM      4    E AAD E   4      36.000  20.000  29.000                       E
ATOM      5    E AAE E   5      36.000  21.000  29.000                       E
```

In addition to these points, the PDB files contain remarks OF `PointsInclusionSphere` parameter inputs that would encompass this pocket plus any `sphere_padding`.
For example, the header of `pocket1.pdb` might look like this:

```text
REMARK CHAIN A: PointsInclusionSphere 35.0 19.0 30.0 5.0
REMARK CHAIN B: PointsInclusionSphere 34.0 19.5 29.0 5.5
REMARK CHAIN C: PointsInclusionSphere 35.75 20.75 27.5 5.94
REMARK CHAIN D: PointsInclusionSphere 36.0 20.0 28.0 5.0
REMARK CHAIN E: PointsInclusionSphere 35.4 20.4 29.2 5.98
```

This means we would use five inclusion spheres when computing volumes and would specify them in a `povme.yml` file like so.

```yaml
points_inclusion_sphere:
  - center: [35.0, 19.0, 30.0]
    radius: 5.0
  - center: [34.0, 19.5, 29.0]
    radius: 5.5
  - center: [35.75, 20.75, 27.5]
    radius: 5.94
  - center: [36.0, 20.0, 28.0]
    radius: 5.0
  - center: [35.4, 20.4, 29.2]
    radius: 5.98
```

## Visualization

The output files can be visualized using molecular visualization tools such as PyMOL or VMD.
These tools allow you to inspect the detected pockets and verify their biological relevance.

<div id="pocket-detect-rel1" class="mol-container"></div>
<script>
document.addEventListener('DOMContentLoaded', (event) => {
    const viewer = molstar.Viewer.create('pocket-detect-rel1', {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: true,
        viewportShowAnimation: false,
        pdbProvider: 'rcsb',
    }).then(viewer => {
      //   viewer.loadStructureFromUrl("/files/example-povme-pockets.pdb", "pdb");
        viewer.loadSnapshotFromUrl("/files/example-povme-pockets.molx", "molx");
    });
});
</script>
