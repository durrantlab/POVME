from loguru import logger
from pydantic import BaseModel

from .io import YamlIO


class POVMEConfig(BaseModel, YamlIO):
    grid_spacing: float = 1.0
    """The distance, in Angstroms, between adjacent points.

    Smaller values improve accuracy but increase computational cost.
    """

    load_points_path: str | None = None
    """Load points from `npy` file."""

    points_inclusion_sphere: list[list[float]] = []
    """A list of spheres to include in the pocket-encompassing region. 

    Each tuple contains (x, y, z, radius).
    """

    points_inclusion_box: list[list[float]] = []
    """A list of rectangular prisms ('boxes') to include in the pocket-encompassing region. 

    Each tuple contains (center_x, center_y, center_z, size_x, size_y, size_z).
    """

    points_exclusion_sphere: list[list[float]] = []
    """A list of spheres to exclude from the pocket-encompassing region. 

    Each tuple contains (x, y, z, radius).
    """

    points_exclusion_box: list[list[float]] = []
    """A list of rectangular prisms ('boxes') to exclude from the pocket-encompassing region. 

    Each tuple contains (center_x, center_y, center_z, size_x, size_y, size_z).
    """

    save_points: bool = False
    """Whether to save the point field to a PDB file for visualization or reuse."""

    load_points_filename: str | None = None
    """Optional filename to load a previously saved point field. 

    Should use the `.pdb.npy` format.
    """

    distance_cutoff: float = 1.09
    """The distance from a receptor atom's van der Waals surface below which points are excluded.

    Default is 1.09 Angstroms, the van der Waals radius of a hydrogen atom.
    """

    convex_hull_exclusion: bool = True
    """Whether to calculate the convex hull of receptor atoms near the pocket and exclude points outside it."""

    contiguous_pocket_seed_sphere: list[list[float]] = []
    """Seed regions for contiguous pocket detection. 

    Each tuple contains (x, y, z, radius).
    """

    contiguous_pocket_seed_box: list[list[float]] = []
    """Seed boxes for contiguous pocket detection.

    Each tuple contains (center_x, center_y, center_z, size_x, size_y, size_z).
    """

    contiguous_points_criteria: int = 4
    """The minimum number of neighboring points required to consider two pocket volumes contiguous."""

    num_processors: int = 4
    """Number of processors to use for the calculation on Unix-based systems."""

    use_disk_not_memory: bool = False
    """Whether to prioritize disk space over memory for large trajectory files."""

    save_individual_pocket_volumes: bool = False
    """Whether to save the pocket-volume points for each frame to separate PDB files."""

    save_pocket_volumes_trajectory: bool = False
    """Whether to save all pocket-volume points for each frame to a single PDB trajectory file."""

    output_equal_num_points_per_frame: bool = False
    """Whether to add extra points at the origin (0.0, 0.0, 0.0) to ensure the same number of points in each frame."""

    save_volumetric_density_map: bool = False
    """Whether to save a volumetric density map in DX format."""

    compress_output: bool = False
    """Whether to compress all output files using gz compression to save disk space."""

    def log(self) -> None:
        """Log all configuration parameters."""
        logger.info("Logging POVME Configuration Settings:")
        config_dict = self.model_dump()
        for key, value in config_dict.items():
            logger.info(f"{key}: {value}")
