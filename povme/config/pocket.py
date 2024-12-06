from loguru import logger
from pydantic import BaseModel

from povme.config.io import YamlIO


class PocketIDConfig(BaseModel, YamlIO):
    pocket_detection_resolution: float = 4.0
    """The distance between probe points used to initially find the pockets"""

    pocket_measuring_resolution: float = 1.0
    """The distance between probe points used to measure identified pockets in greater detail. Should divide --pocket_detection_resolution evenly."""

    clashing_cutoff: float = 3.0
    """In measuring the pockets, any points closer than this cutoff to receptor atoms will be removed."""

    n_neighbors: int = 4
    """In measuring the pockets, any points with fewer than this number of neighbors will be deleted. These are usually just stray points that don't belong to any real pocket."""

    n_spheres: int = 5
    """The number of inclusion spheres to generate for each pocket."""

    sphere_padding: float = 5.0
    """How much larger the radius of the inclusion spheres should be, beyond what is required to encompass the identified pockets."""

    n_cores: int = 1
    """The number of processors to use."""

    def log(self) -> None:
        """Log all configuration parameters."""
        logger.info("Logging Pocket Detection Configuration Settings:")
        config_dict = self.model_dump()
        for key, value in config_dict.items():
            logger.info(f"{key}: {value}")
