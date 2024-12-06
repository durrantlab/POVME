from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from scipy.spatial.distance import cdist


def snap_points(points: npt.NDArray[np.float64], res: float) -> npt.NDArray[np.float64]:
    """
    Aligns a set of 3D points to the nearest points on a uniform grid.

    This function "snaps" each point to the nearest grid point based on a specified
    grid resolution. It is useful for discretizing continuous 3D coordinates into a
    fixed grid.

    Args:
        points:
            A numpy array of shape (n, 3), where `n` is the number of points.
            Each row represents a point in 3D space with [x, y, z] coordinates.
        res:
            The resolution of the grid, defining the distance between adjacent grid points.

    Returns:
        A numpy array of the same shape as `points`, where each point is
        adjusted to lie on the nearest grid point defined by `res`.

    Example:
        ```python
        >>> import numpy as np
        >>> points = np.array([[1.2, 3.5, 7.8], [2.7, 4.1, 6.3]])
        >>> snap_points(points, res=1.0)
        array([[1.0, 4.0, 8.0],
               [3.0, 4.0, 6.0]])
        ```
    """
    points = points + 1e-10  # Adjust for rounding errors
    return np.around(points / res) * res


def build_mesh_grid(
    xs: npt.NDArray[np.float64],
    ys: npt.NDArray[np.float64],
    zs: npt.NDArray[np.float64],
    res: float,
) -> npt.NDArray[np.float64]:
    """
    Creates a 3D grid of points using given `x`, `y`, and `z` coordinate vectors.

    This function generates all possible combinations of `xs`, `ys`, and `zs` to create
    a rectangular 3D grid, then snaps the resulting points to the nearest fixed grid
    using the specified resolution.

    Args:
        xs:
            A 1D numpy array representing the x-coordinates of the grid.
        ys:
            A 1D numpy array representing the y-coordinates of the grid.
        zs:
            A 1D numpy array representing the z-coordinates of the grid.
        res:
            The resolution of the final snapped grid.

    Returns:
        A 2D numpy array of shape (n, 3), where each row is a point on the
        3D grid after snapping to the fixed resolution.

    Example:
        ```python
        >>> import numpy as np
        >>> xs = np.array([0.0, 1.0, 2.0])
        >>> ys = np.array([0.0, 1.0])
        >>> zs = np.array([0.0, 1.0])
        >>> build_mesh_grid(xs, ys, zs, res=0.5)
        array([[0.0, 0.0, 0.0],
               [0.0, 0.0, 0.5],
               [0.0, 0.5, 0.0],
               ...
               [2.0, 1.0, 1.0]])
        ```
    """
    points = np.array(np.meshgrid(xs, ys, zs, indexing="ij")).T.reshape(-1, 3)
    points = snap_points(points, res)
    return points


class Region(ABC):
    """
    Abstract base class representing a 3D region.

    This class provides a blueprint for defining specific types of 3D regions.
    Subclasses must implement methods to provide a string representation of the
    region and to generate a grid of points filling the region.
    """

    def __init__(self):
        pass

    @abstractmethod
    def __str__(self):
        """Returns a string representation of the region."""
        pass

    @abstractmethod
    def get_points(self, res: float) -> npt.NDArray[np.float64]:
        """
        Generates a grid of equally spaced points filling the region.

        Args:
            res:
                The resolution of the grid, defining the spacing between points.

        Returns:
            A numpy array of shape (n, 3), where each row is a [x, y, z]
            coordinate representing a point in the region.
        """
        pass


class SphericalRegion(Region):
    """
    A 3D spherical region defined by a center and a radius.

    This class represents a sphere in 3D space, capable of generating a grid
    of points that fill the sphere. Points outside the sphere boundary are excluded.
    """

    def __init__(self, center: list[float], radius: float) -> None:
        """
        Initialize a spherical region.

        Args:
            center:
                The x, y, and z coordinates of the center of the sphere.
            radius:
                The radius of the sphere in Angstroms.
        """
        super().__init__()
        self.center = np.array(center)
        """`x`, `y`, and `z` coordinates in Angstroms of the center."""
        self.radius = radius
        """Sphere radius in Angstroms."""

    def __str__(self):
        """Returns a string representation of the spherical region."""
        return (
            f"sphere at ({self.center[0]}, {self.center[1]}, {self.center[2]}), "
            f"radius = {self.radius}"
        )

    def get_points(self, res: float) -> npt.NDArray[np.float64]:
        """
        Generates a grid of points filling the spherical region.

        Args:
            res:
                The resolution of the grid, defining the spacing between points.

        Returns:
            A numpy array of shape (n, 3), where each row is a [x, y, z]
            coordinate representing a point within the spherical region.
        """
        xs = np.arange(
            self.center[0] - self.radius,
            self.center[0] + self.radius + res,
            res,
        )
        ys = np.arange(
            self.center[1] - self.radius,
            self.center[1] + self.radius + res,
            res,
        )
        zs = np.arange(
            self.center[2] - self.radius,
            self.center[2] + self.radius + res,
            res,
        )

        result = build_mesh_grid(xs, ys, zs, res)
        # Remove points outside the sphere
        distances = cdist(result, self.center.reshape(1, -1)).flatten()
        inside_sphere = distances < self.radius
        return result[inside_sphere]


class RectangularRegion(Region):
    """
    A 3D rectangular region defined by its center and dimensions.

    This class represents a rectangular box in 3D space, capable of generating a
    grid of points that fill the box.
    """

    def __init__(self, center: list[float], lengths: list[float]) -> None:
        """
        Initialize a rectangular region.

        Args:
            center:
                The x, y, and z coordinates of the center of the box.
            lengths:
                The lengths of the box along the x, y, and z axes.
        """
        super().__init__()
        self.center = np.array(center)
        """A list of three floats representing the x, y, and z coordinates of
            the box's center, in Angstroms."""
        self.lengths = np.array(lengths)
        """A list of three floats representing the lengths of the box along the
            x, y, and z axes, in Angstroms."""

    def __str__(self):
        """Returns a string representation of the rectangular region."""
        return (
            f"box centered at ({self.center[0]}, {self.center[1]}, {self.center[2]}) "
            f"with dimensions ({self.lengths[0]}, {self.lengths[1]}, {self.lengths[2]})"
        )

    def get_points(self, res):
        """
        Generates a grid of points filling the rectangular region.

        Args:
            res:
                The resolution of the grid, defining the spacing between points.

        Returns:
            A numpy array of shape (n, 3), where each row is a [x, y, z]
            coordinate representing a point within the rectangular region.
        """
        xs = np.arange(
            self.center[0] - self.lengths[0] / 2,
            self.center[0] + self.lengths[0] / 2 + res,
            res,
        )
        ys = np.arange(
            self.center[1] - self.lengths[1] / 2,
            self.center[1] + self.lengths[1] / 2 + res,
            res,
        )
        zs = np.arange(
            self.center[2] - self.lengths[2] / 2,
            self.center[2] + self.lengths[2] / 2 + res,
            res,
        )

        return build_mesh_grid(xs, ys, zs, res)


def collect_regions(
    spherical_configs: list[dict[str, list[float] | float]] | None,
    rectangular_configs: list[dict[str, list[float]]] | None,
) -> list[SphericalRegion | RectangularRegion]:
    """
    Creates a collection of 3D regions from spherical and rectangular configurations.

    This function takes configurations for spherical and rectangular regions, instantiates
    the respective region objects, and returns a combined list of these regions. Empty
    or invalid configurations are skipped.

    Args:
        spherical_configs:
            A list of dictionaries, each specifying a spherical region.
            Each dictionary must have:

              - `"center"` (`list[float]`): The x, y, and z coordinates of the
                sphere's center.
              - `"radius"` (`float`): The radius of the sphere.

            If `None`, no spherical regions are added.
        rectangular_configs:
            A list of dictionaries, each specifying a rectangular region.
            Each dictionary must have:

              - `"center"` (`list[float]`): The x, y, and z coordinates of the
                box's center.
              - `"lengths"` (`list[float]`): The lengths of the box along the
                x, y, and z axes.

            If `None`, no rectangular regions are added.

    Returns:
        A list containing `SphericalRegion` and `RectangularRegion` objects
            created based on the provided configurations.

    Example:
        ```python
        >>> spherical_configs = [
        ...     {"center": [0.0, 0.0, 0.0], "radius": 5.0},
        ...     {"center": [10.0, 10.0, 10.0], "radius": 3.0},
        ... ]
        >>> rectangular_configs = [
        ...     {"center": [0.0, 0.0, 0.0], "lengths": [10.0, 5.0, 2.0]}
        ... ]
        >>> regions = collect_regions(spherical_configs, rectangular_configs)
        >>> for region in regions:
        ...     print(region)
        sphere at (0.0, 0.0, 0.0), radius = 5.0
        sphere at (10.0, 10.0, 10.0), radius = 3.0
        box centered at (0.0, 0.0, 0.0) with dimensions (10.0, 5.0, 2.0)
        ```
    """
    regions: list[SphericalRegion | RectangularRegion] = []

    if spherical_configs is not None:
        for config in spherical_configs:  # type: ignore
            if len(config) == 0:
                continue
            regions.append(SphericalRegion(config["center"], config["radius"]))  # type: ignore

    if rectangular_configs is not None:
        for config in rectangular_configs:  # type: ignore
            if len(config) == 0:
                continue
            regions.append(RectangularRegion(config["center"], config["lengths"]))  # type: ignore
    return regions
