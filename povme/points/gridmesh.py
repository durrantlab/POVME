"""Grid mesh construction and manipulation for pocket detection.

This module provides the [`GridMesh`][points.gridmesh.GridMesh] class, which represents a box of
equidistant 3D points used during POVME's pocket-detection phase. The class
supports:

- Generating a regular grid that encompasses a protein's bounding box.
- Removing points outside a convex hull.
- Removing points that clash with protein atoms.
- Filtering isolated points that lack sufficient neighbors.
- Expanding the grid to higher resolution around surviving points.
- Separating surviving points into distinct pockets.
"""

import math
from typing import Any

import numpy as np
import numpy.typing as npt
from loguru import logger
from scipy.ndimage import label as ndimage_label
from scipy.spatial import KDTree

from povme.config import PocketVolumeConfig
from povme.io import write_pdbs
from povme.parallel import RayManager, RayTaskGeneral


class TaskRemovePointsOutsideHull(RayTaskGeneral):
    """A class to remove points outside a convex hull using multiple processors."""

    def process_item(
        self, item: tuple[Any, npt.NDArray[np.float64]]
    ) -> npt.NDArray[np.float64]:
        """Removes points outside the convex hull.

        Args:
            item: A tuple containing:
                - hull: The convex hull object.
                - some_points: A numpy array of points to be tested.

        Returns:
            A numpy array of points that are inside the convex hull.
        """
        try:
            hull, some_points = item
            if len(some_points) == 0:
                return np.array([]).reshape(0, 3)

            # Process all points against each triangle simultaneously.
            inside = np.ones(len(some_points), dtype=bool)
            epsilon = 1.0e-5

            for triangle in hull.hull:
                # Only test points still considered "inside"
                candidate_idx = np.where(inside)[0]
                if len(candidate_idx) == 0:
                    break

                candidate_pts = some_points[candidate_idx]

                # Vector from triangle vertex 0 to each candidate point
                rel_points = candidate_pts - triangle[0]  # (K, 3)

                # Triangle edge vectors
                vec1 = triangle[1] - triangle[0]  # (3,)
                vec2 = triangle[2] - triangle[1]  # (3,)

                # Outward-facing normal of this triangle face
                cross = np.cross(vec1, vec2)  # (3,)

                # Dot product: positive means the point is on the
                # "outside" of this face
                dots = rel_points @ cross  # (K,)

                # Mark points outside this face as outside the hull
                outside_mask = dots > epsilon
                inside[candidate_idx[outside_mask]] = False

            return some_points[inside]
        except Exception as e:
            logger.exception(f"Error in Removing Points Outside Hull: {e}")
            return np.array([])  # Return empty array on error


class GridMesh:
    """A class representing a box of equidistant points."""

    def __init__(self, box: npt.NDArray[np.float64], res: float) -> None:
        """Initialize the class.

        Args:
            box: A numpy array representing two 3D points, (min_x, min_y,
                min_z) and (max_x, max_y, max_z), that define a box.
            res: The space between the points of the box, in the X, Y, and
                Z direction.
        """

        self.write_pdbs = write_pdbs()

        min_x = self.__snap_float(box[0][0], res)
        min_y = self.__snap_float(box[0][1], res)
        min_z = self.__snap_float(box[0][2], res)
        max_x = self.__snap_float(box[1][0], res) + 1.1 * res
        max_y = self.__snap_float(box[1][1], res) + 1.1 * res
        max_z = self.__snap_float(box[1][2], res) + 1.1 * res

        x, y, z = np.mgrid[min_x:max_x:res, min_y:max_y:res, min_z:max_z:res]
        self.points = np.array(list(zip(x.ravel(), y.ravel(), z.ravel())))

    def __snap_float(
        self, val: npt.NDArray[np.float64], res: float
    ) -> npt.NDArray[np.float64]:
        """Snaps an arbitrary point to the nearest grid point.

        Args:
            val: A numpy array corresponding to a 3D point.
            res: The resolution (distance in the X, Y, and Z directions
                between adjacent points) of the grid.

        Returns:
            A numpy array corresponding to a 3D point near val that is on a
                nearby grid point.
        """

        return np.floor(val / res) * res

    def remove_points_outside_convex_hull(self, hull, config):
        """Removes box points that are outside a convex hull.

        Args:
            hull: The convex hull.
            config: Configuration object containing `n_cores`.
        """
        # Prepare input as list of tuples: (hull, some_points)
        chunks = [(hull, t) for t in np.array_split(self.points, config.n_cores)]

        # Initialize RayManager with the appropriate task class
        ray_manager = RayManager(
            task_class=TaskRemovePointsOutsideHull,
            n_cores=config.n_cores,
            use_ray=config.use_ray,
        )
        ray_manager.submit_tasks(items=chunks)
        processed_chunks = ray_manager.get_results()

        # Each element in processed_chunks is either a numpy array of new_pts or an error tuple
        valid_points = []
        for result in processed_chunks:
            if isinstance(result, tuple) and result[0] == "error":
                logger.error(f"Error removing points outside hull: {result[1]}")
            else:
                if isinstance(result, np.ndarray) and result.size > 0:
                    valid_points.append(result)
                elif not isinstance(result, np.ndarray):
                    logger.warning(f"Unexpected result type: {type(result)}")

        if valid_points:
            self.points = np.vstack(valid_points)
        else:
            self.points = np.array([])

    def remove_all_points_close_to_other_points(
        self,
        other_points: npt.NDArray[np.float64],
        dist_cutoff: float,
        config: PocketVolumeConfig,
    ) -> None:
        """Remove grid points that are within a cutoff of protein atoms.

        Args:
            other_points: An `(m, 3)` array of protein-atom coordinates.
            dist_cutoff: Grid points closer than this distance to any
                atom are removed.
            config: Configuration object (`n_cores`, `use_ray`).
        """
        if len(self.points) == 0 or len(other_points) == 0:
            return

        # Build a single KDTree on the grid points
        grid_tree = KDTree(self.points)

        # For each atom, find all grid points within dist_cutoff.
        # query_ball_point on the atom array is efficient: one bulk call.
        atom_tree = KDTree(other_points)
        # pairs[i] is a list of grid-point indices close to atom i
        close_pairs = atom_tree.query_ball_tree(grid_tree, r=dist_cutoff)

        # Collect all unique grid-point indices that are too close
        close_set = set()
        for idx_list in close_pairs:
            close_set.update(idx_list)

        if close_set:
            self.points = np.delete(self.points, sorted(close_set), axis=0)

    def to_pdb(self, let="X"):
        """Converts the points in this box into a PDB representation.

        Args:
            let: An optional string, the chain ID to use. "X" by default.

        Returns:
            A PDB-formatted string.

        """
        return self.write_pdbs.numpy_to_pdb(self.points, let)

    def expand_around_existing_points(self, num_pts, reso):
        """Add points to the current box that surround existing points,
        essentially increasing the resolution of the box.

        For each surviving grid point, new points are placed at all
        integer multiples of `reso` within a cube of half-width
        `num_pts x reso` centred on the original point. Duplicates are removed.

        Args:
            num_pts: An int, the number of points to place on each side of
                the existing points, in the X, Y, and Z directions.
            res: The distance between adjacent added points.
        """
        i = np.arange(-num_pts * reso, num_pts * reso + reso * 0.01, reso)

        # Generate all (K, 3) offset vectors at once
        gx, gy, gz = np.meshgrid(i, i, i, indexing="ij")
        offsets = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

        # Broadcast: (1, N, 3) + (K, 1, 3) -> (K, N, 3) -> (K*N, 3)
        all_new = (self.points[np.newaxis, :, :] + offsets[:, np.newaxis, :]).reshape(
            -1, 3
        )

        self.points = all_new
        self.__unique_points()

    def __unique_points(self):
        """Identifies unique points (rows) in an array of points.

        Args:
            a: A `n x 3` np.array representing 3D points.

        Returns:
            A `n x 2` np.array containing the 3D points that are unique.
        """

        b = np.ascontiguousarray(self.points).view(
            np.dtype((np.void, self.points.dtype.itemsize * self.points.shape[1]))
        )
        unique_points = (
            np.unique(b).view(self.points.dtype).reshape(-1, self.points.shape[1])
        )

        self.points = unique_points

    def filter_isolated_points_until_no_change(self, reso, number_of_neighbors):
        """Iteratively remove points with too few neighbors.

        Points on the fringe of a pocket often have fewer grid neighbors
        than points in the pocket interior.  This method repeatedly
        removes any point with fewer than `number_of_neighbors` neighbors
        (counted within the diagonal distance of one grid cell) until the
        point set stabilizes.

        Args:
            res: The grid spacing. The neighbor cutoff is derived
                as `reso x sqrt(3) x 1.1` to include diagonal (kitty-corner)
                neighbors.
            number_of_neighbors: The minimum number of permissible neighbors.
        """
        cutoff = reso * math.sqrt(3.0) * 1.1

        num_pts = 0
        while num_pts != len(self.points):
            num_pts = len(self.points)
            tree = KDTree(self.points)

            # Count neighbors for each point (subtract 1 to exclude self).
            # return_length=True returns just the count, not the full
            # neighbor lists, saving memory and time.
            counts = (
                np.asarray(
                    tree.query_ball_point(self.points, r=cutoff, return_length=True),
                    dtype=np.int64,
                )
                - 1
            )

            keep_mask = counts >= number_of_neighbors
            self.points = self.points[keep_mask]

    def separate_out_pockets(self) -> list[npt.NDArray[np.float64]]:
        """Partition surviving points into distinct pockets.

        Two points belong to the same pocket if they are connected through
        a chain of grid neighbors (26-connectivity, i.e., including
        diagonal/kitty-corner neighbors).

        1. Points are mapped to integer ``(i, j, k)`` indices by
            subtracting the grid minimum and dividing by the grid spacing.
        2. A 3D boolean volume is constructed and filled at the
            appropriate indices.
        3. [`scipy.ndimage.label`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html)
            performs connected-component
            labelling with 26-connectivity.
        4. Each connected component is extracted as a separate pocket.

        Returns:
            A list of `(n_i, 3)` arrays, one per pocket, sorted in
                descending order of point count (largest pocket first).
        """

        # Infer the grid spacing from the minimum nonzero pairwise distance
        # along any single axis.  Since points are on a regular grid, the
        # smallest positive coordinate difference equals the spacing.
        diffs_x = np.diff(np.unique(np.round(self.points[:, 0], 8)))
        reso = float(diffs_x[diffs_x > 1e-9].min()) if len(diffs_x) > 0 else 1.0

        # Map to integer grid indices
        grid_min = np.min(self.points, axis=0)
        indices = np.round((self.points - grid_min) / reso).astype(int)
        shape = tuple(indices.max(axis=0) + 1)

        # Build 3-D boolean volume
        volume = np.zeros(shape, dtype=bool)
        volume[indices[:, 0], indices[:, 1], indices[:, 2]] = True

        # Connected-component labelling with 26-connectivity
        struct = np.ones((3, 3, 3), dtype=int)
        labels_3d, n_components = ndimage_label(volume, structure=struct)

        # Map labels back to the original point array
        point_labels = labels_3d[indices[:, 0], indices[:, 1], indices[:, 2]]

        pockets: list[npt.NDArray[np.float64]] = []
        for label_id in range(1, n_components + 1):
            pocket_pts = self.points[point_labels == label_id]
            pockets.append(pocket_pts)

        # Sort by size, largest first (matches original behavior)
        pockets.sort(key=lambda p: -len(p))

        return pockets
