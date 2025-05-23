from typing import Any

import numpy as np
import numpy.typing as npt
from loguru import logger
from scipy import spatial

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
            new_pts = [pt for pt in some_points if hull.inside_hull(pt)]
            return np.array(new_pts)
        except Exception as e:
            logger.exception(f"Error in Removing Points Outside Hull: {e}")
            return np.array([])  # Return empty array on error


class TaskGetClosePoints(RayTaskGeneral):
    """A class to identify box points that are near other, user-specified points."""

    def process_item(
        self, item: tuple[spatial.KDTree, float, npt.NDArray[np.float64]]
    ) -> npt.NDArray[np.float64]:
        """Identifies indices of box points close to other points.

        Args:
            item: A tuple containing:
                - box_of_pts_distance_tree: KDTree of box points.
                - dist_cutoff: The cutoff distance.
                - other_points: Numpy array of other points.

        Returns:
            A numpy array of unique indices of box points that are within dist_cutoff.
        """
        try:
            box_of_pts_distance_tree, dist_cutoff, other_points = item

            # Create KDTree for other_points
            other_points_distance_tree = spatial.KDTree(other_points)

            # Find all box points within dist_cutoff of any other point
            sparce_distance_matrix = other_points_distance_tree.sparse_distance_matrix(
                box_of_pts_distance_tree, dist_cutoff
            )

            # Extract unique indices of box points that are close to other points
            indices_of_box_pts_close_to_molecule_points = np.unique(
                sparce_distance_matrix.tocsr().indices
            )

            return indices_of_box_pts_close_to_molecule_points
        except Exception as e:
            logger.exception(f"Error in Getting Close Points: {e}")
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

        x, y, z = np.mgrid[min_x:max_x:res, min_y:max_y:res, min_z:max_z:res]  # type: ignore
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
        """Removes all points in this box that come within the points specified
        in a numpy array

        Args:
            other_points: A numpy array containing the other points.
            dist_cutoff: A float, the cutoff distance to use in determining
                whether or not box points will be removed.
            config: Configuration object containing `n_cores`.
        """

        # note, in newer versions of scipy use cKDTree
        box_of_pts_distance_tree = spatial.KDTree(self.points)

        # Prepare input as list of tuples: (box_of_pts_distance_tree, dist_cutoff, t)
        chunks = [
            (box_of_pts_distance_tree, dist_cutoff, t)
            for t in np.array_split(other_points, config.n_cores)
        ]

        # Initialize RayManager with the appropriate task class
        ray_manager = RayManager(
            task_class=TaskGetClosePoints,
            n_cores=config.n_cores,
            use_ray=config.use_ray,
        )
        ray_manager.submit_tasks(items=chunks)
        processed_chunks = ray_manager.get_results()

        # Each element in processed_chunks is either a numpy array of indices or an error tuple
        valid_indices = []
        for result in processed_chunks:
            if isinstance(result, tuple) and result[0] == "error":
                logger.error(f"Error getting close points: {result[1]}")
            else:
                valid_indices.append(result)

        if valid_indices:
            indices_of_box_pts_close_to_molecule_points = np.unique(
                np.hstack(valid_indices)
            )
            self.points = np.delete(
                self.points, indices_of_box_pts_close_to_molecule_points, axis=0
            )  # remove the ones that are too close to molecule atoms

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

        Args:
            num_pts: An int, the number of points to place on each side of
                the existing points, in the X, Y, and Z directions.
            res: The distance between adjacent added points.

        """

        new_pts = []

        i = np.arange(-num_pts * reso, num_pts * reso + reso * 0.01, reso)
        for xi in i:
            for yi in i:
                for zi in i:
                    vec = np.array([xi, yi, zi])
                    new_pts.append(self.points + vec)
        self.points = np.vstack(new_pts)

        self.__unique_points()

    def __unique_points(self):
        """Identifies unique points (rows) in an array of points.

        Args:
            a: A nx3 np.array representing 3D points.

        Returns:
            A nx2 np.array containing the 3D points that are unique.

        """

        b = np.ascontiguousarray(self.points).view(
            np.dtype((np.void, self.points.dtype.itemsize * self.points.shape[1]))
        )
        unique_points = (
            np.unique(b).view(self.points.dtype).reshape(-1, self.points.shape[1])
        )

        self.points = unique_points

    def filter_isolated_points_until_no_change(self, reso, number_of_neighbors):
        """Keep removing points that don't have enough neighbors, until no
        such points exist.

        Args:
            res: The distance between adjacent points.
            number_of_neighbors: The minimum number of permissible neighbors.

        """

        # calculate the pairwise distances between all box points note, in
        # newer versions of scipy use cKDTree
        box_of_pts_distance_tree = spatial.KDTree(self.points)

        # so kiddy-corner counted as a neighbor
        self.dist_matrix = box_of_pts_distance_tree.sparse_distance_matrix(
            box_of_pts_distance_tree, reso * np.sqrt(3.0) * 1.1
        ).todense()

        # note that the diagnol of self.dist_matrix is zero, as expected, but
        # ones with dist > reso * np.sqrt(3.0) * 1.1 are also 0. Pretty
        # convenient.

        num_pts = 0
        # keep running the pass until there are no changes (points are stable)
        while num_pts != len(self.points):
            num_pts = len(self.points)

            # identify the points that have enough neighbors
            columns_nonzero_count = np.array((self.dist_matrix != 0).sum(0))[0]
            columns_nonzero_count_match_criteria = (
                columns_nonzero_count >= number_of_neighbors
            )
            columns_nonzero_count_match_criteria_index = np.nonzero(
                columns_nonzero_count_match_criteria
            )

            self.__keep_limited_points(columns_nonzero_count_match_criteria_index)

    def __keep_limited_points(self, pt_indices):
        """A support function"""

        # keep only those points
        self.points = self.points[pt_indices]

        # update the distance matrix so it doesn't need to be recalculated
        self.dist_matrix = self.dist_matrix[pt_indices, :][0]
        self.dist_matrix = self.dist_matrix.T
        self.dist_matrix = self.dist_matrix[pt_indices, :][0]
        # self.dist_matrix = self.dist_matrix.T # not necessary because it's a symetrical matrix

    def separate_out_pockets(self) -> list[npt.NDArray[np.float64]]:
        """Separate the points according to the pocket they belong to.
        Determined by looking at patches of contiguous points.

        Returns:
            A list of point arrays, each array corresponding to the points of a
                separate pocket.

        """

        all_pockets = []

        # self.points is an array of 3D points

        # self.dist_matrix is a distance matrix. self.dist_matrix[i,j] is
        # distance between points i and j. But it only contains distances if
        # the points are close (neighbors). Otherwise, 0.

        # Keep going until there are no more points that need to be assigned
        # to a pocket.
        while len(self.points) != 0:
            pocket_indexes = np.array([0])
            num_pts_in_pocket = 0

            # Keep looping into no new unique pockets are added.
            while num_pts_in_pocket != len(pocket_indexes):
                num_pts_in_pocket = len(pocket_indexes)

                # Get all the adjacent points
                indices_of_neighbors = np.nonzero(self.dist_matrix[pocket_indexes, :])[
                    1
                ]

                # Get one of them. Not sure why this was previously in the code...
                # if len(indices_of_neighbors) > 0:
                # In case a point has no neighbors, you need this conditional.
                # one_index_of_neighbor = np.array(indices_of_neighbors)[0]

                # Add that one index to the growing list.
                pocket_indexes = np.hstack(
                    (
                        pocket_indexes,
                        indices_of_neighbors,
                        # one_index_of_neighbor  # Not sure why it used to be this.
                    )
                )

                # Make sure only unique indices ones are retained.
                pocket_indexes = np.unique(pocket_indexes)

            # Save these points (in the pocket) to a list of pockets.
            pocket = self.points[pocket_indexes, :]
            all_pockets.append(pocket)

            # Remove those points from the list of points before trying again.
            self.__delete_limited_points(pocket_indexes)

        # sort the pockets by size
        all_pockets = sorted(all_pockets, key=lambda pts: -len(pts))

        return all_pockets

    def __delete_limited_points(self, pt_indices):
        """A support function"""

        # keep only those points
        self.points = np.delete(self.points, pt_indices, axis=0)

        # update the distance matrix so it doesn't need to be recalculated
        self.dist_matrix = np.delete(self.dist_matrix, pt_indices, axis=0)
        self.dist_matrix = self.dist_matrix.T
        self.dist_matrix = np.delete(self.dist_matrix, pt_indices, axis=0)
