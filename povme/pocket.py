from typing import Any

import np.typing as npt
import numpy as np
from scipy import spatial

from .io import write_pdbs
from .parallel import GeneralTask, MultiThreading

params: dict[str, Any] = {}


class BoxOfPoints:
    """A class representing a box of equidistant points."""

    def __init__(self, box: npt.NDArray[np.float64], reso: int) -> None:
        """Initialize the class.

        Args:
            box: A numpy array representing two 3D points, (min_x, min_y,
                min_z) and (max_x, max_y, max_z), that define a box.
            reso: The space between the points of the box, in the X, Y, and
                Z direction.
        """

        self.write_pdbs = write_pdbs()

        min_x = self.__snap_float(box[0][0], reso)
        min_y = self.__snap_float(box[0][1], reso)
        min_z = self.__snap_float(box[0][2], reso)
        max_x = self.__snap_float(box[1][0], reso) + 1.1 * reso
        max_y = self.__snap_float(box[1][1], reso) + 1.1 * reso
        max_z = self.__snap_float(box[1][2], reso) + 1.1 * reso

        x, y, z = np.mgrid[min_x:max_x:reso, min_y:max_y:reso, min_z:max_z:reso]
        self.points = np.array(list(zip(x.ravel(), y.ravel(), z.ravel())))

    def __snap_float(
        self, val: npt.NDArray[np.float64], reso: int
    ) -> npt.NDArray[np.float64]:
        """Snaps an arbitrary point to the nearest grid point.

        Args:
            val: A numpy array corresponding to a 3D point.
            reso: The resolution (distance in the X, Y, and Z directions
                between adjacent points) of the grid.

        Returns:
            A numpy array corresponding to a 3D point near val that is on a
                nearby grid point.
        """

        return np.floor(val / reso) * reso

    def remove_points_outside_convex_hull(self, hull):
        """Removes box points that are outside a convex hull.

        Args:
            hull: The convex hull.
        """
        chunks = [(hull, t) for t in np.array_split(self.points, params["processors"])]
        tmp = MultiThreading(chunks, params["processors"], self.__MultiIdHullPts)
        self.points = np.vstack(tmp.results)

    class __MultiIdHullPts(GeneralTask):
        """A class to remove points outside a convex hull using multiple
        processors."""

        def value_func(self, items, results_queue):  # so overwriting this function
            """The calculations that will run on a single processor to remove
            points outside a convex hull."""

            hull = items[0]
            some_points = items[1]

            # Note this would be much faster if it were matrix-based instead of
            # point-by-point based. Can preallocate numpy array size because I
            # don't know beforehand how many points will be in the hull
            new_pts = []
            for pt in some_points:
                if hull.inside_hull(pt):
                    new_pts.append(pt)

            if len(new_pts) == 0:
                pass  # here save the results for later compilation
            else:
                self.results.append(np.array(new_pts))

    def remove_all_points_close_to_other_points(
        self, other_points: npt.NDArray[np.float64], dist_cutoff: float
    ) -> None:
        """Removes all points in this box that come within the points specified
        in a numpy array

        Args:
            other_points: A numpy array containing the other points.
            dist_cutoff: A float, the cutoff distance to use in determining
                whether or not box points will be removed.

        """

        # note, in newer versions of scipy use cKDTree
        box_of_pts_distance_tree = spatial.KDTree(self.points)
        chunks = [
            (box_of_pts_distance_tree, dist_cutoff, t)
            for t in np.array_split(other_points, params["processors"])
        ]
        tmp = MultiThreading(chunks, params["processors"], self.__MultiGetClosePoints)
        indices_of_box_pts_close_to_molecule_points = np.unique(np.hstack(tmp.results))

        self.points = np.delete(
            self.points, indices_of_box_pts_close_to_molecule_points, axis=0
        )  # remove the ones that are too close to molecule atoms

    class __MultiGetClosePoints(GeneralTask):
        """A class to remove box points that are near other, user-specified
        points, using multiple processors."""

        def value_func(self, items, results_queue):  # so overwriting this function
            """The calculations that will run on a single processor."""

            box_of_pts_distance_tree = items[0]
            dist_cutoff = items[1]
            other_points = items[2]

            # note, in newer versions of scipy use cKDTree
            other_points_distance_tree = spatial.KDTree(other_points)
            sparce_distance_matrix = other_points_distance_tree.sparse_distance_matrix(
                box_of_pts_distance_tree, dist_cutoff
            )
            indices_of_box_pts_close_to_molecule_points = np.unique(
                sparce_distance_matrix.tocsr().indices
            )

            self.results.append(indices_of_box_pts_close_to_molecule_points)

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
            reso: The distance between adjacent added points.

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
            reso: The distance between adjacent points.
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
