import numpy as np
from scipy.spatial.distance import cdist


class Region:
    """A class for defining regions that will be filled with points."""

    def __init__(self):
        """Initialize some variables."""

        self.center = np.array([9999.9, 9999.9, 9999.9])
        self.radius = 9999.9  # in case the region is a sphere
        # in case the region is a box
        self.box_dimen = np.array([9999.9, 9999.9, 9999.9])

        self.region_type = "SPHERE"  # could also be BOX

    def __str__(self):
        """Returns a string representation of the region."""

        if self.region_type == "SPHERE":
            return (
                "sphere at ("
                + str(self.center[0])
                + ", "
                + str(self.center[1])
                + ", "
                + str(self.center[2])
                + "), radius = "
                + str(self.radius)
            )
        if self.region_type == "BOX":
            return (
                "box centered at ("
                + str(self.center[0])
                + ", "
                + str(self.center[1])
                + ", "
                + str(self.center[2])
                + ") with x,y,z dimensions of ("
                + str(self.box_dimen[0])
                + ", "
                + str(self.box_dimen[1])
                + ", "
                + str(self.box_dimen[2])
                + ")"
            )
        return ""

    def __snap(self, pts, reso):
        """Snaps a set of points to a fixed grid.

        Args:
        pts: A nx3 np.array representing 3D points.
        reso: A float, the resolution of the grid.

        Returns:
        A nx3 np.array with the 3D points snapped to the nearest grid
            point.

        """

        # unfortunately, np.around rounds evenly, so 0.5 rounds to 0.0 and
        # 1.5 rounds to 2.0. very annoying, I'll just add a tiny amount to 0.5
        # => 0.500001 this should work, since user is unlikely to select region
        # center or radius with such precision

        pts = pts + 1e-10
        return np.around(pts / reso) * reso

    def points_set(self, reso):
        """Generates a point field by filling the region with equally spaced
        points.

        Args:
        reso: A float, the resolution of the grid on which the points will
            be placed.

        Returns:
        A nx3 np.array with the 3D points filling the region.

        """

        total_pts = None

        if self.region_type == "BOX":
            xs = np.arange(
                self.center[0] - self.box_dimen[0] / 2,
                self.center[0] + self.box_dimen[0] / 2,
                reso,
            )
            ys = np.arange(
                self.center[1] - self.box_dimen[1] / 2,
                self.center[1] + self.box_dimen[1] / 2,
                reso,
            )
            zs = np.arange(
                self.center[2] - self.box_dimen[2] / 2,
                self.center[2] + self.box_dimen[2] / 2,
                reso,
            )

            total_pts = self._convert_xyz_to_numpy_arr(xs, ys, zs, reso)
        elif self.region_type == "SPHERE":
            total_pts = self._make_sphere_pts(reso)
        return total_pts

    def _make_sphere_pts(self, reso):
        xs = np.arange(self.center[0] - self.radius, self.center[0] + self.radius, reso)
        ys = np.arange(self.center[1] - self.radius, self.center[1] + self.radius, reso)
        zs = np.arange(self.center[2] - self.radius, self.center[2] + self.radius, reso)

        result = self._convert_xyz_to_numpy_arr(xs, ys, zs, reso)
        # now remove all the points outside of this sphere
        index_inside_sphere = np.nonzero(
            cdist(result, np.array([self.center])) < self.radius
        )[0]
        result = result[index_inside_sphere]

        return result

    def _convert_xyz_to_numpy_arr(self, xs, ys, zs, reso):
        result = np.empty((len(xs) * len(ys) * len(zs), 3))

        i = 0
        for x in xs:
            for y in ys:
                for z in zs:
                    result[i][0] = x
                    result[i][1] = y
                    result[i][2] = z

                    i = i + 1

        result = self.__snap(result, reso)

        return result
