import math
import time
from functools import reduce

import numpy as np
import numpy.typing as npt
from loguru import logger
from pymolecule import Molecule
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.spatial import cKDTree

from povme.io import gzopenfile, numpy_to_pdb, openfile, write_to_file
from povme.parallel import RayTaskGeneral


def unique_rows(a):
    """Identifies unique points (rows) in an array of points.

    Args:
        a: A nx3 np.array representing 3D points.

    Returns:
        A nx2 np.array containing the 3D points that are unique.

    """

    a[a == -0.0] = 0.0
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])  # unique_a


class ConvexHull:
    """A class to handle convex-hull calculations.

    `seg_dict` is a dictionary object that contains information about
    segments within the convex hull.

    -   The keys are 2x3 tuples, which represent two ends of a segment in space.
    -   The values of `seg_dict` are the number of times a segment has been part of
        a triangle, either `1` or `2`. `0` times would mean that the segment doesn't
        exist in the dictionary yet.
    """

    def __init__(self, pts: npt.NDArray[np.float64]) -> None:
        """Initializes the ConvexHull class."""

        akl_toussaint_pts = self.akl_toussaint(pts)  # quickly reduces input size
        self.hull = self.gift_wrapping_3d(
            akl_toussaint_pts
        )  # calculate convex hull using gift wrapping algorithm

    def get_seg_dict_num(
        self,
        seg_dict: dict[tuple[tuple[np.float64, ...], ...], int],
        seg_index: tuple[tuple[np.float64, ...], ...],
    ) -> int:
        """This function looks up and returns the value of a seg_index from `seg_dict`.

        Args:
            seg_dict: The dictionary of segment 2x3 tuples as keys, integers as
                values.
            seg_index: the key of the dictionary member we are going to retrieve

        Returns:
            If `seg_index` exists in the keys of `seg_dict`, return the value.
                Otherwise, return 0.

        """
        # we want the index with the greater x-value, so we don't get
        # identical segments in the dictionary more than once
        index = seg_index if seg_index[0][0] > seg_index[1][0] else seg_index[::-1]
        if index in seg_dict:
            return seg_dict[index]
        else:
            return 0

    def increment_seg_dict(
        self,
        seg_dict: dict[tuple[tuple[np.float64, ...], ...], int],
        seg_index: tuple[tuple[np.float64, ...], ...],
    ) -> None:
        """This function increments the values within seg_dict, or initiates them if
        they dont exist yet.

        Args:
            seg_dict: the dictionary of segment 2x3 tuples as keys, integers as values.
            seg_index: the key of the dictionary member we are going to increment.
        """
        # we want the index with the greater x-value, so we don't get
        # identical segments in the dictionary more than once
        if seg_index[0][0] > seg_index[1][0]:
            index = seg_index
        else:
            index = seg_index[::-1]

        # "putting index:", index, "into seg_dict because", index[0][0], ">", index[1][0]

        if index in seg_dict:  # if the entry already exists in seg_dict
            seg_dict[index] += 1  # increment
        else:
            # initiate with a value of 1 because it now exists on a triangle
            seg_dict[index] = 1

    def gift_wrapping_3d(
        self, raw_points: npt.NDArray[np.float64]
    ) -> list[npt.NDArray[np.float64]]:
        """Gift wrapping for 3d convex hull.

        Args:
            raw_points: A nx3 array of points, where each row corresponds to an
                x,y,z point coordinate.

        Returns:
            A convex hull represented by a list of triangles. Each triangle is a
                3x3 array, where each row is an x,y,z coordinate in space. The 3
                rows describe the location of the 3 corners of the triangle. Each
                of the 3 points are arranged so that a cross product will point
                outwards from the hull.

        """

        n = np.shape(raw_points)[0]  # number of points
        point1 = raw_points[0]  # take the first point
        # create a ref vector pointing along x axis
        xaxis = np.array([1, 0, 0])
        maxx = raw_points[0][0]  # initiate highest x value
        points = []  # a list of tuples for easy dictionary lookup
        seg_dict: dict[tuple[tuple[np.float64, ...], ...], int] = (
            {}
        )  # a dictionary that contains the number of triangles a seg is in

        for i in range(n):  # find the n with the largest x value
            point = tuple(raw_points[i])
            points.append(point)
            if point[0] > maxx:
                maxx = point[0]
                point1 = raw_points[i]

        best_dot = -1.0  # initiate dot relative to x-axis
        point2 = np.array(raw_points[1])  # initiate best segment

        # find first/best segment
        for i in range(n):
            pointi = raw_points[i]
            if np.array_equal(pointi, point1):
                continue
            diff_vec = pointi - point1
            diff_len = np.linalg.norm(diff_vec)

            test_dot = np.dot(diff_vec / diff_len, xaxis)
            if test_dot > best_dot:
                best_dot = test_dot
                point2 = pointi

        point1 = tuple(point1)
        point2 = tuple(point2)
        ref_vec = xaxis

        # now find the best triangle
        triangles = []

        seg_list = {(point1, point2)}
        norm_dict = {(point1, point2): xaxis}
        self.increment_seg_dict(seg_dict, (point1, point2))

        counter = 0
        first_time = True

        while (
            seg_list
        ):  # as long as there are unexplored edges of triangles in the hull...

            counter += 1
            seg = seg_list.pop()  # take a segment out of the seg_list
            tuple1 = seg[0]  # the two ends of the segment
            tuple2 = seg[1]
            point1 = np.array(seg[0])
            point2 = np.array(seg[1])
            result = self.get_seg_dict_num(seg_dict, (seg[0], seg[1]))

            if result >= 2:  # then we already have 2 triangles on this segment
                continue  # forget about drawing a triangle for this seg

            # get the norm for a triangle that the segment is part of
            ref_vec = norm_dict[(seg[0], seg[1])]

            best_dot_cross = -1.0
            best_point = None

            for i in range(n):  # look at each point

                pointi = raw_points[i]
                # if np.array_equal(pointi, point1) or np.array_equal(pointi, point2): continue # if we are trying one of the points that are point1 or point2
                diff_vec1 = point2 - point1
                # diff_len1 = np.linalg.norm(diff_vec1)
                diff_vec2 = pointi - point2
                # diff_len2 = np.linalg.norm(diff_vec2)

                # test_cross = np.cross(diff_vec1/diff_len1,diff_vec2/diff_len2)
                # test_cross = np.cross(diff_vec1,diff_vec2)
                test_cross = np.array(
                    [
                        diff_vec1[1] * diff_vec2[2] - diff_vec1[2] * diff_vec2[1],
                        diff_vec1[2] * diff_vec2[0] - diff_vec1[0] * diff_vec2[2],
                        diff_vec1[0] * diff_vec2[1] - diff_vec1[1] * diff_vec2[0],
                    ]
                )  # cross product

                # np.linalg.norm(test_cross) # get the norm of the cross product
                test_cross_len = np.sqrt(
                    test_cross[0] * test_cross[0]
                    + test_cross[1] * test_cross[1]
                    + test_cross[2] * test_cross[2]
                )

                if test_cross_len <= 0.0:
                    continue
                # test_cross_len_inv = 1 / test_cross_len
                test_cross = test_cross / test_cross_len
                dot_cross = np.dot(test_cross, ref_vec)
                # dot_cross = test_cross[0]*ref_vec[0] + test_cross[1]*ref_vec[1] + test_cross[2]*ref_vec[2]
                if dot_cross > best_dot_cross:
                    best_cross = test_cross
                    best_dot_cross = dot_cross
                    best_point = pointi
                    tuple3 = points[i]

            point3 = best_point

            if self.get_seg_dict_num(seg_dict, (tuple2, tuple1)) > 2:
                continue
            if self.get_seg_dict_num(seg_dict, (tuple3, tuple2)) > 2:
                continue
            if self.get_seg_dict_num(seg_dict, (tuple1, tuple3)) > 2:
                continue

            # now we have a triangle from point1 -> point2 -> point3 must test
            # each edge
            if first_time:
                self.increment_seg_dict(seg_dict, (tuple2, tuple1))
                seg_list.add((tuple2, tuple1))
                norm_dict[(tuple2, tuple1)] = best_cross

            self.increment_seg_dict(seg_dict, (tuple3, tuple2))
            seg_list.add((tuple3, tuple2))
            norm_dict[(tuple3, tuple2)] = best_cross

            self.increment_seg_dict(seg_dict, (tuple1, tuple3))
            seg_list.add((tuple1, tuple3))
            norm_dict[(tuple1, tuple3)] = best_cross

            triangles.append((np.array(tuple1), np.array(tuple2), np.array(tuple3)))

            first_time = False

        # print "find all triangles:", time.time() - begintime

        # print "section1:", section1
        # print "section2:", section2
        # print "section3:", section3
        return triangles

    def akl_toussaint(self, points: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """The Akl-Toussaint Heuristic.

        Given a set of points, this definition will create an octahedron whose corners
        are the extremes in x, y, and z directions. Every point within this octahedron
        will be removed because they are not part of the convex hull.
        This causes any expected running time for a convex hull algorithm to be reduced
        to linear time.

        Args:
            points: An nx3 array of x,y,z coordinates.

        Returns:
            All members of original set of points that fall outside the
                Akl-Toussaint octahedron.

        """

        x_high = (-1e99, 0, 0)
        x_low = (1e99, 0, 0)
        y_high = (0, -1e99, 0)
        y_low = (0, 1e99, 0)
        z_high = (0, 0, -1e99)
        z_low = (0, 0, 1e99)

        for point in points:  # find the corners of the octahedron
            if point[0] > x_high[0]:
                x_high = point
            if point[0] < x_low[0]:
                x_low = point
            if point[1] > y_high[1]:
                y_high = point
            if point[1] < y_low[1]:
                y_low = point
            if point[2] > z_high[2]:
                z_high = point
            if point[2] < z_low[2]:
                z_low = point

        octahedron = [  # define the triangles of the surfaces of the octahedron
            np.array((x_high, y_high, z_high)),
            np.array((x_high, z_low, y_high)),
            np.array((x_high, y_low, z_low)),
            np.array((x_high, z_high, y_low)),
            np.array((x_low, y_low, z_high)),
            np.array((x_low, z_low, y_low)),
            np.array((x_low, y_high, z_low)),
            np.array((x_low, z_high, y_high)),
        ]
        new_points = []  # everything outside of the octahedron
        for (
            point
        ) in points:  # now check to see if a point is inside or outside the octahedron
            outside = self.outside_hull(point, octahedron, epsilon=-1.0e-5)
            if outside:
                new_points.append(point)

        return np.array(new_points)  # convert back to an array

    def outside_hull(
        self,
        point: list[float],
        triangles: list[npt.NDArray[np.float64]],
        epsilon: float = 1.0e-5,
    ) -> bool:
        """Given the hull as defined by a list of triangles, this definition
        will return whether a point is within these or not.

        Args:
            point: an x,y,z array that is being tested to see whether it
                exists inside the hull or not.
            triangles: a list of triangles that define the hull.
            epsilon: needed for imprecisions in the floating-point operations.

        Returns:
            True if our_point exists outside of the hull, False otherwise

        """

        our_point = np.array(point)  # convert it to an np.array
        for triangle in triangles:
            # vector from triangle corner 0 to point
            rel_point = our_point - triangle[0]
            # vector from triangle corner 0 to corner 1
            vec1 = triangle[1] - triangle[0]
            # vector from triangle corner 1 to corner 2
            vec2 = triangle[2] - triangle[1]
            # cross product between vec1 and vec2
            our_cross = np.cross(vec1, vec2)
            # dot product to determine whether cross is point inward or
            # outward
            our_dot = np.dot(rel_point, our_cross)
            # if the dot is greater than 0, then its outside
            if np.dot(rel_point, our_cross) > epsilon:
                return True

        return False

    def inside_hull(self, our_point):
        """Determines if a point is inside the hull

        Args:
            our_point: An x,y,z array

        Returns:
            A boolean, True if the point is inside the hull, False otherwise.

        """

        return not self.outside_hull(our_point, self.hull)

    @staticmethod
    def volume(frame_indx, pdb, pts, regions_contig, output_prefix, config):
        # if the user wants to save empty points (points that are removed),
        # then we need a copy of the original
        if config.output_equal_num_points_per_frame:
            pts_deleted = pts.copy()

        # you may need to load it from disk if the user so specified
        if config.use_disk_not_memory:  # so you need to load it from disk
            pym_filename = pdb
            pdb = Molecule()
            pdb.io.load_pym_into(pym_filename)

        # remove the points that are far from the points region anyway
        min_pts = np.min(pts, 0) - config.distance_cutoff - 1
        max_pts = np.max(pts, 0) + config.distance_cutoff + 1

        # identify atoms that are so far away from points that they can be
        coords = pdb.information.coordinates
        x_keep = (coords[:, 0] >= min_pts[0]) & (coords[:, 0] <= max_pts[0])
        y_keep = (coords[:, 1] >= min_pts[1]) & (coords[:, 1] <= max_pts[1])
        z_keep = (coords[:, 2] >= min_pts[2]) & (coords[:, 2] <= max_pts[2])
        keep_indices = np.where(x_keep & y_keep & z_keep)[0]

        # keep only relevant atoms
        if len(keep_indices) > 0:
            pdb = pdb.selections.create_molecule_from_selection(keep_indices)
            coords = pdb.information.coordinates
        else:
            # No atoms to consider, so all points remain if no hull is needed
            # Just calculate the volume and exit early if needed
            volume = len(pts) * (config.grid_spacing**3)
            return (frame_indx, volume, {})

        # get the vdw radii of each protein atom
        vdw = np.ones(len(coords))
        element_stripped = pdb.information.atom_information["element_stripped"]
        vdw[element_stripped == b"H"] = 1.2
        vdw[element_stripped == b"C"] = 1.7
        vdw[element_stripped == b"N"] = 1.55
        vdw[element_stripped == b"O"] = 1.52
        vdw[element_stripped == b"F"] = 1.47
        vdw[element_stripped == b"P"] = 1.8
        vdw[element_stripped == b"S"] = 1.8

        # Create a KD-tree for atom coordinates
        atom_tree = cKDTree(coords)

        # Compute the maximum cutoff for searching
        # Since each atom might have a different vdw radius, we must take the max
        max_cutoff = np.max(vdw) + config.distance_cutoff

        # Query nearby atoms for each point using the KD-tree
        # This returns a list of arrays, where each entry corresponds to a point
        # and contains the indices of atoms within max_cutoff radius.
        candidate_indices = atom_tree.query_ball_point(pts, r=max_cutoff)

        # Now filter points that are too close to at least one atom
        close_pt_indices = []
        for i, atom_inds in enumerate(candidate_indices):
            if len(atom_inds) == 0:
                # No nearby atoms, so this point is safe
                continue
            # Distances to these candidate atoms
            candidate_coords = coords[atom_inds]
            dist_subset = np.sqrt(np.sum((candidate_coords - pts[i])**2, axis=1))

            # Each atom has a different cutoff: vdw[atom_ind] + config.distance_cutoff
            atom_cutoffs = vdw[atom_inds] + config.distance_cutoff
            if np.any(dist_subset < atom_cutoffs):
                close_pt_indices.append(i)

        close_pt_indices = np.array(close_pt_indices, dtype=int)

        # now keep the appropriate points
        pts = np.delete(pts, close_pt_indices, axis=0)

        # exclude points outside convex hull
        if config.convex_hull_exclusion:
            convex_hull_3d = ConvexHull(pts)

            # get the coordinates of the non-hydrogen atoms (faster to discard
            # hydrogens)
            hydros = pdb.selections.select_atoms({"element_stripped": [b"H"]})
            not_hydros = pdb.selections.invert_selection(hydros)
            not_hydros_coors = pdb.information.coordinates[not_hydros]

            # not_hydros = pdb.selections.select_atoms({'name_stripped':['CA']})
            # not_hydros_coors = pdb.information.coordinates[not_hydros]

            # modify pts here.
            # note that the atoms of the pdb frame are in pdb.information.coordinates
            # begintime = time.time() # measure execution time
            akl_toussaint_pts = convex_hull_3d.akl_toussaint(
                not_hydros_coors
            )  # quickly reduces input size
            # print "akl Toussaint:", time.time() - begintime
            begintime = time.time()  # measure execution time
            # calculate convex hull using gift wrapping algorithm
            hull = convex_hull_3d.gift_wrapping_3d(akl_toussaint_pts)
            # print "gift_wrapping:", time.time() - begintime

            # we will need to regenerate the pts list, disregarding those
            # outside the hull
            old_pts = pts
            pts = []
            for pt in old_pts:
                pt_outside = convex_hull_3d.outside_hull(
                    pt, hull
                )  # check if pt is outside hull
                if not pt_outside:
                    # if its not outside the hull, then include it in the
                    # volume measurement
                    pts.append(pt)
            pts = np.array(pts)

        # Now, enforce contiguity if needed
        if len(regions_contig) > 0 and len(pts) > 0:
            # first, for each point, determine how many neighbors it has to
            # count kiddy-corner points too
            cutoff_dist = config.grid_spacing * 1.01 * math.sqrt(3)
            pts_dists = squareform(pdist(pts))
            # minus 1 because an atom shouldn't be considered its own neighor
            neighbor_counts = np.sum(pts_dists < cutoff_dist, axis=0) - 1

            # remove all the points that don't have enough neighbors
            pts = pts[
                np.nonzero(neighbor_counts >= config.contiguous_points_criteria)[0]
            ]

            # get all the points in the defined parameters['ContiguousPocket']
            # seed regions
            contig_pts = regions_contig[0].get_points(config.grid_spacing)
            for Contig in regions_contig[1:]:
                contig_pts = np.vstack(
                    (contig_pts, Contig.get_points(config.grid_spacing))
                )
            contig_pts = unique_rows(contig_pts)

            try:  # error here if there are no points of contiguous seed region outside of protein volume.
                # now just get the ones that are not near the protein
                contig_pts = pts[np.nonzero(cdist(contig_pts, pts) < 1e-7)[1]]

                last_size_of_contig_pts = 0
                while last_size_of_contig_pts != len(contig_pts):
                    last_size_of_contig_pts = len(contig_pts)

                    # now get the indices of all points that are close to the
                    # contig_pts
                    all_pts_close_to_contig_pts_boolean = (
                        cdist(pts, contig_pts) < cutoff_dist
                    )
                    index_all_pts_close_to_contig_pts = np.unique(
                        np.nonzero(all_pts_close_to_contig_pts_boolean)[0]
                    )
                    contig_pts = pts[index_all_pts_close_to_contig_pts]

                pts = contig_pts
            except Exception:
                logger.exception(
                    "Frame "
                    + str(frame_indx)
                    + ": None of the points in the contiguous-pocket seed region\n\t\tare outside the volume of the protein! Assuming a pocket\n\t\tvolume of 0.0 A."
                )
                pts = np.array([])

        # now write the pdb and calculate the volume
        volume = len(pts) * math.pow(config.grid_spacing, 3)

        logger.info("\tFrame " + str(frame_indx) + ": " + repr(volume) + " A^3")
        if config.save_individual_pocket_volumes:
            frame_text = f"REMARK Frame {str(frame_indx)}" + "\n"
            frame_text += f"REMARK Volume = {repr(volume)}" + " Cubic Angstroms\n"
            frame_text += numpy_to_pdb(pts, "X")

            if config.output_equal_num_points_per_frame:
                # you need to find the points that are in pts_deleted but not
                # in pts
                tmp = reduce(
                    lambda x, y: x | np.all(pts_deleted == y, axis=-1),
                    pts,
                    np.zeros(pts_deleted.shape[:1], dtype=np.bool),
                )
                indices = np.where(tmp)[0]
                pts_deleted = np.delete(pts_deleted, indices, axis=0)

                # So extra points will always be at the origin. These can be
                # easily hidden with your visualization software.
                pts_deleted = np.zeros(pts_deleted.shape)
                frame_text = frame_text + numpy_to_pdb(pts_deleted, "X", "XXX")

            frame_text = frame_text + "END\n"

            if config.compress_output:
                fl = gzopenfile(
                    output_prefix + "frame_" + str(frame_indx) + ".pdb.gz",
                    "wb",
                )
            else:
                fl = openfile(
                    output_prefix + "frame_" + str(frame_indx) + ".pdb",
                    "w",
                )
            write_to_file(fl, frame_text, encode=config.compress_output)
            fl.close()

        extra_data_to_add = {}
        if config.save_volumetric_density_map:
            extra_data_to_add["SaveVolumetricDensityMap"] = pts

        return (frame_indx, volume, extra_data_to_add)


class TaskCalcVolume(RayTaskGeneral):
    """A class for calculating the volume."""

    def process_item(self, item):
        """Calculate the volume.

        Args:
            item: A list or tuple containing necessary input data.

        Returns:
            A tuple with frame index, calculated volume, and any extra data.
        """
        try:

            return ConvexHull.volume(*item)
        except Exception as e:
            logger.exception(f"Error in frame {item[0]}: {e}")
            return ("error", str(e))
