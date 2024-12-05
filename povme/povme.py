from typing import Any

import os
import shutil
import sys
import time
from io import StringIO

import numpy as np
from loguru import logger
from pymolecule import Molecule
from scipy.spatial.distance import cdist

from .config import POVMEConfig
from .hull import MultithreadingCalcVolumeTask
from .io import dx_freq, gzopenfile, numpy_to_pdb, openfile, write_to_file
from .parallel import MultiThreading, MultithreadingTaskGeneral
from .points.regions import collect_regions


def get_unique_rows(a):
    """Identifies unique points (rows) in an array of points.

    Arguments:
        a: A nx3 numpy.array representing 3D points.

    Returns:
        A nx2 numpy.array containing the 3D points that are unique.

    """

    a[a == -0.0] = 0.0
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])  # unique_a


class MultithreadingStringToMoleculeTask(MultithreadingTaskGeneral):
    """A class for loading PDB frames (as strings) into pymolecule.Molecule
    objects."""

    def value_func(self, item, results_queue):
        """Convert a PDB string into a pymolecule.Molecule object.

        Args:
            item: A list or tuple, the input data required for the calculation.
            results_queue: A multiprocessing.Queue() object for storing the
                calculation output.

        """

        pdb_string = item[0]
        index = item[1]
        config = item[2]

        # make the pdb object
        str_obj = StringIO(pdb_string)
        tmp = Molecule()
        tmp.io.load_pdb_into_using_file_object(str_obj, False, False, False)

        logger.debug("\tFurther processing frame " + str(index))

        # so load the whole trajectory into memory
        if not config.use_disk_not_memory:
            self.results.append((index, tmp))
        else:  # save to disk, record filename
            pym_filename = f"./.povme_tmp/frame_{str(index)}.pym"
            tmp.io.save_pym(pym_filename, False, False, False, False, False)
            self.results.append((index, pym_filename))


class POVME:
    """The main class to run POVME."""

    def __init__(
        self,
        path_config: str | None = None,
    ) -> None:
        """Initialize POVME.

        Args:
            path_config: Path to a configuration YAML file.
        """

        self.config = POVMEConfig()
        if path_config is not None:
            self.config.from_yaml(path_config)

    def load_multi_frame_pdb(self, filename, config):
        """Load a multi-frame PDB into memory or into separate files
        (depending on user specifications).

        Args:
            filename: A string, the filename of the multi-frame PDB.
            config: POVMEConfig object.

        Returns:
            If the user has requested that the disk be used to save memory, this
                function returns a list of tuples, where the first item in each
                tuple is the frame index, and the second is a filename containing
                the individual frame. If memory is to be used instead of the disk,
                this function returns a list of tuples, where the first item in
                each tuple is the frame index, and the second is a
                pymolecule.Molecule object representing the frame.

        """
        pdb_strings = []
        growing_string = ""

        logger.info("Reading frames from " + filename)

        f = open(filename, "rb")
        while True:

            line = f.readline()

            if len(line) == 0:
                pdb_strings.append(growing_string)
                break
            if line[:3] == b"END":
                pdb_strings.append(growing_string)
                growing_string = ""
            else:
                growing_string = growing_string + line.decode()

        f.close()

        while "" in pdb_strings:
            pdb_strings.remove("")

        # now convert each pdb string into a pymolecule.Molecule object
        molecules = MultiThreading(
            [(pdb_strings[idx], idx + 1, config) for idx in range(len(pdb_strings))],
            config.num_processors,
            MultithreadingStringToMoleculeTask,
        )
        molecules = molecules.results

        return molecules

    def gen_points(self, config):
        logger.info("Generating the pocket-encompassing point field")

        # get all the points of the inclusion regions
        regions_include = collect_regions(
            config.points_inclusion_sphere, config.points_inclusion_box
        )
        pts = regions_include[0].get_points(config.grid_spacing)
        for Included in regions_include[1:]:
            pts = np.vstack((pts, Included.get_points(config.grid_spacing)))
        pts = get_unique_rows(pts)

        # get all the points of the exclusion regions
        regions_exclude = collect_regions(
            config.points_exclusion_sphere, config.points_exclusion_box
        )
        if len(regions_exclude) > 0:
            pts_exclusion = regions_exclude[0].get_points(config.grid_spacing)
            for Excluded in regions_exclude[1:]:
                pts_exclusion = np.vstack(
                    (pts_exclusion, Excluded.get_points(config.grid_spacing))
                )
            pts_exclusion = get_unique_rows(pts_exclusion)

            # remove the exclusion points from the inclusion points I
            # think there ought to be a set-based way of doing this, but
            # I'm going to go for the pairwise comparison. consider
            # rewriting later
            index_to_remove = np.nonzero(cdist(pts, pts_exclusion) < 1e-7)[0]
            pts = np.delete(pts, index_to_remove, axis=0)

        return pts

    @staticmethod
    def write_points(pts, output_prefix, config):
        logger.info("Writing points to PDB file")
        points_filename = output_prefix + "point_field.pdb"

        if config.compress_output:
            afile = gzopenfile(points_filename + ".gz", "wb")
        else:
            afile = openfile(points_filename, "w")

        write_to_file(afile, numpy_to_pdb(pts, "X"), encode=config.compress_output)
        afile.close()

        # save the points as npy
        logger.info("Writing points to NPY file")
        np.save(points_filename + ".npy", pts)

        logger.info(
            "Point field saved to " + points_filename + " to permit visualization"
        )
        logger.info(
            "Point field saved to "
            + points_filename
            + ".npy to optionally load for the volume calculation"
        )

    def write_points_contig(self, regions_contig, output_prefix, config):

        # get all the contiguous points
        contig_pts = regions_contig[0].get_points(config.grid_spacing)
        for Contig in regions_contig[1:]:
            contig_pts = np.vstack((contig_pts, Contig.get_points(config.grid_spacing)))
        contig_pts = get_unique_rows(contig_pts)

        logger.info("\nSaving the contiguous-pocket seed points as a PDB file")

        points_filename = output_prefix + "contiguous_pocket_seed_points.pdb"

        if config.compress_output:
            afile = gzopenfile(points_filename + ".gz", "wb")
        else:
            afile = openfile(points_filename, "w")

        write_to_file(
            afile,
            numpy_to_pdb(contig_pts, "X"),
            encode=config.compress_output,
        )
        afile.close()

        logger.info(
            "Contiguous-pocket seed points saved to "
            + points_filename
            + " to permit visualization"
        )

    def compute_volume(self, path_pdb, pts, regions_contig, output_prefix, config):
        # load the PDB frames
        index_and_pdbs = self.load_multi_frame_pdb(path_pdb, config)

        # calculate all the volumes
        logger.info("Calculating the pocket volume of each frame")
        tmp = MultiThreading(
            [
                (index, pdb_object, pts, regions_contig, output_prefix, config)
                for index, pdb_object in index_and_pdbs
            ],
            config.num_processors,
            MultithreadingCalcVolumeTask,
        )

        # delete the temp swap directory if necessary
        if config.use_disk_not_memory:
            if os.path.exists("./.povme_tmp"):
                shutil.rmtree("./.povme_tmp")
        logger.info("Execution time = " + str(time.time() - self.t_start) + " sec")
        return tmp.results

    @staticmethod
    def write_vol_csv(results_vol, output_prefix, config):
        if config.compress_output:
            f = gzopenfile(
                output_prefix + "volumes.csv.gz",
                "wb",
            )
        else:
            f = openfile(output_prefix + "volumes.csv", "w")

        write_to_file(
            f,
            "frame_idx,volume\n",
            encode=config.compress_output,
        )
        for i in sorted(results_vol.keys()):
            write_to_file(
                f,
                str(i) + "," + str(results_vol[i]) + "\n",
                encode=config.compress_output,
            )
        f.close()

    @staticmethod
    def write_vol_traj(results_vol, output_prefix, config):
        if config.compress_output:
            traj_file = gzopenfile(
                output_prefix + "volume_trajectory.pdb.gz",
                "wb",
            )
        else:
            traj_file = openfile(
                output_prefix + "volume_trajectory.pdb",
                "w",
            )

        for frame_index in range(1, len(list(results_vol.keys())) + 1):
            if config.compress_output:
                frame_file = gzopenfile(
                    output_prefix + "frame_" + str(frame_index) + ".pdb.gz",
                    "rb",
                )
            else:
                frame_file = openfile(
                    output_prefix + "frame_" + str(frame_index) + ".pdb",
                    "r",
                )

            traj_file.write(frame_file.read())
            frame_file.close()

        traj_file.close()

    @staticmethod
    def write_vol_dens(results, output_prefix, config):
        unique_points: dict[str, Any] = {}

        overall_min = np.ones(3) * 1e100
        overall_max = np.ones(3) * -1e100

        for result in results:
            pts = result[2]["SaveVolumetricDensityMap"]

            if len(pts) > 0:
                amin = np.min(pts, axis=0)
                amax = np.max(pts, axis=0)

                overall_min = np.min(np.vstack((overall_min, amin)), axis=0)
                overall_max = np.max(np.vstack((overall_max, amax)), axis=0)

                for pt in pts:
                    pt_key = str(pt[0]) + ";" + str(pt[1]) + ";" + str(pt[2])
                    try:
                        unique_points[pt_key] = unique_points[pt_key] + 1
                    except Exception:
                        unique_points[pt_key] = 1
        if overall_min[0] == 1e100:
            logger.info(
                "ERROR! Cannot save volumetric density file because no volumes present in any frame.",
            )
        else:
            xpts = np.arange(
                overall_min[0],
                overall_max[0] + config.grid_spacing,
                config.grid_spacing,
            )
            ypts = np.arange(
                overall_min[1],
                overall_max[1] + config.grid_spacing,
                config.grid_spacing,
            )
            zpts = np.arange(
                overall_min[2],
                overall_max[2] + config.grid_spacing,
                config.grid_spacing,
            )

            all_pts = np.zeros((len(xpts) * len(ypts) * len(zpts), 4))

            i = 0
            for x in xpts:
                for y in ypts:
                    for z in zpts:
                        key = str(x) + ";" + str(y) + ";" + str(z)
                        all_pts[i][0] = x
                        all_pts[i][1] = y
                        all_pts[i][2] = z

                        try:
                            all_pts[i][3] = unique_points[key]
                        except Exception:
                            pass

                        i = i + 1

            # convert the counts in the fourth column into frequencies
            all_pts[:, 3] = all_pts[:, 3] / len(results)
            dx_freq(all_pts, output_prefix, config)  # save the dx file

    def run(
        self,
        path_pdb: str,
        output_prefix: str | None = None,
    ) -> dict[str, Any]:
        """Start POVME

        Args:
            path_pdb: Path to PDB file. This will overwrite the configuration file.
            output_prefix: Path to output directory including directories.
        """
        self.t_start = time.time()
        config = self.config
        config.log()

        if output_prefix is None:
            output_prefix = (
                "POVME_output."
                + time.strftime("%m-%d-%y")
                + "."
                + time.strftime("%H-%M-%S")
                + "/"
            )

        # If the output prefix includes a directory, create that directory if
        # necessary
        if "/" in output_prefix:
            output_dirname = os.path.dirname(output_prefix)
            os.makedirs(output_dirname, exist_ok=True)

        # create temp swap directory if needed
        if config.use_disk_not_memory:
            if os.path.exists("./.povme_tmp"):
                shutil.rmtree("./.povme_tmp")
            os.mkdir("./.povme_tmp")

        # User specified regions to include in the PDB structure.
        # Thus, we compute these points.
        pts = None
        if config.load_points_path is not None:
            if not os.path.exists(config.load_points_path):
                logger.error(
                    f"points file at {config.load_points_path} does not exist!"
                )
            pts = np.load(config.load_points_path)
        elif (len(config.points_inclusion_box) > 0) or (
            len(config.points_inclusion_sphere) > 0
        ):
            pts = self.gen_points(config)
            if config.save_points:
                self.write_points(pts, output_prefix, config)
        if pts is None:
            logger.error("No points are specified!")
            logger.error("Please specify inclusions or load points from a NumPy file.")
            sys.exit(0)

        # Handle contig TODO:
        regions_contig = collect_regions(
            config.contiguous_pocket_seed_sphere, config.contiguous_pocket_seed_box
        )
        if len(regions_contig) > 0:
            self.write_points_contig(regions_contig, output_prefix, config)

        # Compute volumes of frames in PDB file.
        if not os.path.exists(path_pdb):
            logger.error(f"PDB file {path_pdb} does not exits!")
            sys.exit(0)
        results = self.compute_volume(
            path_pdb, pts, regions_contig, output_prefix, config
        )
        results_vol = {}
        for result in results:
            results_vol[result[0]] = result[1]

        # Save volumes to CSV file
        self.write_vol_csv(results_vol, output_prefix, config)

        # if the user wanted a single trajectory containing all the
        # volumes, generate that here.
        if config.save_pocket_volumes_trajectory:
            self.write_vol_traj(results_vol, output_prefix, config)

        # if the user requested a volumetric density map, then generate it here
        if config.save_volumetric_density_map:
            self.write_vol_dens(results, output_prefix, config)

        return results_vol
