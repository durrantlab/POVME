from typing import Any, Generator

import csv
import os
import shutil
import sys
import time
from io import StringIO

import numpy as np
from loguru import logger
from pymolecule import Molecule
from scipy.spatial.distance import cdist

from povme.config import POVMEConfig
from povme.io import dx_freq, gzopenfile, numpy_to_pdb, openfile, write_to_file
from povme.parallel import RayManager, RayTaskGeneral
from povme.pocket.savers import init_vol_csv, write_vol_csv
from povme.points.hull import ConvexHull
from povme.points.regions import collect_regions


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


def load_multi_frame_pdb_generator(
    filename: str,
) -> Generator[tuple[str, int], None, None]:
    """Generator to yield PDB frames one by one.

    Args:
        filename: Path to the multi-frame PDB file.

    Yields:
        A tuple containing the PDB string and its frame index.
    """
    with open(filename, "r") as f:
        frame_idx = 1
        pdb_lines: list[str] = []
        for line in f:
            if line.startswith("END"):
                if pdb_lines:
                    pdb_string = "".join(pdb_lines)
                    yield (pdb_string, frame_idx)
                    frame_idx += 1
                    pdb_lines = []
            else:
                pdb_lines.append(line)
        # Yield the last frame if the file doesn't end with 'END'
        if pdb_lines:
            pdb_string = "".join(pdb_lines)
            yield (pdb_string, frame_idx)


def collect_pdb_frames_in_chunks(
    filename: str, chunk_size: int
) -> Generator[list[tuple[int, str]], None, None]:
    """
    Read a multi-frame PDB and yield frames in chunks.

    Each yielded chunk is a list of (frame_index, pdb_frame_string).
    """
    frame_buffer: list[str] = []
    frame_index = 0
    chunk = []

    with open(filename, "rb") as f:
        for line in f:
            if line.startswith(b"END"):
                # A frame ended
                if frame_buffer:
                    frame_index += 1
                    chunk.append((frame_index, "".join(frame_buffer)))
                    frame_buffer = []

                    # If we have reached the chunk_size, yield it
                    if len(chunk) == chunk_size:
                        yield chunk
                        chunk = []
            else:
                frame_buffer.append(line.decode())

        # If file does not end with END and we still have a frame collected
        if frame_buffer:
            frame_index += 1
            chunk.append((frame_index, "".join(frame_buffer)))

        # Yield any remaining frames if they don't fill an entire chunk
        if chunk:
            yield chunk


class TaskComputeVolumeFromPDBLines(RayTaskGeneral):

    def process_item(self, item: tuple[Any, ...]) -> tuple[Any, ...]:
        frame_index, pdb_string, config, pts, regions_contig, output_prefix = item

        # Load the PDB from lines
        str_obj = StringIO(pdb_string)
        pdb = Molecule()
        pdb.io.load_pdb_into_using_file_object(str_obj, False, False, False)

        # From here, do the volume calculation steps (adapted from TaskCalcVolume)
        try:
            volumes = ConvexHull.volume(
                frame_index, pdb, pts, regions_contig, output_prefix, config
            )
            return volumes
        except Exception as e:
            logger.exception(f"Error in frame {frame_index}: {e}")
            return ("error", str(e))


class PocketVolume:
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

        logger.info("Saving the contiguous-pocket seed points as a PDB file")

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
        self, path_pdb: str, output_prefix: str | None = None, chunk_size: int = 10
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

        # Initialize RayManager
        ray_manager = RayManager(
            task_class=TaskComputeVolumeFromPDBLines, n_cores=config.n_cores, use_ray=config.use_ray
        )
        init_vol_csv(output_prefix)

        # Collect frames in chunks and submit tasks to RayManager
        for chunk in collect_pdb_frames_in_chunks(path_pdb, chunk_size):
            # Each chunk is [(frame_index, pdb_string), ...]
            tasks = []
            for frame_index, pdb_string in chunk:
                tasks.append(
                    (
                        frame_index,
                        pdb_string,
                        config,
                        pts,
                        regions_contig,
                        output_prefix,
                    )
                )
            ray_manager.submit_tasks(
                tasks,
                chunk_size=len(tasks),  # submit the chunk at once
                save_func=write_vol_csv,  # save intermediate results to CSV
                save_kwargs={"output_prefix": output_prefix},
                save_interval=chunk_size,  # save after every chunk
            )

        results = ray_manager.get_results()
        if len(results) == 0:
            raise RuntimeError("No volume results obtained.")

        # Process final results
        results_vol = {r[0]: r[1] for r in results if r[0] != "error"}

        # if the user wanted a single trajectory containing all the
        # volumes, generate that here.
        if config.save_pocket_volumes_trajectory:
            self.write_vol_traj(results_vol, output_prefix, config)

        # if the user requested a volumetric density map, then generate it here
        if config.save_volumetric_density_map:
            self.write_vol_dens(results, output_prefix, config)

        return results_vol
