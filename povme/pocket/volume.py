"""Pocket volume calculation for single structures and MD trajectories.

This module provides the [`PocketVolume`][pocket.volume.PocketVolume] class, which orchestrates
the end-to-end POVME volume-measurement workflow:

1. Generate (or load) the pocket-encompassing point field from
   user-defined inclusion/exclusion regions.
2. Distribute per-frame volume calculations across workers using
    [`RayManager`][parallel.RayManager].
3. Collect results and optionally write trajectory PDBs and volumetric
    density maps.
"""

import os
import sys
import time
from io import StringIO
from typing import Any, Generator

import numpy as np
from loguru import logger
from pymolecule import Molecule

from povme.config import PocketVolumeConfig
from povme.io import dx_freq, gzopenfile, numpy_to_pdb, openfile, write_to_file
from povme.parallel import RayManager, RayTaskGeneral
from povme.pocket.savers import init_vol_csv, write_vol_csv
from povme.points.hull import ConvexHull
from povme.points.regions import collect_regions


def get_unique_rows(a):
    """Return the unique rows of a 2-D array.

    Each row is treated as an opaque byte string so that
    [`numpy.unique`](https://numpy.org/doc/stable/reference/generated/numpy.unique.html)
    can identify duplicates without floating-point
    tolerance issues (the input is assumed to be snapped to a grid).

    Args:
        a: An array of shape `(n, d)` (typically `d = 3`).

    Returns:
        An array of shape `(m, d)` (`m <= n`) containing only the
        unique rows, in sorted order.
    """
    a[a == -0.0] = 0.0
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])  # unique_a


def remove_exclusion_points(
    pts: np.ndarray,
    pts_exclusion: np.ndarray,
) -> np.ndarray:
    """Remove inclusion points that coincide with exclusion points.

    Both point sets are assumed to lie on the same regular grid (i.e., their
    coordinates are exact multiples of the grid spacing after snapping).
    Instead of computing a full pairwise distance matrix with
    [`cdist`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html),
    this function performs a set-difference using NumPy's structured-array view trick.

    1. Each `(x, y, z)` row is reinterpreted as a single opaque
        [`np.void`](https://numpy.org/doc/stable/reference/arrays.scalars.html) element.
    2. [`numpy.isin`](https://numpy.org/doc/stable/reference/generated/numpy.isin.html)
        performs a hash-based membership test.

    Args:
        pts: An `(N, 3)` array of inclusion grid points.
        pts_exclusion: An `(M, 3)` array of exclusion grid points.

    Returns:
        An `(K, 3)` array (`K <= N`) containing only those rows of
            `pts` that do not appear in *pts_exclusion*.
    """

    def _to_void(a: np.ndarray) -> np.ndarray:
        """View each row of a 2-D array as a single void element."""
        a = np.ascontiguousarray(a)
        return a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))

    pts_void = _to_void(pts)
    exc_void = _to_void(pts_exclusion)
    mask = ~np.isin(pts_void, exc_void).ravel()
    return pts[mask]


def collect_pdb_frames_in_chunks(
    filename: str, chunk_size: int
) -> Generator[list[tuple[int, str]], None, None]:
    """Read a multi-frame PDB and yield frames in chunks.

    Frames are delimited by lines starting with `END`. Frames are
    accumulated into chunks of `chunk_size` before being yielded, reducing
    the overhead of task submission when using parallel workers.

    Args:
        filename: Path to the multi-frame PDB file.
        chunk_size: Maximum number of frames per yielded chunk.

    Yields:
        Lists of `(frame_index, pdb_frame_string)` tuples. The
            `frame_index` is 1-based. The final yielded chunk may contain
            fewer than *chunk_size* frames.
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
        logger.debug(f"Computing volume from PDB lines on frame {frame_index}")
        str_obj = StringIO(pdb_string)
        pdb = Molecule()
        pdb.io.load_pdb_into_using_file_object(str_obj, False, False, False)
        coord_shape = pdb.information.get_coordinates().shape
        if len(coord_shape) < 3:
            logger.debug("Loaded 1 structure")
        else:
            logger.debug(f"Loaded {coord_shape[0]} structure")
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
        config: PocketVolumeConfig | None = None,
    ) -> None:
        """Initialize POVME.

        Args:
            config: Pocket volume calculation configuration.
        """
        if config is None:
            config = PocketVolumeConfig()
        self.config = config

    def gen_points(self, config):
        """Generate the pocket-encompassing point field.

        Constructs a regular grid of 3D points by unioning all
        inclusion-region grids and then subtracting all exclusion-region
        grids.

        Args:
            config: The volume-calculation configuration, which specifies
                inclusion/exclusion spheres and boxes and the grid spacing.

        Returns:
            An `(K, 3)` array of unique grid points that lie inside at
                least one inclusion region and outside all exclusion regions.
        """
        logger.info("Generating the pocket-encompassing point field")

        # Collect inclusion points
        regions_include = collect_regions(
            config.points_inclusion_sphere, config.points_inclusion_box
        )
        pts = regions_include[0].get_points(config.grid_spacing)
        for Included in regions_include[1:]:
            pts = np.vstack((pts, Included.get_points(config.grid_spacing)))
        pts = get_unique_rows(pts)

        # Collect exclusion points and subtract them
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

            pts = remove_exclusion_points(pts, pts_exclusion)

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
        logger.info("Starting pocket volume calculator")
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

        # User specified regions to include in the PDB structure.
        # Thus, we compute these points.
        pts = None
        if config.load_points_path is not None:
            logger.debug("Loading points")
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
        logger.info("Initializing volume calculator manager")
        ray_manager = RayManager(
            task_class=TaskComputeVolumeFromPDBLines,
            n_cores=config.n_cores,
            use_ray=config.use_ray,
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
            logger.debug(f"Submitted {len(tasks)} tasks in chunk")
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
