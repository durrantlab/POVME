import os

import numpy as np
from loguru import logger
from pymolecule import Molecule
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist

from .config import PocketIDConfig
from .hull import ConvexHull
from .io import openfile, write_pdbs
from .pocket import BoxOfPoints


class PocketDetector:
    """The main class to run pocket detection."""

    def __init__(
        self,
        path_config: str | None = None,
    ) -> None:
        """Initialize Pocket detector.

        Args:
            path_config: Path to a configuration YAML file.
        """

        self.config = PocketIDConfig()
        if path_config is not None:
            self.config.from_yaml(path_config)

    def run(self, path_pdb: str, output_prefix: str = "") -> None:

        config = self.config

        # If the output prefix includes a directory, create that directory if
        # necessary
        if "/" in output_prefix:
            output_dirname = os.path.dirname(output_prefix)
            os.makedirs(output_dirname, exist_ok=True)

        # Step 1: Load in the protein

        logger.info("Step 1. Loading the PDB file " + path_pdb + "...")
        molecule = Molecule()
        molecule.io.load_pdb_into(path_pdb)

        # Step 2: Get rid of hydogen atoms. They just slow stuff down.

        print("Step 2. Removing hydrogen atoms...")
        sel = molecule.selections.select_atoms({"element_stripped": b"H"})
        sel = molecule.selections.invert_selection(sel)
        molecule = molecule.selections.get_molecule_from_selection(sel)

        # Step 3: Calculate the convex hull of the protein alpha carbons.
        print("Step 3. Calculating the convex hull of the PDB file...")

        # Get a version of the protein with just the alpha carbons. In my
        # experience, that's better for convex hull identification. Otherwise the
        # program identifies shallow contors in the protein surface as pockets.
        molecule_alpha_carbons = molecule.selections.get_molecule_from_selection(
            molecule.selections.select_atoms({"name_stripped": b"CA"})
        )
        convex_hull_3d = ConvexHull(
            molecule_alpha_carbons.information.get_coordinates()
        )

        # Step 4. Get a box of equispaced points that surround the protein,
        # snapped to reso. I'm putting a whole bunch of other functions in this
        # class as well to manipulate the points of this box.

        logger.info(
            "Step 4. Making a box of points spaced "
            + str(config.pocket_detection_resolution)
            + " A apart that entirely encompasses the protein..."
        )

        # note that the initial box is low resolution (* 4) so convex hull will be
        # very fast
        box_pts = BoxOfPoints(
            molecule.information.get_bounding_box(),
            config.pocket_detection_resolution * 4,
        )

        # Step 5. Remove points outside the convex hull. Gradually fill in
        # protein-occupying region with denser point fields. Faster this way, I
        # think.
        logger.info(
            "Step 5. Removing points that fall outside the protein's convex hull..."
        )
        box_pts.remove_points_outside_convex_hull(convex_hull_3d, config)
        box_pts.expand_around_existing_points(2, config.pocket_detection_resolution * 2)
        box_pts.remove_points_outside_convex_hull(convex_hull_3d, config)
        box_pts.expand_around_existing_points(2, config.pocket_detection_resolution)
        box_pts.remove_points_outside_convex_hull(convex_hull_3d, config)

        # Step 6. Remove the points in this box that are too close to protein
        # atoms. For simplicity's sake, don't worry about atomic radii. Just a
        # simple cutoff.
        logger.info(
            "Step 6. Removing points that come within "
            + str(config.clashing_cutoff)
            + " A of any protein atom..."
        )
        box_pts.remove_all_points_close_to_other_points(
            molecule.information.get_coordinates(), config.clashing_cutoff, config
        )

        # Step 7. Now surround each of these points with higher density points
        # that in the same regions. This is for getting a more detailed view of
        # the identified pockets.
        if config.pocket_measuring_resolution != config.pocket_detection_resolution:
            logger.info(
                "Step 7. Flooding the identified pockets with points spaced "
                + str(config.pocket_measuring_resolution)
                + " A apart for a more detailed measurement of the pocket volume..."
            )
            print("\tAdding points...")
            box_pts.expand_around_existing_points(
                config.pocket_detection_resolution / config.pocket_measuring_resolution,
                config.pocket_measuring_resolution,
            )
            logger.info("\tRemoving points that fall outside the convex hull...")
            box_pts.remove_points_outside_convex_hull(convex_hull_3d, config)
            logger.info(
                "\tRemoving points within "
                + str(config.clashing_cutoff)
                + " A of any protein atom..."
            )
            box_pts.remove_all_points_close_to_other_points(
                molecule.information.get_coordinates(), config.clashing_cutoff, config
            )

        # Step 8. Now start doing a repeated pass filter (keep repeating until no
        # change). Don't know if this is a high pass or low pass filter. I've
        # heard these terms, though, and they sound cool.
        logger.info(
            "Step 8. Removing points until all points have at least "
            + str(config.n_neighbors)
            + " neighbors..."
        )
        box_pts.filter_isolated_points_until_no_change(
            config.pocket_measuring_resolution, config.n_neighbors
        )

        # Step 9. Separate out the pockets so they can be considered in isolation.
        logger.info("Step 9. Partitioning the remaining points by pocket...")
        all_pockets = box_pts.separate_out_pockets()

        # Step 10. Get povme spheres that encompass each pocket, write pockets to
        # separate pdb files
        logger.info("Step 10. Saving the points of each pocket...")
        let_ids = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "U",
            "V",
            "W",
            "X",
            "Y",
            "Z",
        ]
        write_some_pdbs = write_pdbs()

        for i, pts in enumerate(all_pockets):
            path_pdb = output_prefix + "pocket" + str(i + 1) + ".pdb"
            logger.info("Saving " + path_pdb + "...")
            f = openfile(path_pdb, "w")
            f.write("REMARK Pocket #" + str(i + 1) + "\n")

            # do I need to whiten stuff here? not sure what whitening is.

            centroids, idx = kmeans2(pts, config.n_spheres)

            pts_string = ""
            for cluster_num in range(config.n_spheres):
                indexes_for_this_cluster = np.nonzero(idx == cluster_num)[0]
                cluster_pts = pts[indexes_for_this_cluster]
                cluster_center = np.mean(cluster_pts, axis=0)
                try:
                    cluster_radius = np.max(
                        cdist(np.array([cluster_center]), cluster_pts)
                    )
                    f.write(
                        "REMARK CHAIN "
                        + let_ids[cluster_num]
                        + ": PointsInclusionSphere "
                        + str(np.round(cluster_center[0], 2))
                        + " "
                        + str(np.round(cluster_center[1], 2))
                        + " "
                        + str(np.round(cluster_center[2], 2))
                        + " "
                        + str(np.round(cluster_radius + config.sphere_padding, 2))
                        + "\n"
                    )
                    pts_string = pts_string + write_some_pdbs.numpy_to_pdb(
                        cluster_pts, let_ids[cluster_num]
                    )
                except:
                    logger.info(
                        "There was an error, but I don't think it was catastrophic. Could be that one of the pocket clusters was empty."
                    )

            f.write(pts_string)
            f.close()

        logger.info(
            "Done. See the pocket{n}.pdb files. Using a visualization program like VMD, identify which of these files includes the pocket you wish to measure. POVME Pocket ID has divided each pocket volume into "
            + str(config.n_spheres)
            + " sections (i.e., PDB chains). In some cases, the pocket you're interested in might be included in a larger identified pocket, so feel free to use only certain sections of a given pocket as well."
        )
        logger.info(
            "The POVME PointsInclusionSphere commands are located in the header of each pocket{n}.pdb file. A text editor can be used to copy and paste these commands into a POVME input file."
        )
