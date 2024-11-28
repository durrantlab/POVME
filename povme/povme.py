from typing import Any

import os
import shutil
import time
from io import StringIO

import numpy as np
from scipy.spatial.distance import cdist

from . import pymolecule
from .config import ConfigFile
from .hull import MultithreadingCalcVolumeTask
from .io import dx_freq, fix_filename, gzopenfile, numpy_to_pdb, openfile, write_to_file
from .logger import log
from .parallel import MultiThreading, MultithreadingTaskGeneral
from .region import Region


def unique_rows(a):
    """Identifies unique points (rows) in an array of points.

    Arguments:
    a -- A nx3 numpy.array representing 3D points.

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
        parameters = item[2]

        # make the pdb object
        str_obj = StringIO(pdb_string)
        tmp = pymolecule.Molecule()
        tmp.fileio.load_pdb_into_using_file_object(str_obj, False, False, False)

        log("\tFurther processing frame " + str(index), parameters)

        # so load the whole trajectory into memory
        if not parameters["UseDiskNotMemory"]:
            self.results.append((index, tmp))
        else:  # save to disk, record filename
            pym_filename = f"./.povme_tmp/frame_{str(index)}.pym"
            tmp.fileio.save_pym(pym_filename, False, False, False, False, False)
            self.results.append((index, pym_filename))


class RunPOVME:
    """The main class to run POVME."""

    def reference(self, parameters, before=""):
        """Print out a message regarding terms of use."""

        log("", parameters)
        log(
            before
            + "If you use POVME in your research, please cite the following reference:",
            parameters,
        )
        log(
            before
            + '  Durrant, J. D., C. A. de Oliveira, et al. (2011). "POVME: An algorithm',
            parameters,
        )
        log(
            before
            + '  for measuring binding-pocket volumes." J Mol Graph Model 29(5): 773-776.',
            parameters,
        )

    def load_multi_frame_pdb(self, filename, parameters):
        """Load a multi-frame PDB into memory or into separate files
        (depending on user specifications).

        Args:
            filename: A string, the filename of the multi-frame PDB.
            parameters: A python dictionary, where the keys are the user-defined
                parameter names and the values are the corresponding parameter
                values.

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

        log("", parameters)
        log("Reading frames from " + filename, parameters)

        f = open(filename, "rb")
        while True:

            if parameters["NumFrames"] != -1:
                if len(pdb_strings) >= parameters["NumFrames"]:
                    break

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
            [
                (pdb_strings[idx], idx + 1, parameters)
                for idx in range(len(pdb_strings))
            ],
            parameters["NumProcessors"],
            MultithreadingStringToMoleculeTask,
        )
        molecules = molecules.results

        return molecules

    def __init__(
        self,
        path_config: str,
        path_pdb: str | None = None,
        output_prefix: str | None = None,
    ) -> None:
        """Start POVME

        Args:
            path_config: Path to configuration file.
            path_pdb: Path to PDB file. This will overwrite the configuration file.
            output_prefix: Path to output directory including directories.
        """
        start_time = time.time()

        config = ConfigFile(path_config)

        # Process the config file
        parameters: dict[str, Any] = {}

        parameters["GridSpacing"] = 1.0  # default
        parameters["PointsIncludeRegions"] = []
        parameters["PointsExcludeRegions"] = []
        parameters["SavePoints"] = False  # default
        parameters["LoadPointsFilename"] = ""  # default
        parameters["PDBFileName"] = ""  # default
        # default is VDW radius of hydrogen
        parameters["DistanceCutoff"] = 1.09
        parameters["ConvexHullExclusion"] = True
        parameters["ContiguousPocketSeedRegions"] = []
        parameters["ContiguousPointsCriteria"] = 4
        parameters["NumProcessors"] = 4
        parameters["UseDiskNotMemory"] = False
        parameters["OutputFilenamePrefix"] = (
            "POVME_output."
            + time.strftime("%m-%d-%y")
            + "."
            + time.strftime("%H-%M-%S")
            + "/"
        )
        parameters["SaveIndividualPocketVolumes"] = False
        parameters["SavePocketVolumesTrajectory"] = False
        parameters["OutputEqualNumPointsPerFrame"] = False
        parameters["SaveTabbedVolumeFile"] = False
        parameters["SaveVolumetricDensityMap"] = False
        parameters["CompressOutput"] = False
        # This is a parameter for debugging purposes only.
        parameters["NumFrames"] = -1

        float_parameters = ["GridSpacing", "DistanceCutoff"]
        boolean_parameters = [
            "SavePoints",
            "ConvexHullExclusion",
            "CompressOutput",
            "UseDiskNotMemory",
            "SaveVolumetricDensityMap",
            "OutputEqualNumPointsPerFrame",
            "SaveIndividualPocketVolumes",
            "SaveTabbedVolumeFile",
            "SavePocketVolumesTrajectory",
        ]
        int_parameters = ["NumFrames", "ContiguousPointsCriteria", "NumProcessors"]
        filename_parameters = [
            "OutputFilenamePrefix",
            "PDBFileName",
            "LoadPointsFilename",
        ]

        for entity in config.entities:
            if entity[0] == "PDBFILENAME":
                if path_pdb is not None:
                    entity[1] = path_pdb

            if entity[0] == "OUTPUTFILENAMEPREFIX":
                if output_prefix is not None:
                    if output_prefix[-1] != "/":
                        output_prefix += "/"
                    entity[1] = output_prefix

            try:
                index = [p.upper() for p in float_parameters].index(entity[0])
                parameters[float_parameters[index]] = float(entity[1])
            except:
                pass

            try:
                index = [p.upper() for p in boolean_parameters].index(entity[0])

                if entity[1].upper() in ["YES", "TRUE"]:
                    parameters[boolean_parameters[index]] = True
                else:
                    parameters[boolean_parameters[index]] = False
            except:
                pass

            try:
                index = [p.upper() for p in int_parameters].index(entity[0])
                parameters[int_parameters[index]] = int(entity[1])
            except:
                pass

            try:
                index = [p.upper() for p in filename_parameters].index(entity[0])
                parameters[filename_parameters[index]] = fix_filename(
                    entity[1].strip(), False
                ).as_posix()  # so a string
            except:
                pass

            # Regions are handled separately for each parameter...
            if entity[0] == "POINTSINCLUSIONSPHERE":
                Include = Region()
                items = entity[1].split(" ")
                Include.center[0] = float(items[0])
                Include.center[1] = float(items[1])
                Include.center[2] = float(items[2])
                Include.radius = float(items[3])
                Include.region_type = "SPHERE"
                parameters["PointsIncludeRegions"].append(Include)
            elif entity[0] == "POINTSINCLUSIONBOX":
                Include = Region()
                items = entity[1].split(" ")
                Include.center[0] = float(items[0])
                Include.center[1] = float(items[1])
                Include.center[2] = float(items[2])
                Include.box_dimen[0] = float(items[3])
                Include.box_dimen[1] = float(items[4])
                Include.box_dimen[2] = float(items[5])
                Include.region_type = "BOX"
                parameters["PointsIncludeRegions"].append(Include)
            if entity[0] == "CONTIGUOUSPOCKETSEEDSPHERE":
                Contig = Region()
                items = entity[1].split(" ")
                Contig.center[0] = float(items[0])
                Contig.center[1] = float(items[1])
                Contig.center[2] = float(items[2])
                Contig.radius = float(items[3])
                Contig.region_type = "SPHERE"
                parameters["ContiguousPocketSeedRegions"].append(Contig)
            elif entity[0] == "CONTIGUOUSPOCKETSEEDBOX":
                Contig = Region()
                items = entity[1].split(" ")
                Contig.center[0] = float(items[0])
                Contig.center[1] = float(items[1])
                Contig.center[2] = float(items[2])
                Contig.box_dimen[0] = float(items[3])
                Contig.box_dimen[1] = float(items[4])
                Contig.box_dimen[2] = float(items[5])
                Contig.region_type = "BOX"
                parameters["ContiguousPocketSeedRegions"].append(Contig)
            elif entity[0] == "POINTSEXCLUSIONSPHERE":
                Exclude = Region()
                items = entity[1].split(" ")
                Exclude.center[0] = float(items[0])
                Exclude.center[1] = float(items[1])
                Exclude.center[2] = float(items[2])
                Exclude.radius = float(items[3])
                Exclude.region_type = "SPHERE"
                parameters["PointsExcludeRegions"].append(Exclude)
            elif entity[0] == "POINTSEXCLUSIONBOX":
                Exclude = Region()
                items = entity[1].split(" ")
                Exclude.center[0] = float(items[0])
                Exclude.center[1] = float(items[1])
                Exclude.center[2] = float(items[2])
                Exclude.box_dimen[0] = float(items[3])
                Exclude.box_dimen[1] = float(items[4])
                Exclude.box_dimen[2] = float(items[5])
                Exclude.region_type = "BOX"
                parameters["PointsExcludeRegions"].append(Exclude)

        # If there are no ContiguousPocketSeedRegions, don't print out the
        # ContiguousPointsCriteria parameter.
        if len(parameters["ContiguousPocketSeedRegions"]) == 0:
            del parameters["ContiguousPointsCriteria"]

        # If the output prefix includes a directory, create that directory if
        # necessary
        if "/" in parameters["OutputFilenamePrefix"]:
            output_dirname = os.path.dirname(parameters["OutputFilenamePrefix"])
            # if os.path.exists(output_dirname): shutil.rmtree(output_dirname) # So delete the directory if it already exists.
            try:
                os.mkdir(output_dirname)
            except:
                pass

        # print out the header
        self.reference(parameters, "")
        log("", parameters)

        # create temp swap directory if needed
        if parameters["UseDiskNotMemory"]:
            if os.path.exists("./.povme_tmp"):
                shutil.rmtree("./.povme_tmp")
            os.mkdir("./.povme_tmp")

        # print out parameters
        log("Parameters:", parameters)
        for i in list(parameters.keys()):

            if i == "NumFrames" and parameters["NumFrames"] == -1:
                # So only show this parameter if it's value is not the
                # default.
                continue

            if type(parameters[i]) is list:
                for i2 in parameters[i]:
                    if i2 != "":
                        log("\t" + str(i) + ": " + str(i2), parameters)
            else:
                if parameters[i] != "":
                    log("\t" + str(i) + ": " + str(parameters[i]), parameters)

        pts = None
        if len(parameters["PointsIncludeRegions"]) > 0:  # so create the point file

            log("\nGenerating the pocket-encompassing point field", parameters)

            # get all the points of the inclusion regions
            pts = parameters["PointsIncludeRegions"][0].points_set(
                parameters["GridSpacing"]
            )
            for Included in parameters["PointsIncludeRegions"][1:]:
                pts = np.vstack((pts, Included.points_set(parameters["GridSpacing"])))
            pts = unique_rows(pts)

            # get all the points of the exclusion regions
            if len(parameters["PointsExcludeRegions"]) > 0:
                pts_exclusion = parameters["PointsExcludeRegions"][0].points_set(
                    parameters["GridSpacing"]
                )
                for Excluded in parameters["PointsExcludeRegions"][1:]:
                    pts_exclusion = np.vstack(
                        (pts_exclusion, Excluded.points_set(parameters["GridSpacing"]))
                    )
                pts_exclusion = unique_rows(pts_exclusion)

                # remove the exclusion points from the inclusion points I
                # think there ought to be a set-based way of doing this, but
                # I'm going to go for the pairwise comparison. consider
                # rewriting later
                index_to_remove = np.nonzero(cdist(pts, pts_exclusion) < 1e-7)[0]
                pts = np.delete(pts, index_to_remove, axis=0)

            # save the points as PDB
            if parameters["SavePoints"]:

                # First, save the point field itself

                log("\nSaving the point field as a PDB and NPY file", parameters)

                points_filename = parameters["OutputFilenamePrefix"] + "point_field.pdb"

                if parameters["CompressOutput"]:
                    afile = gzopenfile(points_filename + ".gz", "wb")
                else:
                    afile = openfile(points_filename, "w")

                write_to_file(
                    afile, numpy_to_pdb(pts, "X"), encode=parameters["CompressOutput"]
                )
                afile.close()

                # save the points as npy

                np.save(points_filename + ".npy", pts)

                log(
                    "\tPoint field saved to "
                    + points_filename
                    + " to permit visualization",
                    parameters,
                )
                log(
                    "\tPoint field saved to "
                    + points_filename
                    + ".npy to optionally load for the volume calculation",
                    parameters,
                )
                log("", parameters)

                # Now, save the contiguous seed points as well, if specified.
                if len(parameters["ContiguousPocketSeedRegions"]) > 0:
                    # get all the contiguous points
                    contig_pts = parameters["ContiguousPocketSeedRegions"][
                        0
                    ].points_set(parameters["GridSpacing"])
                    for Contig in parameters["ContiguousPocketSeedRegions"][1:]:
                        contig_pts = np.vstack(
                            (contig_pts, Contig.points_set(parameters["GridSpacing"]))
                        )
                    contig_pts = unique_rows(contig_pts)

                    log(
                        "\nSaving the contiguous-pocket seed points as a PDB file",
                        parameters,
                    )

                    points_filename = (
                        parameters["OutputFilenamePrefix"]
                        + "contiguous_pocket_seed_points.pdb"
                    )

                    if parameters["CompressOutput"]:
                        afile = gzopenfile(points_filename + ".gz", "wb")
                    else:
                        afile = openfile(points_filename, "w")

                    write_to_file(
                        afile,
                        numpy_to_pdb(contig_pts, "X"),
                        encode=parameters["CompressOutput"],
                    )
                    afile.close()

                    log(
                        "\tContiguous-pocket seed points saved to "
                        + points_filename
                        + " to permit visualization",
                        parameters,
                    )
                    log("", parameters)

        # so there's a PDB point specified for calculating the volume.
        if parameters["PDBFileName"] != "":

            # load the points in they aren't already present
            if pts is None:
                log("\nLoading the point-field NPY file...", parameters)
                parameters["pts_orig"] = np.load(parameters["LoadPointsFilename"])
            else:
                parameters["pts_orig"] = pts

            # load the PDB frames
            index_and_pdbs = self.load_multi_frame_pdb(
                parameters["PDBFileName"], parameters
            )

            # calculate all the volumes
            log("", parameters)
            log("Calculating the pocket volume of each frame", parameters)
            tmp = MultiThreading(
                [
                    (index, pdb_object, parameters)
                    for index, pdb_object in index_and_pdbs
                ],
                parameters["NumProcessors"],
                MultithreadingCalcVolumeTask,
            )

            # delete the temp swap directory if necessary
            if parameters["UseDiskNotMemory"]:
                if os.path.exists("./.povme_tmp"):
                    shutil.rmtree("./.povme_tmp")

            # display the results
            results_dic = {}
            for result in tmp.results:
                results_dic[result[0]] = result[1]
            log("", parameters)
            log("FRAME        | VOLUME (A^3)", parameters)
            log("-------------+-------------", parameters)
            for i in sorted(results_dic.keys()):
                log(str(i).ljust(13) + "| " + str(results_dic[i]), parameters)

            log("", parameters)
            log(
                "Execution time = " + str(time.time() - start_time) + " sec", parameters
            )
            log("", parameters)

            # if the user requested a separate volume file, save that as well
            if parameters["SaveTabbedVolumeFile"]:
                if parameters["CompressOutput"]:
                    f = gzopenfile(
                        parameters["OutputFilenamePrefix"] + "volumes.tabbed.txt.gz",
                        "wb",
                    )

                else:
                    f = openfile(
                        parameters["OutputFilenamePrefix"] + "volumes.tabbed.txt", "w"
                    )
                for i in sorted(results_dic.keys()):
                    write_to_file(
                        f,
                        str(i) + "\t" + str(results_dic[i]) + "\n",
                        encode=parameters["CompressOutput"],
                    )

                f.close()

            # if the user wanted a single trajectory containing all the
            # volumes, generate that here.
            if parameters["SavePocketVolumesTrajectory"]:
                if parameters["CompressOutput"]:
                    traj_file = gzopenfile(
                        parameters["OutputFilenamePrefix"] + "volume_trajectory.pdb.gz",
                        "wb",
                    )
                else:
                    traj_file = openfile(
                        parameters["OutputFilenamePrefix"] + "volume_trajectory.pdb",
                        "w",
                    )

                for frame_index in range(1, len(list(results_dic.keys())) + 1):
                    if parameters["CompressOutput"]:
                        frame_file = gzopenfile(
                            parameters["OutputFilenamePrefix"]
                            + "frame_"
                            + str(frame_index)
                            + ".pdb.gz",
                            "rb",
                        )
                    else:
                        frame_file = openfile(
                            parameters["OutputFilenamePrefix"]
                            + "frame_"
                            + str(frame_index)
                            + ".pdb",
                            "r",
                        )

                    traj_file.write(frame_file.read())
                    frame_file.close()

                traj_file.close()

            # if the user requested a volumetric density map, then generate it
            # here
            if parameters["SaveVolumetricDensityMap"]:
                unique_points = {}

                overall_min = np.ones(3) * 1e100
                overall_max = np.ones(3) * -1e100

                for result in tmp.results:
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
                            except:
                                unique_points[pt_key] = 1
                if overall_min[0] == 1e100:
                    log(
                        "ERROR! Cannont save volumetric density file because no volumes present in any frame.",
                        parameters,
                    )
                else:
                    xpts = np.arange(
                        overall_min[0],
                        overall_max[0] + parameters["GridSpacing"],
                        parameters["GridSpacing"],
                    )
                    ypts = np.arange(
                        overall_min[1],
                        overall_max[1] + parameters["GridSpacing"],
                        parameters["GridSpacing"],
                    )
                    zpts = np.arange(
                        overall_min[2],
                        overall_max[2] + parameters["GridSpacing"],
                        parameters["GridSpacing"],
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
                                except:
                                    pass

                                i = i + 1

                    # convert the counts in the fourth column into frequencies
                    all_pts[:, 3] = all_pts[:, 3] / len(tmp.results)
                    dx_freq(all_pts, parameters)  # save the dx file

                    # print "To turn into a DX file:"
                    # print all_pts
                    # import cPickle as pickle
                    # pickle.dump(all_pts, open('dill.pickle', 'w'))

        self.results = results_dic
