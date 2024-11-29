import glob
import os
import shutil

from povme import PocketDetector


def test_rel1(path_rel1_config, path_rel1_output):
    dir_output = os.path.dirname(path_rel1_output)
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    pocket_id = PocketDetector(path_rel1_config)
    pocket_id.run("./tests/files/rel1/rel1.pdb", path_rel1_output)

    with open(path_rel1_output + "pocket1.pdb", "r") as f:
        if "PointsInclusionSphere" not in f.read():
            raise Exception(
                "pocket1.pdb did not contain substring 'PointsInclusionSphere'"
            )

    num_output_files = len(glob.glob(path_rel1_output + "*.pdb")) - 1
    if num_output_files != 8:
        raise Exception("Expected 8 output files, but got " + str(num_output_files))
