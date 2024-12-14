import glob
import os
import shutil

from povme.config import PocketDetectConfig
from povme.pocket.detect import PocketDetector


def test_pocket_detect_rofgp2(path_rogfp2_pdb, path_rogfp2_config, path_rogfp2_output):
    dir_output = os.path.dirname(path_rogfp2_output)
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    config = PocketDetectConfig()
    config.from_yaml(path_rogfp2_config)
    pocket_id = PocketDetector(config)
    pocket_id.run(path_rogfp2_pdb, path_rogfp2_output)

    with open(path_rogfp2_output + "pocket1.pdb", "r") as f:
        if "PointsInclusionSphere" not in f.read():
            raise Exception(
                "pocket1.pdb did not contain substring 'PointsInclusionSphere'"
            )

    num_output_files = len(glob.glob(path_rogfp2_output + "*.pdb")) - 1
    if num_output_files != 7:
        raise Exception("Expected 7 output files, but got " + str(num_output_files))


def test_pocket_detect_rel1(path_rel1_config, path_rel1_output):
    dir_output = os.path.dirname(path_rel1_output)
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    config = PocketDetectConfig()
    config.from_yaml(path_rel1_config)
    pocket_id = PocketDetector(config)
    pocket_id.run("./tests/files/rel1/rel1.pdb", path_rel1_output)

    with open(path_rel1_output + "pocket1.pdb", "r") as f:
        if "PointsInclusionSphere" not in f.read():
            raise Exception(
                "pocket1.pdb did not contain substring 'PointsInclusionSphere'"
            )

    num_output_files = len(glob.glob(path_rel1_output + "*.pdb")) - 1
    if num_output_files != 8:
        raise Exception("Expected 8 output files, but got " + str(num_output_files))
