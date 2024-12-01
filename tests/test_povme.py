import glob
import os
import shutil

from povme import POVME


def test_4nss(path_4nss_config, path_4nss_output):
    dir_output = os.path.dirname(path_4nss_output)
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    povme = POVME(path_4nss_config)
    results = povme.run("./tests/files/4nss/4nss.pdb", path_4nss_output)

    expected_vols = set([1673.0, 1493.0, 1711.0, 1854.0, 2023.0])
    actual_vols = set(results.values())
    assert actual_vols == expected_vols
    num_output_files = len(glob.glob(dir_output + "/*"))
    # Was 12, but we currently do not write a log file.
    if num_output_files != 11:
        raise Exception("Expected 11 output files, but got " + str(num_output_files))
