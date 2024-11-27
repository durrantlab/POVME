import glob
import os
import shutil

from povme import RunPOVME


def test_4nss(path_4nss_input, path_4nss_output):
    dir_output = os.path.dirname(path_4nss_output)
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    run_povme = RunPOVME(
        path_4nss_input,
        path_pdb="./tests/files/4nss/4nss.pdb",
        output_prefix=path_4nss_output,
    )
    results = run_povme.results

    expected_vols = set([1673.0, 1493.0, 1711.0, 1854.0, 2023.0])
    actual_vols = set(results.values())
    if len(actual_vols - expected_vols) > 0:
        raise Exception(
            "Expected volumes to be "
            + str(expected_vols)
            + ", but got "
            + str(actual_vols)
        )
    num_output_files = len(glob.glob(dir_output + "/*"))
    if num_output_files != 12:
        raise Exception("Expected 12 output files, but got " + str(num_output_files))
