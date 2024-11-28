from typing import Any

from .io import gzopenfile, write_to_file


def log(astr: str, parameters: dict[str, Any]) -> None:
    """Output POVME statements, either to the screen or to a file.

    Args:
        astr: The string to output.
        parameters: The user-defined parameters.
    """

    # Print the output to the screen.
    print(astr)

    # Save it to the output file as well.
    try:
        if parameters["CompressOutput"]:
            f = gzopenfile(parameters["OutputFilenamePrefix"] + "output.txt.gz", "ab")
        else:
            f = open(parameters["OutputFilenamePrefix"] + "output.txt", "a")

        write_to_file(f, astr + "\n", encode=parameters["CompressOutput"])

        f.close()
    except Exception:
        pass
