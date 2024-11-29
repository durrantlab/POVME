from typing import Any

import gzip
import platform
from pathlib import Path, PureWindowsPath

import numpy as np
import numpy.typing as npt


def fix_filename(path: str, must_exist: bool = True) -> Path:
    """
    Returns a standardized Path object, handling both Windows and Unix-formatted paths.

    Args:
        path: The file path to standardize.
        must_exist: If True, verifies that the path exists.

    Raises:
        FileNotFoundError: If the path does not exist and must_exist is True.
        ValueError: If the path format is invalid.

    Returns:
        Path: A pathlib.Path object representing the standardized path.
    """
    # Determine if the current OS is Windows
    is_windows = platform.system() == "Windows"

    # Attempt to create a Path object using the appropriate class
    try:
        if is_windows:
            path_obj: Path = Path(path)  # On Windows, Path() returns a WindowsPath
        else:
            # On non-Windows systems, check if the path is Windows-formatted
            pure_win_path = PureWindowsPath(path)
            if pure_win_path.is_absolute() and pure_win_path.drive:
                # It's an absolute Windows path
                path_obj = Path(pure_win_path)
            else:
                # It's a Unix path
                path_obj = Path(path)
    except Exception as e:
        raise ValueError(f"Invalid path format: {path}") from e

    if must_exist:
        if not path_obj.exists():
            raise FileNotFoundError(f"File doesn't exist: {path}")

    return path_obj


def openfile(path, mode="r"):
    """Opens a file, regardless of whether unix or windows formatted path.

    Args:
        path: A string, the path.
        mode: A string, the file mode. Defaults to "r".

    Returns:
        The file object.
    """

    return open(fix_filename(path, "r" in mode), mode)


def gzopenfile(path, mode="r"):
    """Opens a gz compressed file, regardless of whether unix or windows
    formatted path.

    Args:
        path: A string, the path.
        mode: A string, the file mode. Defaults to "r".

    Returns:
        The file object.
    """
    return gzip.open(fix_filename(path, "r" in mode), mode)


def write_to_file(fo, s, encode=True):
    """Writes to a file, encoding if necessary.

    Args:
        fo: The file object.
        s: The string to write.
        encode (bool, optional): Whether to encode the file before writing.
            Defaults to True.

    """

    if encode:
        s = s.encode()

    fo.write(s)


def create_pdb_line(numpy_array, index, resname, letter):
    """Create a string formatted according to the PDB standard.

    Args:
        numpy_array: A 1x3 np.array representing a 3D point.
        letter: A string, the atom name/chain/etc to use for the output.

    Returns:
        A string, formatted according to the PDB standard.

    """

    if len(numpy_array) == 2:
        numpy_array = np.array([numpy_array[0], numpy_array[1], 0.0])
    if numpy_array.shape == (1, 3):
        numpy_array = numpy_array[0]

    output = "ATOM "
    output = (
        output
        + str(index % 999999).rjust(6)
        + letter.rjust(5)
        + resname.rjust(4)
        + letter.rjust(2)
        + str(index % 9999).rjust(4)
    )
    output += ("%.3f" % numpy_array[0]).rjust(12)
    output += ("%.3f" % numpy_array[1]).rjust(8)
    output += ("%.3f" % numpy_array[2]).rjust(8)
    output += letter.rjust(24)

    return output


class write_pdbs:
    """A class for converting numpy arrays into PDB-formatted strings."""

    def __create_pdb_line(self, numpy_array, index, resname, letter):
        """Create a string formatted according to the PDB standard.

        Args:
            numpy_array: A 1x3 np.array representing a 3D point.
            index: An integer, the atom index to use in the string.
            resname: A string, the RESNAME to use.
            letter: A string, the atom name/chain/etc to use for the output.

        Returns:
            A string, formatted according to the PDB standard.

        """

        if len(numpy_array) == 2:
            numpy_array = np.array([numpy_array[0], numpy_array[1], 0.0])
        if numpy_array.shape == (1, 3):
            numpy_array = numpy_array[0]

        output = "ATOM "
        output = (
            output
            + str(index % 999999).rjust(6)
            + letter.rjust(5)
            + resname.rjust(4)
            + letter.rjust(2)
            + str(index % 9999).rjust(4)
        )
        output = output + ("%.3f" % numpy_array[0]).rjust(12)
        output = output + ("%.3f" % numpy_array[1]).rjust(8)
        output = output + ("%.3f" % numpy_array[2]).rjust(8)
        output = output + letter.rjust(24)

        return output

    def numpy_to_pdb(self, narray, letter, resname=""):
        """Create a string formatted according to the PDB standard.

        Args:
        narray: A nx3 np.array representing a 3D point.
        letter: A string, the atom name/chain/etc to use for the output.
        resname: An optional string, the RESNAME to use for the output.

        Returns:
        A string, formatted according to the PDB standard.

        """

        if len(narray.flatten()) == 3:
            return self.__create_pdb_line(narray, 1, "AAA", letter) + "\n"
        else:
            if resname == "":
                letters = [
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
                resnames = []
                for l1 in letters:
                    for l2 in letters:
                        for l3 in letters:
                            resnames.append(l1 + l2 + l3)
                resnames.remove("XXX")  # because this is reserved for empty atoms
            else:
                resnames = [resname]

            t = ""
            for i, item in enumerate(narray):
                t = (
                    t
                    + self.__create_pdb_line(
                        item, i + 1, resnames[i % len(resnames)], letter
                    )
                    + "\n"
                )
            return t


def numpy_to_pdb(narray, letter, resname=""):
    """Create a string formatted according to the PDB standard.

    Args:
        narray: A nx3 np.array representing a 3D point.
        letter: A string, the atom name/chain/etc to use for the output.

    Returns:
        A string, formatted according to the PDB standard.

    """

    if len(narray.flatten()) == 3:
        return create_pdb_line(narray, 1, "AAA", letter) + "\n"
    else:
        if resname == "":
            letters = [
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
            resnames = []
            for l1 in letters:
                for l2 in letters:
                    for l3 in letters:
                        resnames.append(l1 + l2 + l3)
            resnames.remove("XXX")  # because this is reserved for empty atoms
        else:
            resnames = [resname]

        t = ""
        for i, item in enumerate(narray):
            t = (
                t
                + create_pdb_line(item, i + 1, resnames[i % len(resnames)], letter)
                + "\n"
            )
        return t


def dx_freq(freq_mat: npt.NDArray[np.float64], parameters: dict[str, Any]) -> None:
    """Generates a DX file that records the frequency that a volume element is
    open.

    Args:
        freq_mat: a Nx4 matrix, where the first 3 columns are the x,y,z coords
            of the point, and the 4th column is the frequency of emptiness for that
            point in space

    """

    # 1. Sort the points into the proper order for a dx file

    # already sorted correctly

    # 2. Obtain key information about the grid
    N = freq_mat.shape[0]  # number of data points

    minx = min(freq_mat[:, 0])
    miny = min(freq_mat[:, 1])
    minz = min(freq_mat[:, 2])  # find the upper and lower corners of the grid
    maxx = max(freq_mat[:, 0])
    maxy = max(freq_mat[:, 1])
    maxz = max(freq_mat[:, 2])

    widthx = maxx - minx  # find the widths of the grid
    widthy = maxy - miny
    widthz = maxz - minz

    xs = np.unique(freq_mat[:, 0])
    ys = np.unique(freq_mat[:, 1])
    zs = np.unique(freq_mat[:, 2])

    resx = xs[1] - xs[0]
    resy = ys[1] - ys[0]
    resz = zs[1] - zs[0]

    # resx = freq_mat[(widthz+1)*(widthy+1),0] - freq_mat[0,0]
    # resy = freq_mat[widthz+1,1] - freq_mat[0,1] # find the resolution of the grid
    # resz = freq_mat[1,2] - freq_mat[0,2]

    nx = (widthx) / resx + 1  # number of grid points in each dimension
    # need to add one because the subtraction leaves out an entire row
    ny = (widthy) / resy + 1
    nz = (widthz) / resz + 1

    # test to make sure all is well with the size of the grid and its
    # dimensions
    assert (
        nx * ny * nz
    ) == N, "Something is wrong with the freq_mat array: it is not a prismatic shape"

    # 3. write the header and footer
    if config.save_volumetric_density_map:
        if config.compress_output:
            dx_file = gzopenfile(output_prefix + "volumetric_density.dx.gz", "wb")
        else:
            dx_file = openfile(output_prefix + "volumetric_density.dx", "w")

        header_template = """# Data from POVME 2.2.2
#
# FREQUENCY (unitless)
#
object 1 class gridpositions counts %d %d %d
origin %8.6e %8.6e %8.6e
delta %8.6e 0.000000e+00 0.000000e+00
delta 0.000000e+00 %8.6e 0.000000e+00
delta 0.000000e+00 0.000000e+00 %8.6e
object 2 class gridconnections counts %d %d %d
object 3 class array type double rank 0 items %d data follows
"""

        header = header_template % (
            nx,
            ny,
            nz,
            minx,
            miny,
            minz,
            resx,
            resy,
            resz,
            nx,
            ny,
            nz,
            N,
        )  # format the header
        footer_template = """
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3"""

        footer = footer_template  # the footer needs no formatting
        write_to_file(dx_file, header, encode=config.compress_output)
        newline_counter = 1
        for i in range(N):  # write the data to the DX file
            write_to_file(
                dx_file, "%8.6e" % freq_mat[i, 3], encode=config.compress_output
            )
            if newline_counter == 3:
                newline_counter = 0
                write_to_file(dx_file, "\n", encode=config.compress_output)
            else:
                write_to_file(dx_file, " ", encode=config.compress_output)
            newline_counter += 1
        write_to_file(dx_file, footer, encode=config.compress_output)
        dx_file.close
