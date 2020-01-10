# POVME 2.1 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html). If you have any questions, comments,
# or suggestions, please don't hesitate to contact me, Jacob Durrant, at
# durrantj [at] pitt [dot] edu.
#
# If you use POVME in your work, please cite:
#
#    Durrant, J. D., C. A. de Oliveira, et al. (2011). "POVME: An algorithm
#    for measuring binding-pocket volumes." J Mol Graph Model 29(5): 773-776.
#
#    Durrant, J. D., L. Votapka, et al. (2014). "Povme 2.0: An Enhanced Tool
#    for Determining Pocket Shape and Volume Characteristics." J Chem Theory
#    Comput 10.11: 5047-56.

from pathlib import Path
import gzip
import os

try:
    from pathlib import WindowsPath as WinPath

    t = WinPath()
except:
    # This is all that's available on linux. But doesn't work with absolute
    # windows paths. On the other hand, paths starting with "C:\" don't make
    # sense on linux anyway.
    from pathlib import PureWindowsPath as WinPath


def fix_filename(path, must_exist=True):
    """Gets the proper filename, regardless of whether windows or unix
    formatted.

    path -- A string, the path.
    must_exist -- A boolean, whether to use file existance in determining which
        path to use.

    Returns:
    A standardized version of the filename.

    """
    if must_exist:
        # Figure out path format based on whether the file exists.
        path_fixed = Path(path)
        if not path_fixed.exists():
            # The path doesn't exist. Could it be a windows-formatted path?
            path_fixed = Path(WinPath(path))
            if not path_fixed.exists():
                # No, it really doesn't exist.

                raise IOError("File doesn't exist: " + path)
    else:
        # Try to figure it out without testing if file exists.
        if not "/" in path and not "\\" in path:
            # No directory, so it doesn't really matter.
            path_fixed = Path(path)
        elif os.path.dirname(path) == "":
            # It must be a windows path.
            path_fixed = Path(WinPath(path))
        else:
            # It's probably a unix path.
            path_fixed = Path(path)

    return path_fixed


def openfile(path, mode="r"):
    """Opens a file, regardless of whether unix or windows formatted path.

    path -- A string, the path.
    mode -- A string, the file mode. Defaults to "r".

    Returns:
    The file object.
    """

    return open(fix_filename(path, "r" in mode), mode)


def gzopenfile(path, mode="r"):
    """Opens a gz compressed file, regardless of whether unix or windows
    formatted path.

    path -- A string, the path.
    mode -- A string, the file mode. Defaults to "r".

    Returns:
    The file object.
    """
    return gzip.open(fix_filename(path, "r" in mode), mode)
