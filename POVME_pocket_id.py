# POVME Pocket ID 1.0 is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html).

# If you have any questions, comments, or suggestions, please don't hesitate
# to contact me, Jacob Durrant, at durrantj [at] pitt [dot] edu.

# POVME Pocket ID 1.0 is a program for identifying protein pockets and
# generating appropriate pocket-encompassing inclusion spheres. These spheres,
# modified as required, can then be used as POVME input.

from povme import pocket_id
import sys

pocket_id.run_pocket_id(sys.argv)
