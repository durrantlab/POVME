Changes
=======

WIP
---

Minor stylistic changes to code.
NOTE: There is a MP error in python3. Need to fix. Use 1 processor for now.

2.2.2
-----

1. Fixed text-encoding bug.

2.2.1
-----

1. Made minor changes to the `README.md` file to clarify the
   `ContiguousPointsCriteria`, `ContiguousPocketSeedSphere`, and
   `ContiguousPocketSeedBox` parameters.
2. Otherwise identical to POVME 2.2.

2.2
---

1. Updated version to 2.2.
2. Fixed an error that prevented POVME from saving its output to uncompressed
   files. (Previously only writing output to compressed output succeeded.)
3. POVME now throws a warning when using Python 3.5 or earlier.
4. When the user does not specify any contiguous-pocket seed regions, POVME no
   longer outputs the default value of the ContiguousPointsCriteria parameter
   to the screen. In this case, POVME does not use the
   ContiguousPointsCriteria parameter anyway, so outputing its value to the
   screen only causes confusion.
5. Applied [Black formatter](https://black.readthedocs.io/en/stable/) to all
   Python code.

2.1
---

1. Updated version to 2.1.
2. Now Python3 compatible. Python2 no longer officially supported.
3. Refactored code so that it can now be called as a library.
4. Moved GUI wrapper to new `depreciated/` folder. This feature is no longer
   supported.
5. Improved formatting (black) and docstrings.
6. Fixed minor bug that assigned atoms with names like "HG21" to element "HG"
   rather than "H". No longer supporting element "HG".
7. Moved examples directory to `povme/examples/`.
8. Can now read both Windows and POSIX path names seamlessly (no need for
   separate ini files).
9. Easy testing now available: `python POVME2.py --test` and `python
   POVME_pocket_id.py --test`
10. Added project roadmap.

2.0.3
-----

1. Updated citation.
2. Convex hull now calculated based on locations of all heavy atoms. In 2.0.1,
   it was all alpha carbons. Take care when comparing POVME 2.0.1 and POVME
   2.0.2/2.0.3 output.
3. Fixed NumPy warning.
4. Spelling error ("Angstroms").

2.0.1
-----

Note that POVME 2.0.1 and all earlier versions can be found at
[https://sourceforge.net/projects/povme/files/](https://sourceforge.net/projects/povme/files/).
