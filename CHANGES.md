Changes
=======

2.1
---

1. Updated version to 2.1.
2. Now Python 3 compatible.
3. Fixed minor bug that assigned atoms with names like "HG21" to element "HG"
   rather than "H". No longer supporting element "HG".
4. Improved formatting (black) and docstrings.
5. Moved GUI wrapper to new `depreciated/` folder. This feature is no longer
   supported.

2.0.3
-----

1. Updated citation.
2. Convex hull now calculated based on locations of all heavy atoms. In 2.0.1,
   it was all alpha carbons. Take care when comparing POVME 2.0.1 and POVME
   2.0.2/2.0.3 output.
3. Fixed NumPy warning.
4. Spelling error ("Angstroms").
