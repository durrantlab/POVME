POVME 2.0.2
===========

0\. License: GNU General Public License version 3
-------------------------------------------------

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see
[https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

1\. Download POVME 2.0
----------------------

Begin by downloading POVME 2.0. An example input file is included with the
download (in the 'examples' directory). This input file is heavily commented
and may be even more useful that the information provided on this website.

POVME is released under the GNU General Public License. If you have any
questions, comments, or suggestions, please don't hesitate to contact me,
[Jacob Durrant](http://durrantlab.com), at durrantj [at] pitt [dot] edu.

If you use POVME in your work, please cite:

1. Durrant, J. D., C. A. de Oliveira, et al. (2011). "POVME: An algorithm for
   measuring binding-pocket volumes." J Mol Graph Model 29(5): 773-776.
2. Durrant, J.D., L. Votapka, J. Sørensen, and R. E. Amaro (2014). "POVME 2.0:
   An Enhanced Tool for Determining Pocket Shape and Volume Characteristics."
   J. Chem. Theory Comput. 10(11):5047-5056.
 
2\. Align a PDB-formatted trajectory
------------------------------------

POVME accepts a multi-frame PDB (Protein Data Bank) file as input. The
computer program [Visual Molecular Dynamics
(VMD)](http://www.ks.uiuc.edu/Research/vmd/) is useful for aligning
trajectories and converting files to the PDB format. Alignment is necessary
because the POVME algorithm assumes the pocket being measured does not
translate or rotate in space. We note also that single-frame PDB files can
likewise serve as POVME input if the user wishes only to measure the volume of
a single pocket.

3\. Define an inclusion region
------------------------------

The user must define an "inclusion" region. This region is constructed from a
combination of user-specified spheres and rectangular prisms. The required
inclusion region should entirely encompass all the binding-pocket
conformations of the trajectory. Specify the spheres and rectangular prisms of
the inclusion region in a text-based POVME input file (e.g., "POVME.in"):

```
PointsInclusionSphere -7.12 2.60 -4.67 6.0
PointsInclusionSphere -2.0 -2.0 -4.0 5.0
PointsInclusionBox -5.0 -7.0 2.0 10.0 10.0 10.0
PointsInclusionBox -2.0 -9.0 -3.0 10.0 5.0 10.0
```

Note that for the spheres, the first three parameters are coordinates, and the
fourth parameter is a radius. For the rectangular prisms (i.e. "boxes"), the
first three parameters are coordinates, and the last three are the dimensions
of the box in the X, Y, and Z directions.

4\. Define an exclusion region
------------------------------

An optional exclusion region defines portions of the inclusion region that
should be ignored, perhaps because they are not truly associated with the
pocket. It is similarly constructed from spheres and boxes:

```
PointsExclusionSphere -2.0 -2.0 -4.0 5.0
PointsExclusionBox -5.0 -7.0 2.0 10.0 10.0 10.0
```
 
5\. Create a field of equidistant points
----------------------------------------

To generate a field of equidistant points that encompasses all the
binding-pocket conformations of the trajectory, POVME first floods the
user-specified inclusion region with points and then removes any points also
contained in the optional exclusion region. You need to specify the distance
separating each of these equidistant points:

`GridSpacing 1.0`
 
6\. How to choose the inclusion and exclusion regions
-----------------------------------------------------

As you can imagine, identifying just the right set of inclusion and exclusion
spheres and boxes to encompass the binding pocket is challenging. One approach
is to define an initial geometry, visualize that geometry together with the
receptor using a program like VMD, and then iteratively add new inclusion and
exclusion regions as required. You can optionally save the point field to a
file called point_field.pdb for visualization:

`SavePoints true`
 
7\. Specify the location of the receptor PDB file to analyze
------------------------------------------------------------

Once you've properly generated a pocket-encompassing point field, you're ready
to use that point field to calculate pocket volumes. Here's how to specify the
location of the PDB receptor file that has the pocket you wish to analyze:

`PDBFileName my_receptor.pdb`

Note that this file can be a trajectory containing multiple frames.
 
8\. Remove points that are near receptor atoms
----------------------------------------------

As the purpose of POVME is to measure the volume of a binding-pocket cavity,
the program next removes any points that are close to receptor atoms, leaving
only those points that are likely to be located within the binding pocket
itself. To specify how close the points can come to the van der Waal's surface
of the receptor before being removed:

`DistanceCutoff 1.09`

Note that if the receptor PDB file contains multiple frames, this will be done
on a frame-by-frame basis.

9\. Remove points outside the receptor's convex hull
----------------------------------------------------

POVME 2.0 introduces an optional new feature for removing points that lie
entirely outside the binding pocket. Specifically, the gift-wrapping algorithm
is used in combination with the Akl-Toussaint heuristic to define the convex
hull of receptor atoms near the user-defined inclusion region. Any points that
fall outside the convex hull are removed. This feature is particularly useful
when the user defines an inclusion region that protrudes into the surrounding
solvent-occupying space. To activate the convex-hull feature:

`ConvexHullExclusion true`
 
10\. Remove points that are not contiguous with the primary pocket
------------------------------------------------------------------

Like the original POVME program, version 2.0 retains the optional ability to
remove isolated patches of points that are not contiguous with the primary
binding pocket. This feature requires that the user define a third region,
again using spheres and rectangular prisms, that always falls within the
primary binding-pocket region, regardless of the trajectory frame considered:

```
ContiguousPocketSeedSphere 67.0 102.0 57.0 4.0
ContiguousPocketSeedBox 50.0 50.0 50.0 10.0 10.0 10.0
```

Two pocket volumes are considered "contiguous" if they share at least X
neighboring points in common, where X is defined by:

`ContiguousPointsCriteria 3`

Note that points that are "kitty-corner" from each other count as neighbors.
 
All pocket-occupying points within or contiguous to this region are retained,
but isolated patches of points that are not directly connected are deleted.

11\. Additional POVME parameters
--------------------------------

Here are some additional POVME parameters you might find helpful:

```
# Tell POVME how to perform the calculations.

NumProcessors               12                             # POVME can use multiple processors on
                                                           # Unix-based systems.

UseDiskNotMemory            false                          # In some cases, your PDB trajectory may
                                                           # be so large that the resulting POVME
                                                           # analysis cannot be easily stored in
                                                           # your computer's memory. If
                                                           # UseDiskNotMemory is set to true, POVME
                                                           # will rely more on your disk space than
                                                           # on memory/RAM.

# Tell POVME how to save the output

OutputFilenamePrefix          ./POVME_test_run/POVME_      # All the files POVME outputs will start
                                                           # with this prefix. POVME automatically
                                                           # creates any required directory
                                                           # (./POVME_test_run/ in this case).

SaveIndividualPocketVolumes   true                         # If true, POVME saves the pocket-volume
                                                           # points of each frame to a separate PDB
                                                           # file. The file names will be like
                                                           # {PREFIX}frame_X.pdb.

SavePocketVolumesTrajectory   true                         # If true, POVME saves all the pocket-
                                                           # volume points of each frame to a single
                                                           # PDB trajectory file. The individual
                                                           # frames are separated by END cards. The
                                                           # file name will be
                                                           # {PREFIX}volume_trajectory.pdb.

OutputEqualNumPointsPerFrame  true                         # Some visualization programs (e.g. VMD)
                                                           # are only compatible with trajectories
                                                           # that have the same number of atoms in
                                                           # each frame. If EqualNumAtomsPerFrame is
                                                           # true, POVME adds extra points at the
                                                           # origin (0.0, 0.0, 0.0) to satisfy this
                                                           # requirement. This affects files created
                                                           # with both SaveIndividualPocketVolumes
                                                           # and SavePocketVolumesTrajectory but
                                                           # does not alter the volume calculation
                                                           # itself.

SaveTabbedVolumeFile          true                         # If true, POVME saves the calculated
                                                           # volumes to a file in a simple tabular
                                                           # format that can be easily pasted into
                                                           # popular spreadsheet programs like 
                                                           # Microsoft Excel. The file is named 
                                                           # {PREFIX}volumes.tabbed.txt

SaveVolumetricDensityMap      true                         # If true, POVME saves a volumetric
                                                           # density map in the DX format. A
                                                           # volumetric density value is associated
                                                           # with each of the pocket-occupying
                                                           # points by calculating the fraction of
                                                           # all trajectory pocket volumes that
                                                           # include the given point. The file is 
                                                           # named {PREFIX}volumetric_density.dx.

CompressOutput                true                         # If you're short on disk space, POVME
                                                           # can automatically compress all output
                                                           # files using gz compression.
```

12\. POVME output
-----------------

By default, POVME writes a number of files to the disk. The calculated pocket
volumes, as well as user-defined parameters and progress messages, are saved
to a simple text-based log file. POVME can also be instructed to save the
volume measurements to a second file in a simple tabular format that can be
easily pasted into popular spreadsheet programs. Pocket-occupying points are
equidistant (1.0 Å by default), so each point is associated with an identical
cubical volume (e.g. 1.0 Å3). The volume of a whole pocket is calculated by
simply summing the individual volumes associated with each unique point.

POVME also optionally saves the pocket-occupying points of each frame to PDB
file(s) on this disk. The user can instruct the program to save these points
to separate files and/or to a single PDB trajectory. Some visualization
programs (e.g. VMD) are only compatible with trajectories that have the same
number of atoms in each frame. POVME can optionally write extra points to the
origin (0.0, 0.0, 0.0) on a frame-by-frame basis to satisfy this requirement.

Finally, POVME also optionally saves a volumetric density map in the Data
Explorer (DX) format. A volumetric density value is associated with each of
the pocket-occupying points by calculating the fraction of all trajectory
pocket volumes that include the given point. If the density map is displayed
as an isosurface, the value of the isosurface expresses the fraction of time
(e.g. over the course of the simulation) that the pocket included the
displayed volume.

(README.md adapted from [this
website](http://rocce-vm0.ucsd.edu/data/sw/hosted/POVME/).)