# Project Roadmap

This document describes planned updates to the POVME2 codebase.

## Web-Based POVME-Preparation GUI

**Projected Dates**: 6/2020-10/2020

**Description**: We will create a web-based graphical user interface (GUI) so
users can easily define and visualize pocket-encompassing protein regions for
subsequent POVME2 analysis. This process is currently tedious and requires
third-party programs such as VMD.

## POVME MD-Trajectory Handling

**Projected Dates**: 8/2020-12/2020

**Description**: We will update POVME2 to improve MD-trajectory handling.
Before starting POVME2, users must now turn to third-party programs to align
their trajectories and convert them to the multi-frame PDB format. With our
updates, POVME2 will accept a wider range of trajectory-file formats and will
optionally perform trajectory alignment itself.

## Port POVME Itself to the Browser

**Projected Dates**: 10/2020-2/2021

**Description**: We will port POVME2 to the browser. Web POVME2 will eliminate
the need to install Python and other dependencies. Instead, users will visit a
simple webpage. Thanks to WebAssembly (Pyodide), POVME2 calculations will run
directly in the browser, without requiring extensive remote computing
infrastructure.

## Online Tools to Analyze and Display POVME2 Results

**Projected Dates**: 1/2021-3/2021

**Description**: We will create online tools for analyzing and displaying
POVME2 results. Users must currently turn to third-party programs to visualize
POVME2 volume distributions, pocket shapes, etc. We will incorporate graphs,
tables, and 3Dmol.js-powered molecular visualizations directly into the POVME2
web app.

## Online Help System with Tutorials

**Projected Dates**: 4/2021-4/2021

**Description**: For each of our proposed web apps, we will create an online
help system together with new tutorials so users understand how the updated
software works.

## Publication

**Projected Dates**: 5/2021-6/2021

**Description**: To let the community know about our improvements, we will
publish peer-reviewed manuscripts describing POVME2 and BINANA updates.
