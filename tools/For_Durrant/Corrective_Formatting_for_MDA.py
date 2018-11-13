import MDAnalysis
from sys import argv
#Get VMD-annotated PDB from user through sys.argv
script, pdb_file = argv

print "This script will create a PDB /n with the correct formatting for MDA"
#Creates a new Universe with user input
u = MDAnalysis.Universe(pdb_file, topology_format = 'PDB')
#Writes a new trajectory in the format MDAnalysis can work with
protein = u.select_atoms("protein")
with MDAnalysis.Writer("Another.pdb", multiframe=True) as pdb:
    for ts in u.trajectory:
        pdb.write(protein)

