import MDAnalysis as mda
from sys import argv

#Get Multimodel PDB file from user through sys.argv
script, pdb_file = argv

#Get number of the POVME frame from User
frame_num = raw_input("What is the frame number you want? ")
frame_num_int = int(frame_num)
frame_num_str = str(frame_num)

#Defining the universe object with the PDB file
u = mda.Universe(pdb_file)
#This designates the atom group(the protein) and position information
#Must be placed before Writer, not after!
protein = u.select_atoms("protein"), u.trajectory[frame_num_int].positions

#Feeds the necessary information into the PDB writer to create a PDB
with mda.Writer("MD_frame_" + frame_num_str + ".pdb",
				start = frame_num_int,
				multiframe = False) as pdb:
						pdb.write(u)
print "Created new PDB file for specified frame, titled 'MD_frame_%r.pdb'" %frame_num
print "Your POVME frame is titled 'POVME_frame_%r.pdb' " % frame_num