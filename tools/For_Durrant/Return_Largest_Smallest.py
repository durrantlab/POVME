from sys import argv
script, filename = argv 

volumes = open(filename)
import collections

#This script creates an ordered dictionary vol_dic of the volumes in the data file.
#Variable names are arbitrary, just used to build the structure of the dictionary
vol_dic = collections.OrderedDict()
for line in volumes:		
	x = line.split(',')     
	y = line.split('\t')
	a = y[0]
	b = y[1]
	c = len(b) - 1
	b = b[0:c]
	vol_dic[a] = float(b)

#Finds maximum value in volumes and stores as max_val
max_val = 0
for val in vol_dic.values():
	if val > max_val:
		max_val = val

#Finds the frame of the maximum value and stores as max_val_frame
for frame, volume in vol_dic.items():
	if volume == max_val:
		max_val_frame = frame

#Finds minimum value in volumes and stores as min_val
min_val = min(vol_dic.values())

#Finds frame of minimum volume and stores as min_val_frame

for frame, volume in vol_dic.items():
	if volume == min_val:
		min_val_frame = frame
		

print """
The maximum pocket volume of this data set is %r, at frame %r,
while the minimum pocket volume of this data set is %r, at frame %r.
 		""" % (max_val, max_val_frame, min_val, min_val_frame)
	


