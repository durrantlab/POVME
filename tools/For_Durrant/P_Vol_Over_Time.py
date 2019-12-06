import matplotlib.pyplot as plt
from sys import argv
script, filename = argv 


volumes = open(filename)

#This script creates an ordered dictionary vol_dic of the volumes in the data file.
#Variable names are arbitrary, just used to build the structure of the dictionary
import collections
vol_dic = collections.OrderedDict()
for line in volumes:		
	x = line.split(',')     
	y = line.split('\t')
	a = y[0]
	b = y[1]
	c = len(b) - 1
	b = b[0:c]
	vol_dic[a] = float(b)

#Creates a list of the volumes from vol_dic
vol_list = list(vol_dic.values())  
#Creates a list of the frames from vol_dic
frame_num = list(vol_dic.keys())   
#Creates scatterplot, x = frames, y = volumes
plt.scatter(frame_num, vol_list, marker = '.', s = 50, c = 'orange')
# X axis label
plt.xlabel('Frame Number')
# Y axis label
plt.ylabel('Binding Pocket Volume')
# Plot Title
plt.title('Binding Pocket Volume Over Time')
#Sets parameters of X axis for viewing
plt.xticks(list(range(-1, len(frame_num), (len(frame_num)/10))))
plt.xlim([- 1, len(frame_num)])
# Shows plot
plt.savefig('Scatterplot.png')