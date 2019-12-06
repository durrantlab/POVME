import matplotlib.pyplot as plt
import numpy as np
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


#Create a list of volume values, store in vol_list
vol_list = list(vol_dic.values())
#Sets bin parameters for the histogram
#bins = np.array([n for n in range(0, len(vol_list)#, 5)])

#Creates histogram, x = Volume, y = Group, parameters
plt.hist(vol_list, bins=20, alpha=0.75, histtype = 'bar', rwidth= .95, color = 'orange', align='mid')

# X axis label
plt.xlabel('Binding Pocket Volume')
# Y axis label
plt.ylabel('Frequency')
# Plot Title
plt.title('Frequency of Binding Pocket Volume')
# Show plot
plt.savefig("Histogram.png")