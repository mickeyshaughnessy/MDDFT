#!/usr/bin/env python

#This script takes a toy database and computes the number of connected components. 

d = [0.5,0.1,0.5,0.1,0.7,0.1,0.5,0.5,0.2,0.7]  #d01,d02,d03,d04,d12,d13,d14,d23,d24,d34
clusters = [2,5,1,6,8] #dummy
def pair(i,j):
	if (i==0 and j==1):
		p = 0
	if (i==0 and j==2):
		p = 1
	if (i==0 and j==3):
		p = 2
	if (i==0 and j==4):
		p = 3
	if (i==1 and j==2):
		p = 4
	if (i==1 and j==3):
		p = 5
	if (i==1 and j==4):
		p = 6
	if (i==2 and j==3):
		p = 7
	if (i==2 and j==4):
		p = 8
	if (i==3 and j==4):
		p = 9
	return p
keys = []
dcut = 0.11
for i in range(len(clusters)):
	keys.append(i)
equiv = keys
changed = 1
while (changed == 1):
	changed = 0
	print 'equiv =', equiv
	for i in range(0,len(keys)-1):
		for j in range(i+1,len(keys)):
			print 'i = ', i, ' j= ', j
			if (d[pair(i,j)] < dcut):
				if(equiv[j] != equiv[i]):
					equiv[j] = equiv[i]
					changed = 1
print len(set(equiv))		

