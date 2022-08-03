#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np
import re

def overlap(read_metadata,map_region):
    region = read_metadata.split("!")[1] # This is a pattern lets say "chr1"
    X = int(read_metadata.split("!")[2])
    Y = int(read_metadata.split("!")[3])
    # Extract the co-ordinates assuming region as a pattern
    pattern = region+':\d+-\d+'
    x = re.findall(pattern,map_region)
    X_ = []
    Y_ = []
    for i in range(len(x)):
        X_.append(int(re.findall('\d+-\d+',x[i])[0].split("-")[0]))
        Y_.append(int(re.findall('\d+-\d+',x[i])[0].split("-")[1]))
    ## Compare here
    olp = 0.0
    for i in range(len(x)):
        olp = olp + (max(0, min(Y, Y_[i]) - max(X, X_[i]))/(Y-X))
    threshold = 0.10
    # print("threshold : ", threshold, " olp : ",olp)
    if olp >= threshold:
        return True
    else:
        return False

graphs = ["Linear","10H","20H","40H","60H","80H","95H"]
wrong_count = []
wrong_reg = []
wrong_path = []
read_count = []
for graph in graphs:
    gfa_ = "wrong_mappings_"+graph+".gaf"
    wrong_gaf = []
    with open(gfa_) as file:
        for line in file:
            wrong_gaf.append(line)
    count_ = "count_mappings_"+graph+".gaf"
    read_count_ = 0
    with open(count_) as file:
        for line in file:
            read_count_ = read_count_ + int(line)
    read_count.append(read_count_)
    count = 0
    for i in range(len(wrong_gaf)):
        read_region = wrong_gaf[i].split("\t")[1].split("!")[1]
        map_region = wrong_gaf[i].split("\t")[2]
        read_metadata = wrong_gaf[i].split("\t")[1]
        #if read_region == map_region :  # No match than wrong region
        #    count = count + 1
        if overlap(read_metadata,map_region) == True: # No overlap then wrong region else wrong path
            count = count + 1
    wrong_count.append(len(wrong_gaf))
    wrong_path.append(count)
wrong_count = np.asarray(wrong_count,dtype=float)
wrong_path = np.asarray(wrong_path,dtype=float)
wrong_count = wrong_count/read_count
wrong_path = wrong_path/read_count

plt.bar(graphs,wrong_count,label='Wrong region',zorder=3)
plt.plot(graphs,wrong_count,marker='x',color='r',zorder=4)
plt.bar(graphs,wrong_path,label='Wrong path',zorder=3)
plt.plot(graphs,wrong_path,marker='x',color='r',zorder=4)
plt.grid(zorder=0)
plt.legend()
plt.xlabel('Graphs')
plt.ylabel('Fraction of wrong read mapping')
plt.savefig("bar_plt.png",dpi=1200)
