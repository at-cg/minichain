#!/usr/bin/env python3

import subprocess
import os
import multiprocessing
import re
import time
import sys
import getopt as getopt


par_threads = 6

## Pass arguments
argv = sys.argv[1:]

if(len(argv)==0):
    print("help: Map_Reads.py -h")
    sys.exit(2)

try:
    opts, args = getopt.getopt(argv, 't:')
except:
    print("usage: Map_Reads.py -t <threads>")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print("usage: Map_Reads.py -t <threads>")
        sys.exit()
    elif opt in ("-t"):
        par_threads = int(arg)
    
print("Mapping threads : " + str(par_threads))

total_threads = multiprocessing.cpu_count()
map_threads = total_threads/par_threads
map_threads = int(map_threads)

Reads = ['PacBio_I', 'PacBio_II', 'ONT']
Graph = 'Graphs/MHC-CHM13.0.gfa.gz'
R = ['0', '1000', '10000', '100000', '1000000', '2000000000']

Metadata = []
for read in Reads:
    if read == 'PacBio_I':
        Read = 'Reads/MHC_CHM13_' + read.split('_')[0] + '_filt.fq.gz'
    else:
        Read = 'Reads/MHC_CHM13_' + read.split('_')[0] + '_filt.part_001.fq.gz ' + ' Reads/MHC_CHM13_' + read.split('_')[0] + '_filt.part_002.fq.gz'
    for r in R:
        Metadata.append([read, r, Read])

# check if Mapped_Reads directory exists
if os.path.exists('Mapped_Reads'):
    os.system('rm -rf Mapped_Reads')
    os.system('mkdir Mapped_Reads')
if not os.path.exists('Mapped_Reads'):
    os.system('mkdir Mapped_Reads')

def Map_Reads(Metadata):
    read, r, Read = Metadata
    Output = 'Mapped_Reads/' + read + '_' + r + '.gaf'
    cmd = 'minichain -t' + str(map_threads) + ' -cx lr -b1 -R' + r + ' ' + Graph + ' ' + Read + ' > ' + Output
    out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    out = out.decode("utf-8")

    print(out)
    # add 6th line from out to count_recomb
    R = ""
    if read == 'PacBio_I':
        out = out.split('\n')
        val = re.findall(r'R: (\d+.\d+)', out[5])
        R = val[0]
    else:
        out = out.split('\n')
        val = re.findall(r'R: (\d+.\d+)', out[6])
        R = val[0]

    # write to file
    with open('Mapped_Reads/' + read + '_' + r + '.txt', 'w') as f:
        f.write('Read\tr\tR\n')
        f.write(read + '\t' + r + '\t' + R + '\n')

time_start = time.time()
# run in parallel
pool = multiprocessing.Pool(processes=par_threads)
pool.map(Map_Reads, Metadata)
pool.close()
pool.join()
time_end = time.time()
print('Mapping time: ' + str(time_end - time_start))