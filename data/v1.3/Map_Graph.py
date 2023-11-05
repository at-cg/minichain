#!/usr/bin/env python3
# author: Ghanshyam Chandra

import os
import time
import sys
import multiprocessing
import subprocess
import re
import getopt as getopt


threads = 48

## Pass arguments
argv = sys.argv[1:]

if(len(argv)==0):
    print("help: Map_Graph.py -h")
    sys.exit(2)

try:
    opts, args = getopt.getopt(argv, 't:')
except:
    print("help: Map_Graph.py -h")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print("help: Map_Graph.py -h")
        sys.exit()
    elif opt in ("-t"):
        threads = int(arg)
    
print("Threads : " + str(threads))

# read fasta file as G1, G2, G3, G4
ref = 'MHC-CHM13.0.fa'

mutation_rates = ['0.1', '1', '5']

for m in mutation_rates:
    fasta_files = []
    fa_dir = "Query_" + m
    align_dir = "align_" + m
    # remove align folder
    if os.path.exists(align_dir):
        os.system("rm -rf " + align_dir)
    os.system("mkdir " + align_dir)
    # read fasta files from the directory
    for file in os.listdir(fa_dir + '/'):
        if file.endswith('.fa'):
            if file == ref:
                continue
            else:
                fasta_files.append(fa_dir + "/" + file)


    # add timer for execution
    start_time = time.time()

    # read gfa files
    graph_dir = "Graphs"

    count_recomb = {}
    align_ids = []
    R = ['0', '1000', '10000', '100000', '1000000', '2000000000']
    # R = ['10000']

    h_gfa = []
    for i in fasta_files:
        for r in R:
            h_gfa.append((i, r))
        

    def align_gfa(idx):
        idx_, r = idx
        graph = graph_dir + "/" + 'MHC-CHM13.0.gfa'
        h = idx_

        # read walks from graph
        print("Reading Walks from Graph ...")
        walks = []
        w_lines = {}
        haps = {}
        with open(graph, 'r') as f:
            for line in f:
                if line.startswith('W'):
                    field = line.split('\t')
                    id = field[1]
                    if id not in walks:
                        walks.append(id)
                        w_lines[id] = field[6]
                        walk = field[6].split('>')[1:]
                        for h_ in walk:
                            h_ = h_.split('\n')[0]
                            if h_ not in haps:
                                haps[h_] = []
                            haps[h_].append(id)
                    else:
                        id_1 = id + ".1"
                        id_2 = id + ".2"
                        walks.append(id_1)
                        walks.append(id_2)
                        walks.remove(id)
                        w_lines[id_1] = w_lines[id]
                        w_lines[id_2] = field[6]
                        del w_lines[id]
                        # del haps[id]
                        # haps[id_1] = []
                        # haps[id_2] = []
                        walk_1 = w_lines[id_1].split('>')[1:]
                        walk_2 = w_lines[id_2].split('>')[1:]
                        for h_ in walk_1:
                            h_ = h_.split('\n')[0]
                            if h_ not in haps:
                                haps[h_] = []
                            haps[h_].append(id_1)
                        for h_ in walk_2:
                            h_ = h_.split('\n')[0]
                            if h_ not in haps:
                                haps[h_] = []
                            haps[h_].append(id_2)
        
        # check length of walks is same as w_lines
        for walk in walks:
            if walk not in w_lines:
                print("Error: walk not in w_lines")
                sys.exit(2)
        
        # verify haps contain all the vertices in graphs
        vertices = []
        edges = []
        with open(graph, 'r') as f:
            for line in f:
                if line.startswith('S'):
                    field = line.split('\t')
                    vertices.append(field[1])
                if line.startswith('L'):
                    field = line.split('\t')
                    edges.append((field[1], field[3]))

        for v in vertices:
            if v not in haps:
                print("Error: vertex not in haps")
                sys.exit(2)

        for u, v in edges:
            if u not in haps or v not in haps:
                print("Error: edge not in haps")
                sys.exit(2)

        print("Graph : " + graph + " Number of Verices : " + str(len(haps)))


        print ("Generating Alignment Files ...")

        # count number of contig switches in a walk
        def count_switches(walk):
            walk = walk.split('>')[1:]
            # Initialize C
            C = {}
            for w in walks:
                w = w.split('\n')[0]
                C[w] = []
                for j in range(len(walk)):
                    C[w].append(2e9)

            # fill for j = 0
            for h_ in haps[walk[0]]:
                C[h_][0] = 0

            # Do DP
            for j in range(1, len(walk)):
                for h_ in haps[walk[j]]:
                    for h__ in haps[walk[j-1]]:
                        if h__ == h_:
                            C[h_][j] = min(C[h_][j], C[h__][j-1])
                        else:
                            C[h_][j] = min(C[h_][j], C[h__][j-1] + 1)

            # find min in last column
            min_ = 2e9
            for h_ in haps[walk[-1]]:
                min_ = min(min_, C[h_][-1])

            return min_ # count of minimum recombination

        print("Processing R = " + r + " and " + h)
        cmd = "bin/minichain -t1 -cx lr " + graph + " " + h + " -R" + r + " --vc -o " + align_dir + "/R"+ r + "_" + h.split('/')[1].split('.fa')[0] + ".gaf"
        # cmd = "minigraph -t1 -cx lr " + graph + " " + h + " > align/R"+ r + "_" + h.split('/')[1].split('.fasta')[0] + ".gaf"
        out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        out = out.decode("utf-8")
        print(out)
        # add 5th line from out to count_recomb
        out = out.split('\n')
        min_ = re.findall(r'\bMin: (\d+)\b', out[5])
        max_ = re.findall(r'\bMax: (\d+)\b', out[5])
        mean_ = re.findall(r'\bMean: (\d+)\b', out[5])
        accuracy_ = re.findall(r'\bAccuracy: (\d+\.\d+)\b', out[5])
        recomb_ = 0
        frac_seq = 0.0
        identity = 0.0
        true_recomb_ = 0
        # read gaf file
        with open(align_dir +"/R" + r + "_" + h.split('/')[1].split('.fa')[0] + ".gaf", 'r') as f:
            for line in f:
                field = line.split('\t')
                true_recomb_ = int(field[0].split('!')[1])
                recomb_ = count_switches(field[5])
                frac_seq = float((int(field[3]) - int(field[2]))/int(field[1])) * 100
                identity = float(int(field[9])/int(field[10])) * 100
                break

        # avoid thread collision
        with open(align_dir + "/R" + r + "_" + h.split('/')[1].split('.fa')[0] + ".txt", 'w') as f:
            f.write("File\tMin\tMax\tMean\tRecomb\tFracSeq\tIdentity\tTrueRecomb\tAccuracy\n")
            f.write(h.split('/')[1].split('.fa')[0]  + "\t" + str(min_[0]) + "\t" + str(max_[0]) + "\t" + str(mean_[0]) + "\t" + str(recomb_) + "\t" + str(frac_seq) + "\t" + str(identity) + "\t" + str(true_recomb_) + "\t" + str(accuracy_[0]) + "\n")



    # run in parallel and print the timeß
    print("Generating Alignment Files ...")
    start_time = time.time()
    x = int(threads/len(R)) # number of threads
    pool = multiprocessing.Pool(processes= x * len(R))
    pool.map(align_gfa, h_gfa)
    print("Time taken to generate alignment files : " + str(time.time() - start_time))

