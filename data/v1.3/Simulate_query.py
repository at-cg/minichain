#!/usr/bin/env python3

import multiprocessing
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
import networkx as nx

node_len = {}
node_seq = {}
GFA  = "Graphs/MHC-CHM13.0.gfa"
G = nx.DiGraph()
Walks = {}
haps_ = {}
# read GFA file
with open(GFA, "r") as f:
    for line in f.readlines():
        if line[0] == "S":
            node = line.split("\t")[1]
            len_ = len(line.split("\t")[2])
            node_len[node] = len_
            node_seq[node] = line.split("\t")[2].split("\n")[0]
            G.add_node(node, len=len_)
        elif line[0] == "L":
            node1 = line.split("\t")[1]
            node2 = line.split("\t")[3]
            G.add_edge(node1, node2)
        elif line[0] == "W":
            id = line.split("\t")[1]
            walk = line.split("\t")[6].split(">")[1:]
            if id not in Walks:
                Walks[id] = []
            for w in walk:
                w = w.split("\n")[0]
                Walks[id].append(w)
                if w not in haps_:
                    haps_[w] = []
                haps_[w].append(id)

topo_sort = list(nx.topological_sort(G))
num_walks = len(Walks.keys())




mutation_rates = ['0.1', '1', '5']

for m in mutation_rates:

    querys = {}
    # if Query folder is there delete it and create a new one
    folder = "Query_" + m
    if os.path.isdir(folder):
        os.system("rm -r " + folder)
        os.system("mkdir " + folder)
    else:
        os.system("mkdir " + folder)


    min_len = 1e9
    max_len = -1
    for walk in Walks:
        sum = 0
        for node in Walks[walk]:
            sum += node_len[node]
        if sum > max_len:
            max_len = sum
        if sum < min_len:
            min_len = sum

    print("Min len: " + str(min_len) + " Max len: " + str(max_len))
        


    num_fasta = 45
    # M =[m] # Mutation rate
    idx = []
    j = 0
    for i in range(num_fasta):
        idx.append((float(m), i, j))
        j += 1
        if j == 3:
            j = 0
        

    def sim_query(idx_):
        mutation, i, k = idx_
        lens  = [1e6, 2e6, 3e6]
        len_q = lens[k]
        # len_q = random.choice(lens)
        f_s = []
        walk_id = random.choice(list(Walks.keys()))
        query = ""
        f_s.append(walk_id)
        walk = Walks[walk_id]
        round = 1
        prev_len = -1
        while len(query) - prev_len > 0:
            prev_len = len(query)
            for v in walk:
                query += node_seq[v]
                if len(query) >= round * len_q:
                    break
            # print("Round : " + str(round) + " Q_len : " + str(len(query)))
            if v == walk[len(walk) - 1]:
                break
            round += 1
            prev_walk_id = walk_id
            walk_id = random.choice(haps_[v])
            while walk_id == prev_walk_id:
                walk_id = random.choice(haps_[v])
            f_s.append(walk_id)
            walk = Walks[walk_id]
            # find idx of v in walk
            idx = 0
            for j in range(len(walk)):
                if walk[j] == v:
                    idx = j
                    break
            walk = walk[(idx + 1):]
            

        idxs = []
        for j in range(1, len(query) - 100):
            idxs.append(j)
        random.shuffle(idxs)

        # Deterministic switch
        count_snps = (mutation * len(query))/100
        snps_idxs = random.sample(idxs, int(count_snps))

        haps = ""
        for f in f_s:
            haps += ">" + f
        print("Number of switches: " + str(round - 1) + " for query " + str(i) + " with Mutation Rate = " + str(mutation) + " Haps: " + str(haps))
        # mutate query at snps_idxs positions
        query = list(query)
        for snp in snps_idxs:
            if query[snp] == "A":
                query[snp] = random.choice(["C", "G", "T"])
            elif query[snp] == "C":
                query[snp] = random.choice(["A", "G", "T"])
            elif query[snp] == "G":
                query[snp] = random.choice(["A", "C", "T"])
            elif query[snp] == "T":
                query[snp] = random.choice(["A", "C", "G"])
        query = "".join(query)
        seq_record = SeqIO.SeqRecord(Seq(query), id="query_" + str(mutation) + "_" + str(i) + "!" + str(round - 1) + "!" + str(haps), description=str(len(query)))
        SeqIO.write(seq_record, folder + "/query_" + str(mutation) + "_" + str(i)+".fa", "fasta")
        

    # use multiprocessing to create querys in parallel
    with multiprocessing.Pool() as pool:
        pool.map(sim_query, idx)