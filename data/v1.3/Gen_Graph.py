#!/usr/bin/env python3
# author: Ghanshyam Chandra

import os
import time
import sys
import multiprocessing
import subprocess
import getopt as getopt


threads = 48

## Pass arguments
argv = sys.argv[1:]

if(len(argv)==0):
    print("help: Gen_Graph.py -h")
    sys.exit(2)

try:
    opts, args = getopt.getopt(argv, 't:')
except:
    print("usage: Gen_Graph.py -t <threads>")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print("usage: Gen_Graph.py -t <threads>")
        sys.exit()
    elif opt in ("-t"):
        threads = int(arg)
    
print("Threads : " + str(threads))

# create directory for graphs
graph_dir = "Graphs"
if not os.path.exists(graph_dir):
    os.makedirs(graph_dir)

ref = 'MHC-CHM13.0.fa'
fa_dir = "Genomes"
# read fasta files from "Genomes" folder
fasta_files = []
for filename in os.listdir(fa_dir):
    if filename.endswith(".fa"):
        if filename == ref: # don't add reference
            continue
        fasta_files.append(fa_dir + "/" + filename)

ref_ = ref.split('.fa')[0]
ref = fa_dir + "/" + ref
fasta_files.append(ref)
# add timer for execution
fasta_gfa_2 = []
data_idx = []
fasta_gfa_2.append((ref_, set(fasta_files)))
data_idx.append(0)
print("Number of fasta files for " + ref_ + " : " + str(len(fasta_files)))


start_time_ = time.time()
mash_file = []
h = []
haplotypes = []
fasta_file = []
tuple_list = []
graph = []
for i in range(len(fasta_gfa_2)):
    mash_file.append(" ")
    h.append(" ")
    graph.append(" ")
    haplotypes.append([])
    fasta_file.append([])
    tuple_list.append([])


def gen_mg_gfa(idx):
    # generate graph in parallel
    h[idx] = fasta_gfa_2[idx][0]
    haplotypes[idx] = fasta_gfa_2[idx][1]

    # add timer for execution
    start_time = time.time()
    
    fasta_file[idx] = []
    for hap in haplotypes[idx]: # reference is already added
        fasta_file[idx].append(hap)

    # sort fasta_file based on mash distance
    mash_file[idx] = "mash_"+ h[idx].split('.fa')[0] +".txt"
    if os.path.exists(mash_file[idx]):
        os.remove(mash_file[idx])

    # # sort fasta_file based on mash distance
    # def get_mash_distance(fasta):
    #     cmd = "mash dist " + ref + " " + fasta + " | grep -v '^#' | head -n 1 | cut -f 3"
    #     out = subprocess.check_output(cmd, shell=True)
    #     out = out.decode("utf-8")
    #     # write to file
    #     with open("mash.txt", 'a') as f:
    #         f.write(fasta + "\t" + out)


    # print computating mash distance
    print("Computing Mash Distance ...")
    for fasta in fasta_file[idx]:
        cmd = "bin/mash dist " + ref + " " + fasta + " | grep -v '^#' | head -n 1 | cut -f 3"
        out = subprocess.check_output(cmd, shell=True)
        out = out.decode("utf-8")
        # write to file
        with open(mash_file[idx], 'a') as f:
            f.write(fasta + "\t" + out)
    # run in parallel
    # num_processes = len(fasta_file)  # Use the number of available CPU cores
    # pool = multiprocessing.Pool(processes=num_processes)
    # pool.map(get_mash_distance, fasta_file)

    # read mash.txt and sort fasta_file based on distance
    fasta_file[idx] = []
    tuple_list[idx] = []
    with open(mash_file[idx], 'r') as f:
        for line in f:
            field = line.split('\t')
            fasta = field[0]
            distance = float(field[1].split('\n')[0])
            tuple_list[idx].append((distance, fasta))

    # sort tuple_list by score
    tuple_list[idx].sort(key=lambda tup : tup[0])

    # add fasta_file in sorted order
    print("Order of fasta files based on mash distance")
    for tup in tuple_list[idx]:
        print(tup[1] + "\t" + str(tup[0]))
        fasta_file[idx].append(tup[1])

    # remove mash.txt
    os.remove(mash_file[idx])

    graph[idx] = graph_dir + "/" + h[idx].split('.fa')[0] + ".gfa"

    cmd = "bin/minigraph -t1 -cxggs " + ref
    for file_name in fasta_file[idx]:
        if file_name == ref: # sanity check
            continue
        cmd += " " + file_name 
    cmd += " > " + graph[idx]

    print('Generating Graph from ' + str(len(fasta_file[idx])) + ' files')

    os.system(cmd)

    print("--- Time : %s seconds ---" % (time.time() - start_time) + " for " + h[idx].split('.fa')[0] + ".gfa")


# generate graph in parallel
num_processes = 48  # Use the number of available CPU cores
pool = multiprocessing.Pool(processes=num_processes)
pool.map(gen_mg_gfa, data_idx)

# add w lines in parallel


for h, haplotypes in fasta_gfa_2:

    # add timer for execution
    start_time = time.time()
    fasta_file = []

    for hap in haplotypes: # reference is already added
        fasta_file.append(hap)

    graph = graph_dir + "/" + h.split('.fa')[0] + ".gfa"

    # exctract reference walk with vg
    ref_walk = ""
    ref_id = ""
    cmd = 'bin/vg convert -g -f ' + graph + ' -W | grep   "P" | head -n 1'
    out = subprocess.check_output(cmd, shell=True)
    out = out.decode("utf-8")
    out = out.split('\t')
    ref_id = out[1]
    temp_walk = out[2].split(',')
    for walk in temp_walk:
        ref_walk +=  ">s" + walk.split('+')[0]


    # sort fasat 
    # perform all to all alignment
    def process_fasta(fasta):
        # cmd = "minigraph -t1 -cx lr --vc " + graph + " " + fasta + " > " +  fasta.split('.fasta')[0] + ".gaf"
        cmd = "bin/minigraph -t1 -cx asm --vc " + graph + " " + fasta + " > " +  fasta.split('.fa')[0] + ".gaf"
        # cmd = "GraphAligner -t1 -x vg -g " + graph + " -f " + fasta + " -a " +  fasta.split('.fasta')[0] + ".gaf"
        # print(cmd)
        os.system(cmd)

    num_processes = threads  # Use the number of available CPU cores
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(process_fasta, fasta_file)

    # read all gaf files
    gaf_files = []
    # List all files in the directory
    files = os.listdir(fa_dir)
    for file_name in files:
        if file_name.endswith(".gaf"):
            gaf_files.append(fa_dir + "/" + file_name)

    # generate w lines
    w_lines = []
    walks = []
    i = 0
    for gaf in gaf_files:
        with open(gaf, 'r') as f:
            for line in f:
                field = line.split('\t')
                if i == 0:
                    print("Pushing Reference Walk ...")
                    # id_ = ref_id.split('#')
                    id = ref_id
                    walks.append(id)
                    walk_ = ref_walk
                    w_lines.append("W\t" + id + "\t1\t" + "MHC" + "\t1\t1\t" + walk_ + "\n")
                    i += 1
                if ref_id == field[0]: # if reference then just push exact walk
                    continue
                else: # else push the walk from gaf
                    is_ref = 0
                    id_ = field[0].split('#')
                    # id = id_[0] + "." +id_[1]
                    id = field[0]
                    walks.append(id)
                    walk = field[5]
                    w_lines.append("W\t" + id + "\t1\t" + "MHC" + "\t1\t1\t" + walk + "\n")
                # break # we just want one walk
            # if is_ref != 1:
            #     walk_ = set(walk.split('>')[1:])
            #     temp_walk = ""
            #     for w in walk_:
            #         w = w.split('\n')[0]
            #         temp_walk += ">" + w
            #     w_lines.append("W\t" + id + "\t1\t" + "MHC" + "\t1\t1\t" + temp_walk + "\n")


    # # generate w lines
    # w_lines = []
    # walks = []
    # i = 0
    # for gaf in gaf_files:
    #     with open(gaf, 'r') as f:
    #         for line in f:
    #             field = line.split('\t')
    #             id_ = field[0].split('#')
    #             id = id_[0] + "." +id_[1]
    #             walks.append(id)
    #             walk = field[5]
    #             w_lines.append("W\t" + id + "\t1\t" + "MHC" + "\t1\t1\t" + walk + "\n")
                    


    # check if count gaf is equal to count w lines
    print("Number of GAF Files: " + str(len(gaf_files)))
    print("Number of W Lines: " + str(len(w_lines)))
    if len(gaf_files) == len(w_lines):
        print("Contigous Walks Generated Successfully")

    # append w lines to gfa
    print("Appending W Lines to GFA ...")
    with open(graph, 'a') as f:
        for line in w_lines:
            f.write(line)

    # Correct graph
    s_lines = []
    vertices = []
    edges = []
    l_lines = []
    w_lines = []
    haps = {}
    walks = []
    with open(graph, 'r') as f:
        for line in f:
            if line.startswith('S'):
                field = line.split('\t')
                id = field[1]
                # haps[id] = [] # used to debug and found \n is the cause
                vertices.append(id)
                s_lines.append(line)
            elif line.startswith('L'):
                field = line.split('\t')
                id_1 = field[1]
                id_2 = field[3]
                edges.append((id_1, id_2))
                l_lines.append(line)
            elif line.startswith('W'):
                field = line.split('\t')
                id = field[1]
                walks.append(id)
                w_lines.append(line)
                walk = field[6].split('>')[1:]
                for w in walk:
                    w = w.split('\n')[0]
                    if w not in haps:
                        haps[w] = []
                    haps[w].append(id)

    for v in vertices:
        if v not in haps:
            print("Error: Vertex " + v + " not found in walks and will be removed")

    for i in range(len(w_lines)):
        walk = w_lines[i]
        field = walk.split('\t')
        id = field[1]
        walk = field[6].split('>')[1:]
        for w in walk:
            w = w.split('\n')[0]
            if id not in haps[w]:
                print("Error: Walk " + id + " not found in haps[" + w + "]")
                sys.exit(2)
    
    # print( "haps is correct")


    # remove redundent vertices from s_lines and l_lines
    s_lines_ = []
    l_lines_ = []
    for line in s_lines:
        field = line.split('\t')
        id = field[1]
        if id in haps:
            s_lines_.append(line)
    for line in l_lines:
        field = line.split('\t')
        id_1 = field[1]
        id_2 = field[3]
        if id_1 in haps and id_2 in haps:
            l_lines_.append(line)
        

    # read every walk and check if it is valid
    for line in s_lines_:
        field = line.split('\t')
        id = field[1]
        if id in haps:
            continue
        else:
            print("Error: Vertex " + id + " not found in walks")
            sys.exit(2)
        
    for line in l_lines_:
        field = line.split('\t')
        id_1 = field[1]
        id_2 = field[3]
        if id_1 in haps and id_2 in haps:
            continue
        else:
            print("Error: Edge " + id_1 + "->" + id_2 + " not found in walks")
            sys.exit(2)
    
    print("Correction of Graph Done ...")

    # delete graph and rewrite
    os.system('rm ' + graph)

    with open(graph, 'w') as f:
        for line in s_lines_:
            f.write(line)
        for line in w_lines:
            f.write(line) 
        for line in l_lines_:
            f.write(line)

    
    # check length of walks is same as w_lines
    if len(walks) != len(w_lines):
        print("Error: walks and w_lines length not same")
        sys.exit(2)

    print("Graph : " + graph + " Number of Vertices : " + str(len(haps)))

    # remove all gaf files
    for gaf in gaf_files:
        os.remove(gaf)


    print("Walks for " + graph + " generated successfully")
    print("--- Time : %s seconds ---" % (time.time() - start_time))

# compress all gfa files
os.system("cd Graphs && gzip *.gfa")

print("--- Total Graph Generation Time : %s seconds ---" % (time.time() - start_time_))