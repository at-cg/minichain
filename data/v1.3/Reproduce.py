#!/usr/bin/env python3

import getopt as getopt
import os
import sys

threads = 48

## Pass arguments
argv = sys.argv[1:]

if(len(argv)==0):
    print("help: Reproduce.py -h")
    sys.exit(2)

try:
    opts, args = getopt.getopt(argv, 't:')
except:
    print("help: Reproduce.py -h")
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print("help: Reproduce.py -h")
        sys.exit()
    elif opt in ("-t"):
        threads = int(arg)
    
print("Threads : " + str(threads))

# Download data and executables
# if bin folder exists then delete it or else create it
if os.path.exists("bin"):
    os.system("rm -rf bin")
os.system("mkdir bin")

os.system("wget -O MHC.agc https://zenodo.org/records/6617246/files/MHC-61.agc?download=1")
os.system("wget https://github.com/refresh-bio/agc/releases/download/v3.0/agc-3.0_x64-linux.tar.gz")
os.system("tar -xvf agc-3.0_x64-linux.tar.gz")
os.system("mv agc-3.0_x64-linux/agc bin/")
os.system("rm -rf agc-3.0_x64-linux*")

# Get vg and minigraph and put in bin folder
os.system("wget https://github.com/vgteam/vg/releases/download/v1.52.0/vg")
os.system("chmod +x vg")
os.system("mv vg bin/")

os.system("git clone https://github.com/lh3/minigraph")
os.system("cd minigraph && make -j")
os.system("cd ..")
os.system("mv minigraph/minigraph bin/")
os.system("rm -rf minigraph")

# Get mash 
os.system("wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar")
os.system("tar -xvf mash-Linux64-v2.3.tar")
os.system("mv mash-Linux64-v2.3/mash bin/")
os.system("rm -rf mash-Linux64-v2.3*")

# Get minichain
os.system("git clone https://github.com/at-cg/minichain")
os.system("cd minichain && make -j")
os.system("cd ..")
os.system("mv minichain/minichain bin/")
os.system("cd minichain && git checkout v1.0 && make -j")
os.system("cd ..")
os.system("mv minichain/minichain bin/minichain_10")
os.system("rm -rf minichain")

os.system("chmod +x bin/*")

# if Genomes folder exists then delete it or else create it
if os.path.exists("Genomes"):
    os.system("rm -rf Genomes")
os.system("mkdir Genomes")

# Extract the genomes from MHC.agc
os.system("bin/agc getcol -o Genomes/ MHC.agc")
os.system("rm MHC.agc")
os.system("rm Genomes/MHC-00GRCh38.fa")


# create conda environment named MC and install python packages numpy, scipy, matplotlib and networkx Biopython getopt seaborn pandas
# check if conda environment named MC exists or not
os.system("source ~/.bashrc && conda create --force -n MC -y && conda activate MC && conda install -c conda-forge -y numpy scipy matplotlib networkx biopython seaborn pandas rich pylatexenc")

map_threads = 6
# Generate the graph
os.system("python3 Gen_Graph.py -t " + str(threads))
# Simulate queries
os.system("source ~/.bashrc && conda activate MC && python3 Simulate_query.py -t " + str(threads))
# Map the queries
os.system("source ~/.bashrc && conda activate MC && python3 Map_Graph.py -t " + str(threads))
# Map the reads
os.system("source ~/.bashrc && conda activate MC && python3 Map_Reads.py -t " + str(map_threads))
# Plot the results
os.system("source ~/.bashrc && conda activate MC && python3 Plot.py")
# Plot the results for mapping
os.system("source ~/.bashrc && conda activate MC && python3 Plot_Map.py")

