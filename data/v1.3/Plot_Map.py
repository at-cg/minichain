#!/usr/bin/env python3

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from rich import print as rprint
from rich.console import Console
from rich.table import Table
from pylatexenc.latex2text import LatexNodes2Text


Reads = ['PacBio_I', 'PacBio_II', 'ONT']
R = ['0', '1000', '10000', '100000', '1000000', '2000000000']
R_labels = ['$0$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$\infty$']
R_labels_ = [LatexNodes2Text().latex_to_text(r) for r in R_labels]

R_ = dict()

# Read data from file
for read in Reads:
    R_[read] = []
    for r in R:
        with open('Mapped_Reads/' + read + '_' + r + '.txt', 'r') as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == 'Read':
                    continue
                R_[read].append(float(line[2]))

# Plot a single bar for R, R+NR_R, R+NR_R+NR_NR for each read
print('Plotting')
os.system('rm Mapped_Reads/Stats_table.txt')
for read in Reads:
    x = np.arange(len(R))
    fig, ax = plt.subplots()
    # Increase the figure size
    fig.set_size_inches(6, 4)
    
    table = Table(title="Reads: " + read)
    if read == 'PacBio_I':
        table = Table(title="Reads: PacBio Sequel II (SRR11292121)")
    elif read == 'PacBio_II':
        table = Table(title="Reads: PacBio Sequel II (SRX5633451)")
    elif read == 'ONT':
        table = Table(title="Reads: Oxford Nanopore (SRR23365080)")
    table.add_column("Recombination Penalty", justify="right", style="cyan", no_wrap=True)
    table.add_column("Complete support", justify="right", style="blue", no_wrap=True)
    for i in range(len(R)):
        table.add_row(R_labels_[i], str(R_[read][i]))
    console = Console(record=True)
    console.print(table, justify="center")

    # save table to pdf file
    console.save_svg("Mapped_Reads/" + read + "_table.svg", title = "Haplotype-aware Chaining")
    with open('Mapped_Reads/Stats_table.txt', 'a') as f:
        rprint(table, file=f)

    val_1 = []
    for i in range(len(R)):
        val_1.append(R_[read][i])
    ax.bar(x, val_1, label='Complete support', zorder = 3)

    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.set_ylabel('Chains supported by reads', fontsize=15)
    # ax.set_title(read + ' reads')
    # puy top y limit as 1
    ax.set_ylim(top=1)

    ax.set_xticks(x)
    ax.set_xticklabels(R_labels, rotation=45)
    ax.set_xlabel('Recombination penalty', fontsize=15)
    plt.grid(True)
    # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # put legend on the top of the plot
    # ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', borderaxespad=0., ncol=3, fontsize=14)
    # # reverse the legend order
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.5, 1.2), loc='upper center', borderaxespad=0., ncol=3)

    fig.tight_layout()
    # scale font size by factor of 1.3
    # plt.rcParams.update({'font.size': 12})
    plt.savefig('Mapped_Reads/' + read + '.pdf')
    plt.close()