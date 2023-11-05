#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import stats
import seaborn as sns
import pandas as pd


mutation_rates = ['0.1', '1', '5']


for m in mutation_rates:

    folder = 'align_' + m

    Query_folder = 'Query_' + m

    fasta_ids = []
    # read gfa files
    for file in os.listdir(Query_folder):
        if file.endswith('.fa'):
            if file != 'MHC-CHM13.0.fa':
                fasta_ids.append(file)

    R = ['0', '1000', '10000', '100000', '1000000', '2000000000']
    # R = ['10000']
    count_recomb = {}
    # read files for min, max, mean, recomb and add to count_recomb
    for r in R:
        count_recomb[r] = {}
        for fa_id in fasta_ids:
            count_recomb[r][fa_id.split('.fa')[0]] = []
            with open(folder + "/R" + r + "_" + fa_id.split('.fa')[0] + ".txt", 'r') as f:
                for line in f:
                    if line.startswith("File"):
                        continue
                    field = line.split('\t')
                    align = field[0]
                    count_recomb[r][align].append(int(field[1]))
                    count_recomb[r][align].append(int(field[2]))
                    count_recomb[r][align].append(int(field[3]))
                    count_recomb[r][align].append(int(field[4]))
                    count_recomb[r][align].append(float(field[5]))
                    count_recomb[r][align].append(float(field[6]))
                    count_recomb[r][align].append(int(field[7]))
                    count_recomb[r][align].append(float(field[8].split('\n')[0]))

    print("============================ Substitution rate : " + m + " ============================")
    print("=======================================================================================")
    # Initialize data lists
    recomb_data = []
    frac_seq_data = []
    identity_data = []
    chain_data = []
    len_data = []
    true_recomb_data = []
    accuracy_data = []
    pearson_corr = []
    kendall_corr = []

    # Initialize x-tick labels
    xtick_labels = []
    xtick_labels_ = []

    # Populate data for each R and set x-tick labels
    for r in R:
        recomb_ = []
        frac_seq_ = []
        identity_ = []
        chain_data_ = []
        groud_truth_ = []
        accuracy_ = []
        mean_rec = 0.0
        mean_acc = 0.0
        count = 0
        for fa_id in fasta_ids:
            align = fa_id.split('.fa')[0]
            # print(r, align, count_recomb[r][align])
            if count_recomb[r][align][4] < 00.00: # 0 freqsec not part of an alignment
                continue
            chain_data_.append(count_recomb[r][align][1]) # max recomb
            recomb_.append(count_recomb[r][align][3])
            mean_rec += count_recomb[r][align][1] # max recomb
            count += 1
            frac_seq_.append(count_recomb[r][align][4])
            identity_.append(count_recomb[r][align][5])
            groud_truth_.append(count_recomb[r][align][6])
            accuracy_.append(count_recomb[r][align][7])
            mean_acc += count_recomb[r][align][7]
        chain_data.append(chain_data_)
        recomb_data.append(recomb_)
        frac_seq_data.append(frac_seq_)
        identity_data.append(identity_)
        true_recomb_data.append(groud_truth_)
        accuracy_data.append(accuracy_)

        len_ = str(len(recomb_))
        if r == '2000000000':
            xtick_labels.append("$\infty$")
            xtick_labels_.append("Recombination penalty = $\infty$")
        elif r == '0':
            xtick_labels.append("0")
            xtick_labels_.append("Recombination penalty =  $0$")
        elif r == '1000':
            xtick_labels.append("$10^3$")
            xtick_labels_.append("$\\rho = 10^3$")
        elif r == '10000':
            xtick_labels.append("$10^4$")
            xtick_labels_.append("Recombination penalty = $10^4$")
        elif r == '100000':
            xtick_labels.append("$10^5$")
            xtick_labels_.append("$\\rho = 10^5$")
        elif r == '200000':
            xtick_labels.append("$2\\times10^5$")
            xtick_labels_.append("$\\rho = 2\\times10^5$")
        elif r == '500000':
            xtick_labels.append("$5\\times10^5$")
            xtick_labels_.append("$\\rho = 5\\times10^5$")
        elif r == '1000000':
            xtick_labels.append("$10^6$")
            xtick_labels_.append("$\\rho = 10^6$")
        elif r == '1':
            xtick_labels.append("$1$")
            xtick_labels_.append("$\\rho = 1$")
        elif r == '10':
            xtick_labels.append("$10$")
            xtick_labels_.append("$\\rho = 10$")

        len_data.append(len_)

        # Calculate Pearson correlation b/w ground truth and chain data
        if r == '2000000000':
            pearson_corr.append(0)
            print("R = " + r + " Mean F1 : " + str(np.round(mean_acc/count, 5)) + " Pearson Coeff : " + str(0))
        else:
            pearson_corr.append(np.corrcoef(groud_truth_, chain_data_)[0][1])
            print("R = " + r + " Mean F1 : " + str(np.round(mean_acc/count, 5)) + " Pearson Coeff : " + str(np.round(np.corrcoef(groud_truth_, chain_data_)[0][1], 5)))

        # find kendall correlation b/w ground truth and chain data
        kendall_corr.append(stats.kendalltau(groud_truth_, chain_data_)[0])

    # use helvetica font for the plot
    plt.rc('font', family='sans-serif') 
    plt.rc('font', serif='Helvetica Neue') 
    plt.rc('text', usetex='false') 
    plt.rcParams.update({'font.size': 22})


    # increase font size
    sns.set(font_scale=1.2)

    # Create a custom color palette with translucent colors
    custom_palette = sns.color_palette("Set1", len(R))
    alpha = 0.9  # Adjust the alpha value for transparency

    # Create a bar plot for Kendall correlation vs. R using Seaborn
    fig8, ax8 = plt.subplots(figsize=(4, 3))
    # sns.barplot(x=R, y=kendall_corr, palette=custom_palette, alpha=alpha, zorder=3)
    sns.barplot(x=R, y=pearson_corr, palette=custom_palette, alpha=alpha, zorder=3)

    # add text N/A for R = 2000000000
    ax8.text(5, 0.05, 'N/A', ha='center', va='center', color='black', fontsize=14)

    # Set x-tick labels
    ax8.set_xticklabels(xtick_labels, rotation=45)
    ax8.set_xlabel("Recombination penalty")
    ax8.set_ylabel("Pearson correlation")
    # ax8.set_title("Kendall Correlation vs Recombination Penalty")

    # Customize the style with Seaborn
    sns.set(style="whitegrid")

    plt.tight_layout()
    # add x and y grid lines with whitegrid
    ax8.grid(color='white', linestyle='-', linewidth=1)

    # Save the plot using Seaborn style
    plt.savefig(folder + "/Pearson_Correlation_vs_R.pdf", bbox_inches='tight', dpi=1200, format='pdf')
    # plt.show()


    # increase font size
    sns.set(font_scale=1.2)

    # Create a single scatterplot matrix for all R values
    fig9, ax9 = plt.subplots(figsize=(4, 3))

    # Define line_kws to make the regression line thicker
    line_kws = {'linewidth': 3}
    # Iterate through R values and plot each dataset with a thicker regression line
    # reduce scatter point size and randomize the points a bit
    for i in range(len(R)):
        if R[i] == '2000000000' or R[i] == '0' or R[i] == '10000':
            # randomize the points a bit so that they are not on top of each other
            sns.regplot(x=true_recomb_data[i], y=chain_data[i], label=f'{xtick_labels_[i]}', line_kws=line_kws, scatter = False, ci=None)
            # use integer y ticks as 0 5 10 15
            ax9.set_yticks(np.arange(0, 16, 5))
        else:
            continue

    # Add a 45-degree line for reference
    x_line = np.linspace(min(true_recomb_data[0]), max(true_recomb_data[0]), 100)
    plt.plot(x_line, x_line, linestyle='--', color='blue', label='Ground truth')

    plt.xlabel("True recombinations")
    plt.ylabel("Chaining recombinations")
    # plt.title("Chaining Recombination vs True Recombination")

    # put legend outside the plot
    # ax9.legend(loc='upper left', borderaxespad=0., fontsize=14, ncols = 4, bbox_to_anchor=(0.0, 1.14, 1.0, 0.102), mode="expand")
    # remove borders from the legend
    plt.savefig(folder + "/line_fit.pdf", bbox_inches='tight', dpi=1200, format='pdf')


    # increase font size
    sns.set(font_scale=1.2)

    # Create box plot for accuracy vs. R
    fig10, ax10 = plt.subplots(figsize=(4, 3))
    # use sns to create boxplot
    # sns.boxplot(accuracy_data, showfliers=False, showmeans=True, meanline=True, meanprops={'marker':'o','markeredgecolor':'black','markerfacecolor':'firebrick'})
    ax10.boxplot(accuracy_data, showfliers=False, showmeans=False, meanline=False, medianprops={'color':'lime', 'linewidth':2.2})
    # add data points
    for i in range(len(accuracy_data)):
        y = accuracy_data[i]
        x = np.random.normal(i + 1, 0.04, size=len(y))
        # increase the size of the points
        ax10.plot(x, y, 'r.', alpha=0.3, markersize=10)
    ax10.set_xticks(np.arange(1, len(R) + 1))
    ax10.set_xticklabels(xtick_labels, rotation=45)
    ax10.set_xlabel("Recombination penalty")
    ax10.set_ylabel("F1-score")
    # ax10.set_title("F1 Score vs Recombination Penalty")
    plt.tight_layout()
    plt.savefig(folder + "/F1_Score_vs_R.pdf", bbox_inches='tight', dpi=1200, format='pdf')
    # plt.show()