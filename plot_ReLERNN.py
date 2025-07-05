#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools



### Plot the ReLERNN inferred recombination map for each population.

populations = ["Al", "EGr", "Ice", "Ire_Scot", 
               "Newfound", "Norw_lago", "Py", "Scot", "Sva", "Swe_lago", "Swe_muta", "WGr"]




for p in populations:

    chr=["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_7", "SUPER_8", "SUPER_9",
    "SUPER_10", "SUPER_11", "SUPER_12", "SUPER_13", "SUPER_14", "SUPER_15", "SUPER_16", "SUPER_17",
    "SUPER_18", "SUPER_19", "SUPER_20", "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24", "SUPER_25",
    "SUPER_26", "SUPER_27", "SUPER_28", "SUPER_29"]


    r_val = []
    window_start = []
    window_end = []
    #rmap_index = []
    genomic_pos = []
    chrom = []

    for c in chr:
        path = "/home/prm/Skrivbord/ReLERNN_plotting/results/"
  
        rmap = pd.read_csv(path+p+"_"+c+"_recombination_map_2Kb_windows.csv", sep = ",")
        r_val.append(rmap['r_mean_window'].to_list())
        window_start.append(rmap['window_start'].to_list())
        window_end.append(rmap['window_end'].to_list())
        #rmap_index.append(rmap['i'].to_list())
        genomic_pos.append(rmap['genomic_pos'].to_list())
        chrom.append(rmap['chr'].to_list())



    r_val_flat = list(itertools.chain(*r_val))
    window_start_flat = list(itertools.chain(*window_start))
    window_end_flat = list(itertools.chain(*window_end))
    #rmap_index_flat = list(itertools.chain(*rmap_index))
    genomic_pos_flat = list(itertools.chain(*genomic_pos))
    chrom_flat = list(itertools.chain(*chrom))

    df = pd.DataFrame({'chr':chrom_flat, 'window_start': window_start_flat, 'window_end': window_end_flat, 
                           'genomic_pos': genomic_pos_flat, 'r_val':r_val_flat})

    df['rmap_index'] = df.index
    
    plot = sns.relplot(data=df, x='rmap_index', y='r_val', aspect=3.7, 
                    hue='chr', palette = 'bright', legend=None, kind="line") 




    chrom_df=df.groupby('chr')['rmap_index'].median()
    plot.ax.set_xlabel('chr'); plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.set_xlim(left = 0)
    plot.ax.set_xticklabels(plot.ax.get_xticklabels(), fontsize=4)
    plot.ax.set(ylim=(0.2e-09, 5e-09))
    plot.ax.set_title(p)


    plt.show()


    plot.savefig(path+p+"_recombination_map.pdf")
    

