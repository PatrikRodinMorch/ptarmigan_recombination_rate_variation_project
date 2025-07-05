#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools
import os
import math


populations = ["Al", "EGr", "Ice", "Ire_Scot", 
               "Newfound", "Norw_lago", "Py", "Scot", "Sva", "Swe_lago", "Swe_muta", "WGr"]


# Break up the windows of varying size chosen by pyrho's hyperparam module into recombinate rate/bp.

for p in populations:

    chr=["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_7", "SUPER_8", "SUPER_9",
    "SUPER_10", "SUPER_11", "SUPER_12", "SUPER_13", "SUPER_14", "SUPER_15", "SUPER_16", "SUPER_17",
    "SUPER_18", "SUPER_19", "SUPER_20", "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24", "SUPER_25",
    "SUPER_26", "SUPER_27", "SUPER_28", "SUPER_29"]




    for c in chr:
        
        r_val_full = []
        chrom = []
        start_full = []

        path = "/home/prm/Skrivbord/pyrho_plotting/results/"+p+"/"

        rmap = pd.read_csv(path+p+"_"+c+".rmap", sep = "\t", names=['start', 'stop', 'r'])
        rmap = rmap.dropna()
        rmap['L'] = rmap['stop'] - rmap['start']
        rmap['L'] = rmap['L'].astype(int)


        r_val = rmap['r'].to_list()
        start = rmap['start'].to_list()
        length_val = rmap['L'].to_list()
        #rmap['I'] = rmap.index


        r_val_full = []
        start_full = []

        start_pos = start[0]
                
        for n in range(len(r_val)):
            r_val_full.append(np.repeat(r_val[n], length_val[n]))
            
            window = length_val[n]
            #print(window)
            
            for n in range(window):
                start_full.append(start_pos)
                start_pos += 1
                    
                                    
        r_val_full_flat = list(itertools.chain(*r_val_full))
        
        start = 0
        stop = 1000000

        r_mean_window = []
        window_start = []
        window_end = []
        genomic_pos = []
        
        chr_len = math.ceil(len(r_val_full_flat) / 1000000)
        
        
        for i in range(chr_len):
            recomb = r_val_full_flat[start:stop]
            recomb_window = sum(recomb)/len(recomb)
            r_mean_window.append(recomb_window)
            window_start.append(start)
            window_end.append(stop)
            genomic_pos.append(start_full[start:stop][0])
                       
                        
            
            start += 1000000
            stop += 1000000


        df = pd.DataFrame({'genomic_pos': genomic_pos, 'window_start': window_start, 
                           'window_end': window_end, 'r_mean_window': r_mean_window})
        df['i'] = df.index
        df['chr'] = c
        df['chr'] = df['chr'].str.replace("SUPER_", "Chr")
        
        df.to_csv(path+p+"_"+c+"_recombination_map_1Mb_windows.csv", index=None)





# %%
