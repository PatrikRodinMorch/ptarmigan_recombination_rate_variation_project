#!/usr/bin/python

import pandas as pd
import math
import numpy as np
import itertools


#Identify recombination coldspots in 2Kb windows from recombination rate estimated by pyrho. A recombination coldspot is here defined as 1/10 of the average recombination rate for the chromosome.

populations = ["Sva", "Al", "EGr", "Ice", "Ire_Scot", "Newfound", "Norw_lago", 
               "Py", "Scot", "Swe_lago", "Swe_muta", "WGr"]

for p in populations:



    chr=["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_7", "SUPER_8", "SUPER_9",
        "SUPER_10", "SUPER_11", "SUPER_12", "SUPER_13", "SUPER_14", "SUPER_15", "SUPER_16", "SUPER_17",
        "SUPER_18", "SUPER_19", "SUPER_20", "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24", "SUPER_25",
        "SUPER_26", "SUPER_27", "SUPER_28", "SUPER_29"]

    
    for c in chr:
        path = "/home/prm/Skrivbord/pyrho_plotting/results/"+p+"/"

        rmap = pd.read_csv(path+p+"_"+c+".rmap", sep = "\t", names=['start', 'stop', 'r'])
        rmap = rmap.dropna()
        rmap['L'] = rmap['stop'] - rmap['start']
        rmap['L'] = rmap['L'].astype(int)


        r_val = rmap['r'].to_list()
        start = rmap['start'].to_list()
        length_val = rmap['L'].to_list()


        r_val_full = []
        start_full = []

        start_pos = start[0]
                
        for n in range(len(r_val)):
            r_val_full.append(np.repeat(r_val[n], length_val[n]))
            
            window = length_val[n]
                        
            for n in range(window):
                start_full.append(start_pos)
                start_pos += 1
                    
                                    
        r_val_full_flat = list(itertools.chain(*r_val_full))
 
        chr_len = math.ceil(len(r_val_full_flat) / 1000)
              
        
        start_win = 0
        stop_win = 2000

  
        r_mean_window = []
        window_start = []
        window_end = []
        genomic_pos_start = []
        genomic_pos_end = []
        cold_spot_10 = []
        chrom = []

        chrom_average = sum(r_val_full_flat)/len(r_val_full_flat)

        #print(sum(r_val_full_flat[start_win:stop_win]))
        #print(len(r_val_full_flat[start_win:stop_win]))
        #print(chr_len)

        
        for i in range(chr_len):      
            
            recomb = r_val_full_flat[start_win:stop_win]
            
            if len(recomb) > 0:            
                recomb_window = sum(recomb)/len(recomb)
                r_mean_window.append(recomb_window)
                window_start.append(start_win)
                window_end.append(stop_win)
                genomic_pos_start.append(start_win+start[0])
                genomic_pos_end.append(stop_win+start[0])
                chrom.append(c)
            
               
                if (recomb_window <= ((1/10)*chrom_average)):
                    cold_spot_10.append("T")
                else:
                    cold_spot_10.append("F")
            
            else:
                pass
            
            start_win += 2000
            stop_win += 2000
        

        df = pd.DataFrame({'chrom': chrom,
                    'genomic_pos_start':genomic_pos_start,
                    'genomic_pos_end':genomic_pos_end,
                    'window_start': window_start,
                    'window_end': window_end,
                    'r_mean_window': r_mean_window,
                    'cold_spot_10': cold_spot_10
                    })



        df.to_csv(path+p+'_'+c+'_recombination_coldspots.csv', sep='\t', index=None)


