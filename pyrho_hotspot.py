#!/usr/bin/python


import pandas as pd
import math
import numpy as np
import itertools


#Identify recombination hotspots in 1Kb windows with a 1Kb step size from recombination rate estimated by pyrho. A recombination hotspot is defined as 10X higher recombination rate compared to the average recombination rate in 40Kb flankin regions.


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
        r_mean_flank_left = []
        r_mean_flank_right = []
        window_start = []
        window_end = []
        genomic_pos_start = []
        genomic_pos_end = []
        hot_spot_10 = []
        flanking_less_40Kb = []
        chrom = []

        
        for i in range(chr_len):
            recomb = r_val_full_flat[start_win:stop_win]
            recomb_window = sum(recomb)/len(recomb)
            r_mean_window.append(recomb_window)
            window_start.append(start_win)
            window_end.append(stop_win)
            genomic_pos_start.append(start_win+start[0])
            genomic_pos_end.append(stop_win+start[0])
            chrom.append(c)
            
            
            if (start_win-40000 <= 0):
                recomb_flanking_left = 0
                recomb_flanking_right = (sum(r_val_full_flat[stop_win:stop_win+40000]))/len(r_val_full_flat[stop_win:stop_win+40000])
                flanking_less_40Kb.append("T")
                
                
            elif (stop_win+40000) >= len(r_val_full_flat):
                recomb_flanking_left = (sum(r_val_full_flat[start_win-40000:start_win])/len(r_val_full_flat[start_win-40000:start_win]))
                recomb_flanking_right = 0
                flanking_less_40Kb.append("T")

            else:
                recomb_flanking_left = (sum(r_val_full_flat[start_win-40000:start_win])/len(r_val_full_flat[start_win-40000:start_win]))
                recomb_flanking_right = (sum(r_val_full_flat[stop_win:stop_win+40000]))/len(r_val_full_flat[stop_win:stop_win+40000])
                flanking_less_40Kb.append("F")
                                     
            r_mean_flank_left.append(recomb_flanking_left)
            r_mean_flank_right.append(recomb_flanking_right)

                
            if (recomb_window >= 10*recomb_flanking_left) and (recomb_window >= 10*recomb_flanking_right):
                hot_spot_10.append("T")
            else:
                hot_spot_10.append("F")
                
                

                
            
            start_win += 1000
            stop_win += 1000
        

        df = pd.DataFrame({'chrom': chrom,
                    'genomic_pos_start':genomic_pos_start,
                    'genomic_pos_end':genomic_pos_end,
                    'window_start': window_start,
                    'window_end': window_end,
                    'r_mean_window': r_mean_window,
                    'r_mean_flank_left': r_mean_flank_left,
                    'r_mean_flank_right': r_mean_flank_right,
                    'hot_spot_10': hot_spot_10})



        df.to_csv(path+p+'_'+c+'_recombination_hotspots.csv', sep='\t', index=None)

