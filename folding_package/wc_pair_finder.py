import os
import numpy
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.analysis.nuclinfo

def wc_pair_finder(number_of_bases, universe_name,  ):
    
    """Write the assigned parameter values into the params file near the corresponding parameter identicator"""
    k=0
    pairs = {}

    for i in range(1, number_of_bases+1):
        for j in range(1, number_of_bases+1):
            distance = mda.analysis.nuclinfo.wc_pair(universe_name, i, j, "RNAA", "RNAA")
            if distance < 3.2 and i != j:
                if universe_name.select_atoms('resid '+str(i)).resnames[0]=='URA' or universe_name.select_atoms('resid '+str(i)).resnames[0]=='CYT':
                    pair_1 = universe_name.select_atoms('resid '+str(i)+' and name N3').ids[0]
                if universe_name.select_atoms('resid '+str(i)).resnames[0]=='ADE' or universe_name.select_atoms('resid '+str(i)).resnames[0]=='GUA':
                    pair_1 = universe_name.select_atoms('resid '+str(i)+' and name N1').ids[0]
                
                if universe_name.select_atoms('resid '+str(j)).resnames[0]=='URA' or universe_name.select_atoms('resid '+str(j)).resnames[0]=='CYT':
                    pair_2 = universe_name.select_atoms('resid '+str(j)+' and name N3').ids[0]
                if universe_name.select_atoms('resid '+str(j)).resnames[0]=='ADE' or universe_name.select_atoms('resid '+str(j)).resnames[0]=='GUA':
                    pair_2 = universe_name.select_atoms('resid '+str(j)+' and name N1').ids[0]    
                
                pairs[k]=(pair_1, pair_2)
                k=k+1
            
    for i in range(0, len(pairs)):
        for j in range(0, len(pairs)):
            if i in pairs and j in pairs:
                if pairs[i][0] == pairs[j][1]:
                    if pairs[i][1] == pairs[j][0]:
                        pairs.pop(j)
    
    return pairs