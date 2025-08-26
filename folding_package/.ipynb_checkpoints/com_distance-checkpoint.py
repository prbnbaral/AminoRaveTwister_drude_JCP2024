from collections import defaultdict
import os
import numpy
import pandas as pd
import MDAnalysis as mda

from folding_package.com import com


# Generates the center of mass coordinates for given number of bases in the system



def com_distance(number_of_bases, max_id, min_id, universe_name):
    n=0
    k=0  
    
    # Set column names
    columns=numpy.array(range(number_of_bases*number_of_bases), dtype=object)
    for i in range(1, number_of_bases+1):
        for j in range(1, number_of_bases+1):
            if i != j:
                columns[n]='d'+str(i)+'_'+str(j)
                n=n+1
            
    dist_COMs = pd.DataFrame(index=['0'],columns=columns)
    
    for i in range(1, number_of_bases+1):
        for j in range(1, number_of_bases+1):
        #[nucleics[i].append(ids) for ids in rna.select_atoms('resid '+str(i)).ids if ids <= max_id]
            if i != j:
                dist_COMs.iloc[0][k]=(((com(max_id, min_id, i, universe_name)-com(max_id, min_id, j, universe_name))**2).sum())**0.5
                k=k+1
    dist_COMs =dist_COMs.dropna(axis='columns')    
    return dist_COMs