from collections import defaultdict
import os
import numpy
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm

from folding_package.com import *

# THIS SUBROUTINE LAST UPDATED ON 10TH OF MAY
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


def com_distance_n(number_of_bases, max_id, min_id, universe_name, res_pairs):
    n=0
    k=0  
    
    # Set column names
    #print("Setting column names")
    columns=numpy.array(range(len(res_pairs)), dtype=object)
    #pbar = tqdm(total=int(len(res_pairs)))
    for i in range(0, len(res_pairs)):
        first = str(res_pairs[i][0]).replace(',', '').replace('[', '').replace(']', '').replace(' ', '')
        second = str(res_pairs[i][1]).replace(',', '').replace('[', '').replace(']', '').replace(' ', '')
        columns[i]='d'+first+'_'+second   
        #pbar.update()
    
    #pbar.close() 
    #print("Column names are set")  
    dist_COMs = pd.DataFrame(index=['0'],columns=columns)
    print("Distance calculation started")
    pbar = tqdm(total=int(len(res_pairs)))
    for i in range(0, len(res_pairs)):
        dist_COMs.iloc[0][k]=(((com_n(max_id, min_id, str(res_pairs[i][0]).replace(',', '').replace('[', '').replace(']', ''), universe_name)-com_n(max_id, min_id, str(res_pairs[i][1]).replace(',', '').replace('[', '').replace(']', ''), universe_name))**2).sum())**0.5
        k=k+1
        pbar.update()
    dist_COMs =dist_COMs.dropna(axis='columns')   
    pbar.close() 
    print("Distance calculation ended")
    return dist_COMs


## COM_DISTANCE_STACK HAS BEEN ADDED ON 10TH OF MAY

def com_distance_stack(number_of_bases, max_id, min_id, universe_name):
    n=0
    k=0  
    
    # Set column names
    
    low=1
    high=2
    columns=[]
    while high<=number_of_bases:
        columns.append('d'+str(low)+'_'+str(high))
        
        low+=1
        high+=1        
    dist_COMs = pd.DataFrame(index=['0'],columns=columns)
    
    low=1
    high=2
    while low<number_of_bases:
        vec_1=com(max_id, min_id, low, universe_name)
        vec_2=com(max_id, min_id, high, universe_name)
        distance=((((vec_1-vec_2)**2).sum())**0.5)
        mag_v1=(vec_1**2).sum()**0.5
        mag_v2=(vec_2**2).sum()**0.5
        # angle between two nucleotides
        #angle=math.degrees(np.arccos(np.dot(vec_1, vec_2)/(mag_v1*mag_v2)))
        
        
        dist_COMs.iloc[0][low-1]=distance
        
        low+=1
        high+=1 
        
    dist_COMs =dist_COMs.dropna(axis='columns')    
    return dist_COMs