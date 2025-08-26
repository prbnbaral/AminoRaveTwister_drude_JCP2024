import os
import numpy
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from itertools import combinations

from folding_package.com_distance import *

# THIS SUBROUTINE HAS BEEN LAST UPDATED ON 15TH OF MAY

def traj_gen_com(number_of_bases, max_id, min_id, frame_step, universe_name):
    cols=com_distance(number_of_bases, max_id, min_id, universe_name).columns
    dist = pd.DataFrame(columns=cols)
    pbar = tqdm(total=int(len(universe_name.trajectory))/frame_step)
    for frame in universe_name.trajectory[0:int(len(universe_name.trajectory)):frame_step]:
        frame=universe_name.trajectory.frame
        dist = dist.append(com_distance(number_of_bases, max_id, min_id, universe_name))
        pbar.update()
    pbar.close()
    return dist


def traj_gen_com_n(number_of_bases, max_id, min_id, frame_step, universe_name, N):
    
    #res_pairs = []
    res_array = numpy.arange(1, number_of_bases+1)
    
    def subsequences(iterable, length):
        return (iterable[i: i + length] for i in range(len(iterable) - length + 1))
    
    
    #res_array = np.arange(1, 4+1)
    #comb = combinations(res_array, N)
    comb = subsequences(res_array, N)
    print("First step is completed")
    res_pairs = list(combinations(comb, 2))
    print("Second step is completed")

    #cols=com_distance_n(number_of_bases, max_id, min_id, universe_name, res_pairs).columns
    
    print("Setting column names")
    cols=numpy.array(range(len(res_pairs)), dtype=object)
    pbar = tqdm(total=int(len(res_pairs)))
    for i in range(0, len(res_pairs)):
        first = str(res_pairs[i][0]).replace(',', '').replace('[', '').replace(']', '').replace(' ', '')
        second = str(res_pairs[i][1]).replace(',', '').replace('[', '').replace(']', '').replace(' ', '')
        cols[i]='d'+first+'_'+second   
        pbar.update()
    
    pbar.close() 
    print("Column names are set")
    dist = pd.DataFrame(columns=cols)
    pbar = tqdm(total=int(len(universe_name.trajectory))/frame_step)
    for frame in universe_name.trajectory[0:int(len(universe_name.trajectory)):frame_step]:
        frame=universe_name.trajectory.frame
        dist = dist.append(com_distance_n(number_of_bases, max_id, min_id, universe_name, res_pairs))
        pbar.update()
    pbar.close()
    return dist

# OP_TRAJ_GEN_STACK IS ADDED ON 15TH OF MAY

def OP_traj_gen_stack(number_of_bases, max_id, min_id, frame_step, universe_name):
    pairs = {}
        
    cols=com_distance_stack(number_of_bases, max_id, min_id, universe_name).columns
    dist = pd.DataFrame(columns=cols)
    pbar = tqdm(total=int(len(universe_name.trajectory))/frame_step)
    for frame in universe_name.trajectory[0:int(len(universe_name.trajectory)):frame_step]:
        frame=universe_name.trajectory.frame
        dist = dist.append(com_distance_stack(number_of_bases, max_id, min_id, universe_name))
        pbar.update()
    pbar.close()
    return dist