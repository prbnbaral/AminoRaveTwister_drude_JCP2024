import os
import numpy
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from folding_package.com_distance import  com_distance



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
