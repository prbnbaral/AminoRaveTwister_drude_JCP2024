import os
import numpy
import pandas as pd
import MDAnalysis as mda
import MDAnalysis as mda
from folding_package.distance_cal import *

def traj_gen(sel_1, sel_2, frame_step, universe_name):
    columns=distance_cal(sel_1, sel_2).columns
    dist = pd.DataFrame(columns=columns)
    for frame in universe_name.trajectory[0:int(len(universe_name.trajectory)):frame_step]:
        frame=universe_name.trajectory.frame
        dist = dist.append(distance_cal(sel_1, sel_2))
    return dist

def traj_gen_wc(pairs, frame_step, universe_name):
    columns=distance_cal_wc(pairs, universe_name).columns
    dist = pd.DataFrame(columns=columns)
    for frame in universe_name.trajectory[0:int(len(universe_name.trajectory)):frame_step]:
        frame=universe_name.trajectory.frame
        dist = dist.append(distance_cal_wc(pairs, universe_name))
    return dist
