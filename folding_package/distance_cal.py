import os
import numpy
import pandas as pd
import MDAnalysis as mda
from folding_package.col_ids import *

def distance_cal(sel_1, sel_2):
    k=0
    num_sel_1 = sel_1.positions.shape[0]
    num_sel_2 = sel_2.positions.shape[0]
    columns=col_ids(sel_1, sel_2, num_sel_1, num_sel_2)
    distances = pd.DataFrame(index=['0'],columns=columns)
    for i in range(num_sel_1):
        for j in range(num_sel_2):
            del_x=sel_1.positions[i][0]-sel_2.positions[j][0]
            del_y=sel_1.positions[i][1]-sel_2.positions[j][1]
            del_z=sel_1.positions[i][2]-sel_2.positions[j][2]
            distances.iloc[0][k]=(del_x**2+del_y**2+del_z**2)**0.5
            k=k+1
    return distances

def distance_cal_wc(pairs, universe_name):
    k=0
    columns=col_ids_wc(pairs)
    distances = pd.DataFrame(index=['0'],columns=columns)
    for i in range(len(pairs)):
        sel_1 = universe_name.select_atoms('index '+str(pairs[list(pairs)[i]][0]))
        sel_2 = universe_name.select_atoms('index '+str(pairs[list(pairs)[i]][1]))
        del_x=sel_1.positions[0][0]-sel_2.positions[0][0]
        del_y=sel_1.positions[0][1]-sel_2.positions[0][1]
        del_z=sel_1.positions[0][2]-sel_2.positions[0][2]
        distances.iloc[0][k]=(del_x**2+del_y**2+del_z**2)**0.5
        k=k+1
    return distances