import os
import numpy
import pandas as pd
import MDAnalysis as mda

def col_ids(sel_1, sel_2, num_sel_1, num_sel_2):
    n=0
    columns=numpy.array(range(num_sel_1*num_sel_2), dtype=object)
    for i in range(num_sel_1):
        for j in range(num_sel_2):
            columns[n]='d'+str(sel_1.ids[i])+'-'+str(sel_2.ids[j])
            n=n+1
    return columns

def col_ids_wc(pairs):
    columns=numpy.array(range(len(pairs)), dtype=object)
    for i in range(len(pairs)):
        columns[i]='d'+str(pairs[list(pairs)[i]][0])+'-'+str(pairs[list(pairs)[i]][1])
    return columns