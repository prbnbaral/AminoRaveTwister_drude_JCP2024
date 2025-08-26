import os
import numpy
import pandas as pd
import MDAnalysis as mda

def com(max_id, min_id, resid, universe_name):
    COM = universe_name.select_atoms('resid '+
                     str(resid)+
                     ' and index '+
                     str(min_id)+
                     ':'+
                     str(max_id)).center_of_mass()
    return COM


def com_n(max_id, min_id, resid, universe_name):
    COM = universe_name.select_atoms('resid '+
                     resid+
                     ' and index '+
                     str(min_id)+
                     ':'+
                     str(max_id)).center_of_mass()
    return COM