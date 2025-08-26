import os
import numpy
import pandas as pd
import MDAnalysis as mda

def res_ids(resid, min_atom_id, max_atom_id, universe_name):
    ids = universe_name.select_atoms('resid '+
                     str(resid)+
                     ' and index '+
                     str(min_atom_id)+
                     ':'+
                     str(max_atom_id)).ids
    return ids


def res_ids_n(resid, min_atom_id, max_atom_id, universe_name):
    ids = universe_name.select_atoms('resid '+
                     resid+
                     ' and index '+
                     str(min_atom_id)+
                     ':'+
                     str(max_atom_id)).ids
    return ids
