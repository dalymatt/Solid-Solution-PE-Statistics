# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:47:47 2020

@author: Matt
"""

def save_obj(Output=None, chem=None):
    import pickle
    name=chem
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(Output, f, pickle.HIGHEST_PROTOCOL)
    return

def load_obj(name):
    import pickle
    with open(name, 'rb') as f:
        return pickle.load(f)

def save_obj_ElasDisc(Output=None, chem=None, stepsx=None, stepsy=None):
    import pickle
    name=chem+'_ElasDisc_'+str(stepsx)+'by'+str(stepsy)
    with open('Data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(Output, f, pickle.HIGHEST_PROTOCOL)
    return

def save_obj_ElasDisc2(Output=None, chem=None, stepsx=None, stepsy=None):
    import pickle
    name=chem+'_ElasDisc2_'+str(stepsx)+'by'+str(stepsy)
    with open('Data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(Output, f, pickle.HIGHEST_PROTOCOL)
    return

def save_obj_ElasDisc90(Output=None, chem=None, stepsx=None, stepsy=None):
    import pickle
    name=chem+'_ElasDisc90_'+str(stepsx)+'by'+str(stepsy)
    with open('Data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(Output, f, pickle.HIGHEST_PROTOCOL)
    return
