#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Load in model precip & PET
Created on Thu Dec 17 17:08:55 2020

@author: lizz
"""
import numpy as np

fpath_p_pet = './data/model_P_PET/'
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios

PET_by_model = {m: {} for m in modelnames}
P_by_model = {m: {} for m in modelnames}
for m in modelnames:
    for s in scenarios:
        example_pet = fpath_p_pet+'{}_{}_PET.txt'.format(m,s)
        example_p = fpath_p_pet+'{}_{}_PREC.txt'.format(m,s)
        PET_by_model[m][s] = np.loadtxt(example_pet)
        P_by_model[m][s] = np.loadtxt(example_p)
    
    
PET_by_basin = {b: {} for b in basin_names} # potential evapotranspiration by basin
P_by_basin = {b: {} for b in basin_names}
AI_by_basin = {b: {} for b in basin_names} #aridity index by basin
for i, b in enumerate(basin_names):
    PET_by_basin[b] = {s: {} for s in scenarios}
    P_by_basin[b] = {s: {} for s in scenarios}
    AI_by_basin[b] = {s: {} for s in scenarios}
    for s in scenarios:
        tempdict_pet = {}
        tempdict_p = {}
        tempdict_ai = {}
        for m in modelnames:
            tempdict_pet[m] = PET_by_model[m][s][i] 
            tempdict_p[m] = P_by_model[m][s][i]
            tempdict_ai[m] = np.divide(tempdict_p[m], tempdict_pet[m])
        PET_by_basin[b][s] = pd.DataFrame.from_dict(tempdict_pet)
        P_by_basin[b][s] = pd.DataFrame.from_dict(tempdict_p)
        AI_by_basin[b][s] = pd.DataFrame.from_dict(tempdict_ai)


fig, ax = plt.subplots()
# for m in modelnames:
#     ax.plot(yrs, PET_by_basin['TARIM']['Rcp4p5'][m], color='k')
ax.plot(yrs, AI_by_basin['COPPER']['Rcp4p5']['CCSM4'], color='k')
ax.set(xlabel='Year', ylabel='Aridity index P/PET', title='Aridity index over time for COPPER basin')
plt.show()