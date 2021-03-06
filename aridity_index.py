#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Load in model precip & PET, compute aridity index, compare with drought stats
Created on Thu Dec 17 17:08:55 2020

@author: lizz
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fpath_p_pet = './data/model_P_PET/'
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios

## Settings in filenames
integration_times = np.arange(3, 28, 4) # all SPEI integration times used
modelnames = ['CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GISS-E2-R', 'INMCM4', 'MIROC-ESM', 'NorESM1-M'] # all models used in comparison
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios
## Basins in the order they are written
basin_names = ['INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE-VINCES',
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE']

yrs = np.linspace(1900, 2101, num=2412)
historical_period = (1900, 1980) # limit analysis to this interval


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
historical_avg_AI = [] # list for average basin aridity index over historical period
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
        historical_avg_AI.append(np.median([np.median(AI_by_basin[b]['Rcp4p5'][m][0:959]) for m in modelnames]))
# historical_avg_AI = [np.median([np.median(AI_by_basin[b]['Rcp4p5'][m][0:959]) for m in modelnames]) for b in basin_names]



fig, ax = plt.subplots()
# for m in modelnames:
#     ax.plot(yrs, PET_by_basin['TARIM']['Rcp4p5'][m], color='k')
ax.plot(yrs, AI_by_basin['COPPER']['Rcp4p5']['CCSM4'], color='k')
ax.set(xlabel='Year', ylabel='Aridity index P/PET', title='Aridity index over time for COPPER basin')
plt.show()


fig1, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4))
for b, aridity in zip(basin_names, historical_avg_AI):
    ax1.scatter(aridity, basin_stats[b][0][1]-basin_stats[b][0][0], color='k')
    ax2.scatter(aridity, basin_stats[b][1][1]-basin_stats[b][1][0], color='k')
    ax3.scatter(aridity, basin_stats[b][2][1]-basin_stats[b][2][0], color='k')
ax1.set(xlabel='Historical P/PET', ylabel='Mean number of droughts 1980-2100')
ax2.set(xlabel='Historical P/PET', ylabel='Mean drought duration 1980-2100')
ax3.set(xlabel='Historical P/PET', ylabel='Mean drought severity 1980-2100')
plt.tight_layout()
plt.show()
    