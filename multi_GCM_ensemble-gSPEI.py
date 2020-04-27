#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Multi-model ensemble SPEI analysis using pandas
Created on Wed Apr 15 15:55:58 2020

@author: EHU
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gSPEI import *

fpath = './data/SPEI_Files/'

## Settings in filenames
integration_times = np.arange(3, 28, 4) # all SPEI integration times used
modelnames = ['CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GISS-E2-R', 'INMCM4', 'MIROC-ESM', 'NorESM1-M'] # all models used in comparison
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios
cases = ['NRunoff', 'WRunoff', 'diff'] # inclusion of glacier runoff

## Basins in the order they are written
basin_names = ['INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE-VINCES',
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE']

yrs = np.linspace(1900, 2101, num=2412)

## Read all in to dict by GCM as in other gSPEI scripts
SPEI_by_model = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath+'NRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[1])
    wrunoff_f_m = fpath+'WRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[1])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']

## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = {b: {} for b in basin_names} # create dictionary indexed by basin name
for i, b in enumerate(basin_names):
    SPEI_by_basin[b] = {case: {} for case in cases}
    for case in cases:
        tempdict = {}
        for m in modelnames:
            tempdict[m] = SPEI_by_model[m][case][i] # pull data from SPEI_by_model into this new dict
        SPEI_by_basin[b][case] = pd.DataFrame.from_dict(tempdict)

## Compute multi-GCM ensemble means for each basin and case
def basin_ensemble_mean(dict_by_basin, basin_id, case):
    basin_df = dict_by_basin[basin_id][case]
    em = basin_df.mean(axis=1) #compute mean among all models at each timestep
    return em # a pandas Series object

def basin_firstquartile(dict_by_basin, basin_id, case):
    basin_df = dict_by_basin[basin_id][case]
    q1 = basin_df.quantile(q=0.25, axis=1)
    return q1

def basin_thirdquartile(dict_by_basin, basin_id, case):
    basin_df = dict_by_basin[basin_id][case]
    q2 = basin_df.quantile(q=0.75, axis=1)
    return q2


## Make example figure--eventually wrap this in a nice helper function in gSPEI
r = basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'WRunoff').rolling(window=12*30).mean()
q1 = basin_firstquartile(SPEI_by_basin, 'TARIM', 'WRunoff').rolling(window=12*30).mean()
q2 = basin_thirdquartile(SPEI_by_basin, 'TARIM', 'WRunoff').rolling(window=12*30).mean()

fig, ax = plt.subplots()
ax.plot(yrs, r, 'k', linewidth=3.0)
ax.plot(yrs, q1, 'b')
ax.plot(yrs, q2, 'b')
ax.fill_between(yrs, q1, q2, color='b', alpha=0.2)
ax.tick_params(axis='both', labelsize=12)
ax.set_xticks([1900,1950, 2000, 2050, 2100])
ax.set_xlabel('Years', fontsize=14)
ax.set_ylabel('Rolling mean SPEI', fontsize=14)
