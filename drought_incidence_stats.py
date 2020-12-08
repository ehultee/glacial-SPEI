#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Compute drought frequency, duration, severity from SPEI 

Created on Wed Nov 25 11:19:09 2020

@author: lizz
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gSPEI as gSPEI
import collections

fpath_conduct = './data/SPEI_Files/variable_stom_conduct/'

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
SPEI_by_model_C = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath_conduct+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath_conduct+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model_C[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model_C[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model_C[m]['diff'] = SPEI_by_model_C[m]['WRunoff'] - SPEI_by_model_C[m]['NRunoff']
    
## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = gSPEI.sort_models_to_basins(SPEI_by_model_C)

## Analyse multi-model ensemble mean & quantiles for drought statistics
r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'WRunoff')
r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'NRunoff')

def find_droughts(series, threshold=0):
    """Identify droughts in a timeseries of SPEI

    Parameters
    ----------
    series : pandas Series
        SPEI time series.
    threshold : float, optional
        Cutoff value, below which we define a drought. The default is 0.

    Returns
    -------
    droughts : OrderedDict
        SPEI values sorted into droughts, with dict keys the indices in yrs of drought onset.

    """
    droughts = collections.OrderedDict()
    for i,v in enumerate(series):
        if (v>=threshold or np.isnan(v)): 
            pass
        else:
            if i==0:
                current_drought = [] #create new for start of time
            elif (series[i-1]>=threshold or np.isnan(series[i-1])):
                current_drought = [] #create new for start of drought
            current_drought.append(v)
            if i==len(series)-1:
                droughts[i] = current_drought #write it out if time ending
            elif series[i+1]>=threshold:
                droughts[i] = current_drought #write it out if drought ending
    return droughts
    

keys_w_bymodel = {m: [] for m in modelnames}
keys_n_bymodel = {m: [] for m in modelnames}
drt_dur_w_bymodel = {m: [] for m in modelnames}
drt_dur_n_bymodel = {m: [] for m in modelnames}
drt_sev_w_bymodel = {m: [] for m in modelnames}
drt_sev_n_bymodel = {m: [] for m in modelnames}

for m in modelnames:
    print(m)
    ser_w = SPEI_by_basin['TARIM']['WRunoff'][m]
    ser_n = SPEI_by_basin['TARIM']['NRunoff'][m]
    droughts_w = find_droughts(ser_w, threshold=0)
    droughts_n = find_droughts(ser_n, threshold=0)
    drts_w_trimmed = collections.OrderedDict({k: droughts_w[k] for k in droughts_w if k>=960}) #960 is first index where yrs > 1980 (glacier model on)
    drts_n_trimmed = collections.OrderedDict({k: droughts_n[k] for k in droughts_n if k>=960})
    
    keys_w_bymodel[m] = drts_w_trimmed.keys()
    keys_n_bymodel[m] = drts_n_trimmed.keys()
    drt_dur_w_bymodel[m] = [len(droughts_w[k]) for k in drts_w_trimmed.keys()]
    drt_dur_n_bymodel[m] = [len(droughts_n[k]) for k in drts_n_trimmed.keys()]
    drt_sev_w_bymodel[m] = [sum(droughts_w[k]) for k in drts_w_trimmed.keys()]
    drt_sev_n_bymodel[m] = [sum(droughts_n[k]) for k in drts_n_trimmed.keys()]

## find ensemble versions
droughts_w = find_droughts(r_w, threshold=0)
droughts_n = find_droughts(r_n, threshold=0)
drts_w_trimmed = collections.OrderedDict({k: droughts_w[k] for k in droughts_w if k>=960}) #960 is first index where yrs > 1980 (glacier model on)
drts_n_trimmed = collections.OrderedDict({k: droughts_n[k] for k in droughts_n if k>=960})

drt_dur_w = [len(droughts_w[k]) for k in drts_w_trimmed.keys()]
drt_dur_n = [len(droughts_n[k]) for k in drts_n_trimmed.keys()]
drt_sev_w = [sum(droughts_w[k]) for k in drts_w_trimmed.keys()]
drt_sev_n = [sum(droughts_n[k]) for k in drts_n_trimmed.keys()]

fig1, ax1 = plt.subplots()
for m in modelnames:
    ax1.scatter(yrs[keys_w_bymodel[m]], drt_dur_w_bymodel[m], color='DarkBlue', alpha=0.3)
    ax1.scatter(yrs[keys_n_bymodel[m]], drt_dur_n_bymodel[m], color='Orange', alpha=0.3)
ax1.scatter(yrs[drts_w_trimmed.keys()], drt_dur_w, color='DarkBlue', alpha=0.8, label='WRunoff')
ax1.scatter(yrs[drts_n_trimmed.keys()], drt_dur_n, color='Orange', alpha=0.8, label='NRunoff')
ax1.set(xlim=(1980, 2100), xlabel=('Year'), ylabel=('Drought duration [mo]'))
ax1.legend(loc='best')
plt.show()

fig2, ax2 = plt.subplots()
for m in modelnames:
    ax2.scatter(yrs[keys_w_bymodel[m]], drt_sev_w_bymodel[m], color='DarkBlue', alpha=0.3)
    ax2.scatter(yrs[keys_n_bymodel[m]], drt_sev_n_bymodel[m], color='Orange', alpha=0.3)
ax2.scatter(yrs[drts_w_trimmed.keys()], drt_sev_w, color='DarkBlue', alpha=0.8, label='WRunoff')
ax2.scatter(yrs[drts_n_trimmed.keys()], drt_sev_n, color='Orange', alpha=0.8, label='NRunoff')
ax2.set(xlim=(1980, 2100), xlabel=('Year'), ylabel=('Accumulated SPEI deficit'))
ax2.legend(loc='best')
plt.show()