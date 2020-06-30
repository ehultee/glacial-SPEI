#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Quantify the effect of variable stomatal conductance

Created on Mon Jun 29 13:29:01 2020

@author: EHU
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gSPEI as gSPEI

fpath_default = './data/SPEI_Files/'
fpath_conduct = './data/SPEI_Files/variable_stom_conduct'

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


## Read default and variable-conductance SPEI in to dicts
SPEI_by_model = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath_default+'NRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath_default+'WRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']
SPEI_by_model_C = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath_conduct+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath_conduct+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model_C[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model_C[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model_C[m]['diff'] = SPEI_by_model_C[m]['WRunoff'] - SPEI_by_model_C[m]['NRunoff']

## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = gSPEI.sort_models_to_basins(SPEI_by_model)
SPEI_by_basin_C = gSPEI.sort_models_to_basins(SPEI_by_model_C)

## Compute pairwise difference in glacial effect due to stomatal conductance
vcd = {b: [] for b in basin_names}
vcd_mean = []
vcd_perdiff = [] #compute percent change in glacial effect
for b in basin_names:
    df = SPEI_by_basin[b]['diff']
    df1 = SPEI_by_basin_C[b]['diff']
    df_diff = pd.DataFrame.subtract(df1, df).mean()
    perdiff = pd.Series.divide(df_diff, df.mean()) #per-model percent difference
    vcd[b].append(v for v in df_diff.values) #per-model mean
    vcd_mean.append(df_diff.mean()) #mean difference across models for each basin
    vcd_perdiff.append(perdiff.mean())
  
## Compare with pairwise difference in NRunoff - does stomatal conductance matter there?
vcd_NR = {b: [] for b in basin_names}
vcd_mean_NR = []
vcd_perdiff_NR = [] #compute percent change in glacial effect
for b in basin_names:
    df = SPEI_by_basin[b]['NRunoff']
    df1 = SPEI_by_basin_C[b]['NRunoff']
    df_diff = pd.DataFrame.subtract(df1, df).mean()
    perdiff = pd.Series.divide(df_diff, df.mean()) #per-model percent difference
    vcd_NR[b].append(v for v in df_diff.values) #per-model mean
    vcd_mean_NR.append(df_diff.mean()) #mean difference across models for each basin
    vcd_perdiff_NR.append(perdiff.mean())

fig, ax = plt.subplots(1)
ax.hist(vcd_perdiff, alpha=0.5)
# ax.hist(vcd_perdiff_NR, alpha=0.5)
ax.set(xlabel='Normalized difference due to stomatal conductance', ylabel='Count', 
       title='SPEI glacial effect')
plt.show()