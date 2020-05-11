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
from matplotlib import cm
import gSPEI as gSPEI

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
    norunoff_f_m = fpath+'NRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath+'WRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']

## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = gSPEI.sort_models_to_basins(SPEI_by_model)

## Compute multi-GCM ensemble means and quartiles
r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'WRunoff').rolling(window=12*30).mean()
r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'NRunoff').rolling(window=12*30).mean()
rm = SPEI_by_basin['TARIM']['WRunoff'].rolling(window=12*30, axis=0).mean()
rm_q1 = rm.quantile(q=0.25, axis=1)
rm_q3 = rm.quantile(q=0.75, axis=1)
rm_n = SPEI_by_basin['TARIM']['NRunoff'].rolling(window=12*30, axis=0).mean()
rm_q1_n = rm_n.quantile(q=0.25, axis=1)
rm_q3_n = rm_n.quantile(q=0.75, axis=1)

## Make example figure
single_models_w = [SPEI_by_basin['TARIM']['WRunoff'][m].rolling(window=12*30).mean() for m in modelnames]
single_models_n = [SPEI_by_basin['TARIM']['NRunoff'][m].rolling(window=12*30).mean() for m in modelnames]


colors_w = cm.get_cmap('Blues')(np.linspace(0.2, 1, num=len(modelnames)))
colors_n = cm.get_cmap('Wistia')(np.linspace(0.2, 1, num=len(modelnames)))
fig, ax = plt.subplots()
ax.plot(yrs, r_w, 'k', linewidth=3.0)
ax.plot(yrs, rm_q1, 'k')
ax.plot(yrs, rm_q3, 'k')
ax.plot(yrs, r_n, 'k', linewidth=3.0, ls=':')
ax.plot(yrs, rm_q1_n, 'k', ls=':')
ax.plot(yrs, rm_q3_n, 'k', ls=':')
# for i in range(len(modelnames)):
#     ax.plot(yrs, single_models_w[i], color=colors_w[i])
#     ax.plot(yrs, single_models_n[i], color=colors_n[i])
ax.fill_between(yrs, rm_q1, rm_q3, color='DarkBlue', alpha=0.2)
ax.fill_between(yrs, rm_q1_n, rm_q3_n, color='DarkOrange', alpha=0.2)
ax.tick_params(axis='both', labelsize=12)
ax.set_xticks([1900,1950, 2000, 2050, 2100])
ax.set_xlabel('Years', fontsize=14)
ax.set_ylabel('Rolling mean SPEI', fontsize=14)


## Calculate changes due to glacial effect at end of century, using ensemble approach
bas_glac_meandiff, quantile_spread = gSPEI.ensemble_glacial_meandiff(SPEI_by_basin)
bas_glac_vardiff, var_spread = gSPEI.ensemble_glacial_vardiff(SPEI_by_basin)
## Calculate median of mean shift and full multi-GCM range, for comparison
# bas_glac_meanmed, mean_spread_full = gSPEI.glacial_meandiff(SPEI_by_model)
# bas_glac_varmed, var_spread_full = gSPEI.glacial_vardiff(SPEI_by_model)

fig1, ax1 = plt.subplots()
ax1.errorbar(x=bas_glac_meandiff, y=bas_glac_vardiff, xerr=quantile_spread, yerr=var_spread, ls='', marker='d', elinewidth=2.0, color='DarkBlue')
# ax1.errorbar(x=bas_glac_meanmed, y=bas_glac_varmed, xerr=mean_spread_full, yerr=var_spread_full, ls='', marker='d', elinewidth=1.0, color='Cyan')
ax1.set_xlabel('Difference in mean SPEI', fontsize=16)
ax1.set_ylabel('Difference in SPEI variance', fontsize=16)
ax1.set(ylim=(-1.5, 1.0), xlim=(-0.5, 5))
ax1.tick_params(axis='both', labelsize=12)
plt.show()