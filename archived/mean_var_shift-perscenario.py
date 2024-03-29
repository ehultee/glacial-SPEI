#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plot glacial effect on SPEI mean and variance in RCP 4.5 versus 8.5
Created on Tue May 18 16:01:26 2021

@author: lizz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import gSPEI as gSPEI

fpath = './data/SPEI_Files/nonparametric-var_stom_c/'

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
SPEI_by_model_4p5 = {m: {} for m in modelnames} # create dictionary indexed by model name
SPEI_by_model_8p5 = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_4 = fpath+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    wrunoff_f_4 = fpath+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    norunoff_f_8 = fpath+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[1])
    wrunoff_f_8 = fpath+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[1])
    SPEI_by_model_4p5[m]['NRunoff'] = np.loadtxt(norunoff_f_4)
    SPEI_by_model_4p5[m]['WRunoff'] = np.loadtxt(wrunoff_f_4)
    SPEI_by_model_8p5[m]['NRunoff'] = np.loadtxt(norunoff_f_8)
    SPEI_by_model_8p5[m]['WRunoff'] = np.loadtxt(wrunoff_f_8)
    SPEI_by_model_4p5[m]['diff'] = SPEI_by_model_4p5[m]['WRunoff'] - SPEI_by_model_4p5[m]['NRunoff']
    SPEI_by_model_8p5[m]['diff'] = SPEI_by_model_8p5[m]['WRunoff'] - SPEI_by_model_8p5[m]['NRunoff']

## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin_4p5_raw = gSPEI.sort_models_to_basins(SPEI_by_model_4p5)
SPEI_by_basin_8p5_raw = gSPEI.sort_models_to_basins(SPEI_by_model_8p5)
SPEI_by_basin_4p5 = {b: {} for b in basin_names}
SPEI_by_basin_8p5 = {b: {} for b in basin_names}

for b in basin_names:
    for c in cases:
        # SPEI_by_basin_4p5[b][c] = SPEI_by_basin_4p5_raw[b][c].fillna(method='ffill')
        # SPEI_by_basin_8p5[b][c] = SPEI_by_basin_8p5_raw[b][c].fillna(method='ffill')
        SPEI_by_basin_4p5[b][c] = SPEI_by_basin_4p5_raw[b][c].fillna(-3)
        SPEI_by_basin_8p5[b][c] = SPEI_by_basin_8p5_raw[b][c].fillna(-3)

## Calculate changes due to glacial effect at end of century, using ensemble approach
meandiff_4p5, quantile_spread_4p5 = gSPEI.ensemble_glacial_meandiff(SPEI_by_basin_4p5)
vardiff_4p5, var_spread_4p5 = gSPEI.ensemble_glacial_vardiff(SPEI_by_basin_4p5)
meandiff_8p5, quantile_spread_8p5 = gSPEI.ensemble_glacial_meandiff(SPEI_by_basin_8p5)
vardiff_8p5, var_spread_8p5 = gSPEI.ensemble_glacial_vardiff(SPEI_by_basin_8p5)

## Calculate full multi-GCM range, for comparison
_, mean_spread_full_4p5 = gSPEI.glacial_meandiff(SPEI_by_model_4p5)
_, var_spread_full_4p5 = gSPEI.glacial_vardiff(SPEI_by_model_4p5)
_, mean_spread_full_8p5 = gSPEI.glacial_meandiff(SPEI_by_model_8p5)
_, var_spread_full_8p5 = gSPEI.glacial_vardiff(SPEI_by_model_8p5)

fig1, ax1 = plt.subplots(figsize=(5,4))
ax1.axhline(y=0, ls=':', color='Grey', alpha=0.5)
ax1.axvline(x=0, ls=':', color='Grey', alpha=0.5)
ax1.errorbar(x=meandiff_8p5, y=vardiff_8p5, xerr=quantile_spread_8p5, yerr=var_spread_8p5, 
             ls='', marker='d', elinewidth=2.0, color='LightGray')
ax1.errorbar(x=meandiff_8p5, y=vardiff_8p5, xerr=mean_spread_full_8p5, yerr=var_spread_full_8p5, 
             ls='', marker='d', elinewidth=1.0, color='LightGray', alpha=0.5) #extend whiskers to full range
ax1.errorbar(x=meandiff_4p5, y=vardiff_4p5, xerr=quantile_spread_4p5, yerr=var_spread_4p5, 
             ls='', marker='d', elinewidth=2.0, color='k')
ax1.errorbar(x=meandiff_4p5, y=vardiff_4p5, xerr=mean_spread_full_4p5, yerr=var_spread_full_4p5, 
             ls='', marker='d', elinewidth=1.0, color='k', alpha=0.5) #extend whiskers to full range
ax1.set_xlabel('Difference in mean SPEI', fontsize=16)
ax1.set_ylabel('Difference in SPEI variance', fontsize=16)
ax1.set(ylim=(-1, 2), yticks=(-1, 0, 1, 2),
        xlim=(-0.1, 2), xticks=(0,1,2))
ax1.tick_params(axis='both', labelsize=12)
plt.tight_layout()
plt.show()


# ## Four-basin example for manuscript
# fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, 
#                                               tight_layout=True, figsize=(9,6))
# example_basins = ('COPPER', 'TARIM', 'RHONE', 'MAJES')
# for example_b, ax in zip(example_basins, (ax1,ax2,ax3,ax4)):
#     r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'WRunoff').rolling(window=12*30).mean()
#     r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'NRunoff').rolling(window=12*30).mean()
#     rm = SPEI_by_basin[example_b]['WRunoff'].rolling(window=12*30, axis=0).mean()
#     rm_q1 = rm.quantile(q=0.25, axis=1)
#     rm_q3 = rm.quantile(q=0.75, axis=1)
#     rm_n = SPEI_by_basin[example_b]['NRunoff'].rolling(window=12*30, axis=0).mean()
#     rm_q1_n = rm_n.quantile(q=0.25, axis=1)
#     rm_q3_n = rm_n.quantile(q=0.75, axis=1)
    
#     ax.plot(yrs, r_w, 'k', linewidth=3.0)
#     ax.plot(yrs, rm_q1, 'k')
#     ax.plot(yrs, rm_q3, 'k')
#     ax.plot(yrs, r_n, 'k', linewidth=3.0, ls=':')
#     ax.plot(yrs, rm_q1_n, 'k', ls=':')
#     ax.plot(yrs, rm_q3_n, 'k', ls=':')
#     # for i in range(len(modelnames)):
#     #     ax.plot(yrs, single_models_w[i], color=colors_w[i], alpha=0.5)
#     #     ax.plot(yrs, single_models_n[i], color=colors_n[i], alpha=0.5)
#     ax.fill_between(yrs, rm_q1, rm_q3, color='DarkBlue', alpha=0.2)
#     ax.fill_between(yrs, rm_q1_n, rm_q3_n, color='DarkOrange', alpha=0.2)
#     ax.tick_params(axis='both', labelsize=12)
#     ax.set_xticks([2000, 2050, 2100])
#     ax.set_xlim(1980, 2100)
# for ax in (ax3, ax4):
#     ax.set_xlabel('Year', fontsize=14)
# for ax in (ax1, ax3):
#     ax.set_ylabel('Rolling mean SPEI', fontsize=14)
# fig1.align_ylabels()