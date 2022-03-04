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
from matplotlib.patches import Rectangle
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
SPEI_by_model = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    wrunoff_f_m = fpath+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']

## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin_raw = gSPEI.sort_models_to_basins(SPEI_by_model)
SPEI_by_basin = {b: {} for b in basin_names}
for b in basin_names:
    for c in cases:
        # SPEI_by_basin[b][c] = SPEI_by_basin_raw[b][c].fillna(method='bfill')
        SPEI_by_basin[b][c] = SPEI_by_basin_raw[b][c].fillna(-3)


## Four-basin example for manuscript
# color_fam = cm.get_cmap('tab20b') 
# color_with = color_fam(0)
# color_no = color_fam(10)
color_with='darkblue' ## going with slightly brighter colours on recc of R2
color_no='gold'
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, 
                                              tight_layout=True, figsize=(9,6))
example_basins = ('COPPER', 'TARIM', 'RHONE', 'MAJES')
example_regions = ('NA','AS','EU','SA')
example_pgs = ('20%', '2.3%', '0.9%', '0.3%')
for example_b, example_r, example_pg, ax in zip(example_basins, example_regions, example_pgs, (ax1,ax2,ax3,ax4)):
    r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'WRunoff').rolling(window=12*30).mean()
    r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'NRunoff').rolling(window=12*30).mean()
    rm = SPEI_by_basin[example_b]['WRunoff'].rolling(window=12*30, axis=0).mean()
    rm_q1 = rm.quantile(q=0.25, axis=1, interpolation='lower')
    rm_q3 = rm.quantile(q=0.75, axis=1, interpolation='higher')
    rm_n = SPEI_by_basin[example_b]['NRunoff'].rolling(window=12*30, axis=0).mean()
    rm_q1_n = rm_n.quantile(q=0.25, axis=1, interpolation='lower')
    rm_q3_n = rm_n.quantile(q=0.75, axis=1, interpolation='higher')
    
    ax.plot(yrs, r_w, 'k', linewidth=3.0)
    ax.plot(yrs, rm_q1, 'k')
    ax.plot(yrs, rm_q3, 'k')
    ax.plot(yrs, r_n, 'k', linewidth=3.0, ls=':')
    ax.plot(yrs, rm_q1_n, 'k', ls=':')
    ax.plot(yrs, rm_q3_n, 'k', ls=':')
    # for i in range(len(modelnames)):
    #     ax.plot(yrs, single_models_w[i], color=colors_w[i], alpha=0.5)
    #     ax.plot(yrs, single_models_n[i], color=colors_n[i], alpha=0.5)
    ax.fill_between(yrs, rm_q1, rm_q3, color=color_with, alpha=0.5)
    ax.fill_between(yrs, rm_q1_n, rm_q3_n, color=color_no, alpha=0.5)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xticks([2000, 2050, 2100])
    ax.set_xlim(1980, 2100)
    extra1 = Rectangle((0,0), 0.1, 0.1, fc='w', fill=False, 
                      edgecolor='none', linewidth=0)
    extra2 = Rectangle((0,0), 0.1, 0.1, fc='w', fill=False, 
                  edgecolor='none', linewidth=0)
    leg = ax.legend([extra1, extra2], 
                    ['{}'.format(example_b), '({}, {})'.format(example_r, example_pg)], 
                    loc='best', 
                    handlelength=0, handletextpad=0, fancybox=True,
                    prop={'size':13})
for ax in (ax3, ax4):
    ax.set_xlabel('Year', fontsize=14)
for ax in (ax1, ax3):
    ax.set_ylabel('Rolling mean SPEI', fontsize=14)
fig1.align_ylabels()