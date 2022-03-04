#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Multi-page composite figure of all basins
Created on Tue Jul 13 09:15:09 2021

@author: lizz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
import gSPEI as gSPEI

## Confirm using latest data - nonparametric standardization, 
## accounting for variable stomatal conductance
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

regions = ['AS', 'AS', 'AS', 'AS', 'NA',
           'AS', 'NA', 'NA', 'NA', 'AS',
           'NA', 'SA', 'NA', 'SA', 'AS',
           'AS', 'NA', 'AS', 'SA', 'SA',
           'NA', 'NA', 'NA', 'EU', 'EU',
           'NA', 'EU', 'NA', 'AS', 'EU',
           'AS', 'EU', 'NA', 'EU', 'AS',
           'EU', 'EU', 'AS', 'AS', 'EU',
           'SA', 'SA', 'NA', 'AS', 'SA',
           'NA', 'SA', 'AS', 'SA', 'SA',
           'NZ', 'SA', 'EU', 'SA', 'EU', 'NA'] ## region tag of each basin above

BasinArea=[1139075,1051731,518011,1233148,64959,1024462,829632,28422,49470,423657,51147,30599,
239678,30760,1745094,258475,668561,191032,5880854,390631,17967,1752001,21211,7527,7311,
118114,97485,42944,2701040,5678,787256,793704,1099380,73066,54103,190522,42862,988062,341227,
25127,15689,11882,7961,58935,107215,29513,24108,411516,130062,18612,17118,41993,
17157,261204,17364,57544] # area of each basin in km2

basin_glacier_area = [26893.8, 24645.4, 16606.7, 15176.7, 12998., 11216., 9535.4, 5614.8, 4304., 
                      3945.4, 3467.6, 3027.8, 2495.1, 2372.3, 2317.4, 2295.9, 1878.4, 1677.3,
                      1634.1, 1601.2, 1583.6, 1519.2, 1337.3, 1251.8, 1098.6, 1032.8, 904.2, 742.3,
                      739.5, 683.4, 485.7, 408.4, 374.7, 347.3, 312.7, 285.0, 269.4, 267.9, 248.4,
                      247.2, 238.1, 198.9, 159.5, 146., 134.5, 86.4, 76.2, 71.2, 64.1, 57.3, 46.5, 
                      40.6, 37.9, 33.3, 32.1, 31.9]

PG = np.array(basin_glacier_area)/np.array(BasinArea)

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
        SPEI_by_basin[b][c] = SPEI_by_basin_raw[b][c].fillna(-3)


## Multi-page plotting - 2 pages of 7x4
batch_size = 20
batched_basins = [basin_names[i:i+batch_size] for i in range(0, len(basin_names), batch_size)]
batched_regions = [regions[i:i+batch_size] for i in range(0, len(basin_names), batch_size)]
batched_pg = [PG[i:i+batch_size] for i in range(0, len(basin_names), batch_size)]
# color_fam = cm.get_cmap('tab20b')
color_with='darkblue' ## going with slightly brighter colours on recc of R2
color_no='gold'

for k in range(len(batched_basins)): ## looping over pages
    fig, axs = plt.subplots(nrows=5, ncols=4, sharex=True, sharey=True,
                           figsize=(8,10), tight_layout=True)
    batch = batched_basins[k]
    batch_r = batched_regions[k]
    batch_pg = batched_pg[k]
    
    for i in range(len(batch)):
        example_b = batch[i]
        example_r = batch_r[i]
        example_pg = batch_pg[i] 
        r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'WRunoff').rolling(window=12*30).mean()
        r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, example_b, 'NRunoff').rolling(window=12*30).mean()
        rm = SPEI_by_basin[example_b]['WRunoff'].rolling(window=12*30, axis=0).mean()
        rm_q1 = rm.quantile(q=0.25, axis=1, interpolation='lower')
        rm_q3 = rm.quantile(q=0.75, axis=1, interpolation='higher')
        rm_n = SPEI_by_basin[example_b]['NRunoff'].rolling(window=12*30, axis=0).mean()
        rm_q1_n = rm_n.quantile(q=0.25, axis=1, interpolation='lower')
        rm_q3_n = rm_n.quantile(q=0.75, axis=1, interpolation='higher')
        
        ax = axs.ravel()[i]
        ax.axhline(0, ls=':', lw=0.5, color='k', alpha=0.5)
        ax.plot(yrs, r_w, 'k', linewidth=3.0)
        ax.plot(yrs, rm_q1, 'k')
        ax.plot(yrs, rm_q3, 'k')
        ax.plot(yrs, r_n, 'k', linewidth=3.0, ls=':')
        ax.plot(yrs, rm_q1_n, 'k', ls=':')
        ax.plot(yrs, rm_q3_n, 'k', ls=':')
        ax.fill_between(yrs, rm_q1, rm_q3, color=color_with, alpha=0.4)
        ax.fill_between(yrs, rm_q1_n, rm_q3_n, color=color_no, alpha=0.4)
        ax.tick_params(axis='both', labelsize=12)
        ax.set(xlim=(1980,2100), ylim=(-1, 1.5), xticks=[2000,2050,2100], 
               yticks=(-1, 0, 1))
        # ax.annotate('{:.1%}'.format(example_pg) , (1990, -0.75))
        # ax.text(0.2, 0.1, str(example_b), transform=ax.transAxes,
        #         ha='left', size=12, weight=500, color='k')
        extra1 = Rectangle((0,0), 0.1, 0.1, fc='w', fill=False, 
                          edgecolor='none', linewidth=0)
        extra2 = Rectangle((0,0), 0.1, 0.1, fc='w', fill=False, 
                  edgecolor='none', linewidth=0)
        leg = ax.legend([extra1, extra2], 
                        ['{}'.format(example_b), '({}, {:.1%})'.format(example_r, example_pg)], 
                        loc='best', handlelength=0, handletextpad=0, fancybox=True)
        if i%4==0: ## label leftmost axes
            ax.set_ylabel('Roll. mean SPEI')
        if i>=(min(len(axs.ravel()),len(batch))-4): ## label bottom axes
            ax.set(xlabel='Year', xticks=[2000,2050,2100], xticklabels=[2000,2050, 2100])
            ax.tick_params(axis='x', reset=True)
    if len(axs.ravel())>len(batch):
        for j in range(len(batch), len(axs.ravel())):
            # axs.ravel()[j].axis('off')
            fig.delaxes(axs.ravel()[j])
    fig.align_ylabels()
    fig.tight_layout()
    fig.set_size_inches(8.00, 9.58)
    fig.savefig('/Users/lizz/Desktop/20220125-batched_recolored_basins-p{}'.format(k+1))
    plt.show()

