#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plot a single SPEI time series and fill droughts

Created on Thu Sep  9 12:06:49 2021

@author: lizz
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
        SPEI_by_basin[b][c] = SPEI_by_basin_raw[b][c].fillna(-3)

## Plot a single series
example_basin='TARIM'
example_model='CCSM4'
example_period = (1980,2010)
example_series = SPEI_by_basin[example_basin]['NRunoff'][example_model]
droughts = gSPEI.find_droughts(example_series, threshold=-1, 
                               period=example_period)
drought_mask = np.full(np.shape(yrs), fill_value=False)
for k in droughts.keys():
    for n in range(k-len(droughts[k]), k+1):
        drought_mask[n] = True

example_color=cm.get_cmap('tab20b')(10)
fig, ax = plt.subplots(figsize=(6,4), tight_layout=True)
ax.axhline(y=0, color='k', alpha=0.8)
ax.axhline(y=-1, ls='-.', color='k')
ax.plot(yrs, example_series, color=example_color)
ax.fill_between(yrs, y1=0, y2=example_series, where=drought_mask, 
                alpha=0.5, color=example_color)
ax.annotate('{} droughts found'.format(len(droughts)), xy=(2000,-3),
            fontweight='bold')
ax.set(xlim=example_period)
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel(r'SPEI$_N$', fontsize=14)
plt.show()


## Comparison series for talk
## Plot a single series

example_series_n = SPEI_by_basin[example_basin]['NRunoff'][example_model]
example_series_w = SPEI_by_basin[example_basin]['WRunoff'][example_model]

droughts = gSPEI.find_droughts(example_series, threshold=-1, 
                               period=example_period)
drought_mask = np.full(np.shape(yrs), fill_value=False)
for k in droughts.keys():
    for n in range(k-len(droughts[k]), k+1):
        drought_mask[n] = True

example_color_n=cm.get_cmap('tab20b')(10)
example_color_w=cm.get_cmap('tab20b')(0)

droughts = gSPEI.find_droughts(example_series_n, threshold=-1, 
                               period=example_period)
drought_mask = np.full(np.shape(yrs), fill_value=False)
for k in droughts.keys():
    for n in range(k-len(droughts[k]), k+1):
        drought_mask[n] = True
fig, ax = plt.subplots(figsize=(15,3), tight_layout=True)
ax.axhline(y=0, color='k', alpha=0.8)
ax.axhline(y=-1, ls='-.', color='k')
ax.plot(yrs, example_series_n, color=example_color_n)
ax.fill_between(yrs, y1=0, y2=example_series_n, where=drought_mask, 
                alpha=0.5, color=example_color_n)
ax.annotate('{} droughts found'.format(len(droughts)), xy=(2000,-3),
            fontweight='bold')
ax.set(xlim=example_period, yticks=(-2, 0, 2))
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel(r'SPEI$_N$', fontsize=14)
plt.show()

droughts = gSPEI.find_droughts(example_series_w, threshold=-1, 
                               period=example_period)
drought_mask = np.full(np.shape(yrs), fill_value=False)
for k in droughts.keys():
    for n in range(k-len(droughts[k]), k+1):
        drought_mask[n] = True
fig, ax = plt.subplots(figsize=(15,3), tight_layout=True)
ax.axhline(y=0, color='k', alpha=0.8)
ax.axhline(y=-1, ls='-.', color='k')
ax.plot(yrs, example_series_w, color=example_color_w)
ax.fill_between(yrs, y1=0, y2=example_series_w, where=drought_mask, 
                alpha=0.5, color=example_color_w)
ax.annotate('{} droughts found'.format(len(droughts)), xy=(2000,-3),
            fontweight='bold')
ax.set(xlim=example_period, yticks=(-2,0,2))
ax.tick_params(axis='both', labelsize=12)
ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel(r'SPEI$_W$', fontsize=14)
plt.show()
