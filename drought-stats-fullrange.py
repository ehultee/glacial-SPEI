#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 21:26:17 2021
Plot range of drought stats as well as mean

@author: lizz
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import gSPEI as gSPEI

## Labels: (P)arametric or (NP)nonparametric;
## Standardization (1) lumped or (2) split by starting month
fpath_NP2 = './data/SPEI_Files/nonparametric-var_stom_c/'

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

yrs = np.linspace(1900, 2101, num=2412)
SPEI_by_model_C = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath_NP2+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    wrunoff_f_m = fpath_NP2+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, scenarios[0])
    SPEI_by_model_C[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model_C[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model_C[m]['diff'] = SPEI_by_model_C[m]['WRunoff'] - SPEI_by_model_C[m]['NRunoff']
    
## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = gSPEI.sort_models_to_basins(SPEI_by_model_C)
for b in basin_names:
    for t in ('NRunoff', 'WRunoff', 'diff'):
        SPEI_by_basin[b][t] = SPEI_by_basin[b][t].fillna(-3)

## Analyse multi-model ensemble mean & quantiles for drought statistics
r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'WRunoff')
r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'NRunoff')


basin_stats_bymodel_hist = {m: {b: gSPEI.basin_summary_stats(SPEI_by_basin, basin_name=b, modelnames=[m], period=(1980,2010)) for b in basin_names} 
                    for m in modelnames}
basin_stats_bymodel_midC = {m: {b: gSPEI.basin_summary_stats(SPEI_by_basin, basin_name=b, modelnames=[m], period=(2030,2060)) for b in basin_names} 
                    for m in modelnames}
basin_stats_bymodel_endC = {m: {b: gSPEI.basin_summary_stats(SPEI_by_basin, basin_name=b, modelnames=[m], period=(2070,2100)) for b in basin_names} 
                    for m in modelnames}

## All stats for historical period
fig1, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4), sharex=True)
for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
    pg = ag/a # percent glaciated
    number_b = []
    duration_b = []
    severity_b = []
    for m in modelnames:
        number_b.append(basin_stats_bymodel_hist[m][b][0][1]-basin_stats_bymodel_hist[m][b][0][0])
        duration_b.append(basin_stats_bymodel_hist[m][b][1][1]-basin_stats_bymodel_hist[m][b][1][0])
        severity_b.append(-1*(basin_stats_bymodel_hist[m][b][2][1]-basin_stats_bymodel_hist[m][b][2][0]))
    ax1.errorbar(pg, np.nanmean(number_b), 
                 yerr=((np.nanmean(number_b)-np.nanmin(number_b), np.nanmax(number_b)-np.nanmean(number_b)),), 
                 color='k', marker='o', lw=1.0)
    ax2.errorbar(pg, np.nanmean(duration_b), 
                 yerr=(((np.nanmean(duration_b)-np.nanmin(duration_b), np.nanmax(duration_b)-np.nanmean(duration_b)),)), 
                 color='k', marker='d', lw=1.0)
    ax3.errorbar(pg, np.nanmean(severity_b), 
                 yerr=(((np.nanmean(severity_b)-np.nanmin(severity_b), np.nanmax(severity_b)-np.nanmean(severity_b)),)), 
                 color='k', marker='*', lw=1.0)
ax1.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts 1980-2010', xscale='log')
ax2.set(xlabel='Glacier area fraction', ylabel='Diff. drought duration 1980-2010', xscale='log')
ax3.set(xlabel='Glacier area fraction', ylabel='Diff. drought deficit 1980-2010', xscale='log',
        xlim=(1E-4, 0.22))
for ax in (ax1, ax2, ax3):
    ax.axhline(0, ls=':', lw=1.0, color='k')
plt.tight_layout()
plt.show()

# ## Drought number over time
fig1, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4), sharex=True, sharey=True)
for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
    pg = ag/a # percent glaciated
    number_hist = []
    number_midC = []
    number_endC = []
    for m in modelnames:
        number_hist.append(basin_stats_bymodel_hist[m][b][0][1]-basin_stats_bymodel_hist[m][b][0][0])
        number_midC.append(basin_stats_bymodel_midC[m][b][0][1]-basin_stats_bymodel_midC[m][b][0][0])
        number_endC.append(basin_stats_bymodel_endC[m][b][0][1]-basin_stats_bymodel_endC[m][b][0][0])
    midC_v_hist = np.nanmean(number_midC)-np.nanmean(number_hist)
    if midC_v_hist >0: # buffering increasing
        midC_color='b'
    elif midC_v_hist==0:
        midC_color='k'
    elif midC_v_hist<0:
        midC_color='r'
    endC_v_midC = np.nanmean(number_endC)-np.nanmean(number_midC)
    # if endC_v_midC >0:
    #     endC_color='b'
    # elif endC_v_midC==0:
    #     endC_color='k'
    # elif endC_v_midC<0:
    #     endC_color='r'
    endC_v_hist = np.nanmean(number_endC)-np.nanmean(number_hist)
    if endC_v_hist >0:
        endC_color='b'
    elif endC_v_hist==0:
        endC_color='k'
    elif endC_v_hist<0:
        endC_color='r'
    ax1.errorbar(pg, np.nanmean(number_hist), 
                 yerr=(((np.nanmean(number_hist)-np.nanmin(number_hist), np.nanmax(number_hist)-np.nanmean(number_hist)),)), 
                 color='k', marker='o', lw=1.0)
    ax2.errorbar(pg, np.nanmean(number_midC), 
                 yerr=(((np.nanmean(number_midC)-np.nanmin(number_midC), np.nanmax(number_midC)-np.nanmean(number_midC)),)), 
                 color=midC_color, marker='o', lw=1.0)
    ax3.errorbar(pg, np.nanmean(number_endC), 
                 yerr=(((np.nanmean(number_endC)-np.nanmin(number_endC), np.nanmax(number_endC)-np.nanmean(number_endC)),)), 
                 color=endC_color, marker='o', lw=1.0)
    # ax3.errorbar(pg, np.nanmean(severity_b), color='k')
ax1.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts 1980-2010', xscale='log')
ax2.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts 2030-2060', xscale='log')
ax3.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts 2070-2100', xscale='log',
        xlim=(1E-4, 0.22))
for ax in (ax1, ax2, ax3):
    ax.axhline(0, ls=':', lw=1.0, color='k')
plt.tight_layout()
plt.show()


# ## Difference of the difference - change over time
# fig1, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4), sharex=True, sharey=True)
# for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
#     pg = ag/a # percent glaciated
#     number_hist = []
#     number_midC = []
#     number_endC = []
#     for m in modelnames:
#         number_hist.append(basin_stats_bymodel_hist[m][b][0][1]-basin_stats_bymodel_hist[m][b][0][0])
#         number_midC.append(basin_stats_bymodel_midC[m][b][0][1]-basin_stats_bymodel_midC[m][b][0][0])
#         number_endC.append(basin_stats_bymodel_endC[m][b][0][1]-basin_stats_bymodel_endC[m][b][0][0])
#     ax1.errorbar(pg, np.nanmean(number_hist), 
#                  yerr=(((np.nanmean(number_hist)-np.nanmin(number_hist), np.nanmax(number_hist)-np.nanmean(number_hist)),)), 
#                  color='k', marker='o', lw=1.0)
#     ax2.errorbar(pg, np.nanmean(number_midC)-np.nanmean(number_hist), 
#                  yerr=(((np.nanmean(number_midC)-np.nanmin(number_midC), np.nanmax(number_midC)-np.nanmean(number_midC)),)), 
#                  color='k', marker='o', lw=1.0)
#     ax3.errorbar(pg, np.nanmean(number_endC)-np.nanmean(number_hist), 
#                  yerr=(((np.nanmean(number_endC)-np.nanmin(number_endC), np.nanmax(number_endC)-np.nanmean(number_endC)),)), 
#                  color='k', marker='o', lw=1.0)
#     # ax3.errorbar(pg, np.nanmean(severity_b), color='k')
# ax1.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts 1980-2010', xscale='log')
# ax2.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts, 2030-2060 vs historical', xscale='log')
# ax3.set(xlabel='Glacier area fraction', ylabel='Diff. number of droughts, 2070-2100 vs historical', xscale='log',
#         xlim=(1E-4, 0.22))
# for ax in (ax1, ax2, ax3):
#     ax.axhline(0, ls=':', lw=1.0, color='k')
# plt.tight_layout()
# plt.show()

  
## Composite of all stats over time
color_fam = cm.get_cmap('tab20b')
inc_color=color_fam(5)
dec_color=color_fam(17)

fig3, ((ax1,ax2, ax3), 
       (ax4,ax5,ax6), 
       (ax7,ax8,ax9)) = plt.subplots(3,3, sharex=True, sharey='row', figsize=(10,12))
for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
    pg = ag/a # percent glaciated
    number_b = []
    duration_b = []
    severity_b = []
    number_midC = []
    number_endC = []
    duration_midC = []
    duration_endC = []
    severity_midC = []
    severity_endC = []

    for m in modelnames:
        number_b.append(basin_stats_bymodel_hist[m][b][0][1]-basin_stats_bymodel_hist[m][b][0][0])
        duration_b.append(basin_stats_bymodel_hist[m][b][1][1]-basin_stats_bymodel_hist[m][b][1][0])
        severity_b.append(-1*(basin_stats_bymodel_hist[m][b][2][1]-basin_stats_bymodel_hist[m][b][2][0]))
        number_midC.append(basin_stats_bymodel_midC[m][b][0][1]-basin_stats_bymodel_midC[m][b][0][0])
        number_endC.append(basin_stats_bymodel_endC[m][b][0][1]-basin_stats_bymodel_endC[m][b][0][0])
        duration_midC.append(basin_stats_bymodel_midC[m][b][1][1]-basin_stats_bymodel_midC[m][b][1][0])
        duration_endC.append(basin_stats_bymodel_endC[m][b][1][1]-basin_stats_bymodel_endC[m][b][1][0])
        severity_midC.append(-1*(basin_stats_bymodel_midC[m][b][2][1]-basin_stats_bymodel_midC[m][b][2][0]))
        severity_endC.append(-1*(basin_stats_bymodel_endC[m][b][2][1]-basin_stats_bymodel_endC[m][b][2][0]))
   
    ## Color code changes over time
    midC_v_hist_n = np.nanmean(number_midC)-np.nanmean(number_b)
    if midC_v_hist_n >0: # buffering on number increasing
        midC_color_n=inc_color
        midC_marker_n='^'
    elif midC_v_hist_n==0:
        midC_color_n='k'
        midC_marker_n='o'
    elif midC_v_hist_n<0:
        midC_color_n=dec_color
        midC_marker_n='v'
    endC_v_hist_n = np.nanmean(number_endC)-np.nanmean(number_b)
    if endC_v_hist_n >0:
        endC_color_n=inc_color
        endC_marker_n='^'
    elif endC_v_hist_n==0:
        endC_color_n='k'
        endC_marker_n='o'
    elif endC_v_hist_n<0:
        endC_color_n=dec_color
        endC_marker_n='v'
        
    midC_v_hist_d = np.nanmean(duration_midC)-np.nanmean(duration_b)
    if midC_v_hist_d >0: # buffering on duration increasing
        midC_color_d=inc_color
        midC_marker_d='^'
    elif midC_v_hist_d==0:
        midC_color_d='k'
        midC_marker_d='o'
    elif midC_v_hist_d<0:
        midC_color_d=dec_color
        midC_marker_d='v'
    endC_v_hist_d = np.nanmean(duration_endC)-np.nanmean(duration_b)
    if endC_v_hist_d >0:
        endC_color_d=inc_color
        endC_marker_d='^'
    elif endC_v_hist_d==0:
        endC_color_d='k'
        endC_marker_d='o'
    elif endC_v_hist_d<0:
        endC_color_d=dec_color
        endC_marker_d='v'
        
    midC_v_hist_s = np.nanmean(severity_midC)-np.nanmean(severity_b)
    if midC_v_hist_s >0: # buffering on duration increasing
        midC_color_s=inc_color
        midC_marker_s='^'
    elif midC_v_hist_s==0:
        midC_color_s='k'
        midC_marker_s='o'
    elif midC_v_hist_s<0:
        midC_color_s=dec_color
        midC_marker_s='v'
    endC_v_hist_s = np.nanmean(severity_endC)-np.nanmean(severity_b)
    if endC_v_hist_s >0:
        endC_color_s=inc_color
        endC_marker_s='^'
    elif endC_v_hist_s==0:
        endC_color_s='k'
        endC_marker_s='o'
    elif endC_v_hist_s<0:
        endC_color_s=dec_color
        endC_marker_s='v'
    ## First column: historical
    ax1.errorbar(pg, np.nanmean(number_b), 
                 yerr=((np.nanmean(number_b)-np.nanmin(number_b), np.nanmax(number_b)-np.nanmean(number_b)),), 
                 color='k', marker='o', lw=1.0)
    ax4.errorbar(pg, np.nanmean(duration_b), 
                 yerr=(((np.nanmean(duration_b)-np.nanmin(duration_b), np.nanmax(duration_b)-np.nanmean(duration_b)),)), 
                 color='k', marker='o', lw=1.0, alpha=0.8)
    ax7.errorbar(pg, np.nanmean(severity_b), 
                 yerr=(((np.nanmean(severity_b)-np.nanmin(severity_b), np.nanmax(severity_b)-np.nanmean(severity_b)),)), 
                 color='k', marker='o', lw=1.0)
    ## Second column: mid-c
    ax2.errorbar(pg, np.nanmean(number_midC), 
                 yerr=(((np.nanmean(number_midC)-np.nanmin(number_midC), np.nanmax(number_midC)-np.nanmean(number_midC)),)), 
                 color=midC_color_n, marker=midC_marker_n, lw=1.0)
    ax5.errorbar(pg, np.nanmean(duration_midC), 
                 yerr=(((np.nanmean(duration_midC)-np.nanmin(duration_midC), np.nanmax(duration_midC)-np.nanmean(duration_midC)),)), 
                 color=midC_color_d, marker=midC_marker_d, lw=1.0)
    ax8.errorbar(pg, np.nanmean(severity_midC), 
                 yerr=(((np.nanmean(severity_midC)-np.nanmin(severity_midC), np.nanmax(severity_midC)-np.nanmean(severity_midC)),)), 
                 color=midC_color_s, marker=midC_marker_s, lw=1.0)
    ## Third column: end of century
    ax3.errorbar(pg, np.nanmean(number_endC), 
                 yerr=(((np.nanmean(number_endC)-np.nanmin(number_endC), np.nanmax(number_endC)-np.nanmean(number_endC)),)), 
                 color=endC_color_n, marker=endC_marker_n, lw=1.0)
    ax6.errorbar(pg, np.nanmean(duration_endC), 
                 yerr=(((np.nanmean(duration_endC)-np.nanmin(duration_endC), np.nanmax(duration_endC)-np.nanmean(duration_endC)),)), 
                 color=endC_color_d, marker=endC_marker_d, lw=1.0)
    ax9.errorbar(pg, np.nanmean(severity_endC), 
                 yerr=(((np.nanmean(severity_endC)-np.nanmin(severity_endC), np.nanmax(severity_endC)-np.nanmean(severity_endC)),)), 
                 color=endC_color_s, marker=endC_marker_s, lw=1.0)

ax1.set(ylabel=r'$\Delta$ Number', title='Historical (1980-2010)', xscale='log')
ax2.set(title='Mid-21st Cent. (2030-2060)')
ax3.set(title='End 21st Cent. (2070-2100)')
ax4.set(ylabel=r'$\Delta$ Duration', xscale='log')
ax7.set(ylabel=r'$\Delta$ Severity', 
        xscale='log',xlabel='Glacier area fraction', xlim=(1E-4, 0.23))
ax8.set(xlabel='Glacier area fraction', xlim=(1E-4, 0.23))
ax9.set(xlabel='Glacier area fraction', xlim=(1E-4, 0.23))
for ax in (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9):
    ax.axhline(0, ls=':', lw=1.0, color='k')
    