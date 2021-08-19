#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Calculate rank correlation between basin percent glaciated area and glacial drought buffering effect
Created on Tue May 18 14:03:41 2021

@author: lizz
"""


import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
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

  
## Composite of all stats over time - Spearman vs % glaciated
mean_number_b = []
mean_duration_b = []
mean_severity_b = []
mean_number_midC = []
mean_number_endC = []
mean_duration_midC = []
mean_duration_endC = []
mean_severity_midC = []
mean_severity_endC = []
percent_glac = []
for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
    percent_glac.append(ag/a) # percent glaciated

    mean_number_b.append(np.nanmean([basin_stats_bymodel_hist[m][b][0][1]-
                                basin_stats_bymodel_hist[m][b][0][0] for m in modelnames]))
    mean_duration_b.append(np.nanmean([basin_stats_bymodel_hist[m][b][1][1]-
                                       basin_stats_bymodel_hist[m][b][1][0] for m in modelnames]))
    mean_severity_b.append(-1*np.nanmean([basin_stats_bymodel_hist[m][b][2][1]-
                                           basin_stats_bymodel_hist[m][b][2][0] for m in modelnames]))
    mean_number_midC.append(np.nanmean([basin_stats_bymodel_midC[m][b][0][1]-
                                        basin_stats_bymodel_midC[m][b][0][0] for m in modelnames]))
    mean_number_endC.append(np.nanmean([basin_stats_bymodel_endC[m][b][0][1]-
                                        basin_stats_bymodel_endC[m][b][0][0] for m in modelnames]))
    mean_duration_midC.append(np.nanmean([basin_stats_bymodel_midC[m][b][1][1]-
                                          basin_stats_bymodel_midC[m][b][1][0] for m in modelnames]))
    mean_duration_endC.append(np.nanmean([basin_stats_bymodel_endC[m][b][1][1]-
                                          basin_stats_bymodel_endC[m][b][1][0] for m in modelnames]))
    mean_severity_midC.append(-1*np.nanmean([basin_stats_bymodel_midC[m][b][2][1]-
                                             basin_stats_bymodel_midC[m][b][2][0] for m in modelnames]))
    mean_severity_endC.append(-1*np.nanmean([basin_stats_bymodel_endC[m][b][2][1]-
                                             basin_stats_bymodel_endC[m][b][2][0] for m in modelnames]))
   
## Compute Spearman's corrs for these buffering measures
rho_n_b = stats.spearmanr(percent_glac, mean_number_b)
rho_d_b = stats.spearmanr(percent_glac, mean_duration_b)
rho_s_b = stats.spearmanr(percent_glac, mean_severity_b)
rho_n_m = stats.spearmanr(percent_glac, mean_number_midC)
rho_d_m = stats.spearmanr(percent_glac, mean_duration_midC)
rho_s_m = stats.spearmanr(percent_glac, mean_severity_midC)
rho_n_e = stats.spearmanr(percent_glac, mean_number_endC)
rho_d_e = stats.spearmanr(percent_glac, mean_duration_endC)
rho_s_e = stats.spearmanr(percent_glac, mean_severity_endC)