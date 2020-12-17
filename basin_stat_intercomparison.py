#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Compare drought buffering effect compared with basin stats

Created on Thu Dec 10 15:39:17 2020

@author: lizz
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gSPEI as gSPEI
import collections

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
    norunoff_f_m = fpath_NP2+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath_NP2+'WRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model_C[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model_C[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model_C[m]['diff'] = SPEI_by_model_C[m]['WRunoff'] - SPEI_by_model_C[m]['NRunoff']
    
## Re-structure dictionary and create pandas DataFrames aggregated by basin
SPEI_by_basin = gSPEI.sort_models_to_basins(SPEI_by_model_C)

## Analyse multi-model ensemble mean & quantiles for drought statistics
r_w = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'WRunoff')
r_n = gSPEI.basin_ensemble_mean(SPEI_by_basin, 'TARIM', 'NRunoff')

def find_droughts(series, threshold=-1, period=(1980, 2100), t_array=yrs):
    """Identify droughts in a timeseries of SPEI

    Parameters
    ----------
    series : pandas Series
        SPEI time series.
    threshold : float, optional
        Cutoff value that must be reached for negative SPEI to be called
        a 'drought'. The default is -1.
    period : tuple, optional
        Start and end years between which to compute stats.  Default is (1980,2100).
    t_array : array, optional
        Array of decimal years corresponding to SPEI timeseries.  
        Default is yrs = np.linspace(1900, 2101, num=2412)

    Returns
    -------
    droughts : OrderedDict
        SPEI values sorted into droughts, with dict keys the indices in yrs of drought onset.

    """
    start_idx = np.argwhere(t_array >= period[0])[0]
    end_idx = np.argwhere(t_array <= period[1])[-1]
    
    droughts = collections.OrderedDict()
    for i,v in enumerate(series):
        if (v>=0 or np.isnan(v)): 
            pass
        else:
            if i==0:
                current_drought = [] #create new for start of time
            elif (series[i-1]>=0 or np.isnan(series[i-1])):
                current_drought = [] #create new for start of drought
            current_drought.append(v)
            if i==len(series)-1:
                droughts[i] = current_drought #write it out if time ending
            elif series[i+1]>=0:
                droughts[i] = current_drought #write it out if drought ending
    droughts_thresholded = collections.OrderedDict(
        {k: droughts[k] for k in droughts.keys() if sum(
            np.less_equal(droughts[k], threshold))>0})
    droughts_trimmed_0 = collections.OrderedDict({k: droughts_thresholded[k] 
                                                  for k in droughts_thresholded.keys() if k>=start_idx})
    droughts_trimmed = collections.OrderedDict({k: droughts_trimmed_0[k] 
                                                  for k in droughts_trimmed_0.keys() if k<end_idx}) #not inclusive here
    return droughts_trimmed

def basin_summary_stats(SPEI_dict, basin_name, modelnames, threshold=-1, 
                        period=(1980,2100), t_array=yrs):
    """Find droughts in all models for a given basin, and report incidence,
    duration, and severity for with/without glacial runoff cases.
    

    Parameters
    ----------
    SPEI_dict : dictionary
        SPEI values with keys [<basin>][<scenario>][<model>].
    basin_name: str
        name of basin 
    modelnames : list
        names of GCMs used.
    threshold : float, optional
        SPEI value that must be reached to define a drought. The default is -1.
    period : tuple, optional
        Start and end years between which to compute stats.  Default is (1980,2100).
    t_array : array, optional
        Array of decimal years corresponding to SPEI timeseries.  
        Default is yrs = np.linspace(1900, 2101, num=2412)
    
    Returns
    -------
    [(mean_drtnumber_w, mean_drtnumber_n), (mean_drtdur_w, mean_drtdur_n), 
    (mean_drtsev_w, mean_drtsev_n)]
        Drought incidence, duration, and severity pairs for this basin in the given period.

    """
    
    keys_w_bymodel = {m: [] for m in modelnames}
    keys_n_bymodel = {m: [] for m in modelnames}
    drt_dur_w_bymodel = {m: [] for m in modelnames}
    drt_dur_n_bymodel = {m: [] for m in modelnames}
    drt_sev_w_bymodel = {m: [] for m in modelnames}
    drt_sev_n_bymodel = {m: [] for m in modelnames}
    
    for m in modelnames:
        ser_w = SPEI_dict[basin_name]['WRunoff'][m]
        ser_n = SPEI_dict[basin_name]['NRunoff'][m]
        droughts_w = find_droughts(ser_w, threshold=threshold, period=period, t_array=t_array)
        droughts_n = find_droughts(ser_n, threshold=threshold, period=period, t_array=t_array)
        # drts_w_trimmed = collections.OrderedDict({k: droughts_w[k] for k in droughts_w.keys() if k>=960}) #960 is first index where yrs > 1980 (glacier model on)
        # drts_n_trimmed = collections.OrderedDict({k: droughts_n[k] for k in droughts_n.keys() if k>=960})
        drts_w_trimmed = droughts_w
        drts_n_trimmed = droughts_n
        
        keys_w_bymodel[m] = drts_w_trimmed.keys()
        keys_n_bymodel[m] = drts_n_trimmed.keys()
        drt_dur_w_bymodel[m] = [len(droughts_w[k]) for k in drts_w_trimmed.keys()]
        drt_dur_n_bymodel[m] = [len(droughts_n[k]) for k in drts_n_trimmed.keys()]
        drt_sev_w_bymodel[m] = [sum(droughts_w[k]) for k in drts_w_trimmed.keys()]
        drt_sev_n_bymodel[m] = [sum(droughts_n[k]) for k in drts_n_trimmed.keys()]
    
    ## Report statistics
    mean_drtnumber_w = np.mean([len(keys_w_bymodel[m]) for m in modelnames])
    mean_drtnumber_n = np.mean([len(keys_n_bymodel[m]) for m in modelnames])
    mean_drtdur_w = np.mean([np.mean(drt_dur_w_bymodel[m]) for m in modelnames])
    mean_drtdur_n = np.mean([np.mean(drt_dur_n_bymodel[m]) for m in modelnames])
    mean_drtsev_w = np.mean([np.mean(drt_sev_w_bymodel[m]) for m in modelnames])
    mean_drtsev_n = np.mean([np.mean(drt_sev_n_bymodel[m]) for m in modelnames])
    
    return [(mean_drtnumber_w, mean_drtnumber_n), (mean_drtdur_w, mean_drtdur_n), (mean_drtsev_w, mean_drtsev_n)]

basin_stats = {b: basin_summary_stats(SPEI_by_basin, basin_name=b, modelnames=modelnames) for b in basin_names}

## plot by basin initial glacial area
fig, (ax1,ax2,ax3) = plt.subplots(1,3)
for b, a in zip(basin_names, basin_glacier_area):
    ax1.scatter(a, basin_stats[b][0][0], color='DarkBlue')
    ax1.scatter(a, basin_stats[b][0][1], color='Orange')
    ax2.scatter(a, basin_stats[b][1][0], color='DarkBlue')
    ax2.scatter(a, basin_stats[b][1][1], color='Orange')
    ax3.scatter(a, basin_stats[b][2][0], color='DarkBlue')
    ax3.scatter(a, basin_stats[b][2][1], color='Orange')
ax1.set(xlabel='Basin glacier area [km2]', ylabel='Mean number of droughts 1980-2100', xscale='log')
ax2.set(xlabel='Basin glacier area [km2]', ylabel='Mean drought duration 1980-2100', xscale='log')
ax3.set(xlabel='Basin glacier area [km2]', ylabel='Mean drought severity 1980-2100', xscale='log')
plt.show()


fig2, (ax4,ax5,ax6) = plt.subplots(1,3)
for b, a, ag in zip(basin_names, BasinArea, basin_glacier_area):
    pg = ag/a # percent glaciated
    ax4.scatter(pg, basin_stats[b][0][0], color='DarkBlue')
    ax4.scatter(pg, basin_stats[b][0][1], color='Orange')
    ax5.scatter(pg, basin_stats[b][1][0], color='DarkBlue')
    ax5.scatter(pg, basin_stats[b][1][1], color='Orange')
    ax6.scatter(pg, basin_stats[b][2][0], color='DarkBlue')
    ax6.scatter(pg, basin_stats[b][2][1], color='Orange')
ax4.set(xlabel='Glacier area fraction', ylabel='Mean number of droughts 1980-2100')
ax5.set(xlabel='Glacier area fraction', ylabel='Mean drought duration 1980-2100')
ax6.set(xlabel='Glacier area fraction', ylabel='Mean drought severity 1980-2100')
plt.show()   