#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Calculate rank correlation between basin percent glaciated area and glacial drought buffering effect
Created on Tue May 18 14:03:41 2021

@author: lizz
"""


import numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import gSPEI as gSPEI
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from cycler import cycler
import seaborn as sns
import pandas as pd

#Change default colour cycle
mpl.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20c([0,4,8,12,16]))
mpl.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20b([0,4,8,12,16]))


## Labels: (P)arametric or (NP)nonparametric;
## Standardization (1) lumped or (2) split by starting month
fpath_NP2 = '/home/joncka/glacier_drought/nonparametric-var_stom_c/'

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

BasinArea = np.asarray(BasinArea)
basin_glacier_area = np.asarray(basin_glacier_area)

rcp = scenarios[0] #rcp45

yrs = np.linspace(1900, 2101, num=2412)
SPEI_by_model_C = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath_NP2+'NRunoff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, rcp)
    wrunoff_f_m = fpath_NP2+'WRunOff_{}_{}_{}_Conduct.txt'.format(integration_times[0], m, rcp)
    SPEI_by_model_C[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model_C[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model_C[m]['diff'] = SPEI_by_model_C[m]['WRunoff'] - SPEI_by_model_C[m]['NRunoff']

    
#Load in climte data and establish mean precip and aridity index
y0 = 1980
y1 = 2010
prcp = np.zeros(len(BasinArea))
pet = np.zeros(len(BasinArea))
indx = np.logical_and(yrs >= y0, yrs <= y1)
for m in modelnames:
    P = np.loadtxt('model_P_PET/' + m + '_' + rcp + '_PREC.txt')[:,indx].mean(axis=1)/BasinArea
    PET = np.loadtxt('model_P_PET/' + m + '_' + rcp + '_PET.txt')[:,indx].mean(axis=1)/BasinArea
    prcp += P
    pet += PET
prcp /= len(modelnames)
pet /= len(modelnames)

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


#Put stats into single array for k-means clustering
X = np.asarray([mean_number_b, mean_duration_b, mean_severity_b, mean_number_midC, mean_number_endC, mean_duration_midC, mean_duration_endC, mean_severity_midC, mean_severity_endC]).transpose()

#Scale dataset so that for each stat: mu=0, sigma=1
X_ = StandardScaler().fit_transform(X)


k = 2 #Number of clusters

#Instantiate kmeans class
kmeans = KMeans(init="random",n_clusters=k,n_init=10,max_iter=300,random_state=42)

#Perform kmeans
kmeans.fit(X_)

#Calc silhouette_score
score = silhouette_score(X_, kmeans.labels_)

#Get median stats for each cluster
med_num = np.zeros((k,3))
med_dur = np.zeros((k,3))
med_sev = np.zeros((k,3))

for i in range(k):
    indx = np.where(kmeans.labels_ == i)
    if len(indx[0]) > 1:
        med_num[i,:] = np.median(X[indx,:].squeeze()[:,[0,3,4]], axis=0)
        med_dur[i,:] = np.median(X[indx,:].squeeze()[:,[1,5,6]], axis=0)
        med_sev[i,:] = np.median(X[indx,:].squeeze()[:,[2,7,8]], axis=0)
    else:
        med_num[i,:] = np.median(X[indx,:].squeeze()[[0,3,4]], axis=0)
        med_dur[i,:] = np.median(X[indx,:].squeeze()[[1,5,6]], axis=0)
        med_sev[i,:] = np.median(X[indx,:].squeeze()[[2,7,8]], axis=0)



#Stats for box plots
X_box = pd.DataFrame(np.hstack((np.vstack((np.zeros(len(X)), kmeans.labels_, X[:,0], X[:,1], X[:,2])), np.vstack((np.zeros(len(X))+1, kmeans.labels_, X[:,3], X[:,5], X[:,7])), np.vstack((np.zeros(len(X))+2, kmeans.labels_, X[:,4], X[:,6], X[:,8])))).transpose(), columns=['Period','Cluster','Number','Duration','Severity'])
X_box['Period'].loc[X_box['Period']==0] = 'baseline'
X_box['Period'].loc[X_box['Period']==1] = 'mid'
X_box['Period'].loc[X_box['Period']==2] = 'late'
X_box['Cluster'] = (X_box['Cluster'] + 1).values.astype('int')


def set_all_font_sizes(ax, fontsize):
    ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize)
    ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    
#Plot results
fig = plt.figure(constrained_layout=True, figsize=(7,8))
gs = fig.add_gridspec(3, 2)
ax_scat = fig.add_subplot(gs[1:,:])
ax_num = fig.add_subplot(gs[0,0])
#ax_dur = fig.add_subplot(gs[0,1])
ax_sev = fig.add_subplot(gs[0,1])


#Plot stats for each cluster
bx = sns.boxplot(x = X_box['Period'], y=X_box['Number'], hue=X_box['Cluster'], ax=ax_num, )
ax_num.set(ylabel='$\Delta$ Number', xlabel='')
bx.legend_.remove()
set_all_font_sizes(ax_num,12)

# bx = sns.boxplot(x = X_box['Period'], y=X_box['Duration'], hue=X_box['Cluster'], ax=ax_dur, )
# ax_dur.set(ylabel='$\Delta$ Duration', xlabel='')
# bx.legend_.remove()
# set_all_font_sizes(ax_dur,12)

bx = sns.boxplot(x = X_box['Period'], y=X_box['Severity'], hue=X_box['Cluster'], ax=ax_sev, )
ax_sev.set(ylabel='$\Delta$ Severity', xlabel='')
bx.legend_.remove()
set_all_font_sizes(ax_sev,12)


# #Plot mean stats for each cluster
# lines0 = ax_num.plot(med_num.transpose(), zorder=2, linestyle='--')
# ax_num.set(ylabel='$\Delta$ Number', xticks=[0,1,2], xticklabels=['baseline','mid','late'])
# set_all_font_sizes(ax_num,12)
# ax_dur.plot(med_dur.transpose(), zorder=3, linestyle='--')
# ax_dur.set(ylabel = '$\Delta$ Duration', xticks=[0,1,2], xticklabels=['baseline','mid','late'])
# set_all_font_sizes(ax_dur,12)
# ax_sev.plot(med_sev.transpose(), zorder=4, linestyle='--')
# ax_sev.set(ylabel = '$\Delta$ Severity', xticks=[0,1,2], xticklabels=['baseline','mid','late'])
# set_all_font_sizes(ax_sev,12)

# for i, l in enumerate(lines0):
#     indx = np.where(kmeans.labels_ == i)
#     if len(indx[0]) > 1:
#         dat_y = X[indx,:].squeeze()[:,[0,3,4]]
#         dat_x = np.zeros(dat_y.shape); dat_x[:,1] = 1; dat_x[:,2] = 2
#     else:
#         dat_y = X[indx,:].squeeze()[[0,3,4]]
#         dat_x = np.zeros(dat_y.shape); dat_x[1] = 1; dat_x[2] = 2
#     ax_num.plot(dat_x.transpose(), dat_y.transpose(), marker='.', linewidth=0.5, color=l.get_color(), zorder=1, alpha=0.5)
    
#     if len(indx[0]) > 1:
#         dat_y = X[indx,:].squeeze()[:,[1,5,6]]
#         dat_x = np.zeros(dat_y.shape); dat_x[:,1] = 1; dat_x[:,2] = 2
#     else:
#         dat_y = X[indx,:].squeeze()[[1,5,6]]
#         dat_x = np.zeros(dat_y.shape); dat_x[1] = 1; dat_x[2] = 2
#     ax_dur.plot(dat_x.transpose(), dat_y.transpose(), marker='.', linewidth=0.5, color=l.get_color(), zorder=1, alpha=0.5)
    
#     if len(indx[0]) > 1:
#         dat_y = X[indx,:].squeeze()[:,[2,7,8]]
#         dat_x = np.zeros(dat_y.shape); dat_x[:,1] = 1; dat_x[:,2] = 2
#     else:
#         dat_y = X[indx,:].squeeze()[[2,7,8]]
#         dat_x = np.zeros(dat_y.shape); dat_x[1] = 1; dat_x[2] = 2
#     ax_sev.plot(dat_x.transpose(), dat_y.transpose(), marker='.', linewidth=0.5, color=l.get_color(), zorder=1, alpha=0.5)

#Plot basin characteristics for each cluster
for i in range(k):
    indx = np.where(kmeans.labels_ == i)
    ax_scat.plot(100*basin_glacier_area[indx]/BasinArea[indx], prcp[indx]/pet[indx], linewidth=0, marker='.', markersize=15, label='cluster ' + str(i+1))
    #ax_scat.plot(BasinArea[indx], prcp[indx]/pet[indx], linewidth=0, marker='.', markersize=15, label=labels[i])
#ax_scat.set_xlabel('Basin Area ($m^2$)')
ax_scat.set_xlabel('Glacier coverage (%)')
ax_scat.set_ylabel('Aridity index (-)')
ax_scat.set_xscale('log')
ax_scat.legend()
set_all_font_sizes(ax_scat,12)
fig.savefig('kmeans_k2.png',dpi=300)



    

