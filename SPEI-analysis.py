## Read in, plot, and analyse SPEI data
## 22 Aug 2019
## Code: EHU | Data: SC

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from gSPEI import *

fpath = 'Documents/GitHub/glacial-SPEI/data/SPEI_Files/' 

## Settings in filenames
integration_times = np.arange(3, 28, 4) # all SPEI integration times used
modelnames = ['CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GISS-E2-R', 'INMCM4', 'MIROC-ESM', 'NorESM1-M'] # all models used in comparison
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios

## Basins in the order they are written
basin_names = ['INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE/VINCES',
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE']

yrs = np.linspace(1900, 2101, num=2412)


## Compare effect across models - read in all to dict
SPEI_by_model = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath+'NRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath+'WRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']

### 30-yr running means
plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model)
plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, which='NRunoff', cmap_name='Greys')
plot_runmean_comparison(basin_id=1, permodel_dict=SPEI_by_model)

## 30-yr period means and variance
modelmeans_1950_1980 = {m: [] for m in modelnames} #dictionary indexed by model name
modelvar_1950_1980 = {m: [] for m in modelnames}
modelmeans_2070_2100 = {m: [] for m in modelnames}
modelvar_2070_2100 = {m: [] for m in modelnames}
mean_shifts = {m: [] for m in modelnames}
var_shifts = {m: [] for m in modelnames}

for m in modelnames:
    means_i = [np.nanmean(SPEI_by_model[m]['diff'][j][600:971]) for j in range(len(basin_names))]
    var_i = [np.nanvar(SPEI_by_model[m]['diff'][j][600:971]) for j in range(len(basin_names))]
    means_f = [np.nanmean(SPEI_by_model[m]['diff'][j][2039:2410]) for j in range(len(basin_names))]
    var_f = [np.nanvar(SPEI_by_model[m]['diff'][j][2039:2410]) for j in range(len(basin_names))]
    modelmeans_1950_1980[m] = means_i
    modelvar_1950_1980[m] = var_i
    modelmeans_2070_2100[m] = means_f
    modelvar_2070_2100[m] = var_f
    mean_shifts[m] = np.array(means_f) - np.array(means_i)
    var_shifts[m] = np.array(var_f) - np.array(var_i)

basin_mean_shifts = {b: [] for b in basin_names} #dictionary indexed by basin name
basin_var_shifts = {b: [] for b in basin_names} #dictionary indexed by model name
basin_meanshift_meds = [] #arrays for plotting
basin_varshift_meds = []
basin_meanshift_range = []
basin_varshift_range = []

for i, b in enumerate(basin_names):
    bmeans_i = [np.nanmean(SPEI_by_model[m]['WRunoff'][i][600:971]) for m in modelnames]
    bvar_i = [np.nanvar(SPEI_by_model[m]['WRunoff'][i][600:971]) for m in modelnames]
    bmeans_f = [np.nanmean(SPEI_by_model[m]['WRunoff'][i][2039:2410]) for m in modelnames]
    bvar_f = [np.nanvar(SPEI_by_model[m]['WRunoff'][i][2039:2410]) for m in modelnames]
    basin_mean_shifts[b] = np.array(bmeans_f) - np.array(bmeans_i)
    basin_var_shifts[b] = np.array(bvar_f) - np.array(bvar_i)
    basin_meanshift_meds.append(np.nanmedian(basin_mean_shifts[b]))
    basin_varshift_meds.append(np.nanmedian(basin_var_shifts[b]))
    basin_meanshift_range.append(np.nanmax(basin_mean_shifts[b]) - np.nanmin(basin_mean_shifts[b]))
    basin_varshift_range.append(np.nanmax(basin_var_shifts[b]) - np.nanmin(basin_var_shifts[b]))


plt.figure('Mean and variance shifts, 2070-2100 versus 1950-1980, per basin for WRunoff case')
plt.errorbar(x=basin_meanshift_meds, y=basin_varshift_meds, xerr=basin_meanshift_range, yerr=basin_varshift_range, ls='')
plt.axes().set_xlabel('Difference in 30-yr mean SPEI', fontsize=16)
plt.axes().set_ylabel('Difference in SPEI variance', fontsize=16)
plt.axes().set_ylim(-0.5, 2.5)
plt.axes().set_xlim(-1, 5)
plt.show()


## Calculate changes due to glacial effect at end of century, using gSPEI functions
bas_glac_meanmed, mean_spread = glacial_meandiff(SPEI_by_model)
bas_glac_varmed, var_spread = glacial_vardiff(SPEI_by_model)

#plt.figure('Mean and variance shifts due to glacial effects in 2070-2100')
#plt.errorbar(x=bas_glac_meanmed, y=bas_glac_varmed, xerr=mean_spread, yerr=var_spread, ls='', marker='d', elinewidth=2.0, color='DarkBlue')
#plt.axes().set_xlabel('Difference in mean SPEI', fontsize=16)
#plt.axes().set_ylabel('Difference in SPEI variance', fontsize=16)
#plt.axes().set_ylim(-1.5, 1.0)
#plt.axes().set_xlim(-0.5, 5)
#plt.show()


### Plot running SPEI variance
SPEIvar_by_model = {m: {'NRunoff': [], 'WRunoff': []} for m in modelnames}
varwindow = 360 # number of months to window in rolling variance
for m in modelnames:
    for i in range(len(basin_names)):
        nr = pd.Series(SPEI_by_model[m]['NRunoff'][i])
        wr = pd.Series(SPEI_by_model[m]['WRunoff'][i])
        v_nr = pd.rolling_var(nr, varwindow)
        v_wr = pd.rolling_var(wr, varwindow)
        SPEIvar_by_model[m]['NRunoff'].append(v_nr)
        SPEIvar_by_model[m]['WRunoff'].append(v_wr)
        
colors_w = cm.get_cmap('Blues')(np.linspace(0.2, 1, num=len(modelnames)))
colors_n = cm.get_cmap('Wistia')(np.linspace(0.2, 1, num=len(modelnames)))
b=2
plt.figure('{} year running average trajectories, {} basin'.format(varwindow, basin_names[b]))
plt.axhline(y=0, color='Gainsboro', linewidth=2.0)
for k,m in enumerate(model_names):
    plt.plot(yrs, SPEIvar_by_model[m]['WRunoff'][b], label=m, color=colors_w[k], linewidth=2.0)
    plt.plot(yrs, SPEIvar_by_model[m]['NRunoff'][b], ls='-.', color=colors_n[k], linewidth=2.0)
plt.legend(loc='best')
plt.tight_layout()
plt.show()
