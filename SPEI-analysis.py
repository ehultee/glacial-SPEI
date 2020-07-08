## Read in, plot, and analyse SPEI data
## 22 Aug 2019
## Code: EHU | Data: SC

import numpy as np
import matplotlib.pyplot as plt
from gSPEI import *

fpath = './data/SPEI_Files/'

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

yrs = np.linspace(1900, 2101, num=2412)


## Compare effect across models - read in all to dict
SPEI_by_model = {m: {} for m in modelnames} # create dictionary indexed by model name
for m in modelnames:
    norunoff_f_m = fpath+'NRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    wrunoff_f_m = fpath+'WRunoff_{}_{}_{}.txt'.format(integration_times[3], m, scenarios[0])
    SPEI_by_model[m]['NRunoff'] = np.loadtxt(norunoff_f_m)
    SPEI_by_model[m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
    SPEI_by_model[m]['diff'] = SPEI_by_model[m]['WRunoff'] - SPEI_by_model[m]['NRunoff']

## 30-yr running means and variance
for i in (1, 4, 26, -7): #plot the basins shown in manuscript main text
    plot_runmean_comparison(basin_id=i, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=True)
    plot_basin_runmean(basin_id=i, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=True)
    plot_basin_runvar(basin_id=i, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=True)

# for i in range(len(basin_names)): #plot all
#     plot_basin_runmean(basin_id=i, permodel_dict=SPEI_by_model, save_plot=True, show_plot=False)
#     plot_basin_runvar(basin_id=i, permodel_dict=SPEI_by_model, save_plot=True, show_plot=False)


# ## Compare series with different SPEI integration times
# ## --adding explicit comparison 15 Apr 2020 in response to reviewer comments
# SPEI_by_itime = {t: {} for t in integration_times} # create dictionary indexed by integration time
# for t in integration_times:
#     SPEI_by_itime[t] = {m: {} for m in modelnames} #nest dictionary by model name
#     for m in modelnames:
#         norunoff_f_m = fpath+'NRunoff_{}_{}_{}.txt'.format(t, m, scenarios[1])
#         wrunoff_f_m = fpath+'WRunoff_{}_{}_{}.txt'.format(t, m, scenarios[1])
#         SPEI_by_itime[t][m]['NRunoff'] = np.loadtxt(norunoff_f_m)
#         SPEI_by_itime[t][m]['WRunoff'] = np.loadtxt(wrunoff_f_m)
#         SPEI_by_itime[t][m]['diff'] = SPEI_by_itime[t][m]['WRunoff'] - SPEI_by_itime[t][m]['NRunoff']

# for t in integration_times:
#     plot_basin_runmean(basin_id=4, permodel_dict=SPEI_by_itime[t], show_plot=False, save_plot=True, output_tag='itime_{}mo'.format(t))
#     plot_basin_runmean(basin_id=26, permodel_dict=SPEI_by_itime[t], show_plot=False, save_plot=True, output_tag='itime_{}mo'.format(t))
#     plot_basin_runmean(basin_id=-7, permodel_dict=SPEI_by_itime[t], show_plot=False, save_plot=True, output_tag='itime_{}mo'.format(t))

