## Read in, plot, and analyse SPEI data
## 22 Aug 2019
## Code: EHU | Data: SC

import numpy as np
import matplotlib.pyplot as plt

fpath = 'Documents/6. MIT/Drought buffering/SPEI_files/' 

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

## Read in two files (no runoff, with runoff for same model)
norunoff_fn = fpath+'NRunoff_{}_{}_{}.txt'.format(integration_times[0], modelnames[0], scenarios[1])
wrunoff_fn = fpath+'WRunoff_{}_{}_{}.txt'.format(integration_times[0], modelnames[0], scenarios[1])

norunoff_array = np.loadtxt(norunoff_fn)
wrunoff_array = np.loadtxt(wrunoff_fn)

yrs = np.linspace(1900, 2101, num=2412)

for k in range(len(basin_names))[2::5]:
    plt.figure('Basin {}'.format(basin_names[k]))
    plt.plot(yrs, norunoff_array[k], color='k')
    plt.plot(yrs, wrunoff_array[k], color='b')
    plt.show()

for k in range(len(basin_names))[2::5]:
    plt.figure('{} Basin glacier effect'.format(basin_names[k]))
    plt.plot(yrs, wrunoff_array[k]-norunoff_array[k], color='k')
    plt.show()
    
## Calculate the effect of including glaciers in each basin
glacierdiff = wrunoff_array - norunoff_array
basin_mean = [np.nanmean(glacierdiff[j]) for j in range(len(basin_names))]
normalized_basin_mean = [np.nanmean(glacierdiff[j])/np.nanmean(norunoff_array[j]) for j in range(len(basin_names))]

## 30-yr means