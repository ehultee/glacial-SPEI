## Functions to support SPEI drought index analysis
## Code: EHU | SPEI data: SC
## 12 Sept 2019

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from datetime import date

## Constants associated with this analysis
yrs = np.linspace(1900, 2101, num=2412)
model_names = ['CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GISS-E2-R', 'INMCM4', 'MIROC-ESM', 'NorESM1-M'] # all models used in comparison
scenarios = ['Rcp4p5', 'Rcp8p5'] # climate scenarios
basin_names = ['INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE/VINCES',
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE']


def plot_basin_runmean(basin_id, permodel_dict, window_yrs=30, cmap_name='viridis', show_plot=True, save_plot=False):
    """Make a plot of running mean difference in SPEI for a given basin, comparing across models.
    Arguments:
        basin_id: integer, index of basin in the standard list "basin_names"
        permodel_dict: dictionary storing SPEI per model, with the structure dict[modelname]['diff'][basinname] = basin difference in SPEI for this model
        window_yrs: number of years to consider in running average.  Default 30
        cmap_name: name of matplotlib colormap from which to select line colors. Default 'viridis'
        show_plot: Boolean, whether to show the resulting plot.  Default True
        save_plot: Boolean, whether to save the plot in the working directory.  Default False
    """
    window_size = 12 * window_yrs # size of window given monthly data
    basin_runavg_bymodel = [np.convolve(permodel_dict[m]['diff'][basin_id], np.ones((window_size,))/window_size, mode='valid') for m in model_names] #compute running means
    colors = cm.get_cmap(cmap_name)(np.linspace(0, 1, num=len(model_names)))
    styles = ('-',':')
    plt.figure('{} year running average glacial effect, {} basin'.format(window_yrs, basin_names[basin_id]))
    for k,m in enumerate(model_names):
        plt.plot(yrs[(window_size/2):-(window_size/2 -1)], basin_runavg_bymodel[k], label=m, color=colors[k], ls=styles[np.mod(k, len(styles))], linewidth=2.0)
    plt.legend(loc='best')
    plt.tight_layout()
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(fname='{}yr_runmean-{}_basin-{}.png'.format(window_yrs, basin_names[basin_id], date.today()))

def glacial_meandiff(permodel_dict, years=(2070, 2100), return_range=True):
    """Calculate the difference in 30-yr mean SPEI with vs. without runoff
    Arguments:
        permodel_dict: dictionary storing SPEI per model, with the structure dict[modelname]['WRunoff'/'NRunoff'][basinname] = basin difference in SPEI for this model
        years: what years of scenario to examine. Default (2070, 2100)
        return_range: whether to return low/high multi-model range. Default True
    Outputs:
        median difference in mean(WRunoff)-mean(NRunoff)
        range in mean difference, expressed as 2xN array of (median-low, high-meadian) for errorbar plotting
    """
    idx_i, idx_f = 12*np.array((years[0]-1900, years[1]-1899)) #add 12 months to last year so that calculation goes for all 12 months
    bas_glac_meanmed = []
    bas_glac_lowmeans = []
    bas_glac_highmeans = []
    
    for i in range(len(basin_names)):
        bmeans_n = [np.nanmean(permodel_dict[m]['NRunoff'][i][idx_i:idx_f]) for m in model_names]
        bmeans_g = [np.nanmean(permodel_dict[m]['WRunoff'][i][idx_i:idx_f]) for m in model_names]
        basin_glacial_meanshift = np.array(bmeans_g) - np.array(bmeans_n)
        bas_glac_meanmed.append(np.nanmedian(basin_glacial_meanshift))
        bas_glac_lowmeans.append(np.nanmedian(basin_glacial_meanshift) - np.nanmin(basin_glacial_meanshift))
        bas_glac_highmeans.append(np.nanmax(basin_glacial_meanshift) - np.nanmedian(basin_glacial_meanshift))

    mean_spread = np.stack((bas_glac_lowmeans, bas_glac_highmeans))
    
    if return_range:
        return bas_glac_meanmed, mean_spread
    else:
        return bas_glac_meanmed
        
def glacial_vardiff(permodel_dict, years=(2070, 2100), return_range=True):
    """Calculate the difference in 30-yr SPEI variance with vs. without runoff
    Arguments:
        permodel_dict: dictionary storing SPEI per model, with the structure dict[modelname]['WRunoff'/'NRunoff'][basinname] = basin difference in SPEI for this model
        years: what years of scenario to examine. Default (2070, 2100)
        return_range: whether to return low/high multi-model range. Default True
    Outputs:
        median difference in var(WRunoff)-var(NRunoff)
        range in variance difference, expressed as 2xN array of (median-low, high-meadian) for errorbar plotting
    """
    idx_i, idx_f = 12*np.array((years[0]-1900, years[1]-1899)) #add 12 months to last year so that calculation goes for all 12 months
    bas_glac_varmed = []
    bas_glac_lowvars = []
    bas_glac_highvars = []
    
    for i in range(len(basin_names)):
        bvar_n = [np.nanvar(permodel_dict[m]['NRunoff'][i][idx_i:idx_f]) for m in model_names]
        bvar_g = [np.nanvar(permodel_dict[m]['WRunoff'][i][idx_i:idx_f]) for m in model_names]
        basin_glacial_varshift = np.array(bvar_g) - np.array(bvar_n)
        bas_glac_varmed.append(np.nanmedian(basin_glacial_varshift))
        bas_glac_lowvars.append(np.nanmedian(basin_glacial_varshift) - np.nanmin(basin_glacial_varshift))
        bas_glac_highvars.append(np.nanmax(basin_glacial_varshift) - np.nanmedian(basin_glacial_varshift))
        
    var_spread = np.stack((bas_glac_lowvars, bas_glac_highvars))

    if return_range:
        return bas_glac_varmed, var_spread
    else:
        return bas_glac_varmed



