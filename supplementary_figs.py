#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 18:59:37 2020

@author: lizz

## Script generation of supplementary figures to avoid manual compositing
"""
import numpy as np
import matplotlib.pyplot as plt

## Figure S2: isolated glacial effect in 56 basins, alphabetised

## Alphabetise records
basin_names = ['INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE-VINCES',
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE']
basin_ids_raw = np.arange(len(basin_names))
basin_ids_sorted = sorted(zip(basin_names, basin_ids_raw), key=lambda x: x[0])


## Plot
page_chunks = np.arange(0, len(basin_names)+1, step=8)
for i in range(len(page_chunks)):
    chunk = range(page_chunks[i], page_chunks[i+1])
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, figsize=(6, 16))
    axs = np.asarray(((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8))).ravel()
    for j, k in enumerate(chunk):
        gSPEI.plot_basin_runmean(basin_id=basin_ids_sorted[k][1], permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=axs[j])
        axs[j].tick_params(axis='both', labelsize=10) #small labels appear large when scaled up to page size
    #ax1.legend(loc='best')
    plt.show()

for i in range(len(page_chunks)):
    chunk = range(page_chunks[i], page_chunks[i+1])
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, figsize=(6, 16))
    axs = np.asarray(((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8))).ravel()
    for j, k in enumerate(chunk):
        gSPEI.plot_basin_runvar(basin_id=basin_ids_sorted[k][1], permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=axs[j])
        axs[j].tick_params(axis='both', labelsize=10) #small labels appear large when scaled up to page size
    #ax1.legend(loc='best')
    plt.show()