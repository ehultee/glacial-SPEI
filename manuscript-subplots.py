#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 21:42:36 2020

@author: lizz
"""
import matplotlib.pyplot as plt
import gSPEI as gSPEI

# make nice figures for publication with subplots

## Figure 1: 4 basins with/without runoff
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=(10,6))
gSPEI.plot_runmean_comparison(basin_id=4, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax1)
gSPEI.plot_runmean_comparison(basin_id=1, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax2)
gSPEI.plot_runmean_comparison(basin_id=26, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax3)
gSPEI.plot_runmean_comparison(basin_id=-7, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax4)
#ax2. legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.show()

# Separate legend for compositing
fig2, ax5 = plt.subplots(1) 
gSPEI.plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, which='WRunoff', cmap_name='Blues', show_labels=False, show_plot=False, save_plot=False, ax=ax5)
gSPEI.plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, which='NRunoff', cmap_name='Wistia', show_labels=False, show_plot=False, save_plot=False, ax=ax5)
ax5.legend(loc='upper left', bbox_to_anchor=(1.05, 1), ncol=2)
plt.tight_layout()
plt.show()


## Figure 2: isolated glacial effect in 4 basins
fig3, ((ax6, ax7), (ax8, ax9)) = plt.subplots(2, 2, sharex=True, figsize=(10,6))
gSPEI.plot_basin_runmean(basin_id=4, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax6)
gSPEI.plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax7)
gSPEI.plot_basin_runmean(basin_id=26, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax8)
gSPEI.plot_basin_runmean(basin_id=-7, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax9)
#ax2. legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.show()

## Separate legend for Figures 2/3
fig4, ax10 = plt.subplots(1)
gSPEI.plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax10)
# plt.gcf().set_visible(False)
ax10.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.show()


## Figure 3: isolated glacial effect on variance in 4 basins
fig5, ((ax11, ax12), (ax13, ax14)) = plt.subplots(2, 2, sharex=True, figsize=(10,6))
gSPEI.plot_basin_runvar(basin_id=4, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax11)
gSPEI.plot_basin_runvar(basin_id=1, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax12)
gSPEI.plot_basin_runvar(basin_id=26, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax13)
gSPEI.plot_basin_runvar(basin_id=-7, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax14)
#ax2. legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.show()