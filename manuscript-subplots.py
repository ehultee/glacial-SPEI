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
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
gSPEI.plot_basin_runmean(basin_id=1, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax1)
gSPEI.plot_basin_runmean(basin_id=4, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax2)
gSPEI.plot_basin_runmean(basin_id=26, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax3)
gSPEI.plot_basin_runmean(basin_id=-7, permodel_dict=SPEI_by_model, show_labels=False, show_plot=False, save_plot=False, ax=ax4)
plt.show()

## Figure 2: isolated glacial effect in 4 basins


## Figure 3: isolated glacial effect on variance in 4 basins
