# glacial-SPEI
Code to quantify the effect of including glacial runoff in Standardized Precipitation Evapotranspiration Index (SPEI).
[![DOI](https://zenodo.org/badge/205238173.svg)](https://zenodo.org/badge/latestdoi/205238173)

## Collaborators
- Sloan Coats (University of Hawaii)
- Jonathan Mackay (British Geological Survey)
- Lizz Ultee (Middlebury College)

## What's included
This repository includes MATLAB code, used to calculate the SPEI for each climate model, scenario, and basin; and Python code, used to analyse the SPEI data.  The Jupyter notebook "SPEI-analysis.ipynb" guides users through an analysis of SPEI with/without glacial runoff for any basin of choice.

## Organization
The main repo contains all analysis scripts.  In subfolders are:
- `archived`: older scripts used to develop parts of the analysis that may have been superseded
- `data`: all data required to reproduce our results
- `SPEI_computation`: code used to produce the basin-total SPEI data from output of GCMs and the glacier model
