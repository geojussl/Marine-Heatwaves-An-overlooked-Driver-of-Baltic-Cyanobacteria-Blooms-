
# Marine-Heatwaves-An-overlooked-Driver-of-Baltic-Cyanobacteria-Blooms-
Detection of marine heatwaves and cyanobacterial bloom intensity in the Baltic Sea using satellite remote sensing (Sentinel-2, NOAA OISST) with Google Earth Engine and statistical analysis in R.

This repository contains the complete analysis workflow for investigating the
relationship between marine heatwaves (MHWs) and harmful cyanobacterial blooms (HABs)
in the Baltic Sea using satellite-based remote sensing data. 

This study is Part of the Course: UGM Umweltsystemanalyse MNF (HS/V) WiSe 2025/26 
EnviroChange 
Author: Justin Lingg-Laham and Jonah van den Bos

## Project Overview

Marine heatwaves have increased in frequency and intensity in many ocean regions,
including stratified marginal seas such as the Baltic Sea. At the same time,
harmful cyanobacterial blooms represent a recurring ecological phenomenon in this region.
This project investigates whether marine heatwaves act as a trigger or an
amplifying factor for cyanobacterial bloom intensity during summer months.

## Methods

The workflow consists of three main components:

### Part 1. Phycocyanin-based bloom detection (Google Earth Engine)
- Sentinel-2 surface reflectance data
- Cloud filtering using CloudMask (Wright et al. 2024)+
- Phycocyanin index calculation following Nesheli et al. (2024)
- Monthly aggregation of bloom intensity metrics (mean, median, 90th percentile)

### Part 2. Marine heatwave detection (Google Earth Engine + R)
- Daily sea surface temperature (SST) from NOAA OISST v2.1
- Summer-only analysis (June–August)
- 1991–2020 climatological baseline using a day-of-year approach (Smith et al. 2025).
- 90th percentile SST threshold calculation
- Export of daily SST and threshold time series

### Part 3. Event detection and statistical analysis (R)
- Marine heatwave event detection following Hobday et al. (2016)
- Calculation of event duration and intensity metrics
- Correlation analyses between MHW characteristics and bloom intensity
- Visualization 


## Repository Structure
