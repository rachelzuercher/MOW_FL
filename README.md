## Identifying the drivers of coral reef biomass in the Florida Keys to assess potential conservation actions

This repository includes data and analysis scripts for the following article:

Zuercher, R., R. Brumbaugh, K. Freeman, R. Layko, D. Kochan, A. Harborne. 2021. Identifying the drivers of coral reef biomass in the Florida Keys to assess potential management actions. *In press at Aquatic Conservation*

---
## Collaborators:
This project is a collaboration between Alastair Harborne's Tropical Fish Ecology Lab (tfel) at Florida International University (FIU) and The Nature Conservancy (TNC). 

*Rachel Zuercher (Florida International University; National Socio-Environmental Synthesis Center)*    
*Alastair Harborne (Florida International University)*    
*Rob Brumbaugh (The Nature Conservancy)*    
*Kathleen Freeman (The Nature Conservancy)*    
*David Kochan (Florida International University)*     
*Rachel Layko (The Nature Conservancy)*    

**Contacts**: rachel.zuercher@gmail.com; aharborne@fiu.edu

---
## Abstract:
1. Coral reef ecosystems and the services that we derive from them are in decline globally, underscoring the need for data-driven management to reduce threats and restore populations. Comparing fishery management approaches is aided by a details understand of the key factors controlling species' abundances. 
2. The aims of this study were to assess the relative importance of biophysical factors compared with fishing impacts on the biomass of reef fishes on Florida's Coral Reef and to evaluate the potential effects of common management interventions on fish biomass. 
3. Fishing impact was estimated using a fishery-independent modelling approach and the biomass of the snapper-grouper complex as a proxy for the effects of fishing. Using a separate subset of data from underwtaer fish surveys, estimated fishing impact was then combined with 18 biophysical variables to model the current biomass of all reef fish species, the snapper-grouper complex, grazing species and species collected for aquaria. 
4. Models explained between 51 and 64% of the variance in fish biomass for the fish groups. The strongest predictor of biomass in the snapper-grouper complex was fishing impact (accounting for 25.2% of the explained variance), whereas reef complexity was the strongest predictor for all other groups. 
5. High-resolution maps were producted from the statistical models, including maps of current fish biomass and maps of potential biomass under several management scenarios: a no-take marine reserve, moderate and extensive coral restoration, and the addition of artificial benthic structure. Adding structure had the largest single impact on predicted fish biomass (23-72% increase from current estimated levels). However, beneficial synergies emerged when combining habitat-based management and fishing closures, with some combinations resulting in a reef-wide averaged 89% increase in biomass relative to current estimated levels. 
6. The results suggest that conservation strategies aimed at protecting and increasing structural complexity of reefs should be an important part of fishery management discussions. 

--- 
## In the repository:
Four .csv files and one .R script are needed to replicate these analyses.

`RVC_impact.csv` -- data file containing each Reef Visual Census survey site used in the fishing impact model (rows), the snapper-grouper biomass for each site (column), and all explanatory variables considered for the fishing impact model (columns)

`ReefPoints_impact.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the fishing impact model

`RVC_biomass.csv` -- data file containing each Reef Visual Census survey site used in the fish biomass models (rows), biomass of all species groups (columns), and all explanatory variables considered for the biomass models, including fishing impact estimated by this project and extrapolated to these sites in ArcGIS (columns)

`ReefPoints_biomass.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the biomass models

`MOW_FL_3.2021.R` script -- runs all analyses and creates plots for Zuercher et al. 2021

The analyses can be replicated by changing the working directory in `MOW_FL_3.2021.R` to the location on your computer where you have stored the .R and .csv files. Additional analyses for this project were conducted in ArcGIS Pro. Spatial data layers are housed privately, but can be requested for the purpose of replication or for additional research. Questions about the code and requests for spatial data layers should be directed to Rachel Zuercher (rachel.zuercher@gmail.com).
