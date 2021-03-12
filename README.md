## Identifying the drivers of coral reef biomass in the Florida Keys to assess potential conservation actions

This repository includes data and analysis scripts for the following article:

Zuercher, R., R. Brumbaugh, K. Freeman, R. Layko, D. Kochan, A. Harborne. 2021. Identifying the drivers of coral reef biomass in the Florida Keys to assess potential management actions. *Submitted to Ecological Applications*

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
## Description:
Coral reef ecosystems and the services that we derive from them are in decline globally, underscoring the need for data-driven and cost-effective management approaches. In response, management actions and policies are being considered to reduce current anthropogenic threats, restore aspects of the ecosystem that have been lost, and increase system resilience. Weighing potential management interventions is significantly aided by a detailed understanding of the key factors controlling species’ abundances, the current status of species’ populations, and the likely outcomes of management initiatives. The aims of this study were to assess the relative importance of biophysical and anthropogenic factors on reef fishes, map fish biomass along the reef tract, and evaluate how various common management interventions might impact fish populations along the Florida reef tract. Using nearly three thousand fish surveys from coral reef habitats, we first modelled fishing impact using snapper and grouper species biomass as a proxy for the effects of fishing. Estimated fishing impact was then combined with 18 biophysical variables to model the current biomass of all reef fish species, the snapper-grouper species complex, grazing species, and species collected for aquaria. Generally, the biomass of fished species decreased with increasing fishing impact, and all species groups were also affected by biophysical factors and especially by benthic complexity. The statistical models allowed the production of high-resolution maps of current fish biomass, as well as maps of potential biomass under different management scenarios: a no-take marine reserve, moderate and extensive coral restoration, and the addition of artificial structure. The addition of structure was predicted to have the largest single impact on fish biomass of different species groups (23-73% increase compared to 5-22% following cessation of fishing). Additionally, our results point to beneficial synergies that occurred when combining habitat-based management with limits on fishing, with some combinations resulting in as much as a 90% increase in biomass. We highlight reef areas with large potential for fish biomass increases and suggest that conservation tools aimed at protecting structural complexity of reefs should remain an important part of the fishery management discussion, but that ways of increasing this complexity should be more widely considered.

--- 
## In the repository:
Four .csv files and one .R script are needed to replicate these analyses.

`RVC_impact.csv` -- data file containing each Reef Visual Census survey site used in the fishing impact model (rows), the snapper-grouper biomass for each site (column), and all explanatory variables considered for the fishing impact model (columns)

`ReefPoints_impact.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the fishing impact model

`RVC_biomass.csv` -- data file containing each Reef Visual Census survey site used in the fish biomass models (rows), biomass of all species groups (columns), and all explanatory variables considered for the biomass models, including fishing impact estimated by this project and extrapolated to these sites in ArcGIS (columns)

`ReefPoints_biomass.csv` -- data file contained every 1 ha reef pixel included in the project (rows) and all significant corresponding explanatory variables for the biomass models

`MOW_FL_3.2021.R` script -- runs all analyses and creates plots for Zuercher et al. 2021

The analyses can be replicated by changing the working directory in `MOW_FL_3.2021.R` to the location on your computer where you have stored the .R and .csv files. Additional analyses for this project were conducted in ArcGIS Pro. Spatial data layers are housed privately, but can be requested for the purpose of replication or for additional research. Questions about the code and requests for spatial data layers should be directed to Rachel Zuercher (rachel.zuercher@gmail.com).
