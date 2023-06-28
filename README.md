# IE-HYPE
This is a modified HYPE hydrological and water quality modelling tool (Version 5.5.1) that replaces 
the operator-splitting sequential calculation numerical scheme of the HYPE version 5.5.1 by the robust
fixed-step implicit Euler numerical scheme. It transforms the differential equations used by HYPE to 
solve percolation, lateral soil fluxe, and evapotranspiration into a state-space formulation and solves 
the fluxes generated from all soil layers simultaneously. This version of HYPE is referred to as IE-HYPE.
