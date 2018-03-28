# microarray-Rscripts
R scripts to analyse microarray time-series data  
+  **vbssmmatrix.R** creates a matrix to load for VBSSM modelling, using a custom dataset (not supplied here) and a genelist (not supplied here). It takes a tab-delimited genelist of a single column with header, overlaps the genelist with the timeseries and creates a matrix ready to use for VBSSM modelling.
It also creates a CATMAprobe_to_TAIR mapping file to use as node attributes with Cytoscape.  
+  **Combined_timeseries_phyper_profiles.R** hosts a function with multiple parameters (see the script for a full explanation). It takes an Arabidopsis code, a microarray CATMA probe ID or a list (or lists) of Arabidopsis AGI codes, and calculates the probability of overlaps between the sets of genelists and the main dataset (not supplied here).
It also creates expression profile plots, heatmaps and similarity matrices.
