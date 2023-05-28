# --- Default variables

# --- Loads the connectivity of the network and builds a connectivity file for each cell, based on the population selected
loadFullConn    = False;    saveFullConn    = False;    updateCellTags  = False;    
loadCellTags    = True

# --- View side comment:
loadCellType    = True              # loads information about cell type classification (based on M1 model analysis)
loadWeightNorm  = True              # (default: True) Load weight normalization values to de-normalize the connection weights during histogram building

generate_post_spks      = True      # loads spike times and builds the dictionary of postsynaptic spikes

# --- View side comment:
# ignoreSpikeTimes:     ( True )    runs analysis on the time window, ignoring the timing of each spike (adds a 'fake' spike at the end time and considers all the pre spikes in the interval [start, end] of the window)
#                       ( False )   uses the spike times to build the spike histogram
ignoreSpikeTimes        = True      # ignores spike times and build a histogram of inputs/s
removeSilentCells       = False     # ignores cells with no spikes when building the spike time histogram

# --- Histogram of Cell Connectivity:
# generate_conns:           builds the dictionary of connections
# generate_histogram_data:  builds the histogram of spikes within the time window analyzed
# save_histogram_data:      saves the histogram of spikes within the time window analyzed
generate_conns,generate_histogram_data,save_histogram_data = False,False,False

# --- Analyzing the QUIET vs MOVEMENT states
compareStates           = False
# --- Analysis for each time range
runDataAnalysis         = True

# --- Plotting for each time range
plotSingleFigs  = False
plot_SPTH_traces        = plotSingleFigs
plot_SPTH_bar           = plotSingleFigs
plot_SPTH_boxplot       = plotSingleFigs
plot_SPTH_violin        = plotSingleFigs
plot_spikes_scatter     = plotSingleFigs

# --- Dimensionality reduction analysis
createDataFrame             = True
plotPCA,plotUMAP,plotKMeans = False,True,True

# --- Analysis combining different time ranges
runPostAnalysis         = True
plotMergedBar,plotMergedBar_cellType,plotMergedBar_kMeans = True,True,True

showPlots,savePlots     = True,False

# --- Printing load info
print('\n\n##############################################')
print('           Default variables loaded')
print('##############################################')