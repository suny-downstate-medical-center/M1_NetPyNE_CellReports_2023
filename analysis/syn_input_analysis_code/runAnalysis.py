#! /Users/joao/opt/anaconda3/bin/python

#####################################################################################################################################################################################################################################################################

# --- Importing libraries
import sys
import math
import glob
import json
import pickle
from tkinter import S
from turtle import width
import netpyne
import numpy as np
import pandas as pd
from cProfile import label
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from collections import Counter

#####################################################################################################################################################################################################################################################################

# --- Snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
import sys
from pydoc import source_synopsis
from matplotlib  import  pyplot  as plt
from numpy import size

print("Matplotlib backend (default): %s" %plt.get_backend())
modules = []
for module in sys.modules:
    if module.startswith('matplotlib'): modules.append(module)
for module in modules:
    sys.modules.pop(module)
import matplotlib
matplotlib.use("MacOSX")
from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (dynamic): %s" %plt.get_backend())

#####################################################################################################################################################################################################################################################################

# --- Importing classes and methods to analyze the data
from AnalyzeData    import BarPlot          as Bp
from AnalyzeData    import ColorMap         as Cmp
from AnalyzeData    import Connectivity     as Cnct
from AnalyzeData    import PlotFigures      as Ptf
from AnalyzeData    import CellConnectivity as CellConn
        
# --- Script Configuration

histogram_data_folder = 'histogram_data'
histogram_figure_folder = 'histogram_figures'

# DO a README FILE

#####################################################################################################################################################################################################################################################################

# --- Population variables  |   Run one at a time
# loadPop=['PT5B']
# loadPop=['IT5B']
loadPop=['PV5B']

# --- Populations to plot   |   Can run multiple pops when plotting
plotPop = loadPop

#####################################################################################################################################################################################################################################################################
# --------------------------------------- #
# --- Loading Default Config settings --- #
# --------------------------------------- #

# --- Loading Default Config settings
from defaultCfg import *

# --------------------------------------- #
# --- Connectivity and Histogram Data --- #
#   Ideally, these variables should be (True) only when updating the connectivity and recalculating the histogram data
#   Otherwise, the connectivity and histograms should be stored and loaded for plotting
# --------------------------------------- #
# --- Updates the files that store the connectivity for each cell 
updateConnectivity = False
if updateConnectivity:  loadFullConn,saveFullConn,updateCellTags = True,True,True
# --- Updates the files that store the histogram data for each population
updateHistDict = False
if updateHistDict:      generate_conns,generate_histogram_data,save_histogram_data = True,True,True    

# --- Changing plots generated based on the population of interest
if loadPop[0]=='PT5B':
    # --- Dimensionality reduction analysis
    createDataFrame = True
    plotPCA,plotUMAP,plotKMeans = False,True,True

    # --- Analysis combining different time ranges
    runPostAnalysis = True
    plotMergedBar,plotMergedBar_cellType,plotMergedBar_kMeans = True,True,True

    # --- Showing and Saving plots
    showPlots,savePlots = True,True

elif loadPop[0]=='IT5B':
    # --- Dimensionality reduction analysis
    createDataFrame             = True
    plotPCA,plotUMAP,plotKMeans = False,False,False

    # --- Analysis combining different time ranges
    runPostAnalysis         = True
    plotMergedBar,plotMergedBar_cellType,plotMergedBar_kMeans = True,True,False

    # --- Showing and Saving plots
    showPlots,savePlots     = True,True

elif loadPop[0]=='PV5B':
    # --- Dimensionality reduction analysis
    createDataFrame             = True
    plotPCA,plotUMAP,plotKMeans = False,False,False

    # --- Analysis combining different time ranges
    runPostAnalysis         = True
    plotMergedBar,plotMergedBar_cellType,plotMergedBar_kMeans = True,False,False

    # --- Showing and Saving plots
    showPlots,savePlots     = True,True

#####################################################################################################################################################################################################################################################################

# --- Sorted list of all populations (internal + long range) from: conns['net']['pops'].keys()
all_pops        = [ 'IT2', 'SOM2', 'PV2', 'IT4', 'IT5A', 'SOM5A', 'PV5A', 'PT5B', 'IT5B', 'SOM5B', 'PV5B', 'CT6', 'IT6', 'SOM6', 'PV6', 'S1', 'S2', 'M2', 'OC', 'TVL', 'TPO', 'cM1']

local_pops      = ['IT2', 'SOM2', 'PV2', 'IT4', 'IT5A', 'SOM5A', 'PV5A', 'PT5B', 'IT5B', 'SOM5B', 'PV5B', 'CT6', 'IT6', 'SOM6', 'PV6']
long_range_pops = ['S1', 'S2', 'M2', 'OC', 'TVL', 'TPO', 'cM1']

#####################################################################################################################################################################################################################################################################

if ignoreSpikeTimes:    hist_method = 'windowAnalysis'
else:                   hist_method = 'spkAligned'

#####################################################################################################################################################################################################################################################################

# --- Creating a Colormap for the plots
map_colors=Cmp.colormap(all_pops,c_map='gist_rainbow')

#####################################################################################################################################################################################################################################################################

dataFolder  = '../data/'
# fileLabel   = 'm1_paper_2019_v56_batch23_v56_batch23_0_0_0_0'
fileLabel   = 'model_outputData'
fileFormat  = '.pkl'
dataFile    = dataFolder+fileLabel+fileFormat
with open(dataFile, 'rb') as fileObj: data = pickle.load(fileObj)

netpyne.sim.allSimData = data['simData']

#####################################################################################################################################################################################################################################################################

# --- Connectivity dataset (output from NetPyne simulation)
conn_full_label         = 'model_connectivity'
conn_full_fileFormat    = '.pkl'
conn_full_fileName      = dataFolder+conn_full_label+conn_full_fileFormat

# --- Folder where .pkl files with individual cell connectivity are stored
connFolder = '../data/conn_info/'
popFolder = '../data/pop_info/'

#####################################################################################################################################################################################################################################################################

# --- Cell tags dataset
# cellTags_label      = 'cell_tags'
cellTags_label      = 'model_cellTags'
cellTags_fileFormat = '.pkl'
cellTags_fileName   = dataFolder+cellTags_label+cellTags_fileFormat

#####################################################################################################################################################################################################################################################################

# --- Weight normalization data
weightNorm_folderPath   = '../weightNorm/'
weightNorm_fileName     = 'weightNorm_dict'
weightNorm_fileFormat   = '.pkl'
weightNorm_filePath     = weightNorm_folderPath+weightNorm_fileName+weightNorm_fileFormat

#####################################################################################################################################################################################################################################################################

# --- Cell type dataset (enhanced vs suppressed cells)
if loadCellType:
    # cellTypeFileLabel   = 'v56_batch19_0_0_0_statDataAll_boxplot'
    cellTypeFileLabel   = 'model_cellType'
    cellTypeFileFormat  = '.json'
    cellTypeFile = dataFolder+cellTypeFileLabel+cellTypeFileFormat
    with open(cellTypeFile, 'r') as cellTypeObj: cellType = json.load(cellTypeObj)
    
    # --- L5B cells
    cellType_dict={'enhanced': cellType['includeAll'][0][4], 'suppressed': cellType['includeAll'][0][5],}

    for ct in cellType_dict.keys():
        overlap=0
        ct_list=list(cellType_dict.keys()); ct_list.remove(ct)
        for gid in cellType_dict[ct]:
            if gid in cellType_dict[ct_list[0]]: overlap+=1
        print('overlap: ',ct,' and ',ct_list[0],': ',overlap)
        
#####################################################################################################################################################################################################################################################################


# --- Loads full Connectivity dataset
if loadFullConn:
    # --- Loads full Connectivity dataset
    with open(conn_full_fileName, 'rb') as conn_full_fileObj: conns = pickle.load(conn_full_fileObj)

    # --- Generates a library of .pkl files with the presynaptic connectivity for each individual postsynaptic cell
    target_pops=['PT5B','PV5B', 'IT5B']
    if saveFullConn: CellConn.generateCellConnectivity(conns, connFolder, target_pops=target_pops)

    # --- Generates a .pkl file with the GIDs of the cells in each pop 
    CellConn.getPopGIDs(conns, popFolder, target_pops=target_pops)
    
    # --- Update dictionary of cell tags (3d position, etc) (without connectivity)
    if updateCellTags: CellConn.generate_cell_tags(conns,cellTags_fileName)
        
#####################################################################################################################################################################################################################################################################

# --- Loads the cell connectivity .pkl files | filename: pop_gid.pkl
pop_gids_fileName = popFolder+loadPop[0]+'_GIDs.pkl'
with open(pop_gids_fileName, 'rb') as pop_gids_fileObj: pop_gids = pickle.load(pop_gids_fileObj); pop_gids.sort()

#####################################################################################################################################################################################################################################################################

# --- Loads only the cell tags, including individual cell information (3d position, etc) (without connectivity)
if loadCellTags:
    # --- Loads the cell connectivity .pkl files | filename: pop_gid.pkl
    with open(cellTags_fileName, 'rb') as cellTags_fileObj: cell_tags_dict = pickle.load(cellTags_fileObj)

#####################################################################################################################################################################################################################################################################

# --- Time range for histogram calculation
debug_code=False
if debug_code:  timeRanges = [[1000,2000]]
else:           timeRanges = [[2000,5000],[6000,9000]]

timeRange_names=[]
for timeRange in timeRanges:
    # --- Defines the time range to iterate over
    if timeRange == ['all']:    timeRange_name = '_all'
    else:                       timeRange_name = '_'+str(timeRange[0])+'_'+str(timeRange[1])
    timeRange_names.append(timeRange_name)

# --- Selecting populations to plot in the figure
select_plot_pops = None   # plots all populations

# --- Load weightNorm dictionary
if loadWeightNorm: weightNorm_dict = Cnct.load_weightNorm(weightNorm_filePath)

# --- Loading presynaptic cell connectivity information
loadSingleCell = False               # loads only a single or a few cells to debug
if loadSingleCell:  all_post_cell_gids = [5133,5134,5533]   # PT5B debug cells
else:               all_post_cell_gids = pop_gids           # Full list of cells in the postsynaptic population

# ---------------------------------------------------------------------------------------------------------------------------------- #
# --- Generate connections and loads spike times

# # === POSTSYNAPTIC SPIKES
# # --- Loads all the spike data for the postsynaptic population and stores in a dictionary format
if generate_post_spks: 
    if ignoreSpikeTimes:
        print('\n\n##############################################')
        print('             IGNORING SPIKE TIMES             ')
        print('##############################################')
        spk_dict_post={}
        for pop in loadPop:
            spk_dict_post.update({pop:{}})
            for post_cell_gid in all_post_cell_gids:
                tFinal=[]
                for timeRange in timeRanges:
                    tFinal.append(timeRange[1]-0.00001)
                spk_dict_post[pop].update({post_cell_gid:tFinal})
        # --- Properties of the time window used to create the spike histogram
        time_slot   = timeRange[1]-timeRange[0]
        n_slots     = 1
        time_bins   = list((i+1)*time_slot for i in range(n_slots))

        # --- generates the original spike dictionary to use as a label to color the cells in the post analysis
        true_spk_dict_post, silent_cells = Cnct.generatePostSpkHist(loadPop, all_post_cell_gids)

    else:
        spk_dict_post, silent_cells = Cnct.generatePostSpkHist(loadPop, all_post_cell_gids)
        # --- Properties of the time window used to create the spike histogram
        time_slot   = 5
        n_slots     = 10
        time_bins   = list((i+1)*time_slot for i in range(n_slots))

if generate_conns:
    # === POSTSYNAPTIC CONNECTIVITY
    # --- Load all connectivity from presynaptic to postsynaptic cell beforehand
    conn_dict_post = Cnct.generatePostSpkConn(loadPop, all_post_cell_gids, connFolder)

    # === PRESYNAPTIC SPIKES
    # --- Obtaining all the GIDs of presynaptic cells to load spikes beforehand    
    spk_dict_pre = Cnct.generatePreSpkHist(conn_dict_post, all_post_cell_gids)
        
# ---------------------------------------------------------------------------------------------------------------------------------- #
# --- Generates spike histogram
if generate_histogram_data:
    # --- Generate plots for multiple timeRanges
    for timeRange_ind, timeRange in enumerate(timeRanges):
        # --- Defines the name of the time range for saving
        timeRange_name = timeRange_names[timeRange_ind]
        # --- Creating empty dictionary to store gids of cells that didn't spike
        post_pop_noSpike={}
        # --- Creating dictionary to store final histogram data
        spike_hist_dict={}; silent_spike_hist_dict={}; quiet_spike_hist_dict={}
        # --------------------------------------------------------------------------------- #
        # --- Iterating over different presynaptic pops
        for pop in spk_dict_post.keys():
            # --------------------------------------------------------------------------------- #
            # --- Expanding histogram dictionary for each postsynaptic population
            spike_hist_dict.update(         {pop:{}})
            silent_spike_hist_dict.update(  {pop:{}})
            quiet_spike_hist_dict.update(   {pop:{}})
            post_pop_noSpike.update(        {pop:[]})
            # --- Creating empty dictionaries to store variables
            cell_dict={}
            post_pop_exc=[]; post_pop_inh=[]
            
            # --------------------------------------------------------------------------------- #
            # --- Identifying the cell in the pop+timeRange with highest number of spikes
            cell_spk_count={}
            most_spikes=[]
            for post_cell_gid in all_post_cell_gids:
                spk_hist=np.histogram(spk_dict_post[pop][post_cell_gid],bins=1,range=timeRange)
                spk_num=spk_hist[0][0]
                # --- Dictionary organized by GID (easier to access)
                cell_spk_count.update({post_cell_gid:spk_num})
            # --- Counter for spikes
            count_spikes = Counter(cell_spk_count)
            # --- Finding 3 highest values
            most_spikes = count_spikes.most_common(3)
            print("Cells with most spikes in pop ", pop,': ',most_spikes)
            # --- Creating a list with spike times of the cell with highest number of spikes
            if not removeSilentCells:
                transposed_spike_times=[]
                for spkt in spk_dict_post[pop][most_spikes[0][0]]:
                    if spkt>timeRange[0] and spkt<=timeRange[1]: transposed_spike_times.append(spkt)
            
            # --------------------------------------------------------------------------------- #
            # --- Iterates through the list of cells in the postsynaptic population
            for post_cell_gid_ind,post_cell_gid in enumerate(all_post_cell_gids):
                # --- Expanding histogram dictionary for each cell in the postsynaptic population
                spike_hist_dict[pop].update({post_cell_gid:{}})
                # --- Status update
                if (post_cell_gid)%100==0:
                    print('\n\n##############################################')
                    print('             Processing cell %s Histogram            '%post_cell_gid)
                    # print('             timerange: ', timeRange)
                    print('             %s cells left'%((len(all_post_cell_gids)-post_cell_gid_ind)))
                    print('##############################################')
                # --------------------------------------------------------------------------------- #
                # --- Presynaptic Connectivity information
                conn_data           = conn_dict_post[pop][post_cell_gid]

                # --- GIDs of the presynaptic cells that project to the given postsynaptic cell
                connected_pre_cell_gids_ = []
                for conn in conn_data:
                    connected_pre_cell_gids_.append(conn[0])
                connected_pre_cell_gids=list(set(connected_pre_cell_gids_))

                # --- Generate Connectivity Dictionary
                conn_dict = Cnct.generate_conn_dict(conn_data=conn_data,weightNorm_dict=weightNorm_dict[pop])
                
                # --- Dictionary only with the spike times of the connected presynaptic cells
                connected_spk_dict_pre={}
                for pre_cell_gid in spk_dict_pre.keys():
                    if pre_cell_gid in connected_pre_cell_gids: connected_spk_dict_pre.update({pre_cell_gid:spk_dict_pre[pre_cell_gid]})
                
                # --------------------------------------------------------------------------------- #
                # --- Spike time information
                post_cell_spk_times = spk_dict_post[pop][post_cell_gid]
                # --- Slices the array of spike times to account for the timeRange to be analyzed
                if len(timeRange)>1:
                    new_spkts=[]
                    for new_spkt in post_cell_spk_times:
                        if new_spkt>timeRange[0] and new_spkt<=timeRange[1]: new_spkts.append(new_spkt)
                    del post_cell_spk_times
                    # --- Spike times within the timeRange analyzed
                    timeRange_post_cell_spk_times = new_spkts

                # --------------------------------------------------------------------------------- #
                # --- Regular spiking cell 
                if len(timeRange_post_cell_spk_times)>0:    post_spikes = timeRange_post_cell_spk_times # 'spiking cell'
                # --- Checking cells that didnt spike
                else:
                    if removeSilentCells:   continue
                    else:                   post_spikes = transposed_spike_times # 'silent and quiet cells'
                # --- Returns the spike histogram by [pre_pop][mech_type][pre_gid]
                post_cell_spike_hist_dict = Cnct.generateSpikeHistogram(conn_dict, time_bins, timeRange_post_cell_spk_times, connected_spk_dict_pre)
                # --- Stores spike histogram for each postsynaptic cell
                spike_hist_dict[pop][post_cell_gid].update(post_cell_spike_hist_dict)

            # --- Saves spike histogram dictionary for the population
            if generate_histogram_data and save_histogram_data:
                print('Saving Histogram data')
                if debug_code:
                    print('Data stored in DEBUG folder')
                    spike_hist_dict_filename='../data/'+histogram_data_folder+'/debug/debug.pkl'
                    with open(spike_hist_dict_filename, 'wb') as f:pickle.dump(spike_hist_dict, f)
                elif ignoreSpikeTimes:
                    spike_hist_dict_filename='../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'_fullWindow.pkl'
                    with open(spike_hist_dict_filename, 'wb') as f:pickle.dump(spike_hist_dict, f)
                else:
                    spike_hist_dict_filename='../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'.pkl'
                    with open(spike_hist_dict_filename, 'wb') as f:pickle.dump(spike_hist_dict, f)
                    
#####################################################################################################################################################################################################################################################################

# --- Index in the time_bins at which the histogram data is capped
if ignoreSpikeTimes:
    max_time_index = 0 ; ax_lim=(-5,110) # sum from 0 to +3000 ms
    max_time = time_bins[0]
    # --- Normalizing to inputs per second in the full window analysis
    perSecond = True; timeScaling = 3 # SECONDS (size of the analysis window)
else:
    max_time_index = 1 ; ax_lim=(-0.5,16) # sum from 0 to 10 ms
    max_time = time_bins[max_time_index]
    perSecond = False; timeScaling = None

#####################################################################################################################################################################################################################################################################

# --- Analyzing the QUIET vs MOVEMENT states
if compareStates:

    # === POSTSYNAPTIC SPIKES
    print('\n\n##############################################')
    print('       Re-loading post cell spike times       ')
    print('##############################################')

    # --- timeRanges for each condition
    load_timeRanges = timeRanges

    # --- Loads all the spike data for the postsynaptic population and stores in a dictionary format
    spk_dict_post={}
    for pop in plotPop:
        spk_dict_post.update({pop:{}})
        post_cell_spk_info=netpyne.analysis.tools.getSpktSpkid(cellGids=all_post_cell_gids)
        spk_dict = Cnct.getSpkDict(post_cell_spk_info,all_cell_GIDs=all_post_cell_gids)
        spk_dict_post[pop].update(spk_dict)
    
        # --- Spike histogram for the different conditions
        compare_states_dict={}
        save_spks={}
        highest_spiking_cells={'q':[],'m':[]}
        save_spks.update({'q':{},'m':{}})
        for timeRange_ind, timeRange in enumerate(load_timeRanges):
            if timeRange_ind==0:    timeRange_key='q'   # quiet
            else:                   timeRange_key='m'   # movement

            compare_states_dict.update({timeRange_key:{'on':{},'off':{}}})

            for post_cell_gid in spk_dict_post[pop].keys():
                spk_hist=np.histogram(spk_dict_post[pop][post_cell_gid],bins=1,range=timeRange)
                spk_num=spk_hist[0][0]
                
                if spk_num>0:   cell_state='on'
                else:           cell_state='off'

                # --- Update dictionary with cell firing
                compare_states_dict[timeRange_key][cell_state].update({post_cell_gid:spk_num})
                # --- Dictionary organized by GID (easier to access)
                save_spks[timeRange_key].update({post_cell_gid:spk_num})

            k = Counter(save_spks[timeRange_key])
            # Finding 3 highest values
            high = k.most_common(3)
            for i in high:
                highest_spiking_cells[timeRange_key].append(i[0])
                print(i[0]," :",i[1]," ")            

        # --- Ratios of cells in each condition
        for network_state in compare_states_dict.keys():
            for cell_state in compare_states_dict[network_state].keys():
                cell_count = len(compare_states_dict[network_state][cell_state].keys())
                print(network_state, cell_state, cell_count)
                print('\tratio: ', cell_count/len(all_post_cell_gids))

        # --- Ratios of cells in each condition
        cell_state={ 'silent': [], 'activated': [], 'silenced': [], 'active': [],}
        for post_cell_gid in all_post_cell_gids:
            if (post_cell_gid in compare_states_dict['q']['off'].keys()) and (post_cell_gid in compare_states_dict['m']['off'].keys()): cell_state['silent'].append(    post_cell_gid)
            if (post_cell_gid in compare_states_dict['q']['off'].keys()) and (post_cell_gid in compare_states_dict['m']['on'].keys()):  cell_state['activated'].append( post_cell_gid)
            if (post_cell_gid in compare_states_dict['q']['on'].keys()) and  (post_cell_gid in compare_states_dict['m']['off'].keys()): cell_state['silenced'].append(  post_cell_gid)
            if (post_cell_gid in compare_states_dict['q']['on'].keys()) and  (post_cell_gid in compare_states_dict['m']['on'].keys()):  cell_state['active'].append(    post_cell_gid)
        plt.figure()
        for ind,key in enumerate(cell_state.keys()):
            cell_num=len(cell_state[key])
            plt.bar(ind,cell_num)
            print(cell_num)

        cell_spikes = np.array([[len(cell_state['silent']),len(cell_state['activated'])],[len(cell_state['silenced']),len(cell_state['active'])],])
        cell_spikes_label = np.array([['Silent','Activated'],['Silenced','Active'],])
        fig, ax = plt.subplots()
        im = ax.imshow(cell_spikes,cmap='viridis')

        # --- Adding labels on grid plot
        for i in range(2):
            for j in range(2):
                text = ax.text(j, i, cell_spikes_label[i, j]+'\n'+str(cell_spikes[i, j]),ha="center", va="center", color="w",size=20,)

        # --- Colormap based on GID
        all_post_cells_colormap = Cmp.colormap(all_post_cell_gids,c_map='jet')
        
        # --- Cell Spatial properties (cell tags)
        for post_cell_ind, post_cell_gid in enumerate(all_post_cell_gids):
            cell_position_x = cell_tags_dict[post_cell_gid]['xnorm']
            cell_position_y = cell_tags_dict[post_cell_gid]['ynorm']

        plt.figure()
        # your input data:
        plt.subplot(221)
        befores_activated=[]
        afters_activated=[]
        for post_cell_ind, post_cell_gid in enumerate(cell_state['activated']):
            before = save_spks['q'][post_cell_gid]
            after  = save_spks['m'][post_cell_gid]

            plt.plot( [1,2], [before, after], c='royalblue',linewidth=0.1)
            plt.plot(1,before,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)
            plt.plot(2,after,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)

            befores_activated.append(before)
            afters_activated.append(after)
        
        mean_befores_activated = np.mean(befores_activated)
        mean_afters_activated  = np.mean(afters_activated)

        plt.plot(1,mean_befores_activated,marker='o',color='k')
        plt.plot(2,mean_afters_activated,marker='o',color='k')
        plt.plot( [1,2], [np.mean(befores_activated), np.mean(afters_activated)], c='k')
        plt.title('Activated: %.0f '%len(cell_state['activated'])) 

        plt.subplot(223)
        befores_silenced=[]
        afters_silenced=[]
        for post_cell_ind, post_cell_gid in enumerate(cell_state['silenced']):
            before = save_spks['q'][post_cell_gid]
            after  = save_spks['m'][post_cell_gid]
            
            plt.plot( [1,2], [before, after], c='r',linewidth=0.1)
            plt.plot(1,before,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)
            plt.plot(2,after,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)
        
            befores_silenced.append(before)
            afters_silenced.append(after)

        mean_befores_silenced = np.mean(befores_silenced)
        mean_afters_silenced  = np.mean(afters_silenced)
        
        plt.plot(1,mean_befores_silenced,marker='o',color='k')
        plt.plot(2,mean_afters_silenced,marker='o',color='k')
        plt.plot( [1,2], [np.mean(befores_silenced), np.mean(afters_silenced)], c='k')
        plt.title('Silenced: %.0f '%len(cell_state['silenced']))

        plt.subplot(122)
        Decreased=0
        Increased=0
        Same=0
        befores_increased=[]
        afters_increased=[]
        befores_decreased=[]
        afters_decreased=[]
        spk_percent_threshold = 30
        for post_cell_ind, post_cell_gid in enumerate(cell_state['active']):
            before = save_spks['q'][post_cell_gid]
            after  = save_spks['m'][post_cell_gid]

            if before>after*(1+(spk_percent_threshold/100)):   # >X% difference
                c='r'
                Decreased+=1
                befores_decreased.append(before)
                afters_decreased.append(after)
            elif after>before*(1+(spk_percent_threshold/100)): # >X% difference
                c='royalblue'
                Increased+=1
                befores_increased.append(before)
                afters_increased.append(after)
            else:
                c='k'
                Same+=1

            plt.plot( [1,2], [before, after], c=c,linewidth=0.1)
            plt.plot(1,before,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)
            plt.plot(2,after,marker='.',color=all_post_cells_colormap[post_cell_gid],alpha=0.5)
        
        mean_befores_increased = np.mean(befores_increased)
        mean_afters_increased  = np.mean(afters_increased)

        plt.plot(1,mean_befores_increased,marker='o',color='k')
        plt.plot(2,mean_afters_increased,marker='o',color='k')
        plt.plot( [1,2], [np.mean(befores_increased), np.mean(afters_increased)], c='royalblue')
        
        mean_befores_decreased = np.mean(befores_decreased)
        mean_afters_decreased = np.mean(afters_decreased)

        plt.plot(1,mean_befores_decreased,marker='o',color='k')
        plt.plot(2,mean_afters_decreased,marker='o',color='k')
        plt.plot( [1,2], [np.mean(befores_decreased), np.mean(afters_decreased)], c='r')
        plt.title('Increased: %.0f | Decreased: %.0f | Same: %.0f'%(Increased,Decreased,Same))

#####################################################################################################################################################################################################################################################################
if runDataAnalysis:
    # --- Types of synaptic mechanisms
    mech_types=['exc','inh']
    # --- Store the output of KMeans
    store_cluster_dictionary={}; store_cluster_embeding={}; store_cluster_dataframe={}; store_kmeans={}
    # --- Generate plots for multiple timeRanges
    for timeRange_ind, timeRange in enumerate(timeRanges):
        # --- Defines the name of the time range for saving
        timeRange_name = timeRange_names[timeRange_ind]
        # --- Name for saving the data and plots
        if   timeRange==['all']:                        network_state = 'All states'
        elif timeRange[0]<5000 or timeRange[0]>=9000:   network_state = 'QUIET'
        elif timeRange[0]>=5000 and timeRange[0]<9000:  network_state = 'MOVEMENT'
        else:                                           network_state = 'UNKNOWN'
        # --- Iterates over the postsynaptic pops
        for pop in plotPop:
            print('\n\n##############################################')
            print('       Plotting data for %s          '%pop)
            print('##############################################')
            try:
                # --- Loading spike histogram dictionary for the population
                if debug_code:          spike_hist_dict_filename='../data/'+histogram_data_folder+'/debug/debug.pkl'
                elif ignoreSpikeTimes:  spike_hist_dict_filename='../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'_fullWindow.pkl'
                else:                   spike_hist_dict_filename='../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'.pkl'

                with open(spike_hist_dict_filename, 'rb') as spike_hist_fileObj: spike_hist_dict = pickle.load(spike_hist_fileObj)
                print('Succesfully loaded file: ', spike_hist_dict_filename)
            except:
                print(pop, ' histogram data missing')
                continue

            pop_spike_hist_dict = spike_hist_dict[pop]
            plot_spike_hist_dict, ordered_pre_pops, valid_post_cell_gids = Ptf.formatData(pop_spike_hist_dict=pop_spike_hist_dict, all_pops=all_pops)

            # --- Dict with spike times for each cell
            pop_spk_dict_post = spk_dict_post[pop]

            #####################################################################################################################################################################################################################################################################
            # --- Plot SPTH traces
            if plot_SPTH_traces:
                select_plot_pops=None
                if select_plot_pops is not None:    pre_pops_name='selected_pre_pops'
                else:                               pre_pops_name='all_pre_pops'
                Ptf.plotSPTHtraces(plot_spike_hist_dict, ordered_pre_pops, time_bins, map_colors, divide_plots=False, select_plot_pops=select_plot_pops)
                if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/1_spth_plot/'+pop+'_spth'+'_histogram'+timeRange_name+'_ms_'+pre_pops_name, dpi=1000)

            #####################################################################################################################################################################################################################################################################
            # --- Spike Time Histogram plot
            if plot_SPTH_bar:
                Ptf.barPlot(pop_spike_hist_dict,ordered_pre_pops,valid_post_cell_gids,max_time_index,all_pops,)
                if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/2_bar_plot/'+pop+'_spth'+'_bar'+timeRange_name+'_ms', dpi=1000)

            #####################################################################################################################################################################################################################################################################
            # --- Bar plot
            if plot_SPTH_boxplot:
                Ptf.boxPlot(    pop_spike_hist_dict, ordered_pre_pops, valid_post_cell_gids, max_time_index,
                                select_pops=[['IT2', 'IT4', 'PT5B', 'TPO', 'TVL'],['PV5A', 'PV5B']])
                if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/3_box_plot/'+pop+'_spth'+'_box'+timeRange_name+'ms', dpi=1000)

            #####################################################################################################################################################################################################################################################################
            # --- Violin plot
            if plot_SPTH_violin:
                Ptf.violinPlot( pop_spike_hist_dict, ordered_pre_pops, valid_post_cell_gids, max_time_index,
                                select_pops=[['IT2', 'IT4', 'PT5B', 'TPO', 'TVL'],['PV5A', 'PV5B']])
                if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/4_violin_plot/'+pop+'_spth'+'_violin'+timeRange_name+'ms', dpi=1000)

            #####################################################################################################################################################################################################################################################################
            # --- Scatter plot
            if plot_spikes_scatter:
                select_max_val = 162
                Ptf.scatterPlot(pop_spk_dict_post, timeRange, all_post_cell_gids, cell_tags_dict, select_max_val = select_max_val, select_colormap='Reds', use_x_position=True)
                if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/5_scatter_plot/'+pop+'_spth'+'_scatter'+timeRange_name+'ms', dpi=1000)

            #####################################################################################################################################################################################################################################################################
            # --- Dataframe for multivariate analysis
            if createDataFrame:

                featuredPops = ordered_pre_pops
                target_data = 'ynorm'
                # target_data = 'spk'
                # if ignoreSpikeTimes:
                #     target_data = ('true_spk',true_spk_dict_post[pop])

                DataFrame_dict, colormap_boudaries = Ptf.formatMultivariateData(    pop_spk_dict_post,pop_spike_hist_dict,
                                                                                    featuredPops,valid_post_cell_gids,
                                                                                    timeRange,max_time_index,cell_tags_dict,
                                                                                    target_data = target_data,
                                                                                    isolate_mech=None
                                                                                    )

                spk_cap=75 # spikes
                if ignoreSpikeTimes:
                    target_data2 = ('true_spk',true_spk_dict_post[pop])
                    #               name               spks                    upper lim
                    target_data4 = ('capped_spk_range',true_spk_dict_post[pop],[0,spk_cap])
                else:
                    # target_data2 = 'spk'
                    target_data2 = ('fixed_spk_range',pop_spk_dict_post,[0,161])
                    #               name               spks                    upper lim
                    target_data4 = ('capped_spk_range',pop_spk_dict_post,      [0,spk_cap])
                    
                DataFrame_dict2, colormap_boudaries2 = Ptf.formatMultivariateData(  pop_spk_dict_post,pop_spike_hist_dict,
                                                                                    featuredPops,valid_post_cell_gids,
                                                                                    timeRange,max_time_index,cell_tags_dict,
                                                                                    target_data = target_data2,
                                                                                    isolate_mech=None)
                
                DataFrame_dict4, colormap_boudaries4 = Ptf.formatMultivariateData(  pop_spk_dict_post,pop_spike_hist_dict,
                                                                                    featuredPops,valid_post_cell_gids,
                                                                                    timeRange,max_time_index,cell_tags_dict,
                                                                                    target_data = target_data4,
                                                                                    isolate_mech=None)
                if ignoreSpikeTimes:
                    fixedSpkRange=False
                    if fixedSpkRange:
                        target_data3 = ('fixed_spk_range',true_spk_dict_post[pop],[0,161])
                        DataFrame_dict3, colormap_boudaries3 = Ptf.formatMultivariateData(  pop_spk_dict_post,pop_spike_hist_dict,
                                                                                                featuredPops,valid_post_cell_gids,
                                                                                                timeRange,max_time_index,cell_tags_dict,
                                                                                                target_data = target_data3,
                                                                                                isolate_mech=None)

                # --- PCA
                if plotPCA: 
                    # --- PCA - Soma cortical depth
                    pca = Ptf.plotPCA(DataFrame_dict,pop_spk_dict_post,n_components=2)
                    plt.title('PCA - Soma cortical depth')
                    plt.rcParams.update({'font.size': 20})
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/6_pca/'+pop+'_spth'+'_pca'+timeRange_name+'ms_'+target_data, dpi=1000)
                    
                    # --- PCA - Spike number
                    pca2 = Ptf.plotPCA(DataFrame_dict2,pop_spk_dict_post,n_components=2)
                    plt.title('PCA - Spike number')
                    plt.rcParams.update({'font.size': 20})
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/6_pca/'+pop+'_spth'+'_pca'+timeRange_name+'ms_'+target_data2[0], dpi=1000)

                    # --- PCA - Capped max value
                    pca4 = Ptf.plotPCA(DataFrame_dict4,pop_spk_dict_post,n_components=2)
                    plt.title('PCA - Capped max value')
                    plt.rcParams.update({'font.size': 20})
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/6_pca/'+pop+'_spth'+'_pca'+timeRange_name+'ms_'+target_data4[0]+'_'+str(target_data4[2][-1])+'_spikes', dpi=1000)

                # --- UMAP
                n_neighbors=200
                convert_to_Hz=True
                if convert_to_Hz:   spk_units = ' Hz'
                else:               spk_units = ' spikes'

                if plotUMAP: 
                    # --- UMAP - Soma cortical depth
                    embedding, df = Ptf.plotUMAP(DataFrame_dict,n_neighbors=n_neighbors,cellType_dict=cellType_dict,color_criteria=target_data,colormap_boudaries=colormap_boudaries)
                    # plt.title('UMAP - Soma cortical depth')
                    plt.rcParams.update({'font.size': 20})
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/7_umap/'+pop+'_spth'+'_umap'+timeRange_name+'ms_'+target_data, dpi=1000)
                    
                    # --- UMAP - Spike number
                    if convert_to_Hz:colormap_boudaries2[1]=(colormap_boudaries2[1]*1000/max_time)
                    embedding2, df2 = Ptf.plotUMAP(DataFrame_dict2,n_neighbors=n_neighbors,cellType_dict=cellType_dict,color_criteria=target_data2,colormap_boudaries=colormap_boudaries2)
                    plt.title('UMAP - Spike count in '+spk_units)
                    plt.rcParams.update({'font.size': 20})
                    # if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/7_umap/'+pop+'_spth'+'_umap'+timeRange_name+'ms_'+target_data2[0], dpi=1000)

                    # --- UMAP - Color spikes based on a fixed value
                    if ignoreSpikeTimes and fixedSpkRange:
                        if convert_to_Hz:colormap_boudaries3[1]=(colormap_boudaries3[1]*1000/max_time)
                        embedding3, df3 = Ptf.plotUMAP(DataFrame_dict3,n_neighbors=n_neighbors,cellType_dict=cellType_dict,color_criteria=target_data3,colormap_boudaries=colormap_boudaries3)
                        plt.title('UMAP - Fixed spike number ('+str(target_data3[2][-1])+spk_units+')')
                        plt.rcParams.update({'font.size': 20})
                        # if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/7_umap/'+pop+'_spth'+'_umap'+timeRange_name+'ms_'+target_data3[0], dpi=1000)
                    
                    # --- UMAP - Capped max value
                    if convert_to_Hz:
                        colormap_boudaries4[1]=(colormap_boudaries4[1]*1000/max_time)
                        fig_units='Hz'
                    else:
                        fig_units='spikes'
                    embedding4, df4 = Ptf.plotUMAP(DataFrame_dict4,n_neighbors=n_neighbors,cellType_dict=cellType_dict,color_criteria=target_data4,colormap_boudaries=colormap_boudaries4)
                    # plt.title('UMAP - Capped max value('+str(round(target_data4[2][-1]))+spk_units+')')
                    plt.rcParams.update({'font.size': 20})
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/7_umap/'+pop+'_spth'+'_umap'+timeRange_name+'ms_'+target_data4[0]+'_'+str(round(target_data4[2][-1]))+'_'+fig_units, dpi=1000)

                if plotKMeans:
                    if network_state == 'MOVEMENT': n_clusters = 2
                    else:                           n_clusters = 2

                    kmeans_dataset      = Ptf.applyKMeans(embedding,n_clusters = n_clusters)
                    cluster_dictionary  = Ptf.plotKMeans(kmeans_dataset, embedding, df, customColors=['darkblue','darkred'])
                    
                    if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/8_kmeans/'+pop+'_spth'+'_kmeans'+timeRange_name+'ms', dpi=1000)

                    # --- Store the output of KMeans (cell GIDs in each cluster)
                    store_cluster_dictionary.update({timeRange_ind:cluster_dictionary})
                    store_cluster_embeding.update({timeRange_ind:embedding})
                    store_cluster_dataframe.update({timeRange_ind:df})
                    store_kmeans.update({timeRange_ind:kmeans_dataset})

                    showKmeansSubplots=False
                    if showKmeansSubplots:
                        Ptf.barplotKMeans2(          cluster_dictionary, pop_spike_hist_dict, featuredPops, max_time_index, divide_plots=False)
                        if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/8_kmeans/'+pop+'_spth'+'_kmeansBarPlot'+timeRange_name+'ms', dpi=1000)

                        boxplot_kMeans=Ptf.boxplotKMeans(  cluster_dictionary, pop_spike_hist_dict, featuredPops, max_time_index, divide_plots=False, )
                        if savePlots: plt.savefig('../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/8_kmeans/'+pop+'_spth'+'_kmeansBoxPlot'+timeRange_name+'ms', dpi=1000)

#####################################################################################################################################################################################################################################################################
# --- Run analysis that merges QUIET and MOVEMENT state into a single plot
if runPostAnalysis:
    from AnalyzeData import PostAnalysis as Pa
    selectPops=False
    if selectPops:  popsFlag='selectPops'
    else:           popsFlag='allPops'

    if selectPops:  select_pre_pops = ['IT2', 'IT4', 'IT5A', 'PT5B', 'PV5A', 'PV5B', 'TVL', 'TPO']
    else:           select_pre_pops = None

    mech_types=['exc','inh']
    plot_ordered_pre_pops=['IT2', 'IT4', 'IT5A', 'SOM5A', 'PV5A', 'PT5B', 'IT5B', 'SOM5B', 'PV5B', 'CT6', 'IT6', 'S1', 'S2', 'M2', 'OC', 'TVL', 'TPO', 'cM1']

    # --- Iterates over the postsynaptic pops
    for pop in plotPop:
        # --- Generate plots for multiple timeRanges
        for timeRange_ind, timeRange in enumerate(timeRanges):
            # --- Defines the name of the time range for saving
            timeRange_name = timeRange_names[timeRange_ind]
            # --- Name for saving the data and plots
            if   timeRange==[2000,5000]:    network_state = 'Q' # Quiet
            elif timeRange==[6000,9000]:    network_state = 'M' # Movement
            else:                           network_state = 'U' # Unknown
            try:
                # --- Loading spike histogram dictionary for the population
                if debug_code:          spike_hist_dict_filename = '../data/'+histogram_data_folder+'/debug/debug.pkl'
                elif ignoreSpikeTimes:  spike_hist_dict_filename = '../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'_fullWindow.pkl'
                else:                   spike_hist_dict_filename = '../data/'+histogram_data_folder+'/spk_histogram_data/'+pop+'_spike_histogram'+timeRange_name+'.pkl'
                    
                with open(spike_hist_dict_filename, 'rb') as spike_hist_fileObj: spike_hist_dict = pickle.load(spike_hist_fileObj)
                print('Succesfully loaded file: ', spike_hist_dict_filename)
            except:
                print(pop, ' histogram data missing')
                continue
            
            pop_spike_hist_dict = spike_hist_dict[pop]
            if network_state == 'Q':
                pop_spike_hist_dict_Q = pop_spike_hist_dict
                plot_spike_hist_dict_Q, ordered_pre_pops_Q, valid_post_cell_gids_Q = Ptf.formatData(pop_spike_hist_dict=pop_spike_hist_dict, all_pops=all_pops)
            elif network_state == 'M':
                pop_spike_hist_dict_M = pop_spike_hist_dict
                plot_spike_hist_dict_M, ordered_pre_pops_M, valid_post_cell_gids_M = Ptf.formatData(pop_spike_hist_dict=pop_spike_hist_dict, all_pops=all_pops)
            else:
                print('Error')
                sys.exit()

        verticalPlot=False
        if verticalPlot: plotOrientation='V' 
        else:            plotOrientation='H'

        # --- Merged bar plot
        if plotMergedBar:
            # --- Merged bar plots with all cells
            spk_hist_Q  = pop_spike_hist_dict_Q
            post_GIDs_Q = valid_post_cell_gids_Q
            
            spk_hist_M  = pop_spike_hist_dict_M
            post_GIDs_M = valid_post_cell_gids_M

            figName     = pop+'_spth'+'_1_bar_merged_'+plotOrientation+timeRange_name+'ms_'+popsFlag
            figData     = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+figName+'.json'
            figFullName = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/9_post_analysis/'+figName
            if pop == 'PV5B':   ax_lim = [0,3000]
            else:               ax_lim = [0,2500]
            Pa.mergedBarPlot(spk_hist_Q,spk_hist_M,plot_ordered_pre_pops,post_GIDs_Q,post_GIDs_M,max_time_index,long_range_pops,select_pre_pops=select_pre_pops,ax_lim=ax_lim,verticalPlot=False,perSecond=perSecond,timeScaling=timeScaling,export_filename=figData)
            # plt.title(pop+' | '+'all cells'+' | '+str(len(post_GIDs_Q))+' total cells')
            plt.title(pop)
            plt.rcParams.update({'font.size': 20})
            plt.tight_layout()
            if savePlots: plt.savefig(figFullName, dpi=1000)

        if plotMergedBar_cellType:
            # --- Enhanced or Suppressed cells
            for cell_type in cellType_dict.keys():
                # --- lists to hold gid values
                post_GIDs_Q=[]; post_GIDs_M=[]
                # --- valid Quiet gids
                for cell_gid_Q in valid_post_cell_gids_Q:   
                    if cell_gid_Q in cellType_dict[cell_type]:  post_GIDs_Q.append(cell_gid_Q)
                # --- valid Move gids
                for cell_gid_M in valid_post_cell_gids_M:   
                    if cell_gid_M in cellType_dict[cell_type]:  post_GIDs_M.append(cell_gid_M)

                print('number of '+cell_type+' cells: Q ',len(post_GIDs_Q),' | M ',len(post_GIDs_M))
                spk_hist_Q  = pop_spike_hist_dict_Q
                spk_hist_M  = pop_spike_hist_dict_M

                figName     = pop+'_spth'+'_2_bar_cellType_'+plotOrientation+timeRange_name+'ms_'+popsFlag+'_'+cell_type+'_cells'
                figData     = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+figName+'.json'
                figFullName = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/9_post_analysis/'+figName
                print(pop+' | '+cell_type+' | '+str(len(post_GIDs_Q))+' total cells')

                Pa.mergedBarPlot(spk_hist_Q,spk_hist_M,plot_ordered_pre_pops,post_GIDs_Q,post_GIDs_M,max_time_index,long_range_pops,select_pre_pops=select_pre_pops,ax_lim=[0,2500],verticalPlot=False,perSecond=perSecond,timeScaling=timeScaling,showLabels=False,export_filename=figData)
                if cell_type == 'enhanced': title_flag=' enhanced'
                else:                       title_flag=' suppressed'
                plt.title(pop+title_flag)
                plt.tight_layout()
                if savePlots: plt.savefig(figFullName, dpi=1000)
                
        # --- Merged bar plots with groups of cells
        if plotMergedBar_kMeans:
            for timeRange_ind, timeRange in enumerate(timeRanges):
                post_GIDs_Cluster0 = store_cluster_dictionary[timeRange_ind][0]
                post_GIDs_Cluster1 = store_cluster_dictionary[timeRange_ind][1]
                spk_hist_Cluster0  = pop_spike_hist_dict
                spk_hist_Cluster1  = pop_spike_hist_dict

                print('Debug: ',        timeRange_ind, timeRange, network_state)
                print('Q vs Default: ', pop_spike_hist_dict_Q==pop_spike_hist_dict)
                print('M vs Default: ', pop_spike_hist_dict_M==pop_spike_hist_dict)
                print('Q vs M: ',       pop_spike_hist_dict_Q==pop_spike_hist_dict_M)

                if timeRange_ind == 0:  hist_dict = pop_spike_hist_dict_Q
                else:                   hist_dict = pop_spike_hist_dict_M

                print('Q vs New dict: ', pop_spike_hist_dict_Q==hist_dict)
                print('M vs New dict: ', pop_spike_hist_dict_M==hist_dict)

                figName     = pop+'_spth'+'_3_bar_kmeans_'+plotOrientation+popsFlag+'_timeRange'+timeRange_names[timeRange_ind]+'_ms'
                figData     = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+figName+'.json'
                figFullName = '../figs/'+histogram_figure_folder+'/'+hist_method+'/'+pop+'/hist_window_'+str(max_time)+'_ms'+'/9_post_analysis/'+figName
                print('State: '+timeRange_names[timeRange_ind]+' cells: Cluster 1 ',len(post_GIDs_Cluster0),' | Cluster 2 ',len(post_GIDs_Cluster1))
                Pa.mergedBarPlot(   hist_dict,hist_dict, # same histogram, but the Clusters are filtered through the GIDs
                                    plot_ordered_pre_pops,
                                    post_GIDs_Cluster0,post_GIDs_Cluster1,
                                    max_time_index,long_range_pops,
                                    select_pre_pops=select_pre_pops,ax_lim=[0,2500],
                                    verticalPlot=False,perSecond=perSecond,
                                    timeScaling=timeScaling,
                                    states = ['Cluster1','Cluster2'],
                                    state_colors=['dodgerblue','orangered'],
                                    export_filename=figData)
                if timeRange_ind == 0:  title_flag=' Quiet'
                else:                   title_flag=' Movement'
                plt.title(pop+title_flag)
                plt.tight_layout()
                if savePlots: plt.savefig(figFullName, dpi=1000)
                
        plot_postAnalysis_KMeans=False
        if plot_postAnalysis_KMeans:
            for timeRange_ind, timeRange in enumerate(store_cluster_dictionary.keys()):
                useQuietClusters=True
                if useQuietClusters:    customLabels=store_kmeans[0].labels_ # --- Uses labels from the clustering performed in the quiet condition
                else:                   customLabels=None
                Ptf.plotKMeans(store_kmeans[timeRange_ind], store_cluster_embeding[timeRange_ind], store_cluster_dataframe[timeRange_ind],customLabels=customLabels)

            plotScatter_kMeans=True
            if plotScatter_kMeans:
                for timeRange_ind, timeRange in enumerate(store_cluster_dictionary.keys()):
                    plt.figure(figsize=(5,10))
                    plt.title(timeRange_names[timeRange_ind]) 
                    for quiet_cell in  store_cluster_dictionary[timeRange_ind][0]:
                        cell_position_x = cell_tags_dict[quiet_cell]['xnorm']
                        cell_position_y = cell_tags_dict[quiet_cell]['ynorm']
                        marker_color = 'darkblue'
                        edge_color = None
                        alpha = 1
                        plt.plot(cell_position_x,cell_position_y,marker='o',color=marker_color,markeredgecolor=edge_color,alpha=alpha)
                    for move_cell in   store_cluster_dictionary[timeRange_ind][1]:
                        cell_position_x = cell_tags_dict[move_cell]['xnorm']
                        cell_position_y = cell_tags_dict[move_cell]['ynorm']
                        marker_color = 'darkred'
                        edge_color = None
                        alpha = 1
                        plt.plot(cell_position_x,cell_position_y,marker='o',color=marker_color,markeredgecolor=edge_color,alpha=alpha)
                    plt.gca().invert_yaxis()
                    plt.tight_layout()

#####################################################################################################################################################################################################################################################################

if showPlots: plt.show()

#####################################################################################################################################################################################################################################################################
