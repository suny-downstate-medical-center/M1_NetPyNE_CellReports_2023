#! /Users/joao/opt/anaconda3/bin/python

import colorsys
from matplotlib import pyplot as plt
import matplotlib as mpl
import netpyne
import numpy as np
import pickle
import math
import json

class ColorMap():
    def colormap(data_list,c_map):
        steprange = range(len(data_list))
        color_map = plt.get_cmap(c_map)
        print("Current map color: %s" %color_map.name)
        map_colors = {}
        for j,k in enumerate(data_list):
            map_colors.update({k:color_map(1.*j/float(len(steprange)))})
        return map_colors

class CellConnectivity():
    def getPopGIDs(conns, popFolder, target_pops):
        # pop_GID_dict={}
        # --- Selects one population from the list of populations specified
        for target_pop in target_pops:
            print('target_pop: ',target_pop)
            # --- Gets the gids of the cells in that population
            target_gids = conns['net']['pops'][target_pop]['cellGids']
            pop_GID_filename=popFolder+target_pop+'_GIDs.pkl'
            # pop_GID_dict.update({target_pop:target_gids})
        
            with open(pop_GID_filename, 'wb') as f: pickle.dump(target_gids, f)

        # return pop_GID_dict

    def generateCellConnectivity(conns, connFolder, target_pops):
        conn_dict={}
        allPops = conns['net']['pops'].keys()
        # --- Selects one population from the list of populations specified
        for target_pop in target_pops:
            print('target_pop: ',target_pop)
            # --- Gets the gids of the cells in that population
            target_gids = conns['net']['pops'][target_pop]['cellGids']
            # --- Creates an empty dictionary for that population to store information
            conn_dict.update({target_pop:{}})
            # --- Loops through the gids of the postsynaptic cells
            for cell_gid in target_gids:
                # --- Generates the filename for the .pkl file with the cell connectivity loaded in the updateConns step
                cell_filename=connFolder+target_pop+'_'+str(cell_gid)+'.pkl'
                print('cell_gid: ',cell_gid,'\tcells remaining: ',target_gids[-1]-cell_gid)
                
                conn_dict[target_pop].update({cell_gid:{}})
                save_conns=[]

                # --- Gets information about each conn from presynaptic cells in an array format [pre_gid, post_sec, post_loc, mech, weight, delay]
                for conn in conns['net']['cells'][cell_gid]['conns']:
                    if type(conn) is list:
                        pre_gid = conn[0]
                    elif type(conn) is dict:
                        pre_gid = conn['preGid']
                        pre_sec = conn['sec']
                        pre_loc = conn['loc']
                        pre_mec = conn['synMech']
                        pre_wei = conn['weight']
                        pre_del = conn['delay']
                        del conn
                        conn = [pre_gid,pre_sec,pre_loc,pre_mec,pre_wei,pre_del]
                    else:
                        return 'invalid conn data format'

                    # pre_gid_list.append(pre_gid)
                    
                    # --- Gets the keys of the pops
                    for pop in allPops:
                        # --- Checks from which pop the presynaptic cell comes from
                        if pre_gid in conns['net']['pops'][pop]['cellGids']:
                            conn.append(pop)
                            save_conns.append(conn)
                            
                # --- Saves the cell connectivity .pkl files | filename: pop_gid.pkl
                with open(cell_filename, 'wb') as f: pickle.dump(save_conns, f)
    
    def generate_cell_tags(conns,fileName):
        # --- List of dictionaries      (    format:     [ {index1:{cell_tags1}}, {index2:{cell_tags2}}, ... ]     )
        try:    cell_tags = [{ind:conns['net']['cells'][ind].tags} for ind in range(len(conns['net']['cells']))]
        except: cell_tags = [{ind:conns['net']['cells'][ind]['tags']} for ind in range(len(conns['net']['cells']))]
        # --- Formatting to a better dictionary structure
        cell_tags_dict={}
        for cell_tag in cell_tags:
            cell_tags_dict.update(cell_tag)

        # --- Saves the cell connectivity .pkl files | filename: cell_tags.pkl
        with open(fileName, 'wb') as f: pickle.dump(cell_tags_dict, f)

class Connectivity():
    # --- Creates a list with all the GIDs of the presynaptic cells
    def get_pre_gids(conn_data):
        pre_gids_=[]
        for conn in conn_data:
            pre_gids_.append(conn[0])
        pre_gids = list(set(pre_gids_))
        return pre_gids

    # --- Creates a dictionary that organizes the spikes of the cells using the GIDs as keys
    # --- obs: SpktSpkid = the output of netpyne.analysis.tools.getSpktSpkid(cellGids=list_of_cell_gids)
    def getSpkDict(SpktSpkid,all_cell_GIDs=None):
        # --- Converting spike time and spike id to lists
        spkts   = SpktSpkid[1]
        spkgids = [int(num) for num in SpktSpkid[2]]
        # --- Creates a set of GIDs to generate the dictionary keys
        if all_cell_GIDs:   # in this case, cells that did not spike are still allocated as an empty list
            # --- Removes duplicates
            spkgids_set=list(set(all_cell_GIDs))
            # --- Sort list of GIDs
            spkgids_set.sort()
        else:
            # --- Removes duplicates
            spkgids_set=list(set(spkgids))
            # --- Sort list of GIDs
            spkgids_set.sort()
        # --- Create Spike dictionary
        spk_dict={}
        for spk_gid in spkgids_set:
            spk_dict.update({spk_gid:[]})
        # --- Creates a dictionary where (key=GID, value=spike times)
        for i in range(len(spkts)):
            spk_dict[spkgids[i]].append(spkts[i])
        return spk_dict
    
    def load_weightNorm(weightNorm_filePath):
        with open(weightNorm_filePath, 'rb') as weightNorm_fileObj: weightNorm_dict_ = pickle.load(weightNorm_fileObj)

        # --- Fixing the dictionary to match pop names
        weightNorm_dict={}
        for cell_type in weightNorm_dict_.keys():
            if  (cell_type=='PT5B_full') or \
                (cell_type=='PT5B_reduced'):    pop_name = 'PT5B'
            elif cell_type.startswith('IT5B'):  pop_name = 'IT5B'
            elif cell_type.startswith('PV4'):   pop_name = 'PV4'
            elif cell_type.startswith('PV5A'):  pop_name = 'PV5A'
            elif cell_type.startswith('SOM5A'): pop_name = 'SOM5A'
            elif cell_type.startswith('PV5B'):  pop_name = 'PV5B'
            elif cell_type.startswith('CT6'):   pop_name = 'CT6'
            elif cell_type.startswith('IT6'):   pop_name = 'IT6'
            # --- A1 pops
            elif cell_type.startswith('ITS4'):  pop_name = 'ITS4'
            elif cell_type.startswith('ITP4'):  pop_name = 'ITP4'
            elif cell_type.startswith('IT5A'):  pop_name = 'IT5A'
            else:
                print(cell_type)
                print('Error loading weight normalization')
                import sys; sys.exit()
            
            weightNorm_dict.update({pop_name:weightNorm_dict_[cell_type]})
        print('Weight norm loaded')
        return weightNorm_dict
    
    def generate_conn_dict(conn_data,simplify_mechs=False,weightNorm_dict=None):

        # --- Obtaining a list of presynaptic populations
        pre_cell_pops_=[]
        for conn_ind,conn in enumerate(conn_data):
            pre_cell_pops_.append(conn[-1])
        pre_cell_pops=list(set(pre_cell_pops_))

        # --- Allocating presynaptic pops as dictionary keys
        conn_dict={}
        temp_conn_dict={} # temp dict to store unsorted GIDs
        for pre_cell_pop in pre_cell_pops:
            conn_dict.update({pre_cell_pop:{'exc':{},'inh':{}}})
            temp_conn_dict.update({pre_cell_pop:{'exc':{},'inh':{}}}) # temp dict to store unsorted GIDs

        # --- Splitting the dictionary in 'exc' and 'inh' connections to facilitate processing later
        for conn_ind,conn in enumerate(conn_data):
            pre_cell_gid = conn[0]
            pre_cell_pop = conn[-1]

            if ('GABA' in conn[3]) and (pre_cell_gid not in conn_dict[pre_cell_pop]['inh'].keys()):         mech_type='inh'
            elif ('GABA' not in conn[3]) and (pre_cell_gid not in conn_dict[pre_cell_pop]['exc'].keys()):   mech_type='exc'
            else:
                print('Unknown Mechanism: ', conn[3])
                continue

            # --- Assigns unsorted GIDs to temporary dictionary
            temp_conn_dict[pre_cell_pop][mech_type].update({pre_cell_gid:{'conns':[],'weights':[],'delays':[]}})
        
        # --- Reassigns sorted GIDs to conn dictionary
        for pre_cell_pop in pre_cell_pops:
            for mech_type in ['exc','inh']:
                key_list = list(temp_conn_dict[pre_cell_pop][mech_type].keys())
                key_list.sort()
                for key in key_list:
                    conn_dict[pre_cell_pop][mech_type].update({key:temp_conn_dict[pre_cell_pop][mech_type][key]})

        # --- Processing connectivity data
        for conn_ind,conn in enumerate(conn_data):
            pre_cell_gid    = conn[0]   # gid
            pre_cell_sec    = conn[1]   # sec
            pre_cell_loc    = conn[2]   # loc
            pre_cell_mech   = conn[3]   # mech
            # --- Type of mechanism (exc or inh)
            if 'GABA' in pre_cell_mech: mech_type = 'inh'
            else:                       mech_type = 'exc'

            # --- Weight
            if weightNorm_dict: pre_cell_weight = conn[4]/weightNorm_dict[conn[1]]  # De-normalizing conn weight
            else:               pre_cell_weight = conn[4]                           # Using normalized conn weight 
            
            pre_cell_delay      = conn[5]   # delay
            pre_cell_pop        = conn[-1]  # pop

            # --- String with all conn parameters
            conn_name = pre_cell_sec+'|'+str(pre_cell_loc)+'|'+pre_cell_mech+'|'+str(pre_cell_weight)+'|'+str(pre_cell_delay)
            
            # --- Storing the detailed connection and separate weight and delay values
            conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['conns'].append(conn_name)
            conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['weights'].append(pre_cell_weight)
            conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['delays'].append(pre_cell_delay)

        # --- Adding weights and Averaging delays of the connections from the same presynaptic cell
        for pre_cell_pop in conn_dict.keys():
            for mech_type in conn_dict[pre_cell_pop].keys():
                for pre_cell_gid in conn_dict[pre_cell_pop][mech_type].keys():
                    conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['sum_weight']=sum(conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['weights'])
                    conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['avg_delay']=np.mean(conn_dict[pre_cell_pop][mech_type][pre_cell_gid]['delays'])
        
        return conn_dict

    def generatePostSpkHist(loadPop, all_post_cell_gids):
        # === POSTSYNAPTIC SPIKES
        # --- Loads all the spike data for the postsynaptic population and stores in a dictionary format
        spk_dict_post_={}
        spk_dict_post={}
        for pop in loadPop:
            
            spk_dict_post_.update({pop:{}})
            spk_dict_post.update({pop:{}})

            post_cell_spk_info=netpyne.analysis.tools.getSpktSpkid(cellGids=all_post_cell_gids)
            spk_dict = Connectivity.getSpkDict(post_cell_spk_info)
            spk_dict_post_[pop].update(spk_dict)

            # --- Adding empty lists for the cells that did not spike throughout the whole simulation
            #     P.s.: these cells didnt fire in the whole sim / different of the 'silent' group, which includes these cells + others that fired only outside of the QUIET/MOVE periods, but mght have fired a e.g. 100 ms, which is out of those ranges
            silent_cells = list(set(all_post_cell_gids).difference(spk_dict_post_[pop].keys()))
            silent_cells.sort()
            for silent_cell in silent_cells:
                spk_dict_post_[pop].update({silent_cell:[]})

            # --- Storing the information in a sorted dictionary
            keylist = list(spk_dict_post_[pop].keys())
            keylist.sort()
            for key in keylist:
                spk_dict_post[pop].update({key:spk_dict_post_[pop][key]})

        return spk_dict_post, silent_cells
    
    def generatePostSpkConn(loadPop, all_post_cell_gids, connFolder):
        # === POSTSYNAPTIC CONNECTIVITY
        # --- Load all connectivity from presynaptic to postsynaptic cell beforehand
        conn_dict_post={}
        for pop in loadPop:
            conn_dict_post.update({pop:{}})
            print('\n\n##############################################')
            print('             Generating %s %s connections            '%(pop,len(all_post_cell_gids)))
            print('##############################################')

            for post_cell_gid_ind, post_cell_gid in enumerate(all_post_cell_gids):
                if (post_cell_gid-all_post_cell_gids[0])%100==0:
                    print('\t %s cells left'%((len(all_post_cell_gids)-post_cell_gid_ind)))
                
                conn_dict_post[pop].update({post_cell_gid:[]})
                loadFile=connFolder+pop+'_'+str(post_cell_gid)+'.pkl'
                with open(loadFile, 'rb') as cell_fileObj: conn_data = pickle.load(cell_fileObj)
                conn_dict_post[pop][post_cell_gid]=conn_data
        
        return conn_dict_post
    
    def generatePreSpkHist(conn_dict_post, all_post_cell_gids):
        # === PRESYNAPTIC SPIKES
        # --- Obtaining all the GIDs of presynaptic cells to load spikes beforehand
        pre_gids_=[]
        for pop in conn_dict_post.keys():
            for post_cell_gid in all_post_cell_gids:
                for i in range(len(conn_dict_post[pop][post_cell_gid])):
                    pre_gids_.append(conn_dict_post[pop][post_cell_gid][i][0])
        pre_cell_gids=list(set(pre_gids_))
        pre_cell_gids.sort()

        # --- Dictionary to store all presynaptic cell spikes by GID
        pre_cell_spks=netpyne.analysis.tools.getSpktSpkid(cellGids=pre_cell_gids)
        spk_dict_pre=Connectivity.getSpkDict(pre_cell_spks)

        # --- Adding empty lists to represent cells that did not spike, but are still connected
        for pre_cell_gid in pre_cell_gids:
            pre_spiking_list = list(spk_dict_pre.keys())
            if pre_cell_gid not in pre_spiking_list: spk_dict_pre.update({pre_cell_gid:[]})

        return spk_dict_pre

    def generateSpikeHistogram(conn_dict, time_bins, post_cell_spk_times, connected_spk_dict_pre):
        
        # ------------------------------------------------------------------------------------------------- #
        #                     Generates STPH for each postsynaptic cell
        # ------------------------------------------------------------------------------------------------- #
        # -------------------------------- Input arguments ------------------------------------------------ #
        # conn_dict                 : connection details
        # time_bins                 : time bins of the histogram to build
        # post_cell_spk_times       : post cell spike times
        # connected_spk_dict_pre    : dict{pre_cell_gid:pre cell spike times} (only connected cells)
        # ------------------------------------------------------------------------------------------------- #
        
        # --- Dictionary to store the spike histogram data for each postsynaptic cell
        post_cell_spike_hist_dict={}
        # --- Iterates over the presynaptic pop
        for pre_pop in conn_dict.keys():
            post_cell_spike_hist_dict.update({pre_pop:{}})
            # --- Iterates over the types of synaptic mechanisms ('exc' or 'inh')
            for mech_type in conn_dict[pre_pop].keys():
                post_cell_spike_hist_dict[pre_pop].update({mech_type:[]})
                # --- List to store sum of weighted normalized histograms
                weighted_spike_histograms=[]
                # --- Iterates over the GIDs of the presynaptic cells
                for pre_cell_gid in conn_dict[pre_pop][mech_type].keys():
                    valid_spike_differences = []
                    # --- Iterates over the list of spike times from post cell
                    for post_cell_spkt in post_cell_spk_times: 
                        # --- Iterates over the list of spike times from pre cell
                        for pre_cell_spkt in connected_spk_dict_pre[pre_cell_gid]: 
                            # --- Delay to propagate the spike
                            pre_spike_delay = conn_dict[pre_pop][mech_type][pre_cell_gid]['avg_delay']
                            # --- Effective spike time (takes into account the spike delay to reach the postsynaptic cell)
                            effective_pre_cell_spkt = pre_cell_spkt+pre_spike_delay
                            # --- Pre spike occurred after Post spike
                            if (effective_pre_cell_spkt)>post_cell_spkt: continue
                            # --- Pre spike occurred before the max interval from Post spike (e.g.: 300-(18+5)>200 = True || 300-(150+5)>200 = False)
                            elif post_cell_spkt-(effective_pre_cell_spkt)>time_bins[-1]: continue
                            else:
                                # --- Spike difference (post-(pre+delay))
                                spike_difference = post_cell_spkt-effective_pre_cell_spkt
                                # print(pre_pop,mech_type,pre_cell_gid,post_cell_spkt,pre_cell_spkt,pre_spike_delay,spike_difference)
                                valid_spike_differences.append(spike_difference)
                    valid_spike_differences.sort()
                    # print(valid_spike_differences)

                    # --- Calculating spike histogram
                    spike_histogram,edges = np.histogram(valid_spike_differences,range=[0,time_bins[-1]],bins=len(time_bins))
                    # print('spike histogram: ', spike_histogram)
                    
                    # --- BREAKS THE CODE SO THAT EMPTY HISTOGRAMS ARE NOT STORED
                    if all(spk_count == 0 for spk_count in spike_histogram): continue
                    
                    list_spike_histogram = list(spike_histogram)
                    # --- Normalizing the histogram by the number of postsynaptic spikes, so that the firing frequency of a postsynaptic cell does not affect the results 
                    normalized_spike_histogram = [bin/(len(post_cell_spk_times)) for bin in list_spike_histogram]
                    # print(pre_pop,mech_type,pre_cell_gid, '\t normalized histogram: ',normalized_spike_histogram)

                    sum_weight = conn_dict[pre_pop][mech_type][pre_cell_gid]['sum_weight']
                    # print(list(normalized_spike_histogram),sum_weight)
                    
                    weighted_spike_histogram=[norm_spk*sum_weight for norm_spk in normalized_spike_histogram]
                    # print(pre_pop,mech_type,pre_cell_gid, '\t weighted histogram: ',weighted_spike_histogram,'\n')

                    weighted_spike_histograms.append(weighted_spike_histogram)
                
                # print(pre_pop,mech_type,pre_cell_gid,'weighted_spike_histograms: ',weighted_spike_histograms)

                # --- Sum of weighted normalized spike histograms (empty histograms removed)
                sum_weighted_spike_histograms_=np.sum(weighted_spike_histograms,axis=0)
                try:    sum_weighted_spike_histograms = list(sum_weighted_spike_histograms_)
                except: sum_weighted_spike_histograms = [0]*len(time_bins)
                # print(type(sum_weighted_spike_histograms), '\t', sum_weighted_spike_histograms)
                
                # --- Final data: List with the sum of weighted normalized spike histograms (empty histograms removed) for each time_bin
                post_cell_spike_hist_dict[pre_pop][mech_type]=sum_weighted_spike_histograms

        return post_cell_spike_hist_dict


class PlotFigures():
    
    def formatData(pop_spike_hist_dict,all_pops):
        
        mech_types = ['exc','inh']

        # --- Dictionary to hold plotting variables 
        #     (because data is organized as dict[post_pop][post_cell_gid][pre_pop][mech_type][pre_cell_gid], but plots are easier organized by dict[mech_type][pre_pop])
        #     (also, [post_pop] is not a major factor, because spike_hist_dict only contains a single [post_pop], which also goes in the filename)
        pre_cell_pops_=[]
        for post_cell_gid in pop_spike_hist_dict.keys():
            # print(pop_spike_hist_dict[post_cell_gid].keys())
            for pre_pop in pop_spike_hist_dict[post_cell_gid].keys():
                pre_cell_pops_.append(pre_pop)
        pre_cell_pops=list(set(pre_cell_pops_))
        pre_cell_pops.sort()

        # --- Ordered populations according to network distribution
        ordered_pre_pops=[]
        for o_pop in all_pops:
            if o_pop in pre_cell_pops: ordered_pre_pops.append(o_pop)

        # print('Pre pop info | \npre_cell_pops: ',pre_cell_pops, '\nordered_pre_pops: ', ordered_pre_pops)

        # --- Allocating keys for the plotting dictionary
        plot_spike_hist_dict={}
        for mech_type in mech_types:
            plot_spike_hist_dict.update({mech_type:{}})
            for pre_cell_pop in ordered_pre_pops:
                plot_spike_hist_dict[mech_type].update({pre_cell_pop:{}})

        # --- Filtering out the empty histograms (post_cell has valid spikes, but out of timeRange)
        valid_post_cell_gids=[]
        for post_cell_gid in pop_spike_hist_dict.keys():
            # --- Filters out empty dictionaries (no histogram, because the cell fired out of timeRange)
            if any(pop_spike_hist_dict[post_cell_gid]): valid_post_cell_gids.append(post_cell_gid)
        
        # --- Calculating Mean and Std across postsynaptic cells
        for pre_cell_pop in ordered_pre_pops:
            for mech_type in mech_types:
                post_cell_hists=[]
                # --- Iterates over valid histograms (cells that fired) - cells that are originally silent/quiet are filtered out
                for post_cell_gid in valid_post_cell_gids:
                    if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                        post_cell_hist = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type]
                        post_cell_hists.append(post_cell_hist)
                post_pop_hist_mean  = np.mean(post_cell_hists,axis=0)
                post_pop_hist_std   = np.std(post_cell_hists,axis=0)
                # print(mech_type, pre_cell_pop, '->', pop,' : ', len(post_cell_hists))
                # print(post_pop_hist_mean,post_pop_hist_std,'\n')
                
                plot_spike_hist_dict[mech_type][pre_cell_pop].update({'mean':post_pop_hist_mean,'std':post_pop_hist_std})

        return plot_spike_hist_dict,ordered_pre_pops,valid_post_cell_gids

    def plotSPTHtraces(plot_spike_hist_dict, ordered_pre_pops, time_bins, map_colors, divide_plots=False, select_plot_pops=None):
        print('SPTH traces method')
        
        # --- P.s.: Select_plot_pops should be a dictionary with the following format
        #     select_plot_pops = {'exc':['pre_pop1', 'pre_pop2', ...], 'inh':['pre_pop2', 'pre_pop5', ...]}
        #     e.g.:
        #     select_plot_pops = {'exc':['PT5B','TVL','TPO','IT2','IT5A'],'inh':['PV5A','PV5B',]}

        # --- Figure
        plt.figure(figsize=(15,10))
        # plt.suptitle('spk_histogram_data: '+pop+' histogram | '+network_state+' | '+timeRange_name+' ms')

        # --- Plot variables
        plot_xlim=(0,55)
        plot_xticks=list(range(5,51,5))
        plot_xticks_labels=[str(-1*(x))+' ← '+str(-1*(x-5)) for x in plot_xticks]
        
        # plot_ylim=None
        plot_ylim=(-1,17.5)
        # plot_yticks=None
        plot_yticks=[0,5,10,15]

        plot_errorbar_capsize=3
        plot_errorbar_alpha=0.75

        plot_pre_pop=ordered_pre_pops

        for mech_ind, mech_type in enumerate(plot_spike_hist_dict.keys()):
            
            mech_color=['b','r']
            plot_errorbar_marker=['^','v']

            if divide_plots: plt.subplot(1, 2, mech_ind+1)

            for pre_pop in plot_pre_pop:
                # --- Change colors for the line plots
                c=map_colors[pre_pop]

                x_data  = time_bins
                y_data  = plot_spike_hist_dict[mech_type][pre_pop]['mean']
                y_std   = plot_spike_hist_dict[mech_type][pre_pop]['std']
                
                if all(y_d == 0 for y_d in y_data) and all(y_s == 0 for y_s in y_std): continue

                if select_plot_pops is not None:
                    if pre_pop not in select_plot_pops[mech_type]: continue

                if len(y_data)==len(x_data):
                    line_style='-'
                    if pre_pop.startswith('PT5B'):c='lightgrey'; line_style='--'# changes color of 'PT5B' to grey
                    plt.plot(x_data,y_data,line_style,color=c,linewidth=3,label=pre_pop+'_'+mech_type)
                    for ind_p,p in enumerate(y_std):
                        # --- Plotting traces only
                        # plt.plot(time_bins[ind_p],mean_exc_spks[ind_p],'o',markerfacecolor='blue',markeredgecolor='w')
                        # --- Plotting errorbars. source: https://stackoverflow.com/questions/22481854/plot-mean-and-standard-deviation
                        plt.errorbar(x_data[ind_p],y_data[ind_p], y_std[ind_p],linestyle='None', color=c, marker=plot_errorbar_marker[mech_ind], markeredgecolor=mech_color[mech_ind], capsize=plot_errorbar_capsize, alpha=plot_errorbar_alpha)
                        # plt.errorbar(x_data[ind_p],y_data[ind_p], y_std[ind_p],linestyle='None', color=mech_color[mech_ind], marker=plot_errorbar_marker, capsize=plot_errorbar_capsize, alpha=plot_errorbar_alpha)

                else: print('skipping plot for ', mech_type, pre_pop)

            # --- Plot formatting
            if divide_plots: plt.title(mech_type)
            else:            plt.title('Presynaptic population')
            plt.legend(loc='upper left')
            plt.xlim(plot_xlim)
            plt.xticks(ticks=plot_xticks, labels=plot_xticks_labels)
            plt.xlabel('Time interval before Postsynaptic spike (ms)')
            if plot_ylim is not None:   plt.ylim(plot_ylim)
            if plot_yticks is not None: plt.yticks(ticks=plot_yticks)                    
            plt.ylabel('Number of presynaptic spikes * synaptic strength')
            plt.gca().invert_xaxis()
        
    def barPlot(pop_spike_hist_dict,ordered_pre_pops,valid_post_cell_gids,max_time_index,all_pops,ax_lim=None):

        mech_types = ['exc','inh']

        # --- Allocating keys for the plotting the bar graph
        plot_bar_hist_dict={}
        for mech_type in mech_types:
            plot_bar_hist_dict.update({mech_type:{}})
            for pre_cell_pop in ordered_pre_pops:
                plot_bar_hist_dict[mech_type].update({pre_cell_pop:{}})

        # --- Calculating Mean and Std across postsynaptic cells
        for mech_type in mech_types:
            for pre_cell_pop in ordered_pre_pops:
                post_cell_hists_windowSum=[]
                for post_cell_gid in valid_post_cell_gids:
                    if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                        # --- Selects a time window from the histogram and adds it together to calculate the total weight
                        post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                        # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                        post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                
                post_pop_hist_windowSum_mean  = np.mean(post_cell_hists_windowSum)
                post_pop_hist_windowSum_std   = np.std(post_cell_hists_windowSum)
                
                plot_bar_hist_dict[mech_type][pre_cell_pop].update({'mean':post_pop_hist_windowSum_mean,'std':post_pop_hist_windowSum_std})
        
        # --- Draw Figure - Bar plot v4
        plt.figure(figsize=(5,5))
        # plt.suptitle('Bar plot v4: '+pop+' histogram | '+str(max_time)+' max time | '+network_state+' | '+timeRange_name+' ms')
        # plt.grid(color='lightgrey')
        
        # --- Plot style
        divide_plots=False

        sum_of_means_exc=0
        sum_of_means_inh=0
        list_of_stds_exc=[] # sum of stds: https://study.com/skill/learn/how-to-calculate-the-standard-deviation-of-the-sum-of-two-random-variables-explanation.html
        list_of_stds_inh=[] # sum of stds: https://study.com/skill/learn/how-to-calculate-the-standard-deviation-of-the-sum-of-two-random-variables-explanation.html

        for pre_pop_ind, pre_pop in enumerate(ordered_pre_pops):
        # for pre_pop_ind, pre_pop in enumerate(plot_bar_hist_dict['exc'].keys()):
            for mech_ind,mech_type in enumerate(plot_bar_hist_dict.keys()):
                if divide_plots: plt.subplot(1, 2, mech_ind+1)

                if mech_type == 'exc':
                    c = 'royalblue'
                    if plot_bar_hist_dict['inh'][pre_pop]['mean'] == 0 and plot_bar_hist_dict['inh'][pre_pop]['std'] == 0: displace=0
                    else: displace=-0.2
                elif mech_type == 'inh':
                    c = 'r'
                    if plot_bar_hist_dict['exc'][pre_pop]['mean'] == 0 and plot_bar_hist_dict['exc'][pre_pop]['std'] == 0: displace=0
                    else: displace=0.2
                else:
                    c = 'k'

                # --- skip empty bars
                if plot_bar_hist_dict[mech_type][pre_pop]['mean'] == 0 and plot_bar_hist_dict[mech_type][pre_pop]['std'] == 0: continue

                x_data = pre_pop_ind

                if divide_plots:    x_name=pre_pop
                else:               x_name=pre_pop+'_'+mech_type

                plt.barh(   x_data+displace,
                            plot_bar_hist_dict[mech_type][pre_pop]['mean'],
                            0.4,
                            xerr=plot_bar_hist_dict[mech_type][pre_pop]['std'],
                            capsize=3,
                            color=c,
                            label=x_name)
                
                if mech_type == 'exc':
                    sum_of_means_exc+=plot_bar_hist_dict[mech_type][pre_pop]['mean']
                    list_of_stds_exc.append(plot_bar_hist_dict[mech_type][pre_pop]['std'])
                else:
                    sum_of_means_inh+=plot_bar_hist_dict[mech_type][pre_pop]['mean']
                    list_of_stds_inh.append(plot_bar_hist_dict[mech_type][pre_pop]['std'])

        # --- Square of the stds
        std_exc_sqr=[std**2 for std in list_of_stds_exc]
        std_inh_sqr=[std**2 for std in list_of_stds_inh]
        
        sum_std_exc=math.sqrt(sum(std_exc_sqr))
        sum_std_inh=math.sqrt(sum(std_inh_sqr))

        # --- plot of sum of bars - exc
        plt.barh(   len(ordered_pre_pops)+1-0.2,
                    # len(plot_bar_hist_dict[mech_type].keys())+1-0.2,
                    sum_of_means_exc,
                    0.4,
                    xerr=sum_std_exc,
                    capsize=3,
                    color='royalblue',)
        # --- plot of sum of bars - inh
        plt.barh(   len(ordered_pre_pops)+1+0.2,
                    # len(plot_bar_hist_dict[mech_type].keys())+1+0.2,
                    sum_of_means_inh,
                    0.4,
                    xerr=sum_std_inh,
                    capsize=3,
                    color='r',)

        # --- Format figure
        plot_labels = ordered_pre_pops+['','Sum']
        # plot_labels = list(plot_bar_hist_dict['exc'].keys())+['','Sum']
        plt.yticks(ticks=list(range(len(plot_labels))),labels=plot_labels)
        plt.gca().invert_yaxis()
        if ax_lim is not None: plt.xlim((ax_lim))

        print('bar plot method')

    def boxPlot(pop_spike_hist_dict, ordered_pre_pops, valid_post_cell_gids, max_time_index, ax_lim=None, select_pops=None):

        mech_types = ['exc','inh']
        mech_colors = ['royalblue','r']
        # --- Boxplot - Figure
        fig = plt.figure(figsize=(5,5))
        # fig.suptitle('Spike Histogram Window Sum Boxplot: '+pop+' histogram | '+str(max_time)+' max time | '+network_state+' | '+timeRange_name+' ms')

        boxplot_data=[]
        boxplot_labels=[]
        boxplot_colors=[]
        
        # --- Calculating Mean and Std across postsynaptic cells
        boxplot_hist_dict={}    # --- Dictionary with data for Boxplot
        for mech_type in mech_types:
            boxplot_hist_dict.update({mech_type:{}})
            for pre_cell_pop in ordered_pre_pops:
                boxplot_hist_dict[mech_type].update({pre_cell_pop:[]})
                post_cell_hists_windowSum=[]
                for post_cell_gid in valid_post_cell_gids:
                    if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                        # --- Selects a time window from the histogram and adds it together to calculate the total weight
                        post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                        # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                        post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                # --- Dictionary with data for Boxplot
                boxplot_hist_dict[mech_type][pre_cell_pop]=post_cell_hists_windowSum

        # --- Boxplot Data
        for mech_ind,mech_type in enumerate(boxplot_hist_dict.keys()):
            for pre_pop in boxplot_hist_dict[mech_type].keys():
                if select_pops is not None:
                    if ((mech_type == mech_types[mech_ind]) and (pre_pop in select_pops[mech_ind])):
                        boxplot_data.append(boxplot_hist_dict[mech_type][pre_pop])
                        boxplot_labels.append(pre_pop)
                        boxplot_colors.append(mech_colors[mech_ind])
                    else: continue # --- quits plotting if pre_pop not in the plot list
                else:
                    boxplot_data.append(boxplot_hist_dict[mech_type][pre_pop])
                    boxplot_labels.append(pre_pop)
                    boxplot_colors.append(mech_colors[mech_ind])     
        
        ax = fig.add_subplot(111)
        bp = ax.boxplot(boxplot_data,patch_artist = True,notch ='True')
        bp_colors = [bp_color for bp_color in boxplot_colors]
        for patch, bp_color in zip(bp['boxes'], bp_colors):
            patch.set_facecolor(bp_color)
        # Changing color and linewidth of whiskers
        for whisker in bp['whiskers']:
            whisker.set(color ='grey',linewidth = 1.5,linestyle =":")
        # changing color and linewidth of medians
        for median in bp['medians']:
            median.set(color ='k',linewidth = 1)
        # changing style of fliers
        for flier in bp['fliers']:
            flier.set(marker ='.',color ='k',alpha = 0.5)

        boxplot_ticks = [tick+1 for tick in range(len(boxplot_labels))]
        ax.set_xticks(ticks=boxplot_ticks,rotation=45,labels=boxplot_labels)
        if ax_lim is not None: ax.set_ylim(ax_lim)

    def violinPlot(pop_spike_hist_dict, ordered_pre_pops, valid_post_cell_gids, max_time_index, ax_lim=None, select_pops=None):

        mech_types = ['exc','inh']
        mech_colors = ['royalblue','r']
        # --- Violin - Figure
        fig = plt.figure(figsize=(5,5))
        # fig.suptitle('Spike Histogram Window Sum Violin: '+pop+' histogram | '+str(max_time)+' max time | '+network_state+' | '+timeRange_name+' ms')

        violinplot_data=[]; violinplot_labels=[]; violinplot_colors=[]
        
        # --- Calculating Mean and Std across postsynaptic cells
        violinplot_hist_dict={}    # --- Dictionary with data for Violin
        for mech_type in mech_types:
            violinplot_hist_dict.update({mech_type:{}})
            for pre_cell_pop in ordered_pre_pops:
                violinplot_hist_dict[mech_type].update({pre_cell_pop:[]})
                post_cell_hists_windowSum=[]
                for post_cell_gid in valid_post_cell_gids:
                    if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                        # --- Selects a time window from the histogram and adds it together to calculate the total weight
                        post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                        # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                        post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                # --- Dictionary with data for Violin
                violinplot_hist_dict[mech_type][pre_cell_pop]=post_cell_hists_windowSum

        # --- Violin
        for mech_ind,mech_type in enumerate(violinplot_hist_dict.keys()):
            for pre_pop in violinplot_hist_dict[mech_type].keys():
                if select_pops is not None:
                    if ((mech_type == mech_types[mech_ind]) and (pre_pop in select_pops[mech_ind])):
                        violinplot_data.append(violinplot_hist_dict[mech_type][pre_pop])
                        violinplot_labels.append(pre_pop)
                        violinplot_colors.append(mech_colors[mech_ind])
                    else: continue # --- quits plotting if pre_pop not in the plot list
                else:
                    violinplot_data.append(violinplot_hist_dict[mech_type][pre_pop])
                    violinplot_labels.append(pre_pop)
                    violinplot_colors.append(mech_colors[mech_ind])  
        
        ax = fig.add_subplot(111)

        bp = ax.violinplot(violinplot_data)
        violinplot_ticks = [tick+1 for tick in range(len(violinplot_labels))]
        ax.set_xticks(ticks=violinplot_ticks,rotation=45,labels=violinplot_labels)
        if ax_lim is not None: ax.set_ylim(ax_lim)

    def scatterPlot(pop_spk_dict_post, timeRange, all_post_cell_gids, cell_tags_dict, select_max_val = None, select_colormap = None, use_x_position=True):
        # --- Scatter - Figure
        fig = plt.figure(figsize=(5,5))

        if select_colormap is not None: c_map = select_colormap
        else: c_map = 'jet'

        if select_max_val is not None: reference_value = select_max_val+1
        else:
            spk_count = [len(pop_spk_dict_post[key]) for key in pop_spk_dict_post.keys()]
            reference_value = max(spk_count)
        scatter_map_colors=ColorMap.colormap(list(range(reference_value)),c_map=c_map)

        for post_cell_ind, post_cell_gid in enumerate(all_post_cell_gids):
            cell_position_x = cell_tags_dict[post_cell_gid]['xnorm']
            cell_position_y = cell_tags_dict[post_cell_gid]['ynorm']
            if post_cell_gid not in (pop_spk_dict_post.keys()): cell_firing=0
            else:
                valid_spikes=[]
                for spkt in pop_spk_dict_post[post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_num=len(valid_spikes)
                cell_firing = spk_num
            
            # --- Different colors for cells with no spikes
            if cell_firing>0:
                marker_color = scatter_map_colors[cell_firing]
                edge_color = None
                alpha = 1
            else:
                marker_color = 'w'
                edge_color = 'k'
                alpha = 0.25

            # --- Choose which data to show
            if use_x_position:
                data_1 = cell_position_x
                data_2 = cell_position_y
            else:
                data_1 = cell_firing
                data_2 = cell_position_y

            plt.plot(data_1,data_2,marker='o',color=marker_color,markeredgecolor=edge_color,alpha=alpha)
        plt.gca().invert_yaxis()

    # --- Method to normalize a dataset
    def NormalizeData(data):
        return (data - np.min(data)) / (np.max(data) - np.min(data))

    def formatMultivariateData(pop_spk_dict_post,pop_spike_hist_dict,featuredPops,valid_post_cell_gids,timeRange,max_time_index,cell_tags_dict,target_data = 'spk',isolate_mech=None,c_map='jet'):
        
        # --- PCA TUTORIAL: 
        #       https://github.com/mGalarnyk/Python_Tutorials/blob/master/Sklearn/PCA/PCA_Data_Visualization_Iris_Dataset_Blog.ipynb
        # --- Pandas Dataframe tutorials:
        #       https://www.includehelp.com/python/dataframe-is-it-pass-by-value-or-pass-by-reference.aspx
        #       https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.from_dict.html

        mech_types = ['exc','inh']
        
        # --- Allocating keys for the Multivariate analysis
        Multivariate_hist_dict={}
        # --- Loops through the list of cells with valid spikes
        for post_cell_gid in valid_post_cell_gids:
            Multivariate_hist_dict.update({post_cell_gid:{}})
            for mech_type in mech_types:
                Multivariate_hist_dict[post_cell_gid].update({mech_type:{}})
                for pre_cell_pop in featuredPops:
                    Multivariate_hist_dict[post_cell_gid][mech_type].update({pre_cell_pop:{}})

        print(Multivariate_hist_dict[post_cell_gid][mech_type].keys())

        # --- Calculating Mean and Std across postsynaptic cells
        for post_cell_gid in Multivariate_hist_dict.keys():
            for mech_type in Multivariate_hist_dict[post_cell_gid].keys():
                for pre_cell_pop in Multivariate_hist_dict[post_cell_gid][mech_type].keys():
                    # --- Checks if pre_cell_pop is part of that cell's dictionary (SOM5A and PV5A projections are not present in all cells)
                    if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                        # --- Selects a time window from the histogram and adds it together to calculate the total weight
                        Multivariate_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                        # Multivariate_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                    # else:
                    #     print(pre_cell_pop, ' not in ', pop_spike_hist_dict[post_cell_gid].keys())
                    Multivariate_hist_dict[post_cell_gid][mech_type].update({pre_cell_pop:sum(Multivariate_hist_windowSum)})

        # --- Generate Colormap
        if target_data == 'spk':        # Based on number of spikes
            # --- Calculating the max number of spikes to generate colormap
            spk_nums=[]
            for post_cell_gid in Multivariate_hist_dict.keys():
                valid_spikes=[]
                for spkt in pop_spk_dict_post[post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_nums.append(len(valid_spikes))                
            max_spk_num=max(spk_nums)
            pca_map_colors=ColorMap.colormap(list(range(max_spk_num+1)),c_map=c_map)
            # --- List with Colormap boundaries
            colormap_boudaries=[0,max_spk_num]

        elif target_data == 'ynorm':    # Based on cell position in the y-axis
            ynorm_list=[]
            for post_cell_gid in Multivariate_hist_dict.keys():
                ynorm_list.append(cell_tags_dict[post_cell_gid]['ynorm'])

            ynorm_list.sort()
            steprange = range(len(ynorm_list))
            color_map = plt.get_cmap(c_map)

            # --- List with Colormap boundaries
            colormap_boudaries=[min(ynorm_list),max(ynorm_list)]

            pca_map_reference={}
            for j,k in enumerate(ynorm_list):
                colormap_val=color_map(1.*j/float(len(steprange)))
                pca_map_reference.update({k:colormap_val})

            pca_map_colors={}
            for post_cell_gid in Multivariate_hist_dict.keys():
                pca_map_colors.update({post_cell_gid:pca_map_reference[cell_tags_dict[post_cell_gid]['ynorm']]})

        elif target_data[0] == 'true_spk':
            # --- Calculating the max number of spikes to generate colormap
            spk_nums=[]
            for post_cell_gid in Multivariate_hist_dict.keys():
                valid_spikes=[]
                for spkt in target_data[1][post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_nums.append(len(valid_spikes))                
            max_spk_num=max(spk_nums)
            pca_map_colors=ColorMap.colormap(list(range(max_spk_num+1)),c_map=c_map)
            # --- List with Colormap boundaries
            colormap_boudaries=[0,max_spk_num]

        elif target_data[0] == 'fixed_spk_range':
            # --- Generate colormap for fixed spike range
            pca_map_colors=ColorMap.colormap(list(range(target_data[2][-1]+1)),c_map=c_map)
            # --- List with Colormap boundaries
            colormap_boudaries=target_data[2]

        elif target_data[0] == 'capped_spk_range':
            # --- Generate colormap for fixed spike range
            pca_map_colors=ColorMap.colormap(list(range(target_data[2][-1]+1)),c_map=c_map)
            # --- List with Colormap boundaries
            colormap_boudaries=target_data[2]

        if isolate_mech is not None: print('plotting only a single mech: ', isolate_mech)

        # --- Creating a dictionary to organize the data in the Pandas DataFrame format
        DataFrame_dict={}
        for post_cell_gid in Multivariate_hist_dict.keys():
            DataFrame_dict.update({post_cell_gid:{}})
            for mech_type in Multivariate_hist_dict[post_cell_gid].keys():
                # --- Isolating a single mech_type in the PCA
                if isolate_mech is not None:
                    if mech_type == isolate_mech: continue
                for pre_cell_pop in Multivariate_hist_dict[post_cell_gid][mech_type].keys():
                    DataFrame_dict[post_cell_gid].update({pre_cell_pop+'_'+mech_type:Multivariate_hist_dict[post_cell_gid][mech_type][pre_cell_pop]})

        for post_cell_ind,post_cell_gid in enumerate(Multivariate_hist_dict.keys()):
            if target_data == 'spk':
                # --- Adding the number of spikes as a feature
                valid_spikes=[]
                for spkt in pop_spk_dict_post[post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_num=len(valid_spikes)
                DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[spk_num]})
            elif target_data == 'ynorm':
                DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[post_cell_gid]})
                # DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[ynorm_norm[post_cell_ind]]})
            elif target_data[0] == 'true_spk':
                # --- Adding the number of spikes as a feature
                valid_spikes=[]
                for spkt in target_data[1][post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_num=len(valid_spikes)
                DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[spk_num]})
            elif target_data[0] == 'fixed_spk_range':
                # --- Adding the number of spikes as a feature
                valid_spikes=[]
                for spkt in target_data[1][post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_num=len(valid_spikes)
                DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[spk_num]})
            elif target_data[0] == 'capped_spk_range':
                # --- Adding the number of spikes as a feature
                valid_spikes=[]
                for spkt in target_data[1][post_cell_gid]:
                    if spkt>=timeRange[0] and spkt<timeRange[1]: valid_spikes.append(spkt)
                spk_num=len(valid_spikes)
                # --- capped data to upper limit
                if spk_num>target_data[2][-1]: spk_num=target_data[2][-1]
                DataFrame_dict[post_cell_gid].update({'target':pca_map_colors[spk_num]})

        return DataFrame_dict, colormap_boudaries

    def createDataFrame(DataFrame_dict):
        import pandas as pd
        from sklearn.preprocessing import StandardScaler
        # --- Creating a Pandas DataFrame
        df              = pd.DataFrame.from_dict(DataFrame_dict, orient='index')
        df_features     = list(df.columns)[:-1] # removing last key, because is the 'target', used as classifier variable
        df_values       = df.loc[:, df_features].values

        # --- Creating a new Dataframe for Target values because the GIDs and indexes were being mixed up
        df_target_      = df.loc[:,['target']].values
        df_target       = pd.DataFrame(df_target_)

        # --- Data Normalization
        df_values_Norm  = StandardScaler().fit_transform(df_values)

        return df,df_features,df_values,df_target,df_values_Norm

    def plotPCA(DataFrame_dict,pop_spk_dict_post,n_components=2):
        
        import pandas as pd
        from sklearn.decomposition import PCA
        df,df_features,df_values,df_target,df_values_Norm = PlotFigures.createDataFrame(DataFrame_dict)
        
        pca = PCA(n_components=n_components)
        principalComponents = pca.fit_transform(df_values_Norm)
        principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])

        finalDf = pd.concat([principalDf, df_target], axis = 1)

        fig = plt.figure(figsize = (5,5))
        ax  = fig.add_subplot(1,1,1) 
        ax.set_xlabel('Principal Component 1', fontsize = 15)
        ax.set_ylabel('Principal Component 2', fontsize = 15)
        colors=list(df[['target']].values)

        post_cell_gid_list=list(pop_spk_dict_post.keys())
        off_PCA_value_gids=[]
        for ind,color in enumerate(colors):
            ax.scatter(finalDf.loc[ind, 'PC1'], finalDf.loc[ind, 'PC2'], c = color, s = 50)
            if finalDf.loc[ind, 'PC1']>5: off_PCA_value_gids.append(post_cell_gid_list[ind])

        # --- Removing boxes from plot axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # --- Removing ticks from plot
        ax.set_xticks(ticks=[],labels=None)
        ax.set_yticks(ticks=[],labels=None)

        return pca

    def plotUMAP(DataFrame_dict,n_neighbors=30,color_criteria='gid',cellType_dict=None,c_map='jet',colormap_boudaries=None):
        import umap
        import sklearn

        df,df_features,df_values,df_target,df_values_Norm = PlotFigures.createDataFrame(DataFrame_dict)
        reducer     = umap.UMAP(random_state=1, n_neighbors=n_neighbors, min_dist=0.0, n_components=2)  # 15, 0.1
        embedding   = reducer.fit_transform(df_values_Norm)
        embedding.shape

        if color_criteria == 'gid': colors=ColorMap.colormap(list(df.index),c_map=c_map)
        elif color_criteria == 'cellType':
            colors={}
            for gid in list(df.index):
                if gid in cellType_dict['enhanced']:    c='seagreen'
                elif gid in cellType_dict['suppressed']:c='gold'
                else:                                   c='k'
                colors.update({gid:c})
        # --- The 'target' values here are already passed as a colormap
        elif (color_criteria == 'ynorm') or (color_criteria == 'spk') or (color_criteria[0] == 'true_spk') or (color_criteria[0] == 'fixed_spk_range') or (color_criteria[0] == 'capped_spk_range'):
            colors_list=list(df[['target']].values)
            colors={}
            for ind,gid in enumerate(list(df.index)):
                colors.update({gid:colors_list[ind][0]})

            if colormap_boudaries is not None:
                # c = np.arange(1, len(colors_list.keys())+1)
                norm = mpl.colors.Normalize(vmin=min(colormap_boudaries), vmax=max(colormap_boudaries))
                colorbar_cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
                colorbar_cmap.set_array([])

        fig, ax = plt.subplots(figsize=(5,5))
        # key is the GID of the cells
        for ind,key in enumerate(colors.keys()):
            marker='.'
            ax.plot( embedding[ind, 0], embedding[ind, 1], color=colors[key], marker=marker, markersize=3,)
        ax.set_aspect('equal', 'datalim')

        # --- Removing boxes from plot axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # --- Removing ticks from plot
        ax.set_xticks(ticks=[],labels=None)
        ax.set_yticks(ticks=[],labels=None)

        if colormap_boudaries is not None: fig.colorbar(colorbar_cmap, ticks=colormap_boudaries)

        return embedding, df
    
    def applyKMeans(dataset,n_clusters = 4):
        # apply k-means to UMAP
        from sklearn import cluster
        kmeans_dataset = cluster.KMeans(n_clusters=n_clusters).fit(dataset)
        return kmeans_dataset

    def plotKMeans(kmeans_dataset, dataset, dataframe,customLabels=None, customColors=None):
        cluster_dictionary={}
        dataFrame_gids=list(dataframe.index)
        if customLabels is not None:    kMeans_labels=customLabels
        else:                           kMeans_labels=kmeans_dataset.labels_

        kMeans_list=list(kMeans_labels)
        kMeans_groups=list(set(kMeans_list))

        if customColors is not None: cluster_colors_list=customColors
        else:
            cluster_colors=ColorMap.colormap(list(range(kmeans_dataset.n_clusters)),c_map='jet')
            cluster_colors_list=list(cluster_colors.values())

        fig, ax = plt.subplots(figsize=(5,5))
        for i, c in zip(range(kmeans_dataset.n_clusters), cluster_colors_list):
            ax.scatter(dataset[kMeans_labels==i, 0], dataset[kMeans_labels==i, 1], color=c, marker='.',s=5,label='Cluster '+str(i+1),)
            ax.set_aspect('equal', 'datalim')
        
        lgnd=ax.legend(loc='lower left',framealpha=0)
        for i,lg in enumerate(lgnd.legendHandles):
            # change the marker size manually for both lines
            lgnd.legendHandles[i]._sizes = [500]

        for i in list(kMeans_groups):
            cluster_dictionary.update({i:[]})

        for ind, i in enumerate(kMeans_list):
            cluster_dictionary[i].append(dataFrame_gids[ind])
        
        # --- Removing boxes from plot axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        # --- Removing ticks from plot
        ax.set_xticks(ticks=[],labels=None)
        ax.set_yticks(ticks=[],labels=None)
        
        return cluster_dictionary
    
    def barplotKMeans2(cluster_dictionary, pop_spike_hist_dict, featuredPops, max_time_index, divide_plots=False):
        
        mech_types=['exc','inh']

        print('barplotKMeans')
        # --- Allocating keys for the plotting the bar graph
        plot_bar_kMeans={}
        for kMeans_group in cluster_dictionary.keys():
            plot_bar_kMeans.update({kMeans_group:{}})
            for mech_type in mech_types:
                plot_bar_kMeans[kMeans_group].update({mech_type:{}})
                for pre_cell_pop in featuredPops:
                    plot_bar_kMeans[kMeans_group][mech_type].update({pre_cell_pop:{}})

        for kMeans_group in cluster_dictionary.keys():
            # --- Calculating Mean and Std across postsynaptic cells
            for mech_type in mech_types:
                for pre_cell_pop in featuredPops:
                    post_cell_hists_windowSum=[]
                    for post_cell_gid in cluster_dictionary[kMeans_group]:
                        # --- Checks if pre_cell_pop is part of that cell's dictionary (SOM5A and PV5A projections are not present in all cells)
                        if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                            # --- Selects a time window from the histogram and adds it together to calculate the total weight
                            post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                            # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                            post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                    
                    post_pop_hist_windowSum_mean  = np.mean(post_cell_hists_windowSum)
                    post_pop_hist_windowSum_std   = np.std(post_cell_hists_windowSum)
                    
                    plot_bar_kMeans[kMeans_group][mech_type][pre_cell_pop].update({'mean':post_pop_hist_windowSum_mean,'std':post_pop_hist_windowSum_std})
        
        # --- Barplot - One fig per kMeans cluster
        fig = plt.figure(figsize=(15,10))
        for kMeans_ind,kMeans_group in enumerate(plot_bar_kMeans.keys()):
            # --- Figure - Bar plot v3
            ax = fig.add_subplot(1, len(plot_bar_kMeans.keys()), kMeans_ind+1)
            # plt.subplot(1, len(plot_bar_kMeans.keys()), kMeans_ind+1)
            for mech_ind,mech_type in enumerate(plot_bar_kMeans[kMeans_group].keys()):
                # if divide_plots: plt.subplot(1, 2, mech_ind+1)
                sum_of_means=0
                if mech_type == 'exc':
                    c = 'royalblue'; displace=-0.2
                elif mech_type == 'inh':
                    c = 'r'; displace=+0.2
                else:
                    c = 'k'
                for pre_pop_ind, pre_pop in enumerate(plot_bar_kMeans[kMeans_group][mech_type].keys()):
                    x_data = pre_pop_ind
                    if divide_plots:    x_name=pre_pop
                    else:               x_name=pre_pop+'_'+mech_type

                    if mech_ind==1: flipper = 1
                    else:           flipper = 1

                    # --- skip empty bars
                    if plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['mean'] == 0 and plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['std'] == 0: continue
                
                    # --- Fix 'nan' values from mean and std calculation
                    if math.isnan(plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['mean']): plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['mean'] = 0

                    if math.isnan(plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['std']): plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['std'] = 0

                    ax.barh(    x_data+displace,
                                flipper*plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['mean'],
                                0.4,
                                xerr=plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['std'],
                                capsize=3,
                                color=c,
                                label=x_name)
                    
                    # --- Calculates the sum of bars
                    sum_of_means+=plot_bar_kMeans[kMeans_group][mech_type][pre_pop]['mean']

                # --- plot of sum of bars
                ax.barh(    len(plot_bar_kMeans[kMeans_group][mech_type].keys())+1+displace,
                            flipper*sum_of_means,
                            0.4,
                            capsize=3,
                            color=c,)

            # --- Format figure
            plot_labels = list(plot_bar_kMeans[kMeans_group]['exc'].keys())+['','Sum']
            if kMeans_ind==0: ax.set_yticks(ticks=list(range(len(plot_labels))),labels=plot_labels)
            else:             ax.set_yticks(ticks=[],labels=None)
            ax.invert_yaxis()

    def boxplotKMeans(cluster_dictionary, pop_spike_hist_dict, featuredPops, max_time_index, divide_plots=False, select_pops=None):
        
        mech_types = ['exc','inh']
        mech_colors = ['royalblue','r']

        print('boxplotKMeans')
        # --- Dictionary with data for Boxplot
        boxplot_kMeans={}    
        for kMeans_group in cluster_dictionary.keys():
            boxplot_kMeans.update({kMeans_group:{}})
            # --- Calculating Mean and Std across postsynaptic cells
            for mech_type in mech_types:
                boxplot_kMeans[kMeans_group].update({mech_type:{}})
                for pre_cell_pop in featuredPops:
                    boxplot_kMeans[kMeans_group][mech_type].update({pre_cell_pop:[]})
                    post_cell_hists_windowSum=[]
                    for post_cell_gid in cluster_dictionary[kMeans_group]:
                        if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                            # --- Selects a time window from the histogram and adds it together to calculate the total weight
                            post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                            # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                            post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                    
                    # --- Dictionary with data for Boxplot
                    boxplot_kMeans[kMeans_group][mech_type][pre_cell_pop]=post_cell_hists_windowSum

        # --- Boxplot - Figure
        fig = plt.figure(figsize=(15,10))
        # fig.suptitle('Spike Histogram Window Sum Boxplot - K-Means: '+pop+' histogram | '+str(max_time)+' max time | '+network_state+' | '+timeRange_name+' ms')

        for kMeans_ind,kMeans_group in enumerate(boxplot_kMeans.keys()):
            # --- Plot style
            # plt.subplot(2, 2, kMeans_ind+1)
            # ax = fig.add_subplot(2, 2, kMeans_ind+1)
            ax = fig.add_subplot(1, len(boxplot_kMeans.keys()), kMeans_ind+1)
            boxplot_data=[]; boxplot_labels=[]; boxplot_colors=[]

            for mech_ind,mech_type in enumerate(boxplot_kMeans[kMeans_group].keys()):
                for pre_pop in boxplot_kMeans[kMeans_group][mech_type].keys():
                    # --- Selecting populations to plot in the figure
                    if select_pops is not None:
                        if ((mech_type == mech_types[mech_ind]) and (pre_pop in select_pops[mech_ind])):
                            boxplot_data.append(boxplot_kMeans[kMeans_group][mech_type][pre_pop])
                            boxplot_labels.append(pre_pop)
                            boxplot_colors.append(mech_colors[mech_ind])
                        else: continue # --- quits plotting if pre_pop not in the plot list
                    else:
                        boxplot_data.append(boxplot_kMeans[kMeans_group][mech_type][pre_pop])
                        boxplot_labels.append(pre_pop)
                        boxplot_colors.append(mech_colors[mech_ind])  

            # --- Vertical Plot Orientation
            vert=False
            bp = ax.boxplot(boxplot_data,patch_artist = True,notch ='True',vert=vert)
            bp_colors = [bp_color for bp_color in boxplot_colors]

            for patch, bp_color in zip(bp['boxes'], bp_colors):
                patch.set_facecolor(bp_color)
            
            # changing color and linewidth of whiskers
            for whisker in bp['whiskers']:
                whisker.set(color ='grey', linewidth = 1.5, linestyle =":")
            # changing color and linewidth of medians
            for median in bp['medians']:
                median.set(color ='k', linewidth = 1)
            # changing style of fliers
            for flier in bp['fliers']:
                flier.set(marker ='.', color ='k', alpha = 0.5)

            boxplot_ticks = [tick+1 for tick in range(len(boxplot_labels))]
            ax.set_title("Num of cells: %s"%(len(boxplot_data[0])))
            
            if vert: ax.set_xticks(ticks=boxplot_ticks,rotation=45,labels=boxplot_labels)                
            else:    ax.set_yticks(ticks=boxplot_ticks,rotation=45,labels=boxplot_labels)                
                
        return boxplot_kMeans

# --- Class to process and plot data merging QUIET and MOVEMENT states
class PostAnalysis():
    def mergedBarPlot(  spk_hist_A,spk_hist_B,
                        ordered_pre_pops,
                        post_GIDs_A,post_GIDs_B,
                        max_time_index,
                        long_range_pops,
                        select_pre_pops=None,
                        ax_lim=None,
                        verticalPlot=True,
                        perSecond=False,
                        timeScaling=None,
                        states = ['Quiet','Movement'],
                        state_colors=['royalblue','crimson'],
                        remove_long_range_inh=True,
                        showLabels=True,
                        export_values=True,
                        export_filename='export_dict.json',
                        plotStatistics=False
                        ):
        print('Merged bar plot')
        import scipy
        # A and B are two different states 
        # (e.g.: 
        #           Quiet vs Move; 
        #           Cluster 0 vs Cluster 1; 
        #           etc
        # )

        if select_pre_pops==None:
            select_pre_pops=ordered_pre_pops
            sumLeftover=False
        else: sumLeftover=True

        # --- Adding missing pops if selected but not present in the original list of pops (will appear as a empty bar)
        for s_pop in select_pre_pops:
            if s_pop not in ordered_pre_pops:
                print('adding ', s_pop, ' to ', ordered_pre_pops)
                ordered_pre_pops.append(s_pop)

        mech_types = ['exc','inh']

        # --- Allocating keys for the plotting the bar graph
        plot_bar_hist_dict={}
        for mech_type in mech_types:
            plot_bar_hist_dict.update({mech_type:{}})
            for pre_cell_pop in ordered_pre_pops:
                plot_bar_hist_dict[mech_type].update({pre_cell_pop:{}})
                for state in states:
                    plot_bar_hist_dict[mech_type][pre_cell_pop].update({state:{'mean':0,'std':0}})
            
        if 'PT5B' not in plot_bar_hist_dict['exc'].keys():
            print('NO PT5B', plot_bar_hist_dict['exc'].keys())
            plot_bar_hist_dict['exc'].update({'PT5B':{states[0]:{'mean':0,'std':0},states[1]:{'mean':0,'std':0}}})
            print('---->', plot_bar_hist_dict['exc'].keys())

        # --- Allocating keys for the exporting values
        if export_values:
            export_dict={}
            for state in states:
                export_dict.update({state:{}})
                for mech_type in mech_types:
                    export_dict[state].update({mech_type:{}})
                    for pre_cell_pop in ordered_pre_pops:
                        export_dict[state][mech_type].update({pre_cell_pop:{'mean':0,'std':0,'vals':[]}})
            
        # --- Calculating Mean and Std across postsynaptic cells
        for mech_type in mech_types:
            for pre_cell_pop in ordered_pre_pops:
                for state in states:
                    if state == states[0]:
                        valid_post_cell_gids = post_GIDs_A
                        pop_spike_hist_dict = spk_hist_A
                    elif state == states[1]:
                        valid_post_cell_gids = post_GIDs_B
                        pop_spike_hist_dict = spk_hist_B

                    post_cell_hists_windowSum=[]
                    for post_cell_gid in valid_post_cell_gids:
                        if pre_cell_pop in pop_spike_hist_dict[post_cell_gid].keys():
                            # --- Selects a time window from the histogram and adds it together to calculate the total weight
                            post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][0:max_time_index+1]
                            # post_cell_hist_windowSum = pop_spike_hist_dict[post_cell_gid][pre_cell_pop][mech_type][1:max_time_index+1]
                            post_cell_hists_windowSum.append(sum(post_cell_hist_windowSum))
                    
                    if perSecond:
                        post_pop_hist_windowSum_mean  = np.mean(post_cell_hists_windowSum)/timeScaling
                        post_pop_hist_windowSum_std   = np.std(post_cell_hists_windowSum)/timeScaling
                        post_pop_hist_windowSum_vals  = [x/timeScaling for x in post_cell_hists_windowSum]
                    else:
                        post_pop_hist_windowSum_mean  = np.mean(post_cell_hists_windowSum)
                        post_pop_hist_windowSum_std   = np.std(post_cell_hists_windowSum)
                        post_pop_hist_windowSum_vals  = post_cell_hists_windowSum
                    
                    plot_bar_hist_dict[mech_type][pre_cell_pop][state].update({'mean':post_pop_hist_windowSum_mean,'std':post_pop_hist_windowSum_std,'vals':post_pop_hist_windowSum_vals})
        
        # --- removing 'nan' values:
        for mech_type in mech_types:
            for pre_cell_pop in ordered_pre_pops:
                for state in states:
                    if math.isnan(plot_bar_hist_dict[mech_type][pre_cell_pop][state]['mean']): plot_bar_hist_dict[mech_type][pre_cell_pop][state]['mean']=0
                    if math.isnan(plot_bar_hist_dict[mech_type][pre_cell_pop][state]['std']):  plot_bar_hist_dict[mech_type][pre_cell_pop][state]['std']=0
            
        # --- Draw Figure - Bar plot v4
        if verticalPlot: fig, ax = plt.subplots(figsize=(10,5))
        else:            fig, ax = plt.subplots(figsize=(6,14))
        
        if ax_lim is not None: ax.set_xlim((ax_lim))

        # --- Plot style
        divide_plots=False

        sum_of_means_exc_Q=0;   sum_of_means_exc_M=0;   sum_of_means_inh_Q=0;   sum_of_means_inh_M=0
        list_of_stds_exc_Q=[];  list_of_stds_exc_M=[];  list_of_stds_inh_Q=[];  list_of_stds_inh_M=[]
        
        leftover_sum_of_means_exc_Q=0;  leftover_sum_of_means_exc_M=0;  leftover_sum_of_means_inh_Q=0;  leftover_sum_of_means_inh_M=0
        leftover_list_of_stds_exc_Q=[]; leftover_list_of_stds_exc_M=[]; leftover_list_of_stds_inh_Q=[]; leftover_list_of_stds_inh_M=[]

        errbar_color='k'; errbar_alpha=0.5

        bar_widht_thick = 0.36
        bar_widht_slim  = 0.23

        for pre_pop_ind, pre_pop in enumerate(select_pre_pops):
            if pre_pop in long_range_pops:
                if remove_long_range_inh: bar_width = bar_widht_thick
                else: bar_width = bar_widht_slim
            else: bar_width = bar_widht_thick

            shift_distance = bar_width/2
            capsize=1.5

            for mech_ind,mech_type in enumerate(plot_bar_hist_dict.keys()):

                if   mech_type == 'exc': c = state_colors[0]
                elif mech_type == 'inh': c = state_colors[1]
                else: c = 'k'

                for state_ind,state in enumerate(states):
                    
                    # --- Plots the pops in <select_pre_pops>
                    if state == states[0]: alpha = 0.5
                    else: alpha = 1

                    if divide_plots: ax.subplot(1, 2, mech_ind+1)
                        
                    # --- Creates 4 slots for long range pops
                    if pre_pop in long_range_pops:
                        if remove_long_range_inh:
                            if   state == states[0]: displace=-shift_distance*1
                            elif state == states[1]: displace=shift_distance*1
                        else:    
                            if   mech_type == 'exc' and state == states[0]: displace=-shift_distance*3
                            elif mech_type == 'exc' and state == states[1]: displace=-shift_distance*1
                            elif mech_type == 'inh' and state == states[0]: displace=shift_distance*1
                            elif mech_type == 'inh' and state == states[1]: displace=shift_distance*3
                    # --- Creates 2 slots for long range pops
                    else:
                        if   state == states[0]: displace=-shift_distance*1
                        elif state == states[1]: displace=shift_distance*1
                        
                    # --- skip empty bars
                    if plot_bar_hist_dict[mech_type][pre_pop][state]['mean'] == 0 and plot_bar_hist_dict[mech_type][pre_pop][state]['std'] == 0: continue
                    elif math.isnan(plot_bar_hist_dict[mech_type][pre_pop][state]['mean']) or math.isnan(plot_bar_hist_dict[mech_type][pre_pop][state]['std']): continue
                    # else:
                    #     print('plotting pop ', pre_pop, state)

                    # --- Store data for SUM bar plot
                    if mech_type == 'exc' and state == states[0]:
                        sum_of_means_exc_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                        list_of_stds_exc_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                    elif mech_type == 'exc' and state == states[1]:
                        sum_of_means_exc_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                        list_of_stds_exc_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                    elif mech_type == 'inh' and state == states[0]:
                        sum_of_means_inh_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                        list_of_stds_inh_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                    elif mech_type == 'inh' and state == states[1]:
                        sum_of_means_inh_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                        list_of_stds_inh_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                    
                    # --- skip plotting long range inh 
                    if remove_long_range_inh:
                        if (pre_pop in long_range_pops) and (mech_type=='inh'): continue

                    # --- Plot data for SUM bar plot
                    x_data = pre_pop_ind

                    if divide_plots:    x_name=pre_pop
                    else:               x_name=pre_pop+'_'+mech_type
                    
                    if verticalPlot:
                        ax.bar(     x_data+displace,
                                    plot_bar_hist_dict[mech_type][pre_pop][state]['mean'],
                                    bar_width,
                                    color=c,
                                    # label=x_name,
                                    alpha=alpha,
                                    label=None,)
                        ax.errorbar(x_data+displace,
                                    plot_bar_hist_dict[mech_type][pre_pop][state]['mean'],
                                    yerr=plot_bar_hist_dict[mech_type][pre_pop][state]['std'],
                                    lolims=True,
                                    uplims=False,
                                    capsize=capsize,
                                    color=errbar_color,
                                    alpha=errbar_alpha,
                                    label=None,)
                    else:
                        ax.barh(    x_data+displace,
                                    plot_bar_hist_dict[mech_type][pre_pop][state]['mean'],
                                    bar_width,
                                    color=c,
                                    # label=x_name,
                                    alpha=alpha,
                                    label=None,)
                        ax.errorbar(plot_bar_hist_dict[mech_type][pre_pop][state]['mean'],
                                    x_data+displace,
                                    xerr=plot_bar_hist_dict[mech_type][pre_pop][state]['std'],
                                    xlolims=True,
                                    xuplims=False,
                                    capsize=capsize,
                                    color=errbar_color,
                                    alpha=errbar_alpha,
                                    label=None,)
                    if export_values: 
                        export_dict[state][mech_type][pre_pop]['mean']=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                        export_dict[state][mech_type][pre_pop]['std']=plot_bar_hist_dict[mech_type][pre_pop][state]['std']
                        export_dict[state][mech_type][pre_pop]['vals']=plot_bar_hist_dict[mech_type][pre_pop][state]['vals']

                if plotStatistics:
                    Data1 = plot_bar_hist_dict[mech_type][pre_pop][states[0]]['vals']
                    Data2 = plot_bar_hist_dict[mech_type][pre_pop][states[1]]['vals']
                    # --- skip empty bars
                    if len(Data1) <= 0 or len(Data2) <= 0: continue
                    if remove_long_range_inh:
                        if mech_type=='inh' and pre_pop in long_range_pops: continue
                    marker_position = max([plot_bar_hist_dict[mech_type][pre_pop][states[0]]['mean']+plot_bar_hist_dict[mech_type][pre_pop][states[0]]['std'],
                                        plot_bar_hist_dict[mech_type][pre_pop][states[1]]['mean']+plot_bar_hist_dict[mech_type][pre_pop][states[1]]['std']])+200
                    stats = scipy.stats.mannwhitneyu(Data1, Data2) 
                    if stats.pvalue>0.05:                           stat_symbol = ''    
                    elif stats.pvalue<=0.05 and stats.pvalue>0.01:  stat_symbol = '*'
                    elif stats.pvalue<=0.01 and stats.pvalue>0.001: stat_symbol = '**'
                    elif stats.pvalue<=0.001:                       stat_symbol = '***'
                    print('plotting ',stats.pvalue,'\t',stat_symbol,'\t',pre_pop_ind,marker_position)
                    ax.text(marker_position,pre_pop_ind+0.18,s=stat_symbol,size=10,rotation='vertical')

        if sumLeftover:
            # --- position of the sum bars in the plots
            sum_bar_position = len(select_pre_pops)+3
            leftover_bar_position = len(select_pre_pops)+1
            plot_labels = select_pre_pops+['','Others']+['','Sum']

            leftover_list=list(set(ordered_pre_pops).difference(select_pre_pops))
            for pre_pop_ind, pre_pop in enumerate(leftover_list):
                for mech_ind,mech_type in enumerate(plot_bar_hist_dict.keys()):
                    for state_ind,state in enumerate(states):
                        # --- Store data for leftover SUM bar plot
                        if mech_type == 'exc' and state == states[0]:
                            sum_of_means_exc_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            list_of_stds_exc_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])

                            leftover_sum_of_means_exc_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            leftover_list_of_stds_exc_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])

                        elif mech_type == 'exc' and state == states[1]:
                            sum_of_means_exc_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            list_of_stds_exc_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                            
                            leftover_sum_of_means_exc_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            leftover_list_of_stds_exc_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                        
                        elif mech_type == 'inh' and state == states[0]:
                            sum_of_means_inh_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            list_of_stds_inh_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                            
                            leftover_sum_of_means_inh_Q+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            leftover_list_of_stds_inh_Q.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
                        
                        elif mech_type == 'inh' and state == states[1]:
                            sum_of_means_inh_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            list_of_stds_inh_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
        
                            leftover_sum_of_means_inh_M+=plot_bar_hist_dict[mech_type][pre_pop][state]['mean']
                            leftover_list_of_stds_inh_M.append(plot_bar_hist_dict[mech_type][pre_pop][state]['std'])
        
            # --- Square of the stds - leftover
            leftover_std_exc_sqr_Q=[std**2 for std in leftover_list_of_stds_exc_Q]
            leftover_std_exc_sqr_M=[std**2 for std in leftover_list_of_stds_exc_M]
            leftover_std_inh_sqr_Q=[std**2 for std in leftover_list_of_stds_inh_Q]
            leftover_std_inh_sqr_M=[std**2 for std in leftover_list_of_stds_inh_M]
            
            leftover_sum_std_exc_Q=math.sqrt(sum(leftover_std_exc_sqr_Q))
            leftover_sum_std_exc_M=math.sqrt(sum(leftover_std_exc_sqr_M))
            leftover_sum_std_inh_Q=math.sqrt(sum(leftover_std_inh_sqr_Q))
            leftover_sum_std_inh_M=math.sqrt(sum(leftover_std_inh_sqr_M))

            positions   =[(leftover_bar_position-shift_distance*3),(leftover_bar_position-shift_distance*1),(leftover_bar_position+shift_distance*1),(leftover_bar_position+shift_distance*3)]
            sum_means   =[leftover_sum_of_means_exc_Q,leftover_sum_of_means_exc_M,leftover_sum_of_means_inh_Q,leftover_sum_of_means_inh_M]
            xerrs       =[leftover_sum_std_exc_Q,leftover_sum_std_exc_M,leftover_sum_std_inh_Q,leftover_sum_std_inh_M]
            colors      =[state_colors[0],state_colors[0],state_colors[1],state_colors[1]]
            alphas      =[0.5,1,0.5,1]

            for ind,sum_mean in enumerate(sum_means):
                if verticalPlot:
                    ax.bar(     positions[ind],
                                sum_means[ind],
                                bar_width,
                                color=colors[ind],
                                alpha=alphas[ind],
                                label=None,)
                    ax.errorbar(positions[ind],
                                sum_means[ind],
                                yerr=xerrs[ind],
                                lolims=True,
                                uplims=False,
                                capsize=capsize,
                                color=errbar_color,
                                alpha=errbar_alpha,
                                label=None,)
                else:                    
                    ax.barh(    positions[ind],
                                sum_means[ind],
                                bar_width,
                                color=colors[ind],
                                alpha=alphas[ind],
                                label=None,)
                    ax.errorbar(sum_means[ind],
                                positions[ind],
                                xerr=xerrs[ind],
                                xlolims=True,
                                xuplims=False,
                                capsize=capsize,
                                color=errbar_color,
                                alpha=errbar_alpha,
                                label=None,)
            if export_values: 
                export_dict.update({'leftover':{}})
                ind=0
                for mech_type in mech_types:
                    for state in states:
                        export_dict['leftover'].update({mech_type+state:{'mean':sum_means[ind],'std':xerrs[ind]}})
                        ind+=1
        else:
            sum_bar_position= len(select_pre_pops)+1
            plot_labels = select_pre_pops+['','Sum']
        
        # --- Square of the stds
        std_exc_sqr_Q=[std**2 for std in list_of_stds_exc_Q]
        std_exc_sqr_M=[std**2 for std in list_of_stds_exc_M]
        std_inh_sqr_Q=[std**2 for std in list_of_stds_inh_Q]
        std_inh_sqr_M=[std**2 for std in list_of_stds_inh_M]
        
        sum_std_exc_Q=math.sqrt(sum(std_exc_sqr_Q))
        sum_std_exc_M=math.sqrt(sum(std_exc_sqr_M))
        sum_std_inh_Q=math.sqrt(sum(std_inh_sqr_Q))
        sum_std_inh_M=math.sqrt(sum(std_inh_sqr_M))

        # --- parameters to draw the sum bars
        sum_bar_width = bar_widht_thick
        sum_shift_distance = sum_bar_width/2
        sum_capsize = 1.5

        positions   =[(sum_bar_position-sum_shift_distance*3),(sum_bar_position-sum_shift_distance*1),(sum_bar_position+sum_shift_distance*1),(sum_bar_position+sum_shift_distance*3)]
        sum_means   =[sum_of_means_exc_Q,sum_of_means_exc_M,sum_of_means_inh_Q,sum_of_means_inh_M]
        xerrs       =[sum_std_exc_Q,sum_std_exc_M,sum_std_inh_Q,sum_std_inh_M]
        colors      =[state_colors[0],state_colors[0],state_colors[1],state_colors[1]]
        alphas      =[0.5,1,0.5,1]
        if showLabels:  labels      =['Exc '+states[0],'Exc '+states[1],'Inh '+states[0],'Inh '+states[1]]
        else:           labels      =['','','','']


        for ind,sum_mean in enumerate(sum_means):
            if verticalPlot:
                ax.bar(     positions[ind],
                            sum_means[ind],
                            sum_bar_width,
                            # bar_width,
                            color=colors[ind],
                            alpha=alphas[ind],
                            label=labels[ind],)
                ax.errorbar(positions[ind],
                            sum_means[ind],
                            yerr=xerrs[ind],
                            lolims=True,
                            uplims=False,
                            capsize=sum_capsize,
                            # capsize=capsize,
                            color=errbar_color,
                            alpha=errbar_alpha,
                            label=None,)
            else:                    
                ax.barh(    positions[ind],
                            sum_means[ind],
                            sum_bar_width,
                            # bar_width,
                            color=colors[ind],
                            alpha=alphas[ind],
                            label=labels[ind],)
                ax.errorbar(sum_means[ind],
                            positions[ind],
                            xerr=xerrs[ind],
                            xlolims=True,
                            xuplims=False,
                            capsize=sum_capsize,
                            # capsize=capsize,
                            color=errbar_color,
                            alpha=errbar_alpha,
                            label=None,)
        if export_values: 
            export_dict.update({'sum':{}})
            ind=0
            for mech_type in mech_types:
                for state in states:
                    export_dict['sum'].update({mech_type+state:{'mean':sum_means[ind],'std':xerrs[ind]}})
                    ind+=1

        # --- Format figure
        if verticalPlot:
            ax.set_xticks(ticks=list(range(len(plot_labels))),labels=plot_labels,rotation=45)
            ax.set_xlabel('Presynaptic population')
            ax.set_ylabel('Estimated synaptic drive')   
        else:
            ax.set_yticks(ticks=list(range(len(plot_labels))),labels=plot_labels)
            # plt.gca().invert_yaxis()
            ax.invert_yaxis()
            ax.set_xlabel('Estimated synaptic drive')
            ax.set_ylabel('Presynaptic population')
            # ax.set_xticks([0,1.0,2.0])

        # --- Removing boxes from plot axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        print('bar plot method')
        # --- Changing font
        plt.rcParams.update({'font.size': 20})
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        
        if export_values: 
            # Writing to sample.json
            json_object = json.dumps(export_dict, indent=4)
            with open(export_filename, "w") as outfile: outfile.write(json_object)

class BarPlot():
    # source: https://stackoverflow.com/questions/14270391/python-matplotlib-multiple-bars
    def bar_plot(ax, data, colors=None, total_width=0.8, single_width=1, legend=True, title=None):
        """Draws a bar plot with multiple bars per data point.

        Parameters
        ----------
        ax : matplotlib.pyplot.axis
            The axis we want to draw our plot on.

        data: dictionary
            A dictionary containing the data we want to plot. Keys are the names of the
            data, the items is a list of the values.

            Example:
            data = {
                "x":[1,2,3],
                "y":[1,2,3],
                "z":[1,2,3],
            }
        colors : array-like, optional
            A list of colors which are used for the bars. If None, the colors
            will be the standard matplotlib color cyle. (default: None)

        total_width : float, optional, default: 0.8
            The width of a bar group. 0.8 means that 80% of the x-axis is covered
            by bars and 20% will be spaces between the bars.

        single_width: float, optional, default: 1
            The relative width of a single bar within a group. 1 means the bars
            will touch eachother within a group, values less than 1 will make
            these bars thinner.

        legend: bool, optional, default: True
            If this is set to true, a legend will be added to the axis.
        """
        # Check if colors where provided, otherwhise use the default color cycle
        if colors is None: colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # Number of bars per group
        n_bars = len(data)
        # The width of a single bar
        bar_width = total_width / n_bars
        # List containing handles for the drawn bars, used for the legend
        bars = []
        # Iterate over all data
        for i, (name, values) in enumerate(data.items()):
            x_offset = (i - n_bars / 2) * bar_width + bar_width / 2 # The offset in x direction of that bar
            # Draw a bar for every value of that type
            for x, y in enumerate(values):
                bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=colors[i % len(colors)])
            bars.append(bar[0]) # Add a handle to the last drawn bar, which we'll need for the legend
        # Draw legend if we need
        if legend: ax.legend(bars, data.keys())
        # Draw legend if we need
        if title: ax.set_title(title)