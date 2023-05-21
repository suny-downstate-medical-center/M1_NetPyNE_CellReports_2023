"""
M1 paper Figure 4

contributors: salvadordura@gmail.com
"""

from matplotlib import pyplot as plt
import os, collections

import IPython as ipy
import pickle as pkl
import json
import numpy as np
import scipy.signal as ss
from netpyne import specs,sim

from shared import *


with open('../../sim/cells/popColors.pkl', 'rb') as fileObj: popColors = pkl.load(fileObj)['popColors']

def loadFile(filename, include):

    if filename.endswith('.json'):
        with open(filename, 'rb') as fileObj:
            data = json.load(fileObj, object_pairs_hook=collections.OrderedDict)
    elif filename.endswith('.pkl'):
        with open(filename, 'rb') as fileObj:
            data = pkl.load(fileObj)

    cfg = specs.SimConfig(data['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data, instantiate=False)
    sim.setSimCfg(cfg)
    
    temp = collections.OrderedDict()

    # order pops by keys in popColors
    for k in include:
        if k in sim.net.params.popParams:
            temp[k] = sim.net.params.popParams.pop(k)

    # add remaining pops at the end
    for k in list(sim.net.params.popParams.keys()):
        temp[k] = sim.net.params.popParams.pop(k)
    
    sim.net.params.popParams = temp
    
    try:
        print('Cells created: ',len(sim.net.allCells))
    except:
        #import IPython; IPython.embed()
        sim.net.createPops()     
        sim.net.createCells()
        sim.setupRecording()
        sim.cfg.createNEURONObj = False
        sim.gatherData() 

    sim.allSimData = data['simData']



###########################
######## MAIN CODE ########
###########################

if __name__ == '__main__':

    filenames = ['../../data/v56_batch39/v56_batch39_0_0_0_0_data.pkl', 
                '../../data/v56_batch38/v56_batch38_0_0_0_0_data.pkl']

    timeRanges = [[1000,5000], [5000, 9000]]
    periodLabels = ['quiet', 'move']
    freqRanges = [[0,4], [30,80]]

    individualPlots = False
    combinedPlots = True
    
    allpops = ['IT2', 'SOM2','PV2','IT4','IT5A','SOM5A','PV5A','IT5B','PT5B','SOM5B','PV5B','IT6','CT6','SOM6','PV6']
    excpops = ['IT2','IT4','IT5A','IT5B','PT5B','SOM5B','PV5B','IT6','CT6']

    loadProcessedData = 1 #False

    fontsiz = 16

    for filename, periodLabel, timeRange, freqRange in zip(filenames, periodLabels, timeRanges, freqRanges):

        if not loadProcessedData:
            loadFile(filename, include=allpops)

            allDataLFP = []
            allDataPSD = []

            # LFP 
            if individualPlots:
                # Overall
                out = sim.plotting.plotLFPTimeSeries( 
                        electrodes = [6], #['avg']+list(range(11)), 
                        timeRange = timeRange, 
                        figSize = (16,12), 
                        rcParams = {'font.size': 20},
                        saveData = False, 
                        saveFig = True,
                        fileName = filename[ :-4]+'_LFP_timeSeries_allpops_%d_%d.png'%(timeRange[0], timeRange[1]), 
                        showFig = False)

                
                out = sim.plotting.plotLFPPSD(
                        plots = ['PSD'], 
                        electrodes = [6],
                        timeRange = timeRange, 
                        minFreq = 0.05,
                        maxFreq = 80,
                        stepFreq = 0.05,
                        orderInverse = False, 
                        figSize = (8,6), 
                        rcParams = {'font.size': 20},
                        saveData = False, 
                        saveFig = True,
                        fileName = filename[:-4]+'_LFP_PSD_allpops_%d_%d.png'%(timeRange[0], timeRange[1]), 
                        showFig = False)


                # By population
                for pop in ['IT2','IT4','IT5A','IT5B','PT5B','SOM5B','PV5B','IT6','CT6']:
                    out = sim.plotting.plotLFPTimeSeries( 
                            electrodes = [6], #['avg']+list(range(11)), 
                            pop = pop,
                            timeRange = timeRange, 
                            figSize = (16,12), 

                            rcParams = {'font.size': 20},
                            saveData = False, 
                            saveFig = True,
                            fileName = filename[ :-4]+'_LFP_timeSeries_%s_%d_%d.png'%(pop, timeRange[0], timeRange[1]), 
                            showFig = False)
                    
                    out = sim.plotting.plotLFPPSD(
                            plots = ['PSD'], 
                            electrodes = [6],
                            pop = pop,
                            timeRange = timeRange, 
                            minFreq = 0.05,
                            maxFreq = 80,
                            stepFreq = 0.05,
                            orderInverse = False, 
                            figSize = (8,6), 
                            rcParams = {'font.size': 20},
                            saveData = False, 
                            saveFig = True,
                            fileName = filename[:-4]+'_LFP_PSD_%s_%d_%d.png'%(pop, timeRange[0], timeRange[1]), 
                            showFig = False)

            if combinedPlots:
                # Overall
                dataLFP = sim.analysis.prepareLFP( 
                        electrodes = [6], #['avg']+list(range(11)), 
                        timeRange = timeRange, 
                        filtFreq = 200,
                        figSize = (16,12), 
                        rcParams = {'font.size': 20},
                        saveData = False, 
                        saveFig = True,
                        fileName = filename[ :-4]+'_LFP_timeSeries_allpops_%d_%d.png'%(timeRange[0], timeRange[1]), 
                        showFig = False)
                
                dataPSD = sim.analysis.preparePSD(
                        electrodes = [6],
                        timeRange = timeRange, 
                        minFreq = 0.05,
                        maxFreq = 80,
                        stepFreq = 0.05,
                        orderInverse = False, 
                        figSize = (8,6), 
                        rcParams = {'font.size': 20},
                        saveData = False, 
                        saveFig = True,
                        fileName = filename[:-4]+'_LFP_PSD_allpops_%d_%d.png'%(timeRange[0], timeRange[1]), 
                        showFig = False)

                allDataLFP.append(dataLFP)
                allDataPSD.append(dataPSD)


                # By population
                for pop in allpops:
                    dataLFP = sim.analysis.prepareLFP( 
                            electrodes = [6], #['avg']+list(range(11)), 
                            pop = pop,
                            timeRange = timeRange, 
                            filtFreq = 200,
                            figSize = (16,12), 
                            rcParams = {'font.size': 20},
                            saveData = False, 
                            saveFig = True,
                            fileName = filename[ :-4]+'_LFP_timeSeries_%s_%d_%d.png'%(pop, timeRange[0], timeRange[1]), 
                            showFig = False)
                    
                    dataPSD = sim.analysis.preparePSD(
                            plots = ['PSD'], 
                            electrodes = [6],
                            pop = pop,
                            timeRange = timeRange, 
                            minFreq = 0.05,
                            maxFreq = 80,
                            stepFreq = 0.05,
                            orderInverse = False, 
                            figSize = (8,6), 
                            rcParams = {'font.size': 20},
                            saveData = False, 
                            saveFig = True,
                            fileName = filename[:-4]+'_LFP_PSD_%s_%d_%d.png'%(pop, timeRange[0], timeRange[1]), 
                            showFig = False)

                    allDataLFP.append(dataLFP) 
                    allDataPSD.append(dataPSD)
            
            # save processed data to file
            with open(filename[ :-4]+'_processed_data.pkl', 'wb') as f:
                pickle.dump([allDataLFP, allDataPSD] ,f)

        # load processed data from file
        else:
            
            with open(filename[ :-4]+'_processed_data.pkl', 'rb') as f:
                [allDataLFP, allDataPSD] = pickle.load(f)


        # calculate pops contributing most to psd
        freqScale = 20
        freqPeaks = [np.sum(x['psdSignal'][0][freqRange[0]*freqScale:freqRange[1]*freqScale]) for x in allDataPSD]
        topPopIndices = np.argsort(freqPeaks)[::-1]
        topPopIndices = topPopIndices[:5]
        

        # plotting
        popLabels = ['All']+allpops #allpops
        popColors['All'] = 'black'

        plt.figure(figsize=(8,4))
        plt.ion()
        fs = 1000/0.025

        for i in topPopIndices:
            # combined PSD 
            # using netpyne PSD
            dataNorm = allDataPSD[i]['psdSignal'][0] #/ np.max(allDataPSD[0]['psdSignal'][0])
            f = allDataPSD[i]['psdFreqs'][0]
            plt.plot(f, dataNorm*1000, label=popLabels[i],  color=popColors[popLabels[i]], linewidth=2)


        # log x and log y (for wavelet/netpyne PSD)
        ax = plt.gca()
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.xlim(1,80) 
        plt.ylim([.4e-3,4e0])

        xstep = 20
        xrange = np.arange(xstep, np.max(f)+1, xstep)
        xticks = [1, 4, 8, 12, 30, 80]# [3, 9, 28, 80]
        ax.set_xticks(xticks)
        ax.set_xticklabels(['%.0f' % x for x in xticks])

        ax = plt.gca()
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(),fontsize=fontsiz)
        plt.xlabel('Frequency (Hz)', fontsize=fontsiz)
        plt.ylabel('LFP power ($\mu$V$^2$/Hz)', fontsize=fontsiz) #'PSD [V**2/Hz]'
        plt.legend(fontsize=fontsiz, loc='upper left', bbox_to_anchor=(1.01, 1))
        plt.subplots_adjust(bottom=0.15, top=0.95 , left=0.1, right=0.82 )

        plt.savefig(filename[:-4]+'_LFP_PSD_%s_combined_Welch.png' % (periodLabel), dpi=300)

        plotEachLFP = 0 

        if plotEachLFP:

            for i in topPopIndices:
                # individual LFP timeseries
                plt.figure(figsize=(8,4))
                lw=0.5
                t = allDataLFP[i]['t']
                plt.plot(t, -allDataLFP[i]['electrodes']['lfps'][0], color=popColors[popLabels[i]], linewidth=lw)
                ax = plt.gca()        
                ax.invert_yaxis()
                plt.axis('off')
                plt.xlabel('time (ms)', fontsize=fontsiz)

                meanSignal = np.mean(-allDataLFP[i]['electrodes']['lfps'][0])
                plt.ylim(meanSignal+0.4,meanSignal-0.5)
                plt.subplots_adjust(bottom=0.0, top=0.9, left=0.1, right=0.9)

                # calculate scalebar size and add scalebar
                scaley = 1000.0  # values in mV but want to convert to uV
                sizey = 100/scaley
                labely = '%.3g $\mu$V'%(sizey*scaley)#)[1:]
                sizex=500
                add_scalebar(ax, hidey=True, matchy=False, hidex=True, matchx=True, sizex=None, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, 
                        unitsx='ms', loc='upper right', pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
                
                plt.title('LFP 0-200 Hz', fontsize=fontsiz, fontweight='bold')
                plt.savefig(filename[:-4]+'_LFP_timeSignal_%s_%s.png' % (popLabels[i], periodLabel),dpi=300)
            
        
