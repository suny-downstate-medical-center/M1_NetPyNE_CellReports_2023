"""
paper.py 

Paper figures

Contributors: salvadordura@gmail.com
"""

import utils

import json
import numpy as np
import scipy
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sb
import os
import pickle
from netpyne.support.scalebar import add_scalebar
from matplotlib import cm
from matplotlib.colors import ListedColormap

# ---------------------------------------------------------------------------------------------------------------
# Population params
allpops = ['IT2','SOM2','PV2','IT4','IT5A','SOM5A','PV5A','IT5B','PT5B','SOM5B','PV5B',
    'IT6','CT6','SOM6','PV6']
excpops = ['IT2','IT4','IT5A','IT5B','PT5B','IT6','CT6']
inhpops = ['SOM2','PV2', 'SOM5A','PV5A', 'SOM5B','PV5B',  'SOM6', 'PV6']
excpopsu = ['IT2','IT4','IT5A','PT5B']
ITsubset = ('IT2', 'IT4', 'IT5A', 'IT6')
SOMsubset = ('SOM2', 'SOM5A','SOM5B', 'SOM6')
PVsubset = ('PV2', 'PV5A', 'PV5B', 'PV6')
try:
    with open('../sim/cells/popColors.pkl', 'rb') as fileObj: 
        popColors = pickle.load(fileObj)['popColors']
    popColors['S2'] = [0.90,0.76,0.00]
    popColors['TPO'] = [52/255.0, 138/255.0, 49/255.0] #[232/255.0, 37/255.0, 101/255.0] #'firebrick' #popColors['S2'] #[253/255.0, 102/255.0, 2/255.0] #
    popColors['M2'] = [0.42,0.67,0.9]
    popColors[ITsubset] = popColors['IT5A']
    popColors[SOMsubset] = popColors['SOM L2-6'] = popColors['SOM2-6'] = popColors['SOM5A']
    popColors[PVsubset] = popColors['PV L2-6'] = popColors['PV2-6'] = popColors['PV5A']
except:
    pass

def loadSimData(dataFolder, batchLabel, simLabel):
    ''' load sim file'''
    root = dataFolder+batchLabel+'/'
    sim,data,out = None, None, None
    if isinstance(simLabel, str):
        if os.path.exists(root+simLabel+'.json'):
            filename = root + simLabel + '.json'
            print(filename)
        elif os.path.exists(root+simLabel+'.pkl'):
            filename = root + simLabel + '.pkl'
            print(filename)

        sim, data, out = utils.plotsFromFile(filename, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0,
                                            include = {'raster': allpops},
                                            popColors=popColors)
    
    return sim, data, out, root

def axisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 

def loadPopGids(dataFolder, popLabels):
    # load all data    
    # sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
    from netpyne import sim
    with open(dataFolder + 'v56_batch2_net/v56_batch2_net_0_0.pkl', 'rb') as f:
        netLoad = pickle.load(f)

    allCells = netLoad['net']['cells']

    L5Bmin=0.47
    L5Bmax=0.8
    L5Bmid = L5Bmin + (L5Bmax - L5Bmin) / 2
    
    # set gids for each pop
    popGids =  {popLabel: [c['gid'] for c in allCells if c['tags']['pop'] == popLabel] for popLabel in popLabels if isinstance(popLabel, str)}
    popGids['upperPT5B'] = tuple([c['gid'] for c in allCells \
        if L5Bmin <= c['tags']['ynorm'] <= L5Bmid and c['tags']['pop']=='PT5B'])
    popGids['lowerPT5B'] = tuple([c['gid'] for c in allCells \
        if L5Bmid <= c['tags']['ynorm'] <= L5Bmax and c['tags']['pop']=='PT5B'])
    popGids['L5B'] = popGids['IT5B'] + popGids['PT5B']

    popNumCells = {popLabel: len(popGids[popLabel]) for popLabel in popGids}

    return popGids, popNumCells