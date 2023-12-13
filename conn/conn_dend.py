"""
scracm.py  - RENAME TO SUBCELLULAR -- conn_dend.py?

Code to analyse sCRACM data

Contributors: salvadordura@gmail.com
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import json
import pickle as pkl
from netpyne import sim, specs
from neuron import h


def setCfg():
    return specs.SimConfig()

def setNetParams():

    """
    netParams.py -High-level specifications for M1 network model using NetPyNE
    """
    netParams = specs.NetParams()   # object of class NetParams to store the network parameters
    netParams.version = 26

    # Cell parameters
    ## PT cell params (full)
    netParams.loadCellParamsRule(label='PT_full', fileName=dataFolder+'PT_full_cellParams.pkl')
    netParams.popParams['PT5B'] =   {'cellModel': 'HH_full', 'cellType': 'PT', 'numCells': 1, 'y': [700,750]}

    ## IT5A cell params (full)
    netParams.loadCellParamsRule(label='IT5A_full', fileName=dataFolder+'IT_full_BS1578_cellParams.pkl')
    netParams.popParams['IT5A'] =   {'cellModel': 'HH_full', 'cellType': 'IT5A', 'numCells': 1, 'y': [700,750]}


    ## IT cell params (full)
    netParams.loadCellParamsRule(label='IT5B_full', fileName=dataFolder+'IT_full_BS1579_cellParams.pkl')
    netParams.popParams['IT5B'] =   {'cellModel': 'HH_full', 'cellType': 'IT5B', 'numCells': 1, 'yRange': [700,750]}

    return netParams


def loadSCRACM (lenY=30, spacing=50, scracmData=None, scracmFile=None, densFile=None):
    """ Set sCRACM params to use as 1d density map 
    lenY: grid Y dimensions
    spacing: 50  # um 
    gridY: grid locations 
    map1d: sCRACM values for each grid location
    """
    gridY = range(0, -spacing*lenY, -spacing) # NEURON's axis for cortical depth goes from 0 (pia) to -cfg.sizeY (WM)
    
    # L2/3 IT -> PT; sCRACM from Suter, 2015; use same for L2/3, L4
    
    if scracmData is None: scracmData = np.loadtxt(scracmFile)
    scracm1d = []
    for jj in range(lenY): scracm1d.append(scracmData[jj])
    scracm1d = [x-min(scracm1d) for x in scracm1d] # offset by min value (negative)

    if densFile:
        densData = np.loadtxt(densFile)
        dens1d = []
        for jj in range(lenY): dens1d.append(densData[jj])
    else:
        dens1d = None    
    return gridY, scracm1d, dens1d


def alignGrid (cell, gridY, fixedSomaY=-735):
    """ Align cell soma with grid
    cell: Cell object to perform operations on
    gridY: grid locations 
    fixedSomaY: yfrac of cell used to extract density from sCRACM; use instead of soma yfrac of each network cell 
    """
    somaLabel = next((s for s in cell.secs if s.startswith('soma')),None)
    if somaLabel:
        somaX, somaY, _ = sim.net._posFromLoc(cell.secs[somaLabel]['hSec'], 0.5) # get cell pos move method to Cell!
        gridY = [y+(somaY-fixedSomaY) for y in gridY] # adjust grid so cell soma is at fixedSomaY
    else:
        print('Cannot find soma in cell')
        sys.exit()
    return gridY


def weightNormGrid (cell, gridY, spacing):
    wnormCounter =[0] * len(gridY) 
    wnormGrid = [0.0] * len(gridY)
    for secName,sec in cell.secs.iteritems():
        if 'weightNorm' in sec:
            hSec = sec['hSec']
            for seg in hSec:
                x,y,z = sim.net._posFromLoc(hSec, seg.x)
                wnorm = sec['weightNorm'][int(round(seg.x*hSec.nseg))-1]
                for i,gridy in enumerate(gridY):
                    if y >= gridy and y < gridy+spacing:
                        wnormCounter[i] += 1
                        wnormGrid[i] += wnorm
        else:
            print(secName+ ' has no weightNorm')
    wnormGrid = [wn/wc if wc>0 else 0 for wn,wc in zip(wnormGrid,wnormCounter)]
    return wnormGrid


def dendLengthGrid (cell, gridY, spacing):
    dendLGrid = [0.0] * len(gridY)
    for secName, sec in cell.secs.iteritems():
        if secName not in ['axon']:
            hSec = sec['hSec']
            for seg in hSec:
                x,y,z = sim.net._posFromLoc(hSec, seg.x)
                segL = hSec.L / hSec.nseg
                for i,gridy in enumerate(gridY):
                    if y >= gridy and y < gridy+spacing:
                        dendLGrid[i] += segL
    return dendLGrid
    


def calcSynDensity (cellType=None, fixedSomaY=None, dendLGrid=None, lenY=30, spacing=50, scracmData=None, scracmFile=None, densFile=None, includeWnorm=True, show=True, save=''):
    # create net
    #cfg, netParams = sim.readCmdLineArgs()
    sim.initialize(simConfig = setCfg(), netParams = setNetParams())  # create network object and set cfg and net params
    sim.net.createPops()                        # instantiate network populations
    sim.net.createCells()                       # instantiate network cells based on defined populations

    # load sCRACM data and set grid
    gridY, scracm1d, dens1d = loadSCRACM(lenY, spacing, scracmData, scracmFile, densFile)

    # align grid with cell soma
    cell = next((c for c in sim.net.cells if c.tags['cellType']==cellType),None)
    #print cellType,cell
    if cellType: gridY = alignGrid(cell, gridY, fixedSomaY = fixedSomaY)

    # calculate weightNorm avg and dend length for each grid position
    if includeWnorm: wnormGrid = weightNormGrid (cell, gridY, spacing)
    if dendLGrid == None:
        dendLGrid = dendLengthGrid (cell, gridY, spacing)

    # normalize
    if includeWnorm: wnormGridNorm = [x/max(wnormGrid) for x in wnormGrid]
    dendLGridNorm = [x/max(dendLGrid) for x in dendLGrid]
    scracm1dNorm = [x/max(scracm1d) for x in scracm1d]
    if densFile: dens1dNorm = [x/max(dens1d) for x in dens1d]

    # calculate syn density = scracm / (wnorm* dendL)
    if includeWnorm: 
        dens1dEstim = [scracm / (wnorm+dendL) for scracm,wnorm,dendL in zip(scracm1dNorm, wnormGridNorm, dendLGridNorm)]
    else:
        dens1dEstim = [scracm / np.sqrt(dendL) for scracm,dendL in zip(scracm1dNorm, dendLGridNorm)]

    dens1dEstim = [x if x not in [np.inf, -np.inf] else 0.0 for x in dens1dEstim]
    dens1dEstim = [x if not np.isnan(x) else 0.0 for x in dens1dEstim]
    dens1dEstimNorm = [x/max(dens1dEstim) for x in dens1dEstim]

    # plot dens1d vs dens1dEstimate
    lw = 2.0
    plt.figure(figsize=(12,6))
    plt.plot(scracm1dNorm, 'y', linewidth=lw, label='scracm')
    plt.plot(dendLGridNorm,  linewidth=lw, label='dendL')
    if includeWnorm: plt.plot(wnormGridNorm, 'k', linewidth=lw, label='wnorm')
    if densFile: plt.plot(dens1dNorm, 'g', linewidth=lw, label='syn density')
    plt.plot(dens1dEstimNorm, 'r',  linewidth=lw, label='syn density (aprox)')
    plt.xlabel('normalized cortical depth ')
    plt.ylabel('normalized strength')
    plt.legend()   
    if save: plt.savefig(save) 
    if show: plt.show()
    
    if densFile: return gridY, dens1dNorm
    else: return gridY, dens1dEstimNorm


def extractColor(colorbarFile, figFile, figDim, checkColorRange=True, show=False, save=False):
    from PIL import Image
    cbar_image = Image.open(colorbarFile) #Can be many different formats.
    cbar = cbar_image.load()
    xlen,ylen=cbar_image.size #Get the width and hight of the image for iterating over
    cbarColors = []
    for y in range(ylen):
        cbarColors.append(cbar[0,y])
    cbarlen = len(cbarColors)
    # cbarmin = min([np.mean(c) for c in cbarColors])
    # cbarmax = max([np.mean(c) for c in cbarColors])

    print('Extracting figure values...')
    fig_image = Image.open(figFile) #Can be many different formats.
    fig = fig_image.load()
    xlen,ylen=fig_image.size #Get the width and hight of the image for iterating over
    print(xlen,ylen, figDim[0], figDim[1])
    xstep = float(xlen) / figDim[0]
    ystep = float(ylen) / figDim[1]
    
    figVals = np.zeros(figDim)
    for (i,xf) in enumerate(np.arange((xstep/4.0),xlen,xstep)):
        for (j,yf) in enumerate(np.arange((ystep/4.0),ylen,ystep)):
            x=int(round(xf))
            y=int(round(yf))
            print (i,j,x,y)
            if x > xlen or y > ylen: continue
            color = fig[x, y][0:3] # if have alpha channel will be 4d
            
            dist = [np.linalg.norm(np.array(color)-np.array(cbcol)) for cbcol in cbarColors]
            minDist = min(dist)
            if checkColorRange:
                incx,incy = 0,0
                #while (np.mean(color) <= cbarmin) or (np.mean(color) >= cbarmax) and inc <= ystep/2.5:
                while minDist > 15.0 and incx <= 0.5*xstep:
                    incy += 2
                    if incy >= 0.5*ystep: 
                        incy = 0
                        incx += 2
                    #print (minDist, x + incx, y + incy)
                    try:
                        color = fig[x+incx,y+incy][0:3]  # if have alpha channel will be 4-d
                    except:
                        pass
                    dist = [np.linalg.norm(np.array(color)-np.array(cbcol)) for cbcol in cbarColors]
                    minDist = min([float(minDist) , min(dist)]) 
            
            index = np.argmin(dist) 
            value = 1.0 - float(index)/float(cbarlen)
            print(value)
            figVals[i][j] = value

    if save:
        figVals.dump(figFile[:-4]+'.pkl')
    if show:
        plt.figure()
        plt.imshow(figVals.T, interpolation='none')
        plt.show()

    return figVals


def synDensPT():
    #########################################################
    #L2/3,TVL,S2,cM1,M2 -> PT (Suter, 2015)
    #########################################################
    ylen = 30  # number of y pixels
    spacing = 50  # size of pixels (um) 
    frac = [-0.3, 0.2]  # fraction of cells to use from each scracm map (relative to center eg. if n=10, frac=[-0.15,0.15] -> [4,5,6])
    fixedSomaY = -735
    gridY = {}  # dict to store the Y grid locations
    synDens = {}  # dict to store syn density along dendrite (extracted from sCRACM data)
    scracmParams = {}
    scracmParams['L2_PT'] = {'file': 'L23_PT_scracm.pkl', 'n': 23}  # M1 L2/3-> L5B PT (Suter 2015) -- use also for L4,L5A->PT ? no
    scracmParams['M2_PT']  = {'file': 'M2_PT_scracm.pkl',  'n': 21}   # M2 -> L5B PT (Suter 2015) 
    scracmParams['TVL_PT'] = {'file': 'TVL_PT_scracm.pkl', 'n': 31}  # TVL -> L5B PT (Suter 2015)
    scracmParams['cM1_PT'] = {'file': 'cM1_PT_scracm.pkl', 'n': 24}  # cM1 -> L5B PT (Suter 2015)
    scracmParams['S2_PT']  = {'file': 'S2_PT_scracm.pkl',  'n': 31}   # S2 -> L5B PT (Suter 2015) -- use also for S1->PT ? yes
    
    # L23 -> PT files
    densFile = rootFolder+'data/csp_syn_input_v6/radial_scracm18_BS0284_memb_BS0477_morph.dat'  # using sfn16 poster
    scracmFile = rootFolder+'data/csp_syn_input_v6/exp_radial_scracm18_BS0284_memb_BS0477_morph.dat'  # from raw data (not extracted from fig)
    
    # create and format fig
    fig1 = plt.figure(figsize=(6,10))
    ax1 = fig1.add_subplot(111)

    # calculate synDens from scracm for each projection type
    for i,(label,params) in enumerate(scracmParams.iteritems()):
        n = params['n']
        #data = extractColor(dataFolder+'colorbar_Suter15.png', dataFolder+params['file'], (n, ylen), checkColorRange=True, show=False, save=True)
        with open(dataFolder+params['file'],'r') as f: data=pkl.load(f)
        crange = [int(round((n/2.0)+(x*n))) for x in frac]
        data1d = np.mean(data[:][crange[0]:crange[1]], axis=0) 
        gridY[label], synDens[label] = calcSynDensity(cellType='PT', fixedSomaY=fixedSomaY, lenY=ylen, spacing=spacing, scracmData=data1d, densFile=None, includeWnorm=False, show=False, save=None)#dataFolder+label+'_dens1d.png')
        #ax1.plot(data1d, color=colorlist[i], linewidth=lw, label=label)
        ax1.plot(synDens[label], np.linspace(0, 1, len(synDens[label])), '-', color=colorlist[i],  linewidth=lw, label=label)

    gridY['S1_PT'], synDens['S1_PT'] = gridY['S2_PT'], synDens['S2_PT']  # use S2 data for S1 as well (very similar)

    dataSave = {'synDens': synDens, 'gridY': gridY, 'fixedSomaY': fixedSomaY}
    with open(outFolder+'conn_dend_PT.json','w') as f: json.dump(dataSave, f)
    fontsiz=20
    ax1.set_ylim([0, 0.9])
    ax1.invert_yaxis()
    ax1.set_ylabel('Normalized cortical depth (NCD) aligned to soma', fontsize=fontsiz)
    ax1.set_xlabel('Normalized synaptic density', fontsize=fontsiz)
    ax1.legend([p.replace('_', '->').replace('L2', 'L2/3') for p in scracmParams],fontsize=fontsiz)
    fig1.savefig(dataFolder+'long_PT.png')


def synDensIT():
    #########################################################
    # TPO, TVL, M2, OC  -> E (L2/3, L5A, L5B, L6) (Hooks 2013)
    #########################################################    
    ylen = 26  # number of y pixels
    spacing = 50  # size of pixels (um) 
    fixedSomaY = -700
    gridY = {}  # dict to store the Y grid locations
    synDens = {}  # dict to store syn density along dendrite (extracted from sCRACM data)
    scracm1d = {}

    # load dendritic density for L4 IT cells (Yamawaki 2015 elife - fig 8D)
    with open(dataFolder+'dendDens_L4.json','r') as f: dendDens_L4 = json.load(f)
    somaY_scracm = -fixedSomaY
    somaY_dens = 450
    alignDist = (somaY_scracm-somaY_dens) / spacing
    dendDens_L4 = [0] * alignDist + dendDens_L4  # insert 0s at beginning
    dendDens_L4 = dendDens_L4[:ylen]  # remove extra elements at end

    longPops = ['TPO', 'TVL', 'M2', 'OC']
    
    for longPop in longPops:
        # create and format fig
        fig1 = plt.figure(figsize=(6,10))
        ax1 = fig1.add_subplot(111)

        # extract scracm from fig 
        #data = extractColor(dataFolder+'colorbar_Hooks13.png', dataFolder+longPop+'_E_scracm.png', (4, ylen), checkColorRange=True, show=False, save=True)
        with open(dataFolder+longPop+'_E_scracm.pkl','r') as f: data=pkl.load(f)

        # parameters: ct=cell type, sy: soma y location; dd: dendritic density
        params = [{'label': longPop+'_L2','ct': None, 'sy': fixedSomaY, 'dd': dendDens_L4}, # dendritic density from Yama15 elife (L4)
                  {'label': longPop+'_L5A', 'ct': 'IT5A', 'sy': fixedSomaY, 'dd': None},
                  {'label': longPop+'_L5B', 'ct': 'IT5B', 'sy': fixedSomaY, 'dd': None},
                  {'label': longPop+'_L6', 'ct': 'PT', 'sy': fixedSomaY, 'dd': None}]

        # calculate synDens from scracm for each projection type
        for i,param in enumerate(params):
            label = param['label']
            scracm1d[label] = data[i]
            gridY[label], synDens[label] = calcSynDensity(cellType=param['ct'], fixedSomaY=param['sy'], dendLGrid=param['dd'], 
                                            lenY=ylen, spacing=50, scracmData=scracm1d[label], densFile=None, includeWnorm=False, show=False, save=None)#dataFolder+label+'_dens1d.png')
            #ax1.plot(data[i], color=colorlist[i], linewidth=lw, label=label)
            ax1.plot(synDens[label], np.linspace(0, 1, len(synDens[label])), '-', color=colorlist[i],  linewidth=lw, label=label)

        dataSave = {'synDens': synDens, 'gridY': gridY, 'fixedSomaY': fixedSomaY}
        with open(outFolder+'conn_dend_IT.json','w') as f: json.dump(dataSave, f)
        ax1.invert_yaxis()
        ax1.set_ylabel('Normalized cortical depth (NCD) aligned to soma')
        ax1.set_xlabel('Normalized synaptic density')
        ax1.legend([p['label'].replace('_', '->') for p in params])
        fig1.savefig(dataFolder+longPop+'_E.png')


# main code
if __name__ == '__main__':

    plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    plt.rcParams.update({'font.size': 16})

    # set folder paths for source data and output data
    rootFolder = '../../'  # should point to root repo folder (m1) -- currently in m1/sim/conn  
    dataFolder = rootFolder+'data/conn/'
    outFolder = rootFolder+'sim/conn/'

    # NUM_COLORS = 8
    # colormap = plt.get_cmap('gist_rainbow')
    # colorlist = [colormap(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    # colorlist=[[0.42,0.67,0.84], [0.90,0.76,0.00], [0.42,0.83,0.59], [0.90,0.32,0.00],
    #     [0.34,0.67,0.67], [0.90,0.59,0.00], [0.42,0.82,0.83], [1.00,0.85,0.00],
    #     [0.33,0.67,0.47], [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
    #     [0.71,0.82,0.41], [0.0,0.2,0.5], [0.70,0.32,0.10]]
    
    colorlist = [[0.90,0.32,0.00], [0.42,0.83,0.59], [0.90,0.59,0.00], [0.90,0.76,0.00],[0.34,0.67,0.67], [0.42,0.67,0.9], [0.8, 0.25, 1.0]]

    lw = 2.0

    # L2/3,TVL,S2,cM1,M2 -> PT (Suter, 2015)
    # synDensPT()
    
    # TVL, TPO, M2, OC -> E (L2/3, L5A, L5B, L6) (Hooks 2013)
    #synDensIT()

    # rest of E -> E = spiny/all dends (no scracm)

    # PV -> IT/PT = perisomatic (no scracm) (Naka16)

    # SOM -> IT/PT = apical dendrites (no scracm) (Naka16)

    # E/I -> I = spiny/all dends (no scracm); 

    # ----------------------------------
    # for netpyne paper: fig 3 - VL -> PT (actually fig 6A)
    ylen = 25
    n = 10
    data = extractColor(dataFolder+'colorbar_Suter15.png', dataFolder+'Suter15_fig6A_v2.png', (n, ylen), checkColorRange=True, show=1, save=True)
    with open(dataFolder+'Suter15_fig6A_v2.pkl','rb') as f:
        d=pkl.load(f)
    for i in range(14, 21):
        d[4, i] = d[5, i]
    d[4, 15] = d[5, 16]

    plt.imshow(d.T, interpolation='none')
    plt.savefig(dataFolder+'Suter15_fig6A_extracted.png')
    
    # ----------------------------------
    # for M1 paper: S2 -> PT
    # ylen = 25
    # n = 10
    # data = extractColor(dataFolder+'colorbar_Suter15.png', dataFolder+'Suter15_fig5D_v2.png', (n, ylen), checkColorRange=True, show=1, save=True)
    # with open(dataFolder+'Suter15_fig5D.pkl','rb') as f:
    #     d=pkl.load(f)
    # # for i in range(12, 17):
    # #      d[5, i] = d[4, i]
    # # d[5, 17] = d[5, 16]
    # # d[6, 14] = d[5, 16]
    # # d[6, 15] = d[5, 16]

    # plt.imshow(d.T, interpolation='none')
    # plt.savefig(dataFolder+'Suter15_fig5D_S2_PT_scracm_extracted_v2.png')

