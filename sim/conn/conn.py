'''

Script to combine M1 experimental connectivity data from multiple papers into single conn matrix for model

smat - strength matrix; smat = pmat * imat
pmat - probability matrix
wmat - weight matrix

'''

import numpy as np
from scipy.io import loadmat, savemat
from scipy import interpolate
#from pylab import *
from matplotlib import pyplot as plt
import sys
import sys

def loadData():
    smat = {}  # dict for strength matrices (smat = pmat * wmat)
    pmat = {}  # dict for probability matrices
    wmat = {}  # dict for weight matrices
    bins = {}  # dict with bin bins

    ''' General coordinate system / layer boundaries '''
    bins['layerLabels'] =   ['pia', 'L1-L2 border', 'L2/3-L4 border',   'L4-L5A border',   'L5A-L5B border',    'L5B-L6 border',   'L6-WM border']  
    bins['layers']      =   [0,     0.1,           0.29,               0.37,               0.47,               0.8,               1.0]


    ''' Weiler, 2008
    Glutamate LSPS maps of unlabeled pyramidal neurons in mouse somatic M1
    Weiler, N., L. Wood, J. Yu, S. A. Solla, and G. M. G. Shepherd (2008, February). Top-down laminar organization of the excitatory network in motor cortex. Nature Neuroscience 11 (3), 360-366.
    http://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=114655
    N=102 neurons, ~11 per each of 9 yfrac-bins
    Maps were 16x16 with 100 um spacing, flattened to 16x1, then binned into postsynaptic yfrac bins 10x1 with 140 um spacing
    This version has raw pA values, it has not been normalized. 
    This matrix was recalculated by Ben on 9/28/2015.   
    It has not been verified vs original study, for absolute value. 
    If normalized to peak (35.9944), it is similar but not exactly the same as the above normalized matrix. 
    '''
    smat['W'] = loadmat(dataFolder+'WeilerMat.mat')['weilerMat_raw_pA']
    binSize = 0.10
    bins['W'] = np.array([np.arange(binSize, 1.0, binSize), np.arange(2*binSize, 1.0+binSize, binSize)]).T  # bin edges
    # add 0-0.1 bin with row of 0s
    smat['W'] = np.vstack((np.array([0]*len(smat['W'])), smat['W']))
    smat['W'] = np.array([[0.0]+list(smat['W'][i]) for i in range(len(smat['W']))])
    bins['W'] = np.vstack((np.array([0,0.1]), bins['W']))

    ''' Yamawaki et al, 2015
    Local inputs to M1-L4 (Figure 5)
    Values extracted from figure using same bin centers as Weiler 
    '''
    smat['Y'] = np.array([0.218, 0.286, 0.447, 0.612, 0.380, 0.286, 0.121, 0.072, 0.075, 0.020])
    smat['Y'][0] = 0.0  # set 0-0.1 bin to 0


    ''' Anderson & Sheets, 2011 
    Glutamate LSPS maps onto CSP in mouse somatic M1
    Anderson, C. T., P. L. Sheets, T. Kiritani, and G. M. Shepherd (2010, June). Sublayer-specific microcircuits of corticospinal and corticostriatal neurons in motor cortex. Nature neuroscience 13 (6), 739-744.
    CSP (PT): N=44 neurons in L5B
    CSTR (IT): N=42 neurons in L5A and L5B
    Maps were 16x16 with 100 um spacing. Somehow converted to 16x16 connectivity matrix, not sure how interpolated/binned to get there.
    '''
    smat['AS_IT'] = loadmat(dataFolder+'andersonCSTRMat.mat')['andersonCSTRMat']
    smat['AS_PT'] = loadmat(dataFolder+'andersonCSPMat.mat')['andersonCSPMat']
    binSize = 0.0625
    bins['AS'] = np.array([np.arange(0.0, 1.0, binSize), np.arange(binSize, 1.0+binSize, binSize)]).T   # bin edges


    ''' Load Lefort, 2009
    Excitatory Microcircuits of C2 Barrel Column
    S. Lefort, C. Tomm, J.C. Floyd Sarria, C.C.H. Petersen (2009).The Excitatory Neuronal Network of the C2 Barrel Column in Mouse Primary Somatosensory Cortex. Neuron 61 301--316
    Use to constrain p_con vs i_con
    Note: L2 connecitvity different for L2A vs L2B; plt.figure 5D also shows conn in 50um bins but data not available
    subpial layer boundaries: L1, 128 +- 1 mm; L2, 269 +- 2mm; L3,418 +- 3mm; L4,588 +- 3mm; L5A,708 +- 4mm; L5B, 890 +- 5 mm; L6, 1154 +- 7 mm.
    '''
    pmat['L'] = loadmat(dataFolder+'LefortMat.mat')['pmat']
    wmat['L'] = loadmat(dataFolder+'LefortMat.mat')['wmat']
    smat['L'] = pmat['L'] * wmat['L']
    bins['L'] = [[128.0/1154,269.0/1154], [269.0/1154,418.0/1154], [418.0/1154,588.0/1154], 
                [588.0/1154,708.0/1154], [708.0/1154,890.0/1154], [890.0/1154,1154.0/1154]]

    
    ''' Kiritani, 2012  
    pmat as fraction of cells from from paired recordings; wmat in mV
    Corticospinal-Corticostriatal connectivity
    T. Kiritani, I.R. Wickersham, H. S. Seung, and G. M. G. Shepherd (2012). Hierarchical connectivity and connection-specific dynamics in the corticospinal-corticostriatal microcircuit in mouse motor cortex. J. Neurosci. 32(14), 4992-5001.
    Use to obtain data for L5 diagonals in connectivity matrix

    bin labels: IT L5, PT L5
    '''
    pmat['K'] = loadmat(dataFolder+'KiritaniMat.mat')['pmat']
    wmat['K'] = loadmat(dataFolder+'KiritaniMat.mat')['wmat']
    smat['K'] = pmat['K'] * wmat['K']
    bins['K'] = [[0.3125, 0.75], [0.4375, 0.8125]]

    return smat, pmat, wmat, bins


def getgridfrombins (bins):
    ngrid = []
    for y in bins:
        for x in bins:
            ngrid.append((x,y))
    return ngrid

def getpts (mat):
    pts = []
    for y in range(mat.shape[0]):
        for x in range(mat.shape[1]):
            pts.append(mat[y,x])
    return pts

def interp2mat (ngridIN, ptsIN, ngridOUT, binsOUT, fillval=0.0):
    ptsOUT = interpolate.griddata(ngridIN, ptsIN, ngridOUT, method='linear', fill_value=fillval, rescale=False)
    matOUT = np.zeros((len(binsOUT),len(binsOUT))); idx=0
    for y,b0 in enumerate(binsOUT):
        for x,b1 in enumerate(binsOUT):
            matOUT[y,x] = ptsOUT[idx]
            idx+=1
    return matOUT


def plotMatsEE(): 
    import matplotlib as mpl
    plt.set_cmap('jet')
    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    mpl.rcParams.update({'font.size': 12})
    #matplotlib.rc('text', usetex=True)
    labels = [('W+AS_norm', 'IT', 'L2/3,4'), ('W+AS_norm', 'IT', 'L5A,5B'), ('W+AS_norm', 'PT', 'L5B'), ('W+AS_norm', 'IT', 'L6'), ('W+AS_norm', 'CT', 'L6')]
    labelPostBins = [('W+AS', 'IT', 'L2/3,4'), ('W+AS', 'IT', 'L5A,5B'), ('W+AS', 'PT', 'L5B'), ('W+AS', 'IT', 'L6'), ('W+AS', 'CT', 'L6')]
    labelPreBins = ['W', 'AS', 'AS', 'W', 'W']

    #ion()
    # smat
    maxVal = max([v.max() for k,v in smat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Relative connection strength (s_{con} = p_{con} \times i_{con})', fontsize=16)
    for i,(label, preBin, postBin) in enumerate(zip(labels,labelPreBins, labelPostBins)):
        plt.subplot(len(labels),1,i+1)
        if i != 4:
            plt.setp(plt.gca().get_xticklabels(), visible=False)

        xylims = (bins[preBin][0][0], bins[preBin][-1][-1], bins[postBin][0][0], bins[postBin][-1][-1])
        im_smat = plt.imshow(smat[label], origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label[1]+' '+label[2])
        # if i == 4:
        #     xlabel(' ', fontsize=18); 
        #if i == 2:
        #    ylabel('postsynaptic normalized cortical depth (NCD)', fontsize=14)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    #plt.figtext(0.02, 0.5, 'postsynaptic normalized cortical depth (NCD)', fontsize=18, rotation=90, ha='center', va='center') # y-axis  

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im_smat, cax=cbar_ax)

    fig.savefig(dataFolder+'connEE_smat.png', dpi=1000)

    # pmat
    maxVal = max([v.max() for k,v in pmat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Probability of connection (p_{con})', fontsize=16)
    for i,(label, preBin, postBin) in enumerate(zip(labels,labelPreBins, labelPostBins)):
        plt.subplot(len(labels),1,i+1)
        xylims = (bins[preBin][0][0], bins[preBin][-1][-1], bins[postBin][0][0], bins[postBin][-1][-1])
        im_pmat=plt.imshow(pmat[label], origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label[1]+' '+label[2])
        if i != 4:
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        #plt.setp(plt.gca().get_yticklabels(), visible=False)
        #xlabel('pre ynorm'); ylabel('post ynorm')
        # if i == 4:
        #      xlabel('presynaptic normalized cortical depth (NCD)', fontsize=18); 
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im_pmat, cax=cbar_ax)

    fig.savefig(dataFolder+'connEE_pmat.png', dpi=1000)

    # wmat
    maxVal = max([v.max() for k,v in wmat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Unitary connection EPSP amplitude (weight) (i_{con})', fontsize=16)
    for i,(label, preBin, postBin) in enumerate(zip(labels,labelPreBins, labelPostBins)):
        plt.subplot(len(labels),1,i+1)
        xylims = (bins[preBin][0][0], bins[preBin][-1][-1], bins[postBin][0][0], bins[postBin][-1][-1])
        im_wmat=plt.imshow(wmat[label], origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label[1]+' '+label[2])
        # if i == 4:
        #     xlabel(' . ', fontsize=18); 
        if i != 4:
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        #plt.setp(plt.gca().get_yticklabels(), visible=False)
        #xlabel('pre ynorm'); ylabel('post ynorm')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im_wmat, cax=cbar_ax)

    fig.savefig(dataFolder+'connEE_wmat.png', dpi=1000)

    #show()

def plotMatsEEHoriz(): 
    import matplotlib as mpl
    plt.set_cmap('jet')
    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    mpl.rcParams.update({'font.size': 12})
    labels = [('W+AS_norm', 'IT', 'L2/3,4'), ('W+AS_norm', 'IT', 'L5A,5B'),  ('W+AS_norm', 'IT', 'L6'), ('W+AS_norm', 'PT', 'L5B'), ('W+AS_norm', 'CT', 'L6')]
    labelPostBins = [('W+AS', 'IT', 'L2/3,4'), ('W+AS', 'IT', 'L5A,5B'),  ('W+AS', 'IT', 'L6'), ('W+AS', 'PT', 'L5B'), ('W+AS', 'CT', 'L6')]
    labelPreBins = ['W', 'AS', 'W', 'AS', 'W']

    # smat
    maxVal = max([v.max() for k,v in smat.items() if k in labels])
    fig=plt.figure(figsize=(10,12))
    ax = [0]*5
    fontsiz = 22

    for i,(label, preBin, postBin) in enumerate(zip(labels,labelPreBins, labelPostBins)):

        binsPre = bins[preBin]
        binsPost = bins[postBin]
        smatLabel = smat[label]

        if i==0:
            ax= plt.subplot2grid((4, 2), (0, 0))  # IT L23

        elif i==1:
            ax = plt.subplot2grid((4, 2), (1, 0), rowspan=2) # IT L5
            binsPost=[[y for y in x] for x in binsPost]
            binsPost[-1] = [0.75,0.8125]
            smatLabel=[[y for y in x] for x in smatLabel]

        elif i==2:
            ax = plt.subplot2grid((4, 2), (3, 0)) # IT L6

        elif i==3:
            ax = plt.subplot2grid((4, 2), (1, 1), rowspan=2) # PT 
            binsPost=[[y for y in x] for x in binsPost]
            binsPost.insert(0,list(bins[('W+AS', 'IT', 'L5A,5B')][0]))
            binsPost[-1] = [0.75,0.8125]
            smatLabel=[[y for y in x] for x in smatLabel]
            smatLabel.insert(0,[0]*len(smatLabel[0]))

        elif i==4:
            ax = plt.subplot2grid((4, 2), (3, 1)) # CT

        if i not in  [2,4]:
            plt.setp(plt.gca().get_xticklabels(), visible=False)

        xylims = (binsPre[0][0], binsPre[-1][-1], binsPost[0][0], binsPost[-1][-1])
        plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0])
        im_smat = plt.imshow(smatLabel, origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        if i in [1,3]: plt.ylim([0.37, 0.8])
        plt.gca().invert_yaxis()
        plt.gca().tick_params(labelsize=fontsiz)
        plt.tight_layout()
        plt.subplots_adjust(top=.95, hspace=0.5)

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.48, 0.85, 0.3, 0.03])
    cbar=fig.colorbar(im_smat, cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=fontsiz)

    fig.set_figheight(8)
    fig.set_figwidth(15)
    fig.savefig(dataFolder+'connEE_smat_hor.png', dpi=600)

    #show()


def rescale(mat, oldBins, newStep):
    newBins = np.arange(oldBins[0][0], oldBins[-1][-1], newStep)
    newMat = np.zeros((len(newBins), len(newBins)))

    for i,x in enumerate(newBins):
        for j,y in enumerate(newBins):
            xindex = next((i for i,bin in enumerate(oldBins) if x>=bin[0] and x<bin[1]), None)
            yindex = next((i for i,bin  in enumerate(oldBins) if y>=bin[0] and y<bin[1]), None)
            newMat[i][j] = mat[xindex, yindex]
    # plt.plt.imshow(newMat)
    # plt.show()
    return newMat #, newBins


def plotMatsEI():  
    import matplotlib as mpl
    plt.set_cmap('jet')
    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    mpl.rcParams.update({'font.size': 12})
    labels = ['LTS', 'FS']
 
    #ion()
    # smat
    maxVal = max([v.max() for k,v in smat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Relative strength matrix (smat = pmat * wmat)', fontsize=16)
    for i,label in enumerate(labels):
        plt.subplot(len(labels),1,i+1)
        xylims = (bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1], bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1])
        im=plt.imshow(rescale(smat[label], bins['FS/LTS'], 0.01), origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label)
        if i != 1:
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        # if i == 1:
        #      xlabel(' . ', fontsize=18); 
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    #plt.figtext(0.02, 0.5, 'postsynaptic normalized cortical depth (NCD)', fontsize=18, rotation=90, ha='center', va='center') # y-axis  

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im, cax=cbar_ax)

    fig.savefig(dataFolder+'connEI_smat.png', dpi=600)

    # pmat
    maxVal = max([v.max() for k,v in pmat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Probability matrix (pmat)', fontsize=16)
    for i,label in enumerate(labels):
        plt.subplot(len(labels),1,i+1)
        xylims = (bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1], bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1])
        im=plt.imshow(rescale(pmat[label], bins['FS/LTS'], 0.01), origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label)
        if i != 1:
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        #plt.setp(plt.gca().get_yticklabels(), visible=False)
        # if i == 1:
             # xlabel('presynaptic normalized cortical depth (NCD)', fontsize=18); 
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im, cax=cbar_ax)

    fig.savefig(dataFolder+'connEI_pmat.png', dpi=600)

    # wmat
    maxVal = max([v.max() for k,v in wmat.items() if k in labels])
    fig=plt.figure(figsize=(8,12))
    #fig.sup#title('Weight matrix (wmat)', fontsize=16)
    for i,label in enumerate(labels):
        plt.subplot(len(labels),1,i+1)
        xylims = (bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1], bins['FS/LTS'][0][0], bins['FS/LTS'][-1][-1])
        im=plt.imshow(rescale(wmat[label], bins['FS/LTS'], 0.01), origin='lower', interpolation='None', aspect='auto', extent=xylims, vmin=0, vmax=maxVal)
        #title('E -> '+label)
        # if i == 1:
        #      xlabel(' . ', fontsize=18); 
        if i != 1:
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        #plt.setp(plt.gca().get_yticklabels(), visible=False)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.1)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    fig.colorbar(im, cax=cbar_ax)

    fig.savefig(dataFolder+'connEI_wmat.png', dpi=1000)
    #show()

# --------------------------------------------------------------------------------------------- #
# MAIN SCRIPT
# --------------------------------------------------------------------------------------------- #

# set folder paths for source data and output data
rootFolder = '../../'  # should point to root repo folder (m1) -- currently in m1/sim/conn  
dataFolder = rootFolder+'data/conn/'
outFolder = rootFolder+'sim/conn/'

# load data
smat, pmat, wmat, bins = loadData()

# 0. Normalize W and AS to each other (raw peak value of W is ~2x AS)
# --> Solution 1: Don't normalize since supposed to be calibrated 
# Solution 2: self-normalize to 1 based on peak value of each dataset
# Solution 3: normalize based on overlapping IT region (avg, max?)

# 1. For postsynaptic IT neurons in layers 2, 3 use the Weiler rows 1-3 (yfrac centers 0.15, 0.25, 0.35 respectively).
W_234_bins = slice(1,4)  # 0.1-0.4
smat[('W+AS', 'IT', 'L2/3,4')] = smat['W'][W_234_bins, :] 
bins[('W+AS', 'IT', 'L2/3,4')] = np.array(bins['W'][W_234_bins, :])

## decrease IT5B->IT2,IT4 (so IT intralaminar > interlaminar)
## note:upper IT5B (0.47-0.5) will not be decreased 
ITinterFactor = 1.0
W_5B_bins = slice(4,9)  # 0.5-1.0
smat[('W+AS', 'IT', 'L2/3,4')][0:2,W_5B_bins] *= ITinterFactor

# 2. For postynaptic IT neurons in layers 5A and 5B, use the A&S CSTR matrix 
# (7 rows: 1 row in L5A, yfrac [0.375,0.4375]; 1 row L5A/L5B, yfrac [0.4375,0.5], 5 rows in L5B, yfrac [0.5, 0.81])
# For ITs in the deepest part of L5B (yfrac center 0.78), value will NaN so will be converted to zero.
AS_5AB_bins = slice(6,13) # 0.375-0.8125
smat[('W+AS', 'IT', 'L5A,5B')] = smat['AS_IT'][AS_5AB_bins, :] 
bins[('W+AS', 'IT', 'L5A,5B')] = np.array(bins['AS'][AS_5AB_bins, :])

## decrease IT5B<->IT5A and IT6->IT5B (so IT intralaminar > interlaminar)
AS_5A_bins = slice(6,7)  # 0.375-0.4375
AS_5B_bins = slice(7,13)  # 0.4375-0.8125
AS_6_bins = slice(13,16)  # 0.8125-1.0
WAS_5A_bins = slice(0,1)  # 0.375-0.4375
WAS_5B_bins = slice(1,7)  # 0.4375-0.8125
smat[('W+AS', 'IT', 'L5A,5B')][WAS_5A_bins, AS_5B_bins] *= ITinterFactor  # IT5B -> IT5A
smat[('W+AS', 'IT', 'L5A,5B')][WAS_5B_bins, AS_5A_bins] *= ITinterFactor   # IT5A -> IT5B
smat[('W+AS', 'IT', 'L5A,5B')][WAS_5B_bins, AS_6_bins] *= ITinterFactor   # IT6 -> IT5B

# 3. For postsynaptic PT neurons in layer 5B, use the A&S CSP matrix (6 rows; yfrac 0.4375:0.8125, see spreadsheet).
smat[('W+AS', 'PT', 'L5B')] = smat['AS_PT'][AS_5B_bins, :] 
bins[('W+AS', 'PT', 'L5B')] = np.array(bins['AS'][AS_5B_bins, :])

# Decrease IT5B/PT5B->PT5B recurrent
L5BrecurrentFactor = 1.0
smat[('W+AS', 'PT', 'L5B')][0:6,5:] *= L5BrecurrentFactor

# L5B recurrent row sum
# print(np.sum(smat[('W+AS', 'PT', 'L5B')][0:6,5:], axis=1))
# print(smat['K'])

# 4. For postsynaptic IT, CT neurons in layer 6, use the Weiler rows 8-9 (yfrac centers 0.85, 0.95 respectively).
# IT->CT / IT->IT6: 0.620 (Yamawaki & Shepherd 2015)
W_6_bins = slice(8,10)
smat[('W+AS', 'IT', 'L6')] = smat['W'][W_6_bins, :] 
bins[('W+AS', 'IT', 'L6')] = bins['W'][W_6_bins, :] 

smat[('W+AS', 'CT', 'L6')] = smat['W'][W_6_bins, :] * 0.62
bins[('W+AS', 'CT', 'L6')] = bins['W'][W_6_bins, :]

# 5. In all cases, truncate negative weights and NaN to 0 when constructing the model.
for k,v in smat.items():
    if isinstance(k,tuple) and k[0] == 'W+AS':
        v[v<0] = 0.0
        v[np.isnan(v)] = 0.0

# 6. Rescale diagonals (increase by 20%) since underestimated in LSPS data (include pmat and wmat)
# Only for Weiler since larger bins, and each bin ~1 layer; 
# Not for AS since not clear sublayer diagonal stronger
diagScaleFactor = 1.0
for k,v in smat.items():
    if isinstance(k,tuple) and k[0] == 'W+AS':
        if k[2] in ['L2/3,4', 'L6']: preBin='W'
        elif k[2] in ['L5A,5B', 'L5B']: preBin='AS'
        for ix, x in enumerate(bins[k]):
            for iy, y in enumerate(bins[preBin]):
                if x[0]==y[0]: 
                    v[ix,iy] = v[ix,iy] * diagScaleFactor

# 7. Avoid overlap between W bins (0.3-0.4; 0.8-0.9) and AS bin (0.375-0.4375; 0.75-0.8125)
bins[('W+AS', 'IT', 'L2/3,4')][2] = [0.3, 0.37]
bins[('W+AS', 'IT', 'L5A,5B')][0] = [0.37, 0.4375]
bins[('W+AS', 'IT', 'L5A,5B')][6] = [0.75,0.8]
bins[('W+AS', 'PT', 'L5B')][5] = [0.75,0.8]

## 8. Rescale all W+AS data so L5A->L5A equivalent to Kiritani: scaling factor = 0.04621/3.41 = 0.0135
AS_L5AtoL5A = np.mean(smat[('W+AS', 'IT', 'L5A,5B')][WAS_5A_bins,AS_5A_bins]) # (1x1 L4/L5A, 0.375-0.4375)
W_L5AtoL5A = smat['W'][3,3]
WAS_L5AtoL5A = np.mean([AS_L5AtoL5A,W_L5AtoL5A])
K_L5AtoL5A = smat['K'][0,0] 
strengthFactor = 0.8  # to obtain reasonable prob conns (max ~0.3)
WAStoK_factor = K_L5AtoL5A / (AS_L5AtoL5A) * strengthFactor

smat[('W+AS_norm', 'IT', 'L2/3,4')] = smat[('W+AS', 'IT', 'L2/3,4')] * WAStoK_factor 
smat[('W+AS_norm', 'IT', 'L5A,5B')] = smat[('W+AS', 'IT', 'L5A,5B')] * WAStoK_factor 
smat[('W+AS_norm', 'PT', 'L5B')] = smat[('W+AS', 'PT', 'L5B')] * WAStoK_factor 
smat[('W+AS_norm', 'IT', 'L6')] = smat[('W+AS', 'IT', 'L6')] * WAStoK_factor 
smat[('W+AS_norm', 'CT', 'L6')] = smat[('W+AS', 'CT', 'L6')] * WAStoK_factor 

# 9a. Increase L4->L2/3 strength based on Yamawaki et al, 2015 
# Suggests similar strength to S1, so use smat L4->L2/3 from Lefort
smat[('W+AS_norm', 'IT', 'L2/3,4')][0:2,3] = smat['L'][0:2,2]

# Decrease L2/3->L4 (= ~0.25 * L4->L2/3, Yamawaki,15)
smat[('W+AS_norm', 'IT', 'L2/3,4')][2,1:3] = smat['L'][0,2] * (1.0/3.9)

# 9b. Use Yamawaki et al 2015 Fig 5 for E -> L4 inputs 
# Input strength in desc order: L4, L5/6, L2/3
# Rescale based on L4->L2/3 = 0.25*L2/3->L4
W_L23_L4 = np.max(smat[('W+AS_norm', 'IT', 'L2/3,4')][2,1:3])  # use max value

YamaScaleFactor = W_L23_L4 / smat['Y'][1] 
smat[('W+AS_norm', 'IT', 'L2/3,4')][2,3:] = smat['Y'][3:] * YamaScaleFactor  # set both bins to same low value

# 10. Separation of strength (smat) into probability (pmat) and unitary weight (wmat) 
## 10a. Check Lefort wmat L5A->L5A approx equivalent to Kiritani: 

### Use L3, L5A, L5B and L6 weights from Lefort 
### removed L2 and L4 because weights too low 
wmat['L_norm'] = wmat['L']
wmat['L_norm'][0,:] = wmat['L_norm'][1,:]
wmat['L_norm'][:,0] = wmat['L_norm'][:,1]
wmat['L_norm'][2,:] = wmat['L_norm'][3,:]
wmat['L_norm'][:,2] = wmat['L_norm'][:,3]

### bound Lefort to 0.3 - 1.0 mV (removes outliers due to small sampling size)
maxW = 1.0
minW = 0.3
wmat['L_norm'][wmat['L']>maxW] = maxW
wmat['L_norm'][wmat['L_norm']<minW] = minW

### rescale to new range of weights consistent with Kiritani
### reduces variability in weights, so W+AS smat depends more on pmat (probabilities)
### keeps weight patterns measured in Lefort
# newMinW = 0.4
# newMaxW = 0.8
# wmat['L_norm'] = ((wmat['L_norm']-minW)/(maxW-minW)) * (newMaxW-newMinW) + newMinW

### L5A->L5A Lefort to Kiritani scaling factor = 0.416/0.660 = 0.63 
### if weights too low, leads to unrealistic probs
# L_L5AtoL5A_w = wmat['L'][3,3]
# K_L5AtoL5A_w = wmat['K'][0,0] 
# LtoK_factor = K_L5AtoL5A_w / L_L5AtoL5A_w
# wmat['L_norm'] = wmat['L'] * LtoK_factor

## 10b. Generate wmat based on the rescaled Lefort 
bincenters = {}
ngrid = {}
pts = {}

### fill around Lefort wmat so interpolation works properly
wmat['L_norm_ext'] = np.insert(wmat['L_norm'], 0, [0.5]*6, axis=0)
wmat['L_norm_ext'] = np.insert(wmat['L_norm_ext'], 0, [0.5]*6, axis=0)
wmat['L_norm_ext'] = np.insert(wmat['L_norm_ext'], 0, [0.5]*8, axis=1)
wmat['L_norm_ext'] = np.insert(wmat['L_norm_ext'], 0, [0.5]*8, axis=1)
wmat['L_norm_ext'] = np.insert(wmat['L_norm_ext'], -1, [0.5]*8, axis=0)
wmat['L_norm_ext'] = np.insert(wmat['L_norm_ext'], -1, [0.5]*9, axis=1)
bins['L_ext'] = [[-0.1,0.0]]+[[0.0,128.0/1154]]+bins['L']+[[1.0, 1.1]]

### calculate bin centers
bincenters['W'] = [np.mean(b) for b in bins['W']]
bincenters['AS'] = [np.mean(b) for b in bins['AS']]
bincenters['L_ext'] = [np.mean(b) for b in bins['L_ext']]

### calculate grid and data points
ngrid['W'] = getgridfrombins(bincenters['W'])
ngrid['AS'] = getgridfrombins(bincenters['AS'])
ngrid['L_ext'] = getgridfrombins(bincenters['L_ext'])
pts['L_norm_ext'] = getpts(wmat['L_norm_ext'])

### rescale Lefort wmat to W+AS dimensions (using interpolation)
wmat['W_normL'] = interp2mat(ngrid['L_ext'], pts['L_norm_ext'], ngrid['W'], bincenters['W'])
wmat['AS_normL'] = interp2mat(ngrid['L_ext'], pts['L_norm_ext'], ngrid['AS'], bincenters['AS'])


### Create W+AS wmat 
wmat[('W+AS_norm', 'IT', 'L2/3,4')] = wmat['W_normL'][W_234_bins, :] 
wmat[('W+AS_norm', 'IT', 'L5A,5B')] = wmat['AS_normL'][AS_5AB_bins, :]  # use weights from Kiritani?
wmat[('W+AS_norm', 'PT', 'L5B')] = wmat['AS_normL'][AS_5B_bins, :]  # # use weights from Kiritani?
wmat[('W+AS_norm', 'IT', 'L6')] = wmat['W_normL'][W_6_bins, :] 
wmat[('W+AS_norm', 'CT', 'L6')] = wmat['W_normL'][W_6_bins, :] 


## 10c. Obtain pmat = smat/wmat
labels = [('W+AS_norm', 'IT', 'L2/3,4'), ('W+AS_norm', 'IT', 'L5A,5B'), ('W+AS_norm', 'PT', 'L5B'), ('W+AS_norm', 'IT', 'L6'), ('W+AS_norm', 'CT', 'L6')]
for label in labels:
    pmat[label] = smat[label] / wmat[label] 


# 11. E->I connections (based on Apicella, 2011)
## 11a. Use final matrices for IT cells but group in L2/3+4, L5A/B and L6 for FS and LTS following Apicella
## To downsample W+AS matrices, group bins into these 3 layers, and use max value so overall E->I wil be stronger than E->E
pmat['LTS'] = np.zeros((3,3))
wmat['LTS'] = np.zeros((3,3))
pmat['FS'] = np.zeros((3,3))
wmat['FS'] = np.zeros((3,3))
bins['FS/LTS'] = [[0.1, 0.37], [0.37, 0.8], [0.8,1.0]]

EtoIprobFactor = 1.0 # same

#print(bins['W'])
W_234_bins = slice(1,2)
W_5_bins = slice(4,7)
W_6_bins = slice(8,9)

AS_234_bins = slice(2,7)
AS_5_bins = slice(7,13)
AS_6_bins = slice(13,16)

for celltype in ['LTS', 'FS']:
    ## postsyn L2/3,L4
    pmat[celltype][0,0] = pmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_234_bins].mean()*EtoIprobFactor
    pmat[celltype][0,1] = pmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_5_bins].mean()*EtoIprobFactor
    pmat[celltype][0,2] = pmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_6_bins].mean()*EtoIprobFactor
    wmat[celltype][0,0] = wmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_234_bins].mean()
    wmat[celltype][0,1] = wmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_5_bins].mean()
    wmat[celltype][0,2] = wmat[('W+AS_norm', 'IT', 'L2/3,4')][:,W_6_bins].mean()

    ## postsyn 5A,5B
    pmat[celltype][1,0] = pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_234_bins].mean()*EtoIprobFactor
    pmat[celltype][1,1] = pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_5_bins].mean()*EtoIprobFactor
    pmat[celltype][1,2] = pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_6_bins].mean()*EtoIprobFactor
    wmat[celltype][1,0] = wmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_234_bins].mean()
    wmat[celltype][1,1] = wmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_5_bins].mean()
    wmat[celltype][1,2] = wmat[('W+AS_norm', 'IT', 'L5A,5B')][:,AS_6_bins].mean()

    ## postsyn L6
    pmat[celltype][2,0] = pmat[('W+AS_norm', 'IT', 'L6')][:,W_234_bins].mean()*EtoIprobFactor
    pmat[celltype][2,1] = pmat[('W+AS_norm', 'IT', 'L6')][:,W_5_bins].mean()*EtoIprobFactor
    pmat[celltype][2,2] = pmat[('W+AS_norm', 'IT', 'L6')][:,W_6_bins].mean()*EtoIprobFactor * 1.5 # increase L6 CT/IT->I (Yamawaki & Shepeherd 2015)
    wmat[celltype][2,0] = wmat[('W+AS_norm', 'IT', 'L6')][:,W_234_bins].mean()
    wmat[celltype][2,1] = wmat[('W+AS_norm', 'IT', 'L6')][:,W_5_bins].mean()
    wmat[celltype][2,2] = wmat[('W+AS_norm', 'IT', 'L6')][:,W_6_bins].mean()

## 11b. Make IT L2/3 -> LTS L5A/B = STRONG (compared to IT L2/3 -> FS L5A/B) 
## IT lower L2/3 -> LTS L5A STRONGER THAN IT lower L2/3 -> L5B (not implemented)
pmat['LTS'][1,0] = 1.5*pmat['LTS'].max() #pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,2:4].max() # increase IT L2/3 -> LTS L5A/B prob
pmat['FS'][1,0] = 0.5*pmat['FS'].min() #pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,2:4].min() # decrease IT L2/3 -> FS L5A/B prob
#wmat['FS'][1,0] = wmat[('W+AS_norm', 'IT', 'L5A,5B')][:,2:4].min() # decrease IT L2/3 -> FS L5A/B weight

## 11c. IT/PT L5A/B  -> FS L5A/B = STRONG (compared to IT/PT L5A/B -> LTS L5A/B)  
pmat['FS'][1,1] = 1.5*pmat['FS'].max() #pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,2:11].max() # increase IT/PT L5A/B -> FS L5A/B prob (use L2/3 values)
pmat['LTS'][1,1] = 0.5*pmat['LTS'].min() #pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,5:6].min() # decrease IT/PT L5A/B -> LTS L5A/B prob (>0)
#pmat['LTS'][1,1] = pmat[('W+AS_norm', 'IT', 'L5A,5B')][:,5:11].min()

## 11.d Calculate smat = pmat * wmat
smat['LTS'] = pmat['LTS'] * wmat['LTS']
smat['FS'] = pmat['FS'] * wmat['FS']

# 12. I->I connections
# Local, intralaminar only; all-to-all but distance-based; high weights
# L5A/B->L5A/B (Naka16)
# Although evidence for L2/3,4,6 -> L5A/B, strongest is intralaminar (Naka16)

# save matrices
savePickle = 1
saveMat = 0

if savePickle:
    import pickle
    with open(outFolder+'conn.pkl', 'wb') as fileObj:
        pickle.dump({'smat': smat, 'pmat': pmat, 'wmat':wmat, 'bins': bins}, fileObj)

if saveMat:
    savemat(outFolder+'conn.mat', {'smat': smat, 'pmat': pmat, 'wmat':wmat, 'bins': bins})


'''
    labels = [('W+AS_norm', 'IT', 'L2/3,4'), ('W+AS_norm', 'IT', 'L5A,5B'), ('W+AS_norm', 'PT', 'L5B'), ('W+AS_norm', 'IT', 'L6'), ('W+AS_norm', 'CT', 'L6')]
    labelPostBins = [('W+AS', 'IT', 'L2/3,4'), ('W+AS', 'IT', 'L5A,5B'), ('W+AS', 'PT', 'L5B'), ('W+AS', 'IT', 'L6'), ('W+AS', 'CT', 'L6')]
    labelPreBins = ['W', 'AS', 'AS', 'W', 'W']
'''

# print(bins[('W+AS', 'IT', 'L2/3,4')])
print(bins[('W+AS', 'IT', 'L5A,5B')])
print(bins[('W+AS', 'PT', 'L5B')])
# print(bins[('W+AS', 'IT', 'L6')])
# print(bins[('W+AS', 'CT', 'L6')])

# print(bins['W'])
# print(bins['AS'])

# print(wmat[('W+AS_norm', 'IT', 'L2/3,4')])

# print(wmat[('W+AS_norm', 'IT', 'L6')])
# print(wmat[('W+AS_norm', 'CT', 'L6')])

# print(pmat[('W+AS_norm', 'IT', 'L6')].max())
# print(pmat[('W+AS_norm', 'CT', 'L6')].max())



# plot matrices
plot = 1
if plot: 
    plotMatsEEHoriz()
    #plotMatsEE()
    #plotMatsEI()