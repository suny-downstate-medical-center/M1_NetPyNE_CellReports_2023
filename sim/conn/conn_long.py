'''

Script to combine M1 experimental long-rang connectivity data from multiple papers into single conn matrix for model

smat - strength matrix; smat = pmat * imat
cmat - convergence matrix
wmat - weight matrix

'''

import numpy as np
from scipy.io import loadmat, savemat
from pprint import pprint
from scipy import interpolate
from pylab import *
from pprint import pprint
from collections import OrderedDict
from matplotlib import pyplot as plt

colorList=[[0,0,0],[0.42,0.67,0.84], [0.90,0.76,0.00], [0.42,0.83,0.59], [0.90,0.32,0.00],
        [0.34,0.67,0.67], [0.90,0.59,0.00], [0.42,0.82,0.83], [1.00,0.85,0.00],
        [0.33,0.67,0.47], [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
        [0.71,0.82,0.41], [0.0,0.2,0.5], [0.70,0.32,0.10]]

colorList = [[0,0,0], [0.90,0.32,0.00], [0.42,0.83,0.59], [0.90,0.59,0.00], [0.90,0.76,0.00],[0.34,0.67,0.67], [0.42,0.67,0.9], [0.8, 0.25, 1.0]]

#colorList = ['r', 'g', 'b', 'y', 'k', 'r', 'g', 'b', 'y']


def plotMats():    
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import colorConverter
    from scipy.interpolate import spline,interp1d
    from matplotlib import colors as mcolors
    def cc(arg):
        return mcolors.to_rgba(arg, alpha=0.6)

    longPops = ['TPO', 'TVL', 'S1', 'S2', 'cM1', 'M2', 'OC']
    cellTypes = ['IT', 'PT', 'CT']#, 'PV', 'SOM']
    #cellTypes = ['CT', 'PV', 'SOM']
    #cellTypes = ['IT']
    EorIvals = ['exc', 'inh']

    #ion()
    # cmat
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    matplotlib.rcParams.update({'font.size': 16}) #, 'font.weight':'bold'})
    for EorI in ['exc']:#EorIvals:
        fig=figure(figsize=(20,7))
        #fig.suptitle('Convergence matrices', fontsize=16)
        i=0

        zs = np.arange(0, len(longPops), 1.0)
            
        for ct in cellTypes:
            #ax = plt.subplot(111, projection='3d')
            ax = fig.add_subplot(1,len(cellTypes),i+1, projection='3d')
            
            i=i+1
            title('Long-range ' + EorI + ' inputs -> '+ ct, fontsize=16, fontweight='bold')
            ax.elev=30 # 90 
            ax.azim=-45#-45 # -90
            verts = []
            for longPop in longPops:
                xconv = cmat[(longPop,ct,EorI)]
                ybins = [mean(x) for x in bins[(longPop, ct)]]
                if len(xconv) == 1:
                    ybins = [(1.0+0.0*ix)*x for pair in bins[(longPop, ct)] for ix,x in enumerate(pair)]
                    xconv = [(1.0+0.0*pt)*x for x in cmat[(longPop,ct,EorI)] for pt in (0,1)]

                y_smooth = linspace(min(ybins), max(ybins), 50)
                xfunc = interp1d(ybins, xconv, kind='linear')
                x_smooth = list(xfunc(y_smooth))
                y_smooth=list(y_smooth)

                x_smooth.insert(0,0)
                x_smooth.append(0)
                y_smooth.insert(0,min(y_smooth)-0.03)
                y_smooth.append(max(y_smooth)+0.03)

                verts.append(list(zip(y_smooth,x_smooth)))

            # verts.insert(0,[(0,0),(0,0.1)])
            # verts.append([(0,0.1), (0,0)])

            poly = PolyCollection(verts, linewidths=[2.0]*len(verts),
                                        edgecolor=[cc(c)[0:3]+(1.0,) for c in colorList[1:len(verts)+1]], 
                                        #edgecolor=[0.1, 0.1, 0.1], 
                                        facecolors = [cc(c) for c in colorList[1:len(verts)+1]]) #[0:len(verts)])


            ax.add_collection3d(poly, zs=zs, zdir='y')

            popLabels = [' '] + longPops # + [' ']
            ax.set_xlabel('\npostsynaptic NCD', fontsize=18, fontweight='normal',linespacing=2)
            ax.set_xlim3d(0, 1)
            ax.set_zlabel('\nconvergence', fontsize=18, fontweight='normal', linespacing=0.5)
            zmax=max([p[1] for points in verts for p in points])
            ax.set_zlim3d(0, zmax)
            ax.set_zticks(arange(0, zmax, 40))
            ax.set_ylabel('\nLong range input', fontsize=18, fontweight='normal', linespacing=2)
            ax.set_ylim3d(-1, max(zs))
            #ax.set_yticks(arange(0, max(zs)-1))
            ax.set_yticklabels(popLabels)


            tight_layout()
            subplots_adjust(top=0.9, bottom=0.1)
        if fixedNumSyns:
            filename = dataFolder+'conn_long_cmat_fixedNumSyns_3D_'+cellTypes[0]+EorI+'.png'
        else:
            filename = dataFolder+'conn_long_cmat_varNumSyns_'+EorI+'.png'

        #plt.show()

        #fig.savefig(filename, dpi=600,  bbox_inches='tight')# bbox_extra_artists=(lgd,),)
        fig.savefig(filename, dpi=300)# bbox_inches='tight')

    return ax



def plotMats2D(labels=0, cbar=0):
    from scipy.interpolate import spline,interp1d

    longPops = ['TPO', 'TVL', 'S1', 'S2', 'cM1', 'M2', 'OC']
    cellTypes = ['IT', 'PT', 'CT']#, 'PV', 'SOM']
    EorIvals = ['exc', 'inh']
    # cmat
    ybinsNew = np.arange(0, 1.0, 0.0001)
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    matplotlib.rcParams.update({'font.size': 16}) #, 'font.weight':'bold'})
    plt.set_cmap('jet')
    fontsiz=26
    step = 0.0001
    
    maxVal = max([max(v) for k,v in cmat.items()])

    for EorI in ['exc']:#EorIvals:
        fig=figure(figsize=(15,6))            
        for isubplot, ct in enumerate(cellTypes):
            fig.add_subplot(1, len(cellTypes), isubplot+1)
            if labels: plt.title('Long-range ' + EorI + ' inputs -> '+ ct, fontsize=16, fontweight='bold')
            mat2d = []
            for longPop in longPops:
                #convs = [cmat[(longPop,ct,EorI)][0]]+cmat[(longPop,ct,EorI)]+[cmat[(longPop,ct,EorI)][-1]]
                #ybins = [min(bins[(longPop, ct)][0])]+[mean(x) for x in bins[(longPop, ct)]]+[max(bins[(longPop, ct)][-1])]
                ybinsAll = [bins[(longPop, ct)][0][0]]+[x[1] for x in bins[(longPop, ct)]]

                convsExt, ybinsExt = [], []
                prevConv = 0
                for conv, ybin in zip(cmat[(longPop,ct,EorI)]+[cmat[(longPop,ct,EorI)][-1]], ybinsAll):
                    ybinsExt.append(ybin-step)
                    ybinsExt.append(ybin+step)
                    convsExt.append(prevConv)
                    convsExt.append(conv)
                    prevConv = conv

                f = interp1d(ybinsExt, convsExt, fill_value=(0, 0), bounds_error=False, kind='nearest')
                convsNew = f(ybinsNew)
                mat2d.append(convsNew)

                #from IPython import embed; embed()


            im_conv = plt.imshow(np.array(mat2d).T, origin='lower', interpolation='None', aspect='auto', extent=[0,7, 0,1], vmin=0, vmax=maxVal)
            
            if labels: plt.ylabel('postsynaptic NCD', fontsize=16, fontweight='normal',linespacing=2)
            if labels: plt.xlabel('Long range input', fontsize=16, fontweight='normal', linespacing=2)
            plt.xlim(0, len(longPops))
            yaxis=[x/10.0 for x in np.arange(0, 11,2)]
            #print yaxis
            plt.gca().set_yticks(yaxis)
            plt.gca().set_yticklabels(yaxis,fontsize=fontsiz)
            plt.gca().set_xticks([x+0.5 for x in range(7)])
            plt.gca().set_xticklabels(longPops, fontsize=fontsiz, rotation=45)
            tight_layout()
            plt.gca().invert_yaxis()
            subplots_adjust(top=0.9, bottom=0.15)
        if cbar: 
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=fontsiz)
            subplots_adjust(top=0.9, bottom=0.15, right=0.8)
            if labels:
                cbar.ax.set_ylabel('convergence')


        filename = dataFolder+'conn_long_cmat_2D_'+EorI+'.png'

        #plt.show()

        #fig.savefig(filename, dpi=600,  bbox_inches='tight')# bbox_extra_artists=(lgd,),)
        fig.savefig(filename, dpi=300)# bbox_inches='tight')

    return fig


def plotPies():
    fracs = {}
    fracLabels = {}
    fracLabels['EI'] = ['Exc', 'Inh']
    fracs['EI'] = [fracSyns['M1'][k] for k in ['exc','inh']]
    fracLabels['ELL'] = ['Long (exc)', 'Local (exc)']
    fracs['ELL'] = [fracSyns['M1'][k] for k in ['exc_long','exc_local']]
    fracLabels['ILL'] = ['Long (inh)', 'Local (inh)']
    fracs['ILL'] = [fracSyns['M1'][k] for k in ['inh_long','inh_local']]
    fracLabels['L'] = ['TPO', 'TVL', 'S1', 'S2', 'cM1', 'M2', 'OC', 'Other']
    fracs['L'] = [fracSyns['M1'][k] for k in fracLabels['L']]

    fracsTitles = {'EI':'Exc vs inh inputs', 'ELL': 'Exc long-range vs local inputs', 'ILL': 'Inh long-range vs local inputs', 'L':'Long-range inputs'} 
    
    print (fracs)

    for k in fracs.keys():
        fracLabel,frac = fracLabels[k],fracs[k]
        # make a square figure and axes
        figure(1, figsize=(9,9))
        ax = axes([0.1, 0.1, 0.8, 0.8])

        # The slices will be ordered and plotted counter-clockwise.
        labels = fracLabel

        ## RED/BLUE for E/I!!
        print (k)
        if k=='EI':
             colors = [colorList[4], colorList[1]]
        else:
            colors = colorList[1:len(frac)+1]
        mpl.rcParams['font.size'] = 20.0
        mpl.rcParams['font.weight'] = 'bold'

        # DONT CUT LABELS!!
        pie(frac, labels=labels, autopct='%1.0f%%', pctdistance=0.7, labeldistance=1.05, radius=0.8, shadow=True, startangle=0, colors=colors)
        title(fracsTitles[k])
        subplots_adjust(top=0.9)
        savefig(dataFolder+k+'_frac.png')
        show()



# --------------------------------------------------------------------------------------------- #
# MAIN SCRIPT
# --------------------------------------------------------------------------------------------- #

# set folder paths for source data and output data
rootFolder = '../../'  # should point to root repo folder (m1) -- currently in m1/sim/conn  
dataFolder = rootFolder+'data/conn/'
outFolder = rootFolder+'sim/conn/'

# load data
smat = {}  # dict for strength matrices (smat = cmat * wmat) - not used here
cmat = {}  # dict for convergence matrices
wmat = {}  # dict for weight matrices
bins = {}  # dict for bin intervals
numSyns = {} # dict with number of syns per cell
synsPerConn = {}  # dict with syns per connection
numCells = {}  # dict for the number of artificial spike generators (NetStims) per presyn population
delays = {}  # dict for the delays 
rates = {}  # dict for the spontaneous firing rate of each presyn pop

## cell types
cellTypes = ['IT', 'PT', 'CT', 'SOM', 'PV']

## General coordinate system / layer boundaries '''
bins['layerLabels'] =   ['pia', 'L1-L2 border', 'L2/3-L4 border',   'L4-L5A border',   'L5A-L5B border',    'L5B-L6 border',   'L6-WM border']  
bins['layers']      =   [0,     0.1,           0.29,               0.37,               0.47,               0.8,               1.0]
layers = bins['layers']

# ------------------------------------------------------------------------------------------------------------------
# 1) Set total num of long range syns/cell:
#  a) Start from Ben's estimate (dend length and 1.5 spines/um) for our M1 PT cell = ~16400
#  b) Estimate num for other cell types / layers based on Markram,2015/Meyer 2012 (rat S1), scaling based on PT = 10400
#  c) Make sure consistent with Shuz 1989 (~7500); not sure what to make of DeFelipe 2002 (21983) - ignore since estimate ?
#
# Note: Maybe set same num syns for all pops so the %-based syn input distributions from papers have the desired effect ?! 

synDensity = 1.5  # spines/um
PTlength = 10933  # total spiny dend length (um)
numSyns['M1'] = {}
numSyns['M1']['PT5B'] = int(round(PTlength * synDensity)) # total number of syns (16400)

## Meyer et al, 2010
numSyns['Meyer'] = {'IT2': 2155, 'IT3': 2973, 'IT4': 1777, 'STL4': 1275, 'IT5': 2887, 'PT5B': 5178, 'IT6': 1954}

## scaling Meyer to M1 PT
numSynsFactor = float(numSyns['M1']['PT5B']) / numSyns['Meyer']['PT5B']  # scale factor between mouse M1 PT cell and Meyer 2010 rat S1

## Use variable or fixed (mean) number of syns for each layers
fixedNumSyns = True  # use diff number of syns for each layer; otherwise use mean for all
if fixedNumSyns:
    meanNumSyns = mean(list(numSyns['Meyer'].values()))
    numSyns['Meyer'] = {k: meanNumSyns for k in numSyns['Meyer']}

## tot num syns exc pops
m1EPops = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6']
s1EPops = ['IT2', 'IT4', 'IT5', 'IT5', 'PT5B', 'IT6', 'IT6']

for m1pop, s1pop in zip(m1EPops, s1EPops):
    numSyns['M1'][m1pop] = int(round(numSyns['Meyer'][s1pop] * numSynsFactor))

## tot num syns inh pops
m1IPops = ['PV2', 'SOM2', 'PV5A', 'SOM5A', 'PV5B', 'SOM5B', 'PV6', 'SOM6']
s1IPops = ['IT2', 'IT2', 'IT5', 'IT5', 'IT5', 'IT5', 'IT6', 'IT6']

numSynsFactorInh =  1.0  # Potjans & Diesman used ~0.8; Mark15 mentions ~10:1 of E vs I thalamic input synapses (consistent with proportion of E:I cells)
for m1pop, s1pop in zip(m1IPops, s1IPops):
    numSyns['M1'][m1pop] = int(round(numSyns['Meyer'][s1pop] * numSynsFactor * numSynsFactorInh))

#pprint(numSyns['M1'])
#pprint(mean(numSyns['M1'].values()))

# ------------------------------------------------------------------------------------------------------------------
# 2) Set num of total syns/cell for each input region: 
#  a) Set % of exc vs inh, and long vs local
#  overall long:long = 70-30%
fracSyns = {'M1': {}}  
fracSyns['M1']['exc'] = 0.8  # exc/inh inputs = 84/16% (DeFelipe 2002); 75%/25% (Markram 2015)
fracSyns['M1']['inh'] = 0.2
fracSyns['M1']['exc_long'] = 0.8  # local/long exc = 20-80% (Markram 2015; Stepanyants 2009); 25-75% (Suter 2013) 
fracSyns['M1']['exc_local'] = 0.2
fracSyns['M1']['inh_long'] = 0.3  # local/long inh = 70-30% (Stepanyants 2009)
fracSyns['M1']['inh_local'] = 0.7 

## Thalamus VL and PO
##  mouse TVAL->M1 PT = 18% (Suter, 2013); rat VPM,POm->S1 = 7.4%, 4% (Meyer, 2010b); rat Thal -> L4 S1 rat (Bruno 2006): ~15%; PO=0.7, VAL=3.5 (Oh 2014)
fracSyns['M1']['TVL'] = 0.15  # Thalamus, motor-related, ventro-lateral, core, Cerebellar-relaying
fracSyns['M1']['TPO'] = 0.15  # Thalamus, sensory-related, medial posterior, matrix; less syns thatn VL (Meyer 2011; Zhang 2016; Oh 2014; Bopp 2017)

## S1
## upper-limb=1.16, barrel=0.86 (Oh 2014, norm strength); upper-limb=8.61%, barrel=29% (Zhang 2016, % presyn cells); 
fracSyns['M1']['S1'] = 0.15  

## S2
## 0.83 (Oh 2014, norm strength); ~1/12th long-range syns (Suter 2013 notes)
fracSyns['M1']['S2'] = 0.15

## M2
## 0.72 (Oh 2014, norm strength)
fracSyns['M1']['M2'] = 0.10 

## contralateral M1
## 1.27 (Oh 2014, norm strength)
fracSyns['M1']['cM1'] = 0.15 

## OC
## 0.32+0.36 (Oh 2014, norm strength)
fracSyns['M1']['OC'] = 0.10 

## Other (not modeled)
fracSyns['M1']['Other'] = 0.5


# ------------------------------------------------------------------------------------------------------------------
# 3) Set convergence based on layer/yfrac

#############
## Thalamus 
inputPops = ['TPO', 'TVL']
normInput = {'Yama15elife': {}, ('Yama15', 'IT'): {}, ('Yama15', 'PT'): {}, ('Yama15', 'CT'): {}, ('Yama15', 'SOM'): {}, ('Yama15', 'PV'): {}}

## Yamawaki,2015 elife (forelimb)
normInput['Yama15elife']['TPO'] = [0.26, 0.72, 0.72, 0.58, 0.25, 0.090, 0.044, 0.022]
normInput['Yama15elife']['TVL'] = [0.30, 0.43, 0.64, 0.24, 0.34, 0.77, 0.69, 0.061] 
normInput['Yama15elife']['TVL'][4:7] = [np.mean(normInput['Yama15elife']['TVL'][4:7])]*3  # make all L5B homogeneous as suggested by gordon 
bins['Yama15elife'] = [[0.1,0.2], [0.2,0.3], [0.3,0.4], [0.4,0.5], [0.5,0.6], [0.6,0.7], [0.7,0.8], [0.8,0.9]]

## yfrac based bins but corrected so each bin is within a single layer (L2/3 has 2 bins; L5B has 3 bins; L6 has 2 bins; rest 1 bin)
## corresponds approximately to original bin ranges
bins['Yama15_labels'] = ['L23u', 'L23l', 'L4', 'L5A', 'L5Bu', 'L5Bm', 'L5Bl', 'L6u', 'L6l'] 
bins['Yama15'] = [[layers[1], layers[1]+(layers[2]-layers[1])/2], 
                [layers[1]+(layers[2]-layers[1])/2, layers[2]], 
                [layers[2], layers[3]], [layers[3], layers[4]], \
                [layers[4], layers[4]+(layers[5]-layers[4])/3], 
                [layers[4]+(layers[5]-layers[4])/3, layers[4]+2*(layers[5]-layers[4])/3], 
                [layers[4]+2*(layers[5]-layers[4])/3, layers[5]], 
                [layers[5], layers[5]+(1.0-layers[5])/2], 
                [layers[5]+(1.0-layers[5])/2, layers[6]]]


## Extend based on Yamawaki,2015 jns (TVL)
## CT/PT = 0.05; CT/IT5B = 0.16; CT/IT6 (yfrac<0.85) = 0.33; CT=IT6 (yfrac>0.85) = 0.2
## 3 eqn system: CT/PT5B = 0.05; CT/IT5B = 0.16; (PT5B+IT5B)/2 = 0.6 -> solutions: CT=0.0457, PT=0.914, IT5B=0.286
## To match Yamawaki jns ratios for IT/PT/CT, Yamawaki elife values were adapted as follows: 
## PT = (0.914/0.6)*L5B; IT5B = (0.286/0.6)*L5B; CT=0.0457; IT6 upper = 3*CT; IT6 lower = CT  
## TPO keep same values but adapted to use same structure as VL 
## Naka16 confirmed in Table 1: VL/PO->IT=low,high; VL/PO->PT=high,low 
normInput[('Yama15', 'PT')]['TPO'] = normInput['Yama15elife']['TPO'][4:6+1]  
normInput[('Yama15', 'PT')]['TVL'] = [x*(0.914/0.6) for x in normInput['Yama15elife']['TVL'][4:6+1]]  # PT = (0.914/0.6)*L5B
bins[('TVL', 'PT')] = bins[('TPO', 'PT')] = bins['Yama15'][4:6+1]

normInput[('Yama15', 'CT')]['TPO'] = [normInput['Yama15elife']['TPO'][7]] * 2  # add extra L6 
normInput[('Yama15', 'CT')]['TVL'] = [0.0457, 0.0457] # CT=0.0457 (derived from 3 eqn system above based on yamawaki jns ratios)
bins[('TVL', 'CT')] = bins[('TPO', 'CT')] = bins['Yama15'][7:8+1]

normInput[('Yama15', 'IT')]['TPO'] = normInput['Yama15elife']['TPO'] + [normInput['Yama15elife']['TPO'][-1]]  # add extra L6 
normInput[('Yama15', 'IT')]['TVL'] = normInput['Yama15elife']['TVL'] + [normInput['Yama15elife']['TVL'][-1]]  # add extra L6 bin 
normInput[('Yama15', 'IT')]['TVL'][4:6+1] = [x*0.286/0.6 for x in normInput[('Yama15', 'IT')]['TVL'][4:6+1]] # IT5B = (0.286/0.6)*L5B
normInput[('Yama15', 'IT')]['TVL'][2] = np.mean(normInput[('Yama15', 'IT')]['TVL'][0:2]) * (0.78/0.3) # increase TVL->L4 based on Yama15 Fig 6F: ratio of L2/3:L4 = 0.3:0.78
normInput[('Yama15', 'IT')]['TVL'][7] = normInput[('Yama15', 'CT')]['TVL'][0]*3 # IT6 upper = 3*CT;
normInput[('Yama15', 'IT')]['TVL'][8] = normInput[('Yama15', 'CT')]['TVL'][1] # IT6 lower = CT 
bins[('TVL', 'IT')] = bins[('TPO', 'IT')] = bins['Yama15'][0:8+1]

bins[('TVL', 'PV')] = bins[('TVL', 'SOM')] = bins[('TVL', 'IT')]
bins[('TPO', 'PV')] = bins[('TPO', 'SOM')] =  bins[('TPO', 'IT')]

## For inhib cells assume original Yamawaki elife yfrac-based values
normInput[('Yama15', 'SOM')]['TPO'] = normInput['Yama15elife']['TPO'] + [normInput['Yama15elife']['TPO'][-1]] # add extra L6 
normInput[('Yama15', 'PV')]['TPO'] = normInput['Yama15elife']['TPO'] + [normInput['Yama15elife']['TPO'][-1]]
normInput[('Yama15', 'SOM')]['TVL'] = normInput['Yama15elife']['TVL'] + [normInput['Yama15elife']['TVL'][-1]]
normInput[('Yama15', 'PV')]['TVL'] = normInput['Yama15elife']['TVL'] + [normInput['Yama15elife']['TVL'][-1]]

## normalize so max=1.0
for inputPop in inputPops:
    maxValue = max([x for ct in cellTypes for x in normInput[('Yama15', ct)][inputPop]])
    for ct in cellTypes:
        normInput[('Yama15', ct)][inputPop] = [x/maxValue for x in normInput[('Yama15', ct)][inputPop]]

## syns per conn (= synaptic contacts per connection)
## (Bruno 2006, in vivo rat S1) = 7 syns/con; (Markram 2015, Rat S1)  = 8.1+-4.2 syns/con
synsPerConn['TPO'] = 5
synsPerConn['TVL'] = 5

## calculate convergence
M1BinPops = {}  # pop corresponding to each bin (thalamus) to obtain syns per cell
M1BinPops['IT'] = ['IT2', 'IT2', 'IT4', 'IT5A', 'IT5B',  'IT5B',  'IT5B', 'IT6', 'IT6']
M1BinPops['PT'] = ['PT5B', 'PT5B', 'PT5B']
M1BinPops['CT'] = ['CT6', 'CT6']
M1BinPops['PV'] = ['PV2', 'PV2', 'PV5A', 'PV5A', 'PV5B', 'PV5B', 'PV5B', 'PV6', 'PV6']
M1BinPops['SOM'] = ['SOM2', 'SOM2', 'SOM5A', 'SOM5A', 'SOM5B', 'SOM5B', 'SOM5B', 'SOM6', 'SOM6']

for inputPop in inputPops:
    for ct in cellTypes:
        for EorI in ['exc', 'inh']:
            # num of syns = total * exc * exc_long_range * fraction for this thalamic input 
            nsyns = [int(round(numSyns['M1'][M1BinPop] * fracSyns['M1'][EorI] * fracSyns['M1'][EorI+'_long'] * fracSyns['M1'][inputPop])) for M1BinPop in M1BinPops[ct]]
            # convergence = num syns / syns per conn * normalized input
            cmat[(inputPop, ct, EorI)] = [float(round(nsyn / synsPerConn[inputPop] * ninput)) for nsyn,ninput in zip(nsyns, normInput[('Yama15', ct)][inputPop])]


## Check Yamawaki 2015 jns ratios in terms of the number of syns per cell type
# print cmat 
# print 'Ratio of num of syns for comparison with Yamawaki 2015 jns (note tot syns different for each cell type):'
# print 'CT/PT (Yamawaki=0.05) = ', cmat[('TVL','CT','exc')][0] / mean(cmat[('TVL','PT','exc')])
# print 'CT/IT5B (Yamawaki=0.16) ', cmat[('TVL','CT','exc')][0] / mean(cmat[('TVL','IT','exc')][4:7])
# print 'CT/IT6upper (Yamawaki=0.33) ', cmat[('TVL','CT','exc')][0] / cmat[('TVL','IT','exc')][7]
# print 'CT/IT6lower (Yamawaki=1) ', cmat[('TVL','CT','exc')][0] / cmat[('TVL','IT','exc')][8]


#############
## S1
## L23, 0.7, L5A, 1.0, L5B, 0.14, L6, 0.1 (Mao 2011)
## L4 cortical input = very low so set to 0.1 (Yama15 elife)
normInput['Mao11'] = [0.7, 1.0, 0.14, 0.1]
bins['Mao11'] =  [[layers[1], layers[2]], [layers[2],layers[4]], [layers[4], layers[5]], [layers[5],1.0]]

## Adapt to use in model by setting specific bins for each cell type
normInput[('Mao11', 'IT')] = [0.7, 0.1, 1.0, 0.14, 0.1]
bins[('S1', 'IT')] =  [[layers[1],layers[2]], [layers[2],layers[3]], [layers[3], layers[4]], [layers[4], layers[5]], [layers[5],1.0]]

normInput[('Mao11', 'PT')] = [normInput[('Mao11', 'IT')][3]]
bins[('S1', 'PT')] =  [bins[('S1', 'IT')][3]]

normInput[('Mao11', 'CT')] = [normInput[('Mao11', 'IT')][4]]
bins[('S1', 'CT')] =  [bins[('S1', 'IT')][4]]

normInput[('Mao11', 'PV')] = normInput[('Mao11', 'IT')]
bins[('S1', 'PV')] = bins[('S1', 'IT')] 

normInput[('Mao11', 'SOM')] = normInput[('Mao11', 'IT')] # input? yes (Wall 2016) vs no (Harris & Shepeherd 2015) 
bins[('S1', 'SOM')] =  bins[('S1', 'IT')] 

## syns per conn
synsPerConn['S1'] = 5

## calculate convergence
M1BinPops = {}  # pop corresponding to each bin (S1) to obtain syns per cell
M1BinPops['IT'] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'IT6']
M1BinPops['PT'] = ['PT5B']
M1BinPops['CT'] = ['CT6']
M1BinPops['PV'] = ['PV2', 'PV5A', 'PV5A', 'PV5B',  'PV6']
M1BinPops['SOM'] = ['SOM2', 'SOM5A', 'SOM5A', 'SOM5B', 'SOM6']

inputPop = 'S1'
for ct in cellTypes:
    for EorI in ['exc', 'inh']:
        # num of syns = total * exc * exc_long_range * fraction for this thalamic input 
        nsyns = [int(round(numSyns['M1'][M1BinPop] * fracSyns['M1'][EorI] * fracSyns['M1'][EorI+'_long'] * fracSyns['M1'][inputPop])) for M1BinPop in M1BinPops[ct]]
        
        # convergence = num syns / syns per conn * normalized input
        cmat[(inputPop, ct, EorI)] = [float(round(nsyn / synsPerConn[inputPop] * ninput)) for nsyn,ninput in zip(nsyns, normInput[('Mao11', ct)])]


#############
## S2
## L23=1.0, L5A=0.46, L5B=0.35, L6=0.20; in L5B decreases with slope=-5.6 (Suter 2015)
## divide L5B into 3 bins and calculate strength of each bin based on slope
## bin steps: (layers[5]-layers[4])/3 
## Interpolate L5Bupper and L5Blower from values of L5A (0.46), L5Bmid (0.35) and L6 (0.2)
## Also add layer 4 but set to low input (0.1) (Yama15elife) 
L5Bwidth = layers[5]-layers[4]
#normInput['Suter15'] = [1.0, 0.1, 0.46, (0.46+0.35)/2, 0.35, (0.35+0.20)/2, 0.20]
normInput['Suter15'] = [1.0, 0.1, 0.46, 0.35, 0.35, 0.35, 0.20]
bins['Suter15'] =  [[layers[1],layers[2]], 
                    [layers[2],layers[3]], 
                    [layers[3],layers[4]], 
                    [layers[4], layers[4]+L5Bwidth/3], 
                    [layers[4]+L5Bwidth/3, layers[4]+2*L5Bwidth/3], 
                    [layers[4]+2*L5Bwidth/3, layers[5]], 
                    [layers[5],1.0]]

## Adapt to use in model by setting specific bins for each cell type
normInput[('Suter15', 'IT')] = normInput['Suter15'] 
bins[('S2', 'IT')] =  bins['Suter15'] 

normInput[('Suter15', 'PT')] = normInput['Suter15'][3:5+1]
bins[('S2', 'PT')] =  bins['Suter15'] [3:5+1]

normInput[('Suter15', 'CT')] = [normInput['Suter15'][6]]
bins[('S2', 'CT')] =  [bins['Suter15'][6]]

normInput[('Suter15', 'PV')] = normInput['Suter15'] 
bins[('S2', 'PV')] = bins['Suter15'] 

normInput[('Suter15', 'SOM')] = normInput['Suter15'] # input? yes (Wall 2016) vs no (Harris & Shepeherd 2015) 
bins[('S2', 'SOM')] = bins['Suter15'] 

## syns per conn
synsPerConn['S2'] = 5

## calculate convergence
M1BinPops = {}  # pop corresponding to each bin (S2) to obtain syns per cell
M1BinPops['IT'] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'IT5B', 'IT5B', 'IT6']
M1BinPops['PT'] = ['PT5B', 'PT5B', 'PT5B']
M1BinPops['CT'] = ['CT6']
M1BinPops['PV'] = ['PV2', 'PV5A', 'PV5A', 'PV5B', 'PV5B', 'PV5B', 'PV6']
M1BinPops['SOM'] = ['SOM2', 'SOM5A', 'SOM5A', 'SOM5B', 'SOM5B', 'SOM5B', 'SOM6']

inputPop = 'S2'
for ct in cellTypes:
    for EorI in ['exc', 'inh']:
        # num of syns = total * exc * exc_long_range * fraction for this thalamic input 
        nsyns = [int(round(numSyns['M1'][M1BinPop] * fracSyns['M1'][EorI] * fracSyns['M1'][EorI+'_long'] * fracSyns['M1'][inputPop])) for M1BinPop in M1BinPops[ct]]
        
        # convergence = num syns / syns per conn * normalized input
        cmat[(inputPop, ct, EorI)] = [float(round(nsyn / synsPerConn[inputPop] * ninput)) for nsyn,ninput in zip(nsyns, normInput[('Suter15', ct)])]


#############
## M2 and contralateral M1
## L23: 0.17, L5A: 0.11, L5B: 1.0, L6: 0.44; L5Blower significantly higher than upper (Hooks 2013)
## divided L5B in upper vs lower: L23: 0.17, L5A: 0.11, L5Bu: 1.0, L5Bl: 1.0, L6: 0.44
## but kept same value for upper vs lower LT5B since not clear from graph (also gmgs email)
## Also add layer 4 as subset of L5A
normInput['Hooks13_M2'] = [0.17, 0.11, 0.11, 1.0, 1.0, 0.44] 
bins['Hooks13_M2'] =  [[layers[1], layers[2]], 
                        [layers[2], layers[3]], 
                        [layers[3], layers[4]], 
                        [layers[4], layers[4]+(layers[5]-layers[4])/2], 
                        [layers[4]+(layers[5]-layers[4])/2, layers[5]],
                        [layers[5], 1.0]]


## Adapt to use in model by setting specific bins for each cell type
normInput[('Hooks13_M2', 'IT')] = normInput['Hooks13_M2'] 
bins[('M2', 'IT')] = bins['Hooks13_M2'] 

normInput[('Hooks13_M2', 'PT')] = normInput['Hooks13_M2'][3:4+1]
bins[('M2', 'PT')] = bins['Hooks13_M2'] [3:4+1]

normInput[('Hooks13_M2', 'CT')] = [normInput['Hooks13_M2'][5]]
bins[('M2', 'CT')] = [bins['Hooks13_M2'][5]]

normInput[('Hooks13_M2', 'PV')] = normInput['Hooks13_M2'] 
bins[('M2', 'PV')] = bins['Hooks13_M2'] 

normInput[('Hooks13_M2', 'SOM')] = normInput['Hooks13_M2'] # input? yes (Wall 2016) vs no (Harris & Shepeherd 2015) 
bins[('M2', 'SOM')] = bins['Hooks13_M2'] 

for ct in cellTypes:
    bins[('cM1', ct)] = bins[('M2', ct)]

## syns per conn
synsPerConn['cM1'] = 5
synsPerConn['M2'] = 5

## calculate convergence
M1BinPops = {}  # pop corresponding to each bin (S2) to obtain syns per cell
M1BinPops['IT'] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'IT5B', 'IT6']
M1BinPops['PT'] = ['PT5B', 'PT5B']
M1BinPops['CT'] = ['CT6']
M1BinPops['PV'] = ['PV2', 'PV5A', 'PV5A', 'PV5B', 'PV5B', 'PV6']
M1BinPops['SOM'] = ['SOM2', 'SOM5A', 'SOM5A', 'SOM5B', 'SOM5B', 'SOM6']


for pop in ['cM1', 'M2']:
    for ct in cellTypes:
        for EorI in ['exc', 'inh']:
            # num of syns = total * exc * exc_long_range * fraction for this thalamic input 
            nsyns = [int(round(numSyns['M1'][M1BinPop] * fracSyns['M1'][EorI] * fracSyns['M1'][EorI+'_long'] * fracSyns['M1'][pop])) for M1BinPop in M1BinPops[ct]]
            
            # convergence = num syns / syns per conn * normalized input
            cmat[(pop, ct, EorI)] = [float(round(nsyn / synsPerConn[pop] * ninput)) for nsyn,ninput in zip(nsyns, normInput[('Hooks13_M2', ct)])]


#############
## OC
## L23: 0.19, L5A: 0.30, L5B, 0.29, L6, 1.0 (Hooks 2013)
## L4: low cortical input (Yama15elife)
## Add layer 4
normInput['Hooks13_OC'] = [0.19, 0.1, 0.30, 0.29, 1.0]
bins['Hooks13_OC'] =  [[layers[1],layers[2]], [layers[2],layers[3]], [layers[3],layers[4]], [layers[4], layers[5]], [layers[5],1.0]]


## Adapt to use in model by setting specific bins for each cell type
normInput[('Hooks13_OC', 'IT')] = normInput['Hooks13_OC'] 
bins[('OC', 'IT')] = bins['Hooks13_OC'] 

normInput[('Hooks13_OC', 'PT')] = [normInput['Hooks13_OC'][3]]
bins[('OC', 'PT')] = [bins['Hooks13_OC'][3]]

normInput[('Hooks13_OC', 'CT')] = [normInput['Hooks13_OC'][4]]
bins[('OC', 'CT')] = [bins['Hooks13_OC'][4]]

normInput[('Hooks13_OC', 'PV')] = normInput['Hooks13_OC'] 
bins[('OC', 'PV')] = bins['Hooks13_OC']

normInput[('Hooks13_OC', 'SOM')] = normInput[('Hooks13_OC', 'IT')] # input? yes (Wall 2016) vs no (Harris & Shepeherd 2015) 
bins[('OC', 'SOM')] = bins['Hooks13_OC']

## syns per conn
synsPerConn['OC'] = 5

## calculate convergence
M1BinPops = {}  # pop corresponding to each bin (S1) to obtain syns per cell
M1BinPops['IT'] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'IT6']
M1BinPops['PT'] = ['PT5B']
M1BinPops['CT'] = ['CT6']
M1BinPops['PV'] = ['PV2', 'PV5A', 'PV5A', 'PV5B',  'PV6']
M1BinPops['SOM'] = ['SOM2', 'SOM5A', 'SOM5A', 'SOM5B', 'SOM6']

inputPop = 'OC'
for ct in cellTypes:
    for EorI in ['exc', 'inh']:
        # num of syns = total * exc * exc_long_range * fraction for this thalamic input 
        nsyns = [int(round(numSyns['M1'][M1BinPop] * fracSyns['M1'][EorI] * fracSyns['M1'][EorI+'_long'] * fracSyns['M1'][inputPop])) for M1BinPop in M1BinPops[ct]]
        
        # convergence = num syns / syns per conn * normalized input
        cmat[(inputPop, ct, EorI)] = [float(round(nsyn / synsPerConn[inputPop] * ninput)) for nsyn,ninput in zip(nsyns, normInput[('Hooks13_OC', ct)])]


# pprint(normInput)
# pprint(cmat)#[('OC', 'PT')])
# pprint(bins)

# ------------------------------------------------------------------------------------------------------------------
# 4) Weights 
## Thalamus; Hu et al 2016 ~0.6mV, Constantinople et al 2013 E+FS= ~0.57mV, LTS=1.44mV
wmat['TPO'] = 0.6
wmat['TVL'] = 0.6

## S1, S2, contraM1, M2, OC
wmat['S1'] = 0.5
wmat['S2'] = 0.5
wmat['cM1'] = 0.5
wmat['M2'] = 0.5
wmat['OC'] = 0.5

# ------------------------------------------------------------------------------------------------------------------
# 5) delays 
## (Hu, 2016, mouse S1) = 2.2ms; (Constantinople, 2013, in vivo rat S1) = 11,8,10,10,11 (range 10-30) ms
delays['TPO'] = 5
delays['TVL'] = 5

## S1, S2, contraM1, M2, OC
delays['S1'] = 5
delays['S2'] = 5
delays['cM1'] = 5
delays['M2'] = 5
delays['OC'] = 5

# ------------------------------------------------------------------------------------------------------------------
# 6) Num cells in each regions
## TPO and TVL: L5B PT ~ 44% probability (Constantinople 2013) and 42.5% (Bruno 2006)
## set to ~2x max convergence = 600; tradeoff between reducing correlation vs computation time 
maxConv = 300

numCells['TPO'] = 2*maxConv  
numCells['TVL'] = 2*maxConv

## S1, S2, contraM1, M2, OC
numCells['S1'] = 2*maxConv
numCells['S2'] = 2*maxConv
numCells['cM1'] = 2*maxConv
numCells['M2'] = 2*maxConv
numCells['OC'] = 2*maxConv

# ------------------------------------------------------------------------------------------------------------------
# 7) avg spontaneous firing rates 
## TPO and TVL
## LGN ~= 1-10 Hz (rat); VB ~= 2.9 +- 0.7 Hz (Hirata, 2006, rat); Thal = 8 Hz (Potjans 2013, rat/cat S1/V1);  
rates['TPO'] = [0,5]
rates['TVL'] = [0,5]

## S1, S2, contraM1, M2, OC
## in vivo mouse S1->M1 L23 = 0.1 Hz (Yamashita 2013); M1 IT = 0-10 Hz (Isomura, 2009; Jacob 2012; Li 2016)
rates['S1'] = [0,5]
rates['S2'] = [0,5]
rates['cM1'] = [0,5]
rates['M2'] = [0,5]
rates['OC'] = [0,5]

print ('TPO', bins[('TPO', 'IT')])
print ('TVL', bins[('TVL', 'IT')])
print ('S1', bins[('S1', 'IT')])
print ('S2', bins[('S2', 'IT')])
print ('cM1', bins[('cM1', 'IT')])
print ('M2', bins[('M2', 'IT')])
print ('OC', bins[('OC', 'IT')])


# ------------------------------------------------------------------------------------------------------------------
# save matrices
savePickle = 1
saveMat = 0

data = {'cmat': cmat, 'wmat': wmat, 'bins': bins, 'numSyns': numSyns, 'synsPerConn': synsPerConn, 
            'delays': delays, 'numCells': numCells, 'rates': rates}

if savePickle:
    import pickle
    with open(outFolder+'conn_long.pkl', 'wb') as fileObj:        
        pickle.dump(data, fileObj)

if saveMat:
    savemat(outFolder+'conn.mat', data)

# plot matrices
plotMat = 0
plotMat2D = 1
if plotMat: ax=plotMats()
if plotMat2D: ax=plotMats2D(labels=0, cbar=0)

plotPie = 0
if plotPie: ax=plotPies()


