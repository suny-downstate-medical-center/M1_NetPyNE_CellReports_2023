"""
utils.py 

General functions to analyse simulation data

Contributors: salvadordura@gmail.com
"""
import json
import pickle
import numpy as np
from pylab import *
from itertools import product
import pandas as pd
import seaborn as sb
from scipy import stats
from pprint import pprint
from collections import OrderedDict
try:
  basestring
except NameError:
  basestring = str

from netpyne import specs

def toPandas(params, data):
    if 'simData' in data[data.keys()[0]]:
        rows = [list(d['paramValues'])+[s for s in d['simData'].values()] for d in data.values()]
        cols = [str(d['label']) for d in params]+[s for s in data[data.keys()[0]]['simData'].keys()]
    else:
        rows = [list(d['paramValues'])+[s for s in d.values()] for d in data.values()]
        cols = [str(d['label']) for d in params]+[s for s in data[data.keys()[0]].keys()]
    
    df = pd.DataFrame(rows, columns=cols) 
    df['simLabel'] = data.keys()

    colRename=[]
    for col in list(df.columns):
        if col.startswith("[u'"):
            colName = col.replace(", u'","_'").replace("[u","").replace("'","").replace("]","").replace(", ","_")
            colRename.append(colName)
        else: 
            colRename.append(col)
    print(colRename) 
    df.columns = colRename

    return df


def setPlotFormat(numColors=8):
    # plt.style.use('ggplot')
    #plt.style.use(['dark_background', 'presentation'])
    plt.style.use('seaborn-whitegrid')

    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 'large'

    NUM_COLORS = numColors
    colormap = plt.get_cmap('nipy_spectral')
    colorlist = [colormap(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    
    # add red and blue at the beginning
    # colorlist.insert(0, [1, 0, 0])
    # colorlist.insert(0, [0, 0, 1])

    plt.rc('axes', prop_cycle=(cycler('color', colorlist)))


def compare(source_file, target_file, source_key=None, target_key=None):
    from deepdiff import DeepDiff 
    with open(source_file, 'r') as fileObj:
        if source_file.endswith('.json'):
            source = json.load(fileObj, object_pairs_hook=OrderedDict)
        elif source_file.endswith('.pkl'):
            source = pickle.load(fileObj)
    if source_key: source = source[source_key]

    with open(target_file, 'r') as fileObj:
        if target_file.endswith('.json'):
            target = json.load(fileObj, object_pairs_hook=OrderedDict)
        elif source_file.endswith('.pkl'):
            target = pickle.load(fileObj)
    if target_key: target = target[target_key]
    
    ddiff = DeepDiff(source, target)
    pprint(ddiff)
    return ddiff


def readBatchData(dataFolder, batchLabel, loadAll=False, saveAll=True, vars=None, maxCombs=None, listCombs=None):
    # load from previously saved file with all data
    if loadAll:
        print('\nLoading single file with all data...')
        filename = '%s/%s/%s_allData.json' % (dataFolder, batchLabel, batchLabel)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj, object_pairs_hook=OrderedDict)
        params = dataLoad['params']
        data = dataLoad['data']
        return params, data

    if isinstance(listCombs, basestring):
        filename = str(listCombs)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj)
        listCombs = dataLoad['paramsMatch']

    # read the batch file and cfg
    batchFile = '%s/%s/%s_batch.json' % (dataFolder, batchLabel, batchLabel)
    with open(batchFile, 'r') as fileObj:
        b = json.load(fileObj)['batch']

    # read params labels and ranges
    params = b['params']

    # reorder so grouped params come first
    preorder = [p for p in params if 'group' in p and p['group']]
    for p in params:
        if p not in preorder: preorder.append(p)
    params = preorder

    # read vars from all files - store in dict 
    if b['method'] == 'grid':
        labelList, valuesList = zip(*[(p['label'], p['values']) for p in params])
        valueCombinations = product(*(valuesList))
        indexCombinations = product(*[range(len(x)) for x in valuesList])
        data = {}
        print('Reading data...')
        missing = 0
        for i,(iComb, pComb) in enumerate(zip(indexCombinations, valueCombinations)):
            if (not maxCombs or i<= maxCombs) and (not listCombs or list(pComb) in listCombs):
                print(i, iComb)
                # read output file
                iCombStr = ''.join([''.join('_'+str(i)) for i in iComb])
                simLabel = b['batchLabel']+iCombStr
                outFile = b['saveFolder']+'/'+simLabel+'.json'
                try:
                    with open(outFile, 'r') as fileObj:
                        output = json.load(fileObj, object_pairs_hook=OrderedDict)
                    # save output file in data dict
                    data[iCombStr] = {}  
                    data[iCombStr]['paramValues'] = pComb  # store param values
                    if not vars: vars = output.keys()

                    for key in vars:
                        if isinstance(key, tuple):
                            container = output
                            for ikey in range(len(key)-1):
                                container = container[key[ikey]]
                            data[iCombStr][key[1]] = container[key[-1]]

                        elif isinstance(key, basestring): 
                            data[iCombStr][key] = output[key]

                except:
                    print('... file missing')
                    missing = missing + 1
                    output = {}
            else:
                missing = missing + 1

        print('%d files missing' % (missing))

        # save
        if saveAll:
            print('Saving to single file with all data')
            filename = '%s/%s/%s_allData.json' % (dataFolder, batchLabel, batchLabel)
            dataSave = {'params': params, 'data': data}
            with open(filename, 'w') as fileObj:
                json.dump(dataSave, fileObj)
        
        return params, data



def plotsFromFile(filename, raster=True, stats=True, rates=False, syncs=False, hist=True, psd=True, grang=True, traces=True, plotAll=False, 
                timeRange=None, textTop='', include={'raster':[]}, popColors=None, orderBy='gid'):
    from netpyne import specs
    from collections import OrderedDict

    # try:
    if filename.endswith('.json'):
        with open(filename, 'rb') as fileObj:
            data = json.load(fileObj, object_pairs_hook=OrderedDict)
    elif filename.endswith('.pkl'):
        with open(filename, 'rb') as fileObj:
            data = pickle.load(fileObj)
    # except:
    #     print('Error opening file %s' % (filename))
    #     return 0
    
    sim,data, out = plotsFromData(data, raster=raster, stats=stats, rates=rates, syncs=syncs, hist=hist, psd=psd, traces=traces, grang=grang,
                             plotAll=plotAll, timeRange=timeRange, textTop=textTop, include=include, popColors=popColors, orderBy=orderBy)
    return sim,data, out


def plotsFromData(data, textTop = '', raster=True, stats=False, rates=False, syncs=False, hist=True, psd=True, traces=True, grang=True, 
                plotAll=0, timeRange=None, include=None, popColors=None, orderBy='gid'):
    import matplotlib.pyplot as plt
    from netpyne import specs,sim
    if 'simConfig' in data:
        cfg = specs.SimConfig(data['simConfig'])
    else:
        cfg = specs.SimConfig()
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data, instantiate=False)
    sim.setSimCfg(cfg)
    
    import collections
    temp = collections.OrderedDict()
    if len(include['raster']) == 0:
        include['raster'] = popColors.keys()

    # order pops by keys in popColors
    for k in include['raster']:
        if k in sim.net.params.popParams:
            temp[k] = sim.net.params.popParams.pop(k)

    # add remaining pops at the end
    for k in list(sim.net.params.popParams.keys()):
        temp[k] = sim.net.params.popParams.pop(k)
    
    sim.net.params.popParams = temp
    
    try:
        print('Cells created: ',len(sim.net.allCells))
    except:
        sim.net.createPops()     
        sim.net.createCells()
        sim.setupRecording()
        sim.cfg.createNEURONObj = False
        sim.gatherData() 

    sim.allSimData = data['simData']

    plt.style.use('seaborn-ticks') 

    out = None

    if raster or plotAll:
        fig1 = sim.analysis.plotRaster(include=include['raster'], timeRange=timeRange, labels='overlay', popRates=False, orderInverse=True, lw=0, markerSize=5, 
            marker='.', popColors=popColors, showFig=0, saveFig=0, figSize=(10,8), orderBy=orderBy)
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.ylabel('Neuron number')
        plt.title('')
        plt.savefig(cfg.filename+'_raster_%d_%d_%s.png'%(timeRange[0], timeRange[1], orderBy), dpi=600)

    if stats or plotAll:
        filename = cfg.filename+'_%d_%d'%(timeRange[0], timeRange[1])
        fig1 = sim.analysis.plotSpikeStats(include=include['stats'], figSize=(4,8), timeRange=timeRange, stats = ['rate', 'isicv'], popColors=popColors, showFig=0, saveFig=filename)

    if rates or plotAll:
        midpoint = (timeRange[1]+timeRange[0]) / 2.0
        timeRanges = [[timeRange[0], midpoint], [midpoint, timeRange[1]]]
        colors = [popColors[inc] if inc in popColors else (0,0,0) for inc in include['rates']]
        out = sim.analysis.plotRates(include=include['rates'], timeRanges=timeRanges, figSize=(4,2), timeRangeLabels=['Pre stim', 'Post stim'], colors=colors, showFig=0, saveFig=1)

    if syncs or plotAll:
        midpoint = (timeRange[1]+timeRange[0]) / 2.0
        timeRanges = [[timeRange[0], midpoint], [midpoint, timeRange[1]]]
        colors = [popColors[inc] if inc in popColors else (0,0,0) for inc in include['syncs']]
        fig1 = sim.analysis.plotSyncs(include=include['syncs'], timeRanges=timeRanges, figSize=(5,4), timeRangeLabels=['Pre stim', 'Post stim'], colors=colors, showFig=0, saveFig=1)

    if hist or plotAll:
        fig2 = sim.analysis.plotSpikeHist(include=include['hist'], yaxis='rate', binSize=5, graphType='bar', timeRange=timeRange,  popColors=popColors, showFig=False, saveFig=0, figSize=(12,8))
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.savefig(cfg.filename+'_spikeHist_replot_%d_%d.png'%(timeRange[0], timeRange[1]))

    if psd or plotAll:
        popColors['allCells'] = 'k'
        fig3 = sim.analysis.plotRatePSD(include=include['psd'], timeRange=timeRange, Fs=160, smooth=16 , showFig=0, saveFig=0, popColors=popColors, figSize=(11,5))
        # ylim=[-55, 0]
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.savefig(cfg.filename+'_spikePSD_%d_%d.png'%(timeRange[0], timeRange[1]), dpi=300)

    if traces or plotAll:
        colors = popColors; colors=['r','b']
        fig4 = sim.analysis.plotTraces(include=include['traces'], timeRange=timeRange, overlay = True, oneFigPer = 'trace', rerun = False, ylim=[-90, 30], axis='on', figSize = (12,4), saveData = None, saveFig = 1, showFig = 0)

    return sim, data, out


def getSpksPop(data, timeRange=None):
    from netpyne import specs,sim
    cfg = specs.SimConfig(data['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data)
    sim.setSimCfg(cfg)
    sim.net.createPops()     
    sim.net.createCells()
    sim.gatherData() 

    sim.allSimData = data['simData']
    spkts = sim.allSimData['spkt']
    spkids = sim.allSimData['spkid']

    popLabels = sim.net.allPops.keys()
    gidPops = [cell['tags']['pop'] for cell in sim.net.allCells]
    popNumCells = [float(gidPops.count(pop)) for pop in popLabels]
    numCellsPop = {pop:float(gidPops.count(pop)) for pop in popLabels}
    spktPop = {}
    for pop, popNum in zip(popLabels, popNumCells):
        if timeRange:
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop if timeRange[0]<=spkt<=timeRange[1]]
        else:
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop]

    return spktPop, numCellsPop


def readBatchUpdatePopRates(dataFolder, batchLabel, timeRange=None, saveAll=True, loadAll=False):
    # load from previously saved file with all data
    if loadAll:
        from netpyne import specs
        print('\nLoading single file with all data...')
        filename = '%s/%s/%s_allData_popRates.json' % (dataFolder, batchLabel, batchLabel)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj, object_pairs_hook=OrderedDict)
        params = dataLoad['params']
        data = dataLoad['data']
        return params, data

    # load single sim to get netParams and simConfig
    params, data = readBatchData(dataFolder, batchLabel, maxCombs=1)  

    # recreate net to get cell+pop params
    from netpyne import specs,sim
    cfg = specs.SimConfig(data.values()[0]['simConfig'])
    cfg.createNEURONObj = False
    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data.values()[0])
    sim.setSimCfg(cfg)
    sim.net.createPops()     
    sim.net.createCells()
    sim.gatherData() 

    # caclulate pop stats
    if not timeRange: timeRange = [0, cfg.duration]
    tsecs = (timeRange[1]-timeRange[0])/1e3 
    popLabels = sim.net.allPops.keys()
    gidPops = [cell['tags']['pop'] for cell in sim.net.allCells]
    popNumCells = [float(gidPops.count(pop)) for pop in popLabels]
    #numCellsPop = {pop:float(gidPops.count(pop)) for pop in popLabels}
    popLabels = sim.net.allPops.keys()

    # load spkt and spkid
    params, data = readBatchData(dataFolder, batchLabel, vars=[('simData','spkt'),('simData','spkid'),('simData','popRates')], saveAll=False)
    for d in data.values():
        spkts = d['spkt']
        spkids = d['spkid']

        spktPop = {}
        for pop, popNum in zip(popLabels, popNumCells):
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop if timeRange[0]<=spkt<=timeRange[1]]
            d['popRates'][pop] = len(spktPop[pop]) / popNum / tsecs  # calcualte new popRate
            d['spkt'] = None  # remove to save space
            d['spkid'] = None  

    # save
    if saveAll:
        print('Saving to single file with all data')
        filename = '%s/%s/%s_allData_popRates.json' % (dataFolder, batchLabel, batchLabel)
        dataSave = {'params': params, 'data': data}
        with open(filename, 'w') as fileObj:
            json.dump(dataSave, fileObj)

    return params, data





# ---------
# Function to plot boxplot of model vs experiment
def stats_boxplot(dfStats, x, y, hue, color, overlayLine, quietExp, moveExp, quietModel, moveModel, figSize, fontsize, palette=None, order=None):
    #my_pal = {'Model': color, 'Experiment': '#999999'}
    if palette:
        my_pal = palette
    else:
        my_pal = {'Model': 'royalblue', 'Experiment': 'orange'}

    fig=plt.figure(figsize=figSize)
    ax = sb.boxplot(x=x, y=y, hue=hue, data=dfStats, palette=my_pal, showfliers=False, order=order)#, boxprops=dict(alpha=.75), linewidth=0.5, width=0.3)  #
    
    #plt.title (title, fontsize=fontsiz)
    ax.get_legend().remove()
    plt.tight_layout()
    plt.ylabel('Firing Rate (Hz)', fontsize=fontsize)
    plt.xlabel('')
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    fig.subplots_adjust(left=0.13)

    n = len(my_pal)
    width=0.8/n

    onlyStar = True

    if overlayLine:
        for i,(title,color) in enumerate(my_pal.items()):
            if 'Experiment' in title:
                data = [np.mean(quietExp), np.mean(moveExp)]
            elif isinstance(quietModel, dict):
                data = [np.mean(quietModel[title]), np.mean(moveModel[title])]
            else:
                data = [np.mean(quietModel), np.mean(moveModel)]
            
            if onlyStar:
                plt.scatter(-0.4+(i*width)+width/2, data[0], marker='*', s=80, color='white', edgecolors='black') 
                plt.scatter(1-0.4+(i*width)+width/2, data[1], marker='*', s=80, color='white', edgecolors='black') 
            else:
                plt.plot([-0.4+(i*width)+width/2, 1-0.4+(i*width)+width/2], data, marker='*', markersize=10, linewidth=2, mec='black', mfc=color, color=color) 

    # fixed axis scale
    maxRate = 30
    plt.ylim(0, maxRate)



# ---------
# Function to plot lineplot with errorbar of model vs experiment
#
# Note: individual model line plots are taken from statDataAllQuiet and statDataAllMove, which contain all data points 
# ie. filtering of >0.01 Hz has not been applied to these raw variables, so may not correspond to quietExp and moveExp
def stats_lineplot(title, quietExp, moveExp, quietModel, moveModel, statDataAllQuiet, statDataAllMove, modelInd, figSize, fontsize, printOutput=False, linestyle='-', colors={}):
                    
    fig=plt.figure(figsize=figSize)

    avgFunc = np.mean 
    varFunc = np.std 

    capsize = 8 # set to 0 for legend
    lw = 2
    elw=1.0
    marker = 'o'
    ms = 15  # set to 10 for legend
    mec = 'black'
    labelExp = 'Experiment'
    

    # plot model
    if not isinstance(quietModel, dict) and not isinstance(moveModel, dict)  :  # convert to common dict format
        quietModel = {'Model': quietModel}
        moveModel = {'Model': moveModel}
        colors = {'Model': 'blue'}

    for k in quietModel:
        ax2 = plt.errorbar([0.2,0.8], [avgFunc(quietModel[k]), avgFunc(moveModel[k])], yerr=[varFunc(quietModel[k]), varFunc(moveModel[k])], 
            linewidth=lw, capsize=capsize, color=colors[k], marker=marker, markersize=ms, mec='black', elinewidth=elw,alpha=1.0) #mec='blue', mec='black')#, mfc='blue') 
        for cap in ax2[1]:
            cap.set_markeredgewidth(elw)


    # plot experiment (move to top for legend)
    if linestyle != '-': # either show exp data in lower alpha; or not include exp data at all 
        pass 
        # ax1 = plt.errorbar([0.2,0.8], [avgFunc(quietExp), avgFunc(moveExp)], yerr=[varFunc(quietExp), varFunc(moveExp)], 
        #     linewidth=lw, capsize=capsize, color='orange', alpha=0.5, marker=marker, markersize=ms, mec='black', linestyle=linestyle, elinewidth=elw) #mec='orange',)#, mfc='orange') 
        # ax1[-1][0].set_linestyle(linestyle)

    else:
        ax1 = plt.errorbar([0.2,0.8], [avgFunc(quietExp), avgFunc(moveExp)], yerr=[varFunc(quietExp), varFunc(moveExp)], 
            linewidth=lw, capsize=capsize, color='orange', marker=marker, markersize=ms, mec='black', linestyle=linestyle, elinewidth=elw, alpha=1.0) #mec='orange',)#, mfc='orange') 
    
        for cap in ax1[1]:
            cap.set_markeredgewidth(elw)


    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


    trialsQuiet = []
    trialsMove = []

    for statTrialQuiet, statTrialMove in zip(statDataAllQuiet, statDataAllMove):

        plt.plot([0.2,0.8], [avgFunc(statTrialQuiet[modelInd]), avgFunc(statTrialMove[modelInd])],
        linewidth=0.25, color='royalblue', alpha=0.8)# capsize=capsize, elinewidth=0.25) #, capsize=capsize) 

        trialsQuiet.append(avgFunc(statTrialQuiet[modelInd]))
        trialsMove.append(avgFunc(statTrialMove[modelInd]))

    plt.xticks([0.2,0.8], ['Quiet', 'Movement'], fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.yticks([0,5,10,15,20,25],[0,5,10,15,20,25])
    plt.xlim(0,1)
    plt.ylim(0,25)
    plt.ylabel('Firing Rate (Hz)', fontsize=fontsize)
    plt.subplots_adjust(left=0.18, bottom=0.1)


    if printOutput:
            print(title)

            for k in quietModel:
                if len(k) > 0: print(k)
                print('MODEL: %s +- %s (across cells):   quiet=%.1f+-%.1f , move=%.1f+-%.1f' 
                % (avgFunc.__name__, varFunc.__name__, avgFunc(quietModel[k]), varFunc(quietModel[k]), avgFunc(moveModel[k]), varFunc(moveModel[k])))

            print('EXP: %s +- %s (across cells):   quiet=%.1f+-%.1f , move=%.1f+-%.1f' 
            % (avgFunc.__name__, varFunc.__name__, avgFunc(quietExp), varFunc(quietExp), avgFunc(moveExp), varFunc(moveExp)))
            


# -------------
#  statistical tests 
def stats_tests(title, quietExp, moveExp, quietModel, moveModel, sameLength=True):
    import scipy

    def pvalueStars(prob):
        if prob <= 0.001:
            probstar = '***' 
        elif prob <= 0.01:
            probstar = '**' 
        elif prob <= 0.05:
            probstar = '*' 
        else:
            probstar = 'ns'
        return probstar


    if sameLength:
        if len(quietModel) > len(moveModel):
            quietModel = quietModel[:len(moveModel)]
        else:
            moveModel = moveModel[:len(quietModel)]

    # Mann-Whitney U
    dquiet_mw, probquiet_mw = scipy.stats.mannwhitneyu(quietExp, quietModel) 
    dmove_mw, probmove_mw = scipy.stats.mannwhitneyu(moveExp, moveModel) 
    probquietstar_mw = pvalueStars(probquiet_mw)
    probmovestar_mw = pvalueStars(probmove_mw)

    # confidence interval / significance level, alpha; typically 0.05
    alpha = 0.05
    c_alpha = {0.2: 1.07, 0.15: 1.14, 0.1: 1.22, 0.05: 1.36, 0.025: 1.48, 0.01: 1.63, 0.005: 1.73}
    
    dstarquiet = c_alpha[alpha] * np.sqrt((len(quietExp)+len(quietModel) / len(quietExp)*len(quietModel)))   
    dstarmove = c_alpha[alpha] * np.sqrt((len(moveExp)+len(moveModel) / len(moveExp)*len(moveModel)))   

    # null hypothesis (same distribution) rejected if alpha >= prob
    # i.e. same distribution NOT rejected if alpha < prob
    
    print('\nMann-Whitney U:')
    print('%s QUIET, N_exp=%d, N_model=%d, D=%f, D*=%f, alpha=%f,  prob=%g %s' % (title, len(quietExp), len(quietModel), dquiet_mw, dstarquiet, alpha, probquiet_mw, probquietstar_mw))
    print('%s MOVE, N_exp=%d, N_model=%d, D=%f, D*=%f, alpha=%f,  prob=%g %s' % (title, len(moveExp), len(moveModel), dmove_mw, dstarmove, alpha, probmove_mw, probmovestar_mw))
    print('')
