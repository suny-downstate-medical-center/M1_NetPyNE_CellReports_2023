"""
M1 paper Figure 2 

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis

# ----------------------------------------------------------------
def fig_quiet():
    ''' Figure showing spontaneous activity: 
        - raster plot  
        - spike histogram
        - voltage traces
        - boxplot stats model and experiment'''


    # ---------------------------------------------------------------------------------------------------------------
    # Config

    raster = 0
    histogram = 0
    traces = 1
    stats_boxplot = 0   # boxplot of rates

    fontsiz = 16

    dataFolder = '../data/'
    if  raster or histogram or traces:    
        batchLabel = 'v56_batch18'  
        simLabel = 'v56_batch18_0_0_0'  # same as 'v56_batch7_7_3_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

    elif stats_boxplot: # 50 sec N=1
        
        timeRange = [1000, 5000] 
        include = excpops+[SOMsubset,PVsubset]
        xlim = [0, 40]  
        labelsModel = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'SOM L2-6', 'PV L2-6']
        multipleTrials = True
        loadAll = 1

        # single trial
        if not multipleTrials:
            batchLabel = 'v56_batch7'
            simLabel = 'v56_batch7_3_0_0'  
            sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        # multiple trials
        else:
            batchLabel = 'v56_batch7'
            simLabel =  ['v56_batch7_3_'+str(iseed)+'_'+str(iconn) for iseed in range(5) for iconn in range(5)]
            root = dataFolder+batchLabel+'/'
            
            if not loadAll:
        
                popGids, popNumCells = loadPopGids(dataFolder, labelsModel+['SOM2','SOM4','SOM5A','SOM5B','SOM6','PV2','PV4','PV5A','PV5B','PV6'])
                popGids['SOM L2-6'] = popGids['SOM2']+popGids['SOM4']+popGids['SOM5A']+popGids['SOM5B']+popGids['SOM6']
                popGids['PV L2-6'] = popGids['PV2']+popGids['PV4']+popGids['PV5A']+popGids['PV5B']+popGids['PV6']

                # load data from all sim files                
                statDataAll = [[]]*len(simLabel)
                for isim,simLab in enumerate(simLabel):  # get data from each sim 
                    filename = root+simLab+'.json'
                    print(filename)
                    sim,data,out = utils.plotsFromFile(filename, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, popColors=popColors)

                    statDataAll[isim] = []
                    for isubset, subset in enumerate(labelsModel):
                        cellGids = popGids[subset]
                        dfSpks, _, _ = analysis.utils.getSpktSpkid(cellGids, timeRange, sim)  # get quiet spikes
                        dfRates = dfSpks.groupby("spkid").count().div((timeRange[1]-timeRange[0])/1000.0)   
                        dfRates = dfRates.reindex(cellGids, fill_value=0.0)
                        statDataAll[isim].append(list(dfRates.spkt))
                        print(subset, len(dfRates.spkt), np.mean(dfRates.spkt), np.median(dfRates.spkt), scipy.stats.iqr(dfRates.spkt))

                # save
                statData = {'statDataAll': statDataAll} #,  'gidsData': gidsData, 'ynormsData': ynormsData}
                with open(root+'%s_statDataAll_boxplot.json'%(simLabel[0][:-2]), 'w') as fileObj:
                    json.dump(statData, fileObj)

                # load All        
            else:
                with open(root+'%s_statDataAll_boxplot.json'%(simLabel[0][:-2]), 'r') as fileObj:
                    statData = json.load(fileObj)
                filename = root+simLabel[0]+'.json'


    plt.style.use('seaborn-ticks') 

    # ---------------------------------------------------------------------------------------------------------------
    # raster
    if raster:
        timeRange = [2000, 4000] 
        include = allpops
        orderBy = ['pop', 'y']


        from netpyne.analysis.spikes_legacy import plotRaster
        fig1 = plotRaster(include=include, timeRange=timeRange, labels='overlay', 
            popRates=0, orderInverse=True, lw=0, markerSize=3.5, marker='.', popColors=popColors, 
            showFig=0, saveFig=0, figSize=(8.5*0.9*0.75, 7), orderBy=orderBy)# 
        ax = plt.gca()

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner
        plt.xticks([2000, 3000, 4000], ['2', '3', '4'], fontsize=fontsiz)
        plt.yticks([0, 5000, 10000], [0, 5000, 10000], fontsize=fontsiz)
        
        plt.ylabel('Neuron ID', fontsize=fontsiz) #Neurons (ordered by NCD within each pop)')
        plt.xlabel('Time (s)', fontsize=fontsiz)
        
        plt.title('')
        filename='%s%s_raster_%d_%d_%s_small.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
        plt.savefig(filename, dpi=600)


    # ---------------------------------------------------------------------------------------------------------------
    # histogram
    if histogram:

        timeRange = [2000, 4000] 
        include = allpops
        orderBy = ['pop', 'y']

        binSize=5
        measure = 'count'
        graphType = 'bar'

        from netpyne.analysis.spikes_legacy import plotSpikeHist
        fig_hist = plotSpikeHist(include=excpops, timeRange=timeRange, figSize=(8.5*0.9*0.75, 2.5), 
            popColors=popColors, legend=False, showFig=0, saveFig=0, linewidth=0.5, binSize=binSize, graphType='bar', 
            axis=True, measure=measure, scalebarLoc='upper center')

        ax=plt.gca()
        ax.get_legend().remove()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15, top=1.0, right=0.9, left=0.15)
        plt.yticks(fontsize=fontsiz)
        plt.xticks([2000, 3000, 4000], ['2', '3', '4'], fontsize=fontsiz)

        filename='%s%s_spikehist_%d_%d_bin-%d_%s_%s.png'%(root, simLabel, timeRange[0], timeRange[1], binSize, measure, graphType)
        plt.savefig(filename, dpi=600)


    # ---------------------------------------------------------------------------------------------------------------
    # traces
    if traces:
        # to plot all cells recorded
        batchLabel = 'v56_batch23'  
        simLabel = 'v56_batch23_merged'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        sim.cfg.recordTraces = {'V_soma': {}}

        fontsiz = 20
        timeRange = [2000, 4000]


        # plot subset of cells with the correct format 
        from netpyne.support.scalebar import add_scalebar

        def addScaleBar(timeRange=timeRange, loc=1):
            ax = plt.gca()
            sizex =  (timeRange[1]-timeRange[0])/20.0
            add_scalebar(ax, hidex=False, hidey=True, matchx=False, matchy=True, sizex=sizex, sizey=None, unitsx='ms', unitsy='mV', scalex=1, scaley=1, loc=loc, pad=-1, borderpad=0.5, sep=4, prop=None, barcolor="black", barwidth=3)
            plt.axis('off')

        # IT/PT cells for quiet
        cells = [4504] 
        timeRange = [2000, 4000]
        recordStep = 0.1

        t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep)

        for gid in cells: 
            fullTrace = sim.allSimData['V_soma']['cell_'+str(gid)]
            vtrace = np.array(fullTrace[int(timeRange[0]/recordStep):int(timeRange[1]/recordStep)])
            
            plt.figure(figsize=(15*1.5*0.5, 2.5))
            plt.plot(t[:len(vtrace)], vtrace, linewidth=1, color=popColors['IT5B'])
            baseline = min(vtrace)
            plt.ylim(baseline, baseline+60) # truncate AT 60mV above baseline!
            plt.axis('off')

            ax = plt.gca()
            plt.title('')
            plt.tight_layout()
            plt.savefig('%s/%s_V_soma_cell_%d_%d_%d_x1.5_axis.png' % (root, simLabel, gid, timeRange[0], timeRange[1]), dpi=200)


    # ---------------------------------------------------------------------------------------------------------------
    # stats boxplot
    if stats_boxplot:    

        fontsiz = 20

        # Model stats - single trial (N=1)
        if not multipleTrials:
            fig1,modelData = sim.analysis.plotSpikeStats(include=include, figSize=(8,4), timeRange=timeRange, xlim=xlim,
                stats = ['rate'], legendLabels=labelsModel, fontSize = fontsiz, popColors=popColors, showFig=0, dpi=300, saveFig=False)
        
        # Model stats - multiple trials (N=25)
        else:
            
            minValue = 0.01
            # combine data
            statDataAll = statData['statDataAll']

            statDataCombined = {}
            for ipop in range(len(include)):
                statDataCombined[ipop] = []
                for isim in range(len(simLabel)):
                    nonzeroData = statDataAll[isim][ipop]
                    nonzeroData = [x for x in nonzeroData if x>=minValue]
                    statDataCombined[ipop].extend(nonzeroData)

            modelData = {'statData': statDataCombined}
            statData = modelData['statData']

            # save combiend data
            with open(root+'%s_statData_boxplot.json'%(batchLabel), 'w') as fileObj:
                json.dump(statData, fileObj)

        
        # Combine model stats with experimental data stats 
        expData = {}
        from scipy.io import loadmat
        
        ## Este18 L5 RS and PV baseline  
        matData = loadmat(dataFolder+'/Este18/data.mat')  # [0][0][2][0] = RS; [0][0][2][0] = PV+

        expData['Este18_L5_RS_baseline'] = list(matData['baseline'][0][0][0][0]) + \
                                        list(matData['sound_base'][0][0][0][0]) + \
                                        list(matData['vibration_base'][0][0][0][0]) + \
                                        list(matData['movement_base'][0][0][0][0])
                                    
        expData['Este18_PV_baseline'] = list(matData['baseline'][0][0][2][0]) + \
                                        list(matData['sound_base'][0][0][2][0]) + \
                                        list(matData['vibration_base'][0][0][2][0]) + \
                                        list(matData['movement_base'][0][0][2][0]) 
        
        ## Zagh15 L5 RS and PV baseline 
        with open(dataFolder+'/Zagh15/spike_rates.json', 'r') as f:  #
            jsdonData = json.load(f)
            expData['Zagh15_L5_RS_baseline'] = jsdonData['pre-stimulus']['RSU']
            expData['Zagh15_PV_baseline'] = jsdonData['pre-stimulus']['FSU']

        ## Schi15 
        import schi15
        df = schi15.readExcelFiringRatesAndMetada()
        df.dropna()

        expData['Schi15_L5B_RS_baseline'] = list(df.query("condition == 'main' and code!='23'").quietFR)  # main, quiet (medium VL, medium ih)
        expData['Schi15_IT2_RS_baseline'] = list(df.query("condition=='main' and code=='23'").quietFR)  # main, quiet (medium VL, medium ih)
        expData['Schi15_IT5B_RS_baseline'] = list(df.query("condition=='main' and cell_class=='IT'").quietFR)  # main, quiet (medium VL, medium ih) 
        expData['Schi15_PT5B_RS_baseline']  = list(df.query("condition=='main' and cell_class=='PT'").quietFR)  # main, quiet (medium VL, medium ih) N=3! remove
        

        ## Dacr19 - waiting

        ## LiCh15 (ALM) L5 RS and PV baseline
        matData = loadmat(dataFolder+'/LiCh15/ALM_compiled_all_data.mat')
        expData['LiCh15_L5_RS_baseline'] = []
        expData['LiCh15_PV_baseline'] = []
        for i, (x1, x2) in enumerate(zip(list(matData['spk_count_yes_all'][:, 4]),list(matData['spk_count_no_all'][:, 4]))):
            if not np.isnan(np.mean([x1,x2])):
                if matData['cellType_all'][i, 0] == 1:  # Pyr   
                    expData['LiCh15_L5_RS_baseline'].append(np.mean([x1,x2]))
                elif matData['cellType_all'][i,0] == 2: # FS
                    expData['LiCh15_PV_baseline'].append(np.mean([x1,x2]))


        ## Econ18 (ALM) L5 PT baseline
        matData = loadmat(dataFolder+'/Econ18/Econ18.mat')
        expData['Econ18_PT5B_RS_baseline'] = [x[0] for x in list(matData['ratePTdown'])+list(matData['ratePTup'])]


        # ensure exp data values are >= minValue
        for k,v in expData.items():
            vmin = [x for x in v if x >= minValue]
            expData[k] = vmin

        # combine model and experimental data
        
        ## statData[0] - IT2
        ## statData[1] - IT4
        ## statData[2] - IT5A
        ## statData[3] - IT5B
        ## statData[4] - PT5B


        statData =    [statData[0]] \
                    + [expData['Schi15_IT2_RS_baseline']] \
                    + [statData[1]]+[statData[2]]+[statData[3]] \
                    + [expData['Schi15_IT5B_RS_baseline']] \
                    + [statData[4]] \
                    + [expData['Econ18_PT5B_RS_baseline']] \
                    + [statData[3]+statData[4]] \
                    + [expData['Schi15_L5B_RS_baseline']] \
                    + [statData[2]+statData[3]+statData[4]] \
                    + [expData['Este18_L5_RS_baseline'], expData['Zagh15_L5_RS_baseline'], expData['LiCh15_L5_RS_baseline']] \
                    + [statData[5]] + [statData[6]]                     

        labels = [  'L2/3 IT',
                    'L2 IT \n(Schiemann 2015)', \
                    'L4 IT', 'L5A IT', 'L5B IT', \
                    'L5B IT \n(Schiemann 2015)', \
                    'L5B PT', \
                    'L5B PT \n(Economo 2018)', \
                    'L5B IT,PT', \
                    'L5B IT,PT \n(Schiemann 2015)', \
                    'L5 IT,PT', \
                    'L5 IT,PT \n(Estebanez 2018)', 'L5 IT,PT \n(Zagha 2015)', 'L5 IT,PT \n(Li 2015)',\
                    'L6 IT', 'L6 CT']


        labels_orig = [  'IT2/3',
                    'IT2 \n(Schiemann 2015)', \
                    'IT4', 'IT5A', 'IT5B', \
                    'IT5B \n(Schiemann 2015)', \
                    'PT5B', \
                    'PT5B \n(Economo 2018)', \
                    'IT5B+PT5B', \
                    'IT5B+PT5B \n(Schiemann 2015)', \
                    'IT5A+IT5B+PT5B', \
                    'IT5A+IT5B+PT5B \n(Estebanez 2018)', 'IT5A+IT5B+PT5B \n(Zagha 2015)', 'IT5A+IT5B+PT5B \n(Li 2015)',\
                    'IT6', 'CT6']
        
        popColors['IT5A+IT5B+PT5B'] = '#BC01FF'
        popColors['IT5B+PT5B'] = '#800080'
        popColors['IT2/3'] = popColors['IT2']

        # update pop color labels
        for i,l in enumerate(labels_orig):
            if '(' not in l:
                popColors[labels[i]] = popColors[l]

        colors = [popColors[p] if p in popColors else '#999999' for p in labels]


        # plot boxplot comparing model and experimental data 
        plt.figure(figsize=(10,14))
        meanpointprops = dict(marker = (5, 1, 0), markeredgecolor = 'black', markerfacecolor = 'white')
        
        flierprops = dict(marker='.', markerfacecolor='gray', markersize=2,linestyle='none', markeredgecolor='gray')
        
        bp=plt.boxplot(statData[::-1], labels=labels[::-1], notch=False, flierprops=flierprops, meanprops=meanpointprops, whis=1.5, widths=0.6, vert=False, showmeans=True, patch_artist=True)  #labels[::-1] #positions=np.array(range(len(statData)))+0.4,

        # mean and std in errbar
        capsize = 8
        lw = 2
        elw=3.0
        marker = 'o'
        ms = 15
        mec = 'black'

        x = [np.mean(v) for v in statData[::-1]]
        xerr = [np.std(v) for v in statData[::-1]]
        y = np.array(range(len(x)))[::-1]

        plt.xlabel('Rate (Hz)', fontsize=fontsiz)
        plt.ylabel('', fontsize = fontsiz)
        plt.subplots_adjust(left=0.3,right=0.95, top=0.9, bottom=0.1)

        icolor=0
        borderColor = 'k'
        for i in range(0, len(bp['boxes'])):
            icolor = i
            bp['boxes'][i].set_facecolor(colors[::-1][icolor])
            bp['boxes'][i].set_linewidth(2)
            # we have two whiskers!
            bp['whiskers'][i*2].set_color(borderColor)
            bp['whiskers'][i*2 + 1].set_color(borderColor)
            bp['whiskers'][i*2].set_linewidth(2)
            bp['whiskers'][i*2 + 1].set_linewidth(2)
            bp['medians'][i].set_color(borderColor)
            bp['medians'][i].set_linewidth(3)
            for c in bp['caps']:
                c.set_color(borderColor)
                c.set_linewidth(2)

        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.tick_params(axis='x', length=0)
        #ax.tick_params(axis='y', direction='out')
        ax.grid(axis='x', color="0.9", linestyle='-', linewidth=1)
        ax.set_axisbelow(True)
        if xlim: ax.set_xlim(xlim)
        axisFontSize(ax, fontsiz-4)
        

        filename = root+batchLabel+'_boxplot_%d_%d'%(timeRange[0], timeRange[1])
        #filename = root+batchLabel+'_errbar_%d_%d'%(timeRange[0], timeRange[1])
        plt.savefig(filename, dpi=600)


    # stats and rank-test p-value
    for i,stat in enumerate(statData):
        print('%s: %.2f, %.2f, %.2f, %.2f' % (labels[i], np.mean(stat), np.std(stat), np.median(stat), scipy.stats.iqr(stat)))

    # L5 IT,PT exp1 vs exp2
    scipy.stats.mannwhitneyu(statData[11], statData[12])

    # L5 IT,PT model vs exp2
    scipy.stats.mannwhitneyu(statData[10], statData[12])



# Main code
if __name__ == '__main__': 
    fig_quiet()
