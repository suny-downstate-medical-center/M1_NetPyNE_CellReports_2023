"""
M1 paper Figure 5

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis

# ----------------------------------------------------------------
def fig_narblock():
    ''' Figure NA-R block activity: 
        - raster plot of 1-7 sec 
        - traces 1-7 sec (compare to exp?)
        - stats (boxplot,line,scatter) comparing quite+move to exp
    '''

    # ---------------------------------------------------------------------------------------------------------------
    # Config

    raster = 0
    histogram = 0
    traces = 0
    stats = 1  # boxplot of rates
    
    fontsiz = 26

    dataFolder = '../../data/'

    # -------------------------------------------------------------------------------------------
    # Raster plot
    if  raster:  # 2 sec N=1  
        batchLabel = 'v56_batch22'  
        simLabel = 'v56_batch22_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        timeRange = [4000, 7000] #[2000, 4000]
        include = allpops
        orderBy = ['pop', 'y']

        #8.5*0.9*0.75 
        from netpyne.analysis.spikes_legacy import plotRaster
        fig1 = plotRaster(include=include, timeRange=timeRange, labels='overlay', 
            popRates=0, orderInverse=True, lw=0, markerSize=3.5, marker='.', popColors=popColors, 
            showFig=0, saveFig=0, figSize=(8.5*1.5*4/6, 7), orderBy=orderBy)# 
        ax = plt.gca()

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner
        plt.xticks([4000, 5000, 6000, 7000], ['4', '5', '6', '7'], fontsize=fontsiz)
        plt.yticks([0, 5000, 10000], [0, 5000, 10000], fontsize=fontsiz)
        
        plt.ylabel('Neuron ID', fontsize=fontsiz) #Neurons (ordered by NCD within each pop)')
        plt.xlabel('Time (s)', fontsize=fontsiz)
        
        plt.title('')
        filename='%s%s_raster_%d_%d_%s.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
        plt.savefig(filename, dpi=600)

    # -------------------------------------------------------------------------------------------
    # Histogram
    if  histogram: 

        batchLabel = 'v56_batch22'  
        simLabel = 'v56_batch22_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        fontsiz = 26
        timeRange = [4000, 7000]
        
        binSize=10
        measure = 'count'
        graphType = 'bar'

        from netpyne.analysis.spikes_legacy import plotSpikeHist
        fig_hist = plotSpikeHist(include=excpops, timeRange=timeRange, figSize=(8.5*1.5*(4/6), 2.5), 
            popColors=popColors, legend=False, showFig=0, saveFig=0, linewidth=0.5, binSize=binSize, graphType='bar', 
            axis=True, measure=measure, scalebarLoc='upper center')

        ax=plt.gca()
        ax.get_legend().remove()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.17, top=1.0, right=0.9, left=0.15)
        plt.yticks(fontsize=fontsiz)
        plt.xticks([4000, 5000, 6000, 7000], ['0', '1', '2', '3'], fontsize=fontsiz)

        filename='%s%s_spikehist_%d_%d_bin-%d_%s_%s.png'%(root, simLabel, timeRange[0], timeRange[1], binSize, measure, graphType)
        plt.savefig(filename, dpi=600)



    # -------------------------------------------------------------------------------------------
    # Traces plot
    if traces:
        batchLabel = 'v56_batch22'  
        simLabel = 'v56_batch22_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        fontsiz = 20
        timeRange = [4000, 7000]
        include = [2901, 5709] #[2900-2904; 5709-5728]
        colors = [popColors[p] for p in ['IT5A', 'PT5B']] #, 'CT6', 'S2', 'M2','PT5B', 'CT6', 'S2', 'M2' ]]

        include = [5709] #[2900-2904; 5709-5728]
        colors = [popColors[p] for p in ['PT5B']] #, 'CT6', 'S2', 'M2','PT5B', 'CT6', 'S2', 'M2' ]]


        fig4 = sim.analysis.plotTraces(include=include, timeRange=timeRange, colors=colors, 
            overlay=True, oneFigPer='trace', rerun=False, ylim=[-85, 90], axis='off', 
            figSize=(15*0.7,2.5), saveData=None, saveFig=0, showFig=0)

        ax = plt.gca()
        ax.get_legend().remove()
        plt.title('')
        plt.tight_layout()
        plt.savefig('%s%s_traces_%d_%d_PT-%d_x1.5.png'%(root, simLabel, timeRange[0], timeRange[1], include[0]), dpi=200)



    # ---------------------------------------------------------------------------------------------------------------
    # stats plots (boxplot, line plot, scatter)

    if stats:  # 50 sec N=1
        
        timeRangeQuiet = [1000, 5000]
        timeRangeMove = [5000, 9000]
        include = ['IT5B', 'PT5B', ['IT5B','PT5B']]
        xlim = [0, 40]  #[0,70] # ok to cut off max and flyers (Bill19)
        labelsModel = ['IT5B', 'PT5B', 'L5B', 'L5Benh', 'L5Bsupp']
        multipleTrials = True
        loadAll = 1

        # single trial
        if not multipleTrials:
            batchLabel = 'v56_batch22'
            simLabel = 'v56_batch22_0_0'
            sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
        
            # fix
            fig1, modelData = sim.analysis.plotSpikeStats(include=include, figSize=(8,4), timeRange=timeRangeMove, xlim=xlim,
               stats = ['rate'], includeRate0=True, legendLabels=None, fontSize = fontsiz, popColors=popColors, showFig=0, dpi=300, saveFig=False)
        
            statDataMove = modelData['statData'] 
        
        # multiple trials
        else:
            # ih = 1.0
            batchLabel = 'v56_batch22'
            simLabel = ['v56_batch22_' + str(iseed) + '_' + str(iconn) for iseed in range(5) for iconn in range(5)]

            batchLabelMain = 'v56_batch19'
            simLabelMain = ['v56_batch19_0_0_' + str(iseed) + '_' + str(iconn) for iseed in range(5) for iconn in range(5)]
            
            root = dataFolder + batchLabel + '/'
            rootMain = dataFolder + batchLabelMain + '/'
            
            plt.style.use('seaborn-ticks')
            
            # read data from sim files
            if not loadAll:
                
                # load data from all sim files                
                statDataAllMove = [] 
                statDataAllQuiet = [] 
                includeAll = [] 
                
                popGids, popNumCells = loadPopGids(dataFolder, labelsModel)

                for isim, simLab in enumerate(simLabel):  # get data from each sim 
                    statDataAllMove.append([]) 
                    statDataAllQuiet.append([])
                    includeAll.append([]) 
                    
                    filename = root + simLab + '.json'
                    print('Loading %s ... ' % (filename))
                
                    # recreate network
                    try: 
                        sim, data, out = utils.plotsFromFile(filename, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, popColors=popColors)
                    except:
                        filenameTemp = root + simLabel[isim - 1] + '.json' 
                        print('Missing file %s ; using %s instead for now ...' % (filename, filenameTemp))
                        sim, data, out = utils.plotsFromFile(filenameTemp, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, popColors=popColors)

                    # calculate enhanced vs suppressed to include in stats
                    print('Classifying enhanced vs suppressed neurons ...')
                    spkts, spkids = sim.allSimData['spkt'], sim.allSimData['spkid']
                    L5Bgids = popGids['L5B'] 
                    L5Benh = []
                    L5Bsupp = []

                    calculateEnhSup = True
                    if calculateEnhSup:
                        for gid in L5Bgids:
                            _, spktcell, _ = analysis.utils.getSpktSpkid([gid], [timeRangeQuiet[0], timeRangeMove[1]], sim)
                                            
                            rateQuiet = len([s for s in spktcell if timeRangeQuiet[0] < s < timeRangeQuiet[1]]) \
                                        * 1000 / (timeRangeQuiet[1] - timeRangeQuiet[0])
                            rateMove = len([s for s in spktcell if timeRangeMove[0] < s < timeRangeMove[1]]) \
                                        * 1000 / (timeRangeMove[1] - timeRangeMove[0])
                            if rateMove > 1.05 * rateQuiet:
                                L5Benh.append(gid)
                            elif rateMove < 0.95 * rateQuiet:
                                L5Bsupp.append(gid)
                                        
                    popGids['L5Benh'] = tuple(L5Benh)
                    popGids['L5Bsupp'] = tuple(L5Bsupp)

                    includeAll[-1] = include + [popGids['L5Benh'], popGids['L5Bsupp']]

                    # # calculate stats for different pops/subpops for move period
                    # _, dataMove = sim.analysis.plotSpikeStats(include=include, graphType='none', includeRate0=1, timeRange=timeRangeMove, stats=['rate'], showFig=0, saveFig=0)
                    # statDataAllMove[isim] = dataMove['statData']
                    # dfspks, spkt, spkid =
                    
                    # # calculate stats for different pops/subpops for quiet period
                    # _, dataQuiet = sim.analysis.plotSpikeStats(include=include, graphType='none', includeRate0=1, timeRange=timeRangeQuiet, stats=['rate'], showFig=0, saveFig=0)
                    # statDataAllQuiet[isim] = dataQuiet['statData']

                   
                    for subset in labelsModel:
                        cellGids = popGids[subset]                        
                        
                        dfSpksMove, _, _ = analysis.utils.getSpktSpkid(cellGids, timeRangeMove, sim)  # get move spikes 
                        dfSpksQuiet, _, _ = analysis.utils.getSpktSpkid(cellGids, timeRangeQuiet, sim)  # get quiet spikes
                        
                        dfRatesMove = dfSpksMove.groupby("spkid").count().div((timeRangeMove[1]-timeRangeMove[0])/1000.0)   
                        dfRatesQuiet = dfSpksQuiet.groupby("spkid").count().div((timeRangeQuiet[1]-timeRangeQuiet[0])/1000.0)   

                        # include cells with rate 0Hz so can compare quiet vs move
                        dfRatesMove = dfRatesMove.reindex(cellGids, fill_value=0.0)
                        dfRatesQuiet = dfRatesQuiet.reindex(cellGids, fill_value=0.0)

                        statDataAllMove[-1].append(list(dfRatesMove.spkt))
                        statDataAllQuiet[-1].append(list(dfRatesQuiet.spkt))

                        print(subset, len(dfRatesMove.spkt), len(dfRatesQuiet.spkt))
                    
                # save
                statData = {'statDataAllMove': statDataAllMove,  'statDataAllQuiet': statDataAllQuiet, 'includeAll': includeAll}
                with open(root+'%s_statDataAll_boxplot.json' % (simLabel[0][:-2]), 'w') as fileObj:
                    json.dump(statData, fileObj)

            # load All data from combined data file
            else:
                with open(root+'%s_statDataAll_boxplot.json'%(simLabel[0][:-2]), 'r') as fileObj:
                    statData = json.load(fileObj)
                filename = root + simLabel[0] + '.json'
                
                # also load main condition data to compare
                with open(rootMain+'%s_statDataAll_boxplot.json'%(simLabelMain[0][:-2]), 'r') as fileObj:
                    statDataMain = json.load(fileObj)
                filenameMain = root + simLabel[0] + '.json'


            statDataAllMove = statData['statDataAllMove']
            statDataAllQuiet = statData['statDataAllQuiet']
            includeAll = statData['includeAll']

            # # code to combine IT5B+PT5B into L5B
            # for i in range(1):
            #     statDataAllMove[i][2] = statDataAllMove[i][0]+statDataAllMove[i][1]
            #     statDataAllQuiet[i][2] = statDataAllQuiet[i][0]+statDataAllQuiet[i][1]

            # combine data
            minValueQuiet = 0.0 # 0.01  # min for 4 sec sims is 0.25
            minValueMove = 0.0 #0.01
            print('\nminValueQuiet = %.2f Hz; minValueMove = %.2f Hz' % (minValueQuiet, minValueMove))

            statDataMove = {}
            statDataQuiet = {}

            for ipop in range(len(statDataAllMove[0])):
                statDataMove[ipop] = []
                statDataQuiet[ipop] = []

                for isim in range(len(simLabel)):
                    moveRates = statDataAllMove[isim][ipop]
                    quietRates = statDataAllQuiet[isim][ipop]

                    try:
                        nonzeroDataQuiet, nonzeroDataMove = zip(*[(q,m) for (q,m) in zip(quietRates, moveRates) if q>=minValueQuiet and m>=minValueMove])
                    except:
                        nonzeroDataQuiet, nonzeroDataMove = [0.0], [0.0]

                    statDataQuiet[ipop].extend(nonzeroDataQuiet)
                    statDataMove[ipop].extend(nonzeroDataMove)


                # check this code below - use code from mthinact.py 
                # combine data Main
                
                statDataAllMoveMain = statDataMain['statDataAllMove']
                statDataAllQuietMain = statDataMain['statDataAllQuiet']
                includeAllMain = statDataMain['includeAll']
                
                minValue = 0.0
                statDataMoveMain = {}
                statDataQuietMain = {}

                for ipop in range(len(statDataAllMoveMain[0])):
                    statDataMoveMain[ipop] = []
                    statDataQuietMain[ipop] = []

                    for isim in range(len(simLabel)):
                        nonzeroData = statDataAllMoveMain[isim][ipop]
                        nonzeroData = [x for x in nonzeroData if x>=minValue]
                        statDataMoveMain[ipop].extend(nonzeroData)

                        nonzeroData = statDataAllQuietMain[isim][ipop]
                        nonzeroData = [x for x in nonzeroData if x>=minValue]
                        statDataQuietMain[ipop].extend(nonzeroData)
                
        
        # ------------------------------------------------------------
        # Combine model stats with experimental data stats 
        expData = {}
        from scipy.io import loadmat
    
        ## Schi15 
        import schi15
        dfSchi = schi15.readExcelFiringRatesAndMetada()
        dfSchi.dropna()

        expData['Schi15_L5B_RS_quiet'] = list(dfSchi.query("condition == 'main' and code!='23'").quietFR)  # main, quiet (medium VL, medium ih)
        expData['Schi15_IT2_RS_quiet'] = list(dfSchi.query("condition=='main' and code=='23'").quietFR)  
        expData['Schi15_IT5B_RS_quiet'] = list(dfSchi.query("condition=='main' and cell_class=='IT'").quietFR)  
        expData['Schi15_PT5B_RS_quiet']  = list(dfSchi.query("condition=='main' and cell_class=='PT'").quietFR) 
        expData['Schi15_L5Benh_RS_quiet'] = list(dfSchi.query("condition == 'main' and type=='enh'").quietFR)  
        expData['Schi15_L5Bsupp_RS_quiet'] = list(dfSchi.query("condition == 'main' and type=='supp'").quietFR)  
        expData['Schi15_L5Benh5%_RS_quiet'] = list(dfSchi.query("condition == 'main' and moveFR>1.05*quietFR").quietFR)  
        expData['Schi15_L5Bsupp5%_RS_quiet'] = list(dfSchi.query("condition == 'main' and moveFR<0.95*quietFR").quietFR)
        
        expData['Schi15_L5Benh5%_RS_move'] = list(dfSchi.query("condition == 'main' and moveFR>1.05*quietFR").moveFR) # main, move (high VL, low ih)

        expData['Schi15_L5B_RS_quiet_na-block'] = list(dfSchi.query("condition == 'na-block' and code!='23'").quietFR)  # na-block, quiet (medium VL, high ih)
        expData['Schi15_IT2_RS_quiet_na-block'] = list(dfSchi.query("condition=='na-block' and code=='23'").quietFR)  
        expData['Schi15_IT5B_RS_quiet_na-block'] = list(dfSchi.query("condition=='na-block' and cell_class=='IT'").quietFR)  
        expData['Schi15_PT5B_RS_quiet_na-block']  = list(dfSchi.query("condition=='na-block' and cell_class=='PT'").quietFR) 
        expData['Schi15_L5Benh_RS_quiet_na-block'] = list(dfSchi.query("condition == 'na-block' and type=='enh'").quietFR)  
        expData['Schi15_L5Bsupp_RS_quiet_na-block'] = list(dfSchi.query("condition == 'na-block' and type=='supp'").quietFR)  
        expData['Schi15_L5Benh5%_RS_quiet_na-block'] = list(dfSchi.query("condition == 'na-block' and moveFR>1.05*quietFR").quietFR)  
        expData['Schi15_L5Bsupp5%_RS_quiet_na-block'] = list(dfSchi.query("condition == 'na-block' and moveFR<0.95*quietFR").quietFR)  

        expData['Schi15_L5B_RS_move_na-block'] = list(dfSchi.query("condition == 'na-block' and code!='23'").moveFR)  # na-block, move (high VL, high ih)
        expData['Schi15_IT2_RS_move_na-block'] = list(dfSchi.query("condition=='na-block' and code=='23'").moveFR)  
        expData['Schi15_IT5B_RS_move_na-block'] = list(dfSchi.query("condition=='na-block' and cell_class=='IT'").moveFR)  
        expData['Schi15_PT5B_RS_move_na-block']  = list(dfSchi.query("condition=='na-block' and cell_class=='PT'").moveFR)  
        expData['Schi15_L5Benh_RS_move_na-block'] = list(dfSchi.query("condition == 'na-block' and type=='enh'").moveFR)  
        expData['Schi15_L5Bsupp_RS_move_na-block'] = list(dfSchi.query("condition == 'na-block' and type=='supp'").moveFR)  
        expData['Schi15_L5Benh5%_RS_move_na-block'] = list(dfSchi.query("condition == 'na-block' and moveFR>1.05*quietFR").moveFR)  
        expData['Schi15_L5Bsupp5%_RS_move_na-block'] = list(dfSchi.query("condition == 'na-block' and moveFR<0.95*quietFR").moveFR)  



        # ------------------------------------------------------------
        # plot combined model and experimental data firing rate stats 
        ## statData[0] - IT2
        ## statData[2] - IT5A
        ## statData[3] - IT5B
        ## statData[4] - PT5B

        conds = ['Quiet', 'Move']
        sources = ['Model', 'Experiment']

        titles = ['L5B IT', 'L5B PT', 'L5B', 'L5B enh', 'L5B supp']
        modelInds = [0, 1, 2, 3, 4]
        expLabels = ['Schi15_IT5B_RS', 'Schi15_L5B_RS', 'Schi15_L5B_RS', 'Schi15_L5Benh5%_RS', 'Schi15_L5Bsupp5%_RS']
        
        colors = [popColors['IT5B'], popColors['PT5B'], '#800080', 'orange', 'purple']

        # what to plot
        boxplot = 0
        lineplot = 0
        lineplotMultipleSims = 1
        barplot_diff = 0
        densityplot = 0
        scatterplot = 0
        stattests = 0
        combined = 0

        processedModelData = {}
        saveProcessedModelData = 1 - int(loadAll) #True

        for title, modelInd, expLabel, color in zip(titles, modelInds, expLabels, colors):
            statDataRows = [['Quiet', 'Model', statDataQuiet[modelInd]],
                            ['Quiet', 'Experiment', expData[expLabel+'_quiet_na-block']],
                            ['Move', 'Model', statDataMove[modelInd]],
                            ['Move', 'Experiment', expData[expLabel+'_move_na-block']]]
            statDataCols = ['cond', 'source', 'rates']
            
            dfStats = pd.DataFrame(statDataRows, columns=statDataCols)

            if len(dfStats.iloc[1].rates) == 0 or len(dfStats.iloc[3].rates) == 0:
                continue

            dfStats = utils.explode(dfStats, ['rates']) 

            quietExp = dfStats.query('cond=="Quiet" and source=="Experiment"').rates
            quietModel = dfStats.query('cond=="Quiet" and source=="Model"').rates
            moveExp = dfStats.query('cond=="Move" and source=="Experiment"').rates
            moveModel = dfStats.query('cond=="Move" and source=="Model"').rates

            quietExp, moveExp = zip(*[(q,m) for q,m in zip(quietExp, moveExp) if q>=minValueQuiet and m>=minValueMove])
            quietExp=pd.Series(quietExp)
            moveExp=pd.Series(moveExp)

            processedModelData[title] = {'quietModel': quietModel, 'moveModel': moveModel, 'dfStats': dfStats}

            # ---------
            # boxplot       
            if boxplot:    
                overlayLine = True
                utils.stats_boxplot(dfStats=dfStats, x='cond', y='rates', hue='source', color=color, 
                    overlayLine=overlayLine, quietExp=quietExp, moveExp=moveExp, quietModel=quietModel, moveModel=moveModel, figSize=(8,8), fontsize=fontsiz)  
                filename = root+batchLabel+'_boxplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                plt.savefig(filename, dpi=300)            
            
            # ---------
            # line plot with std
            # Note: individual model line plots are taken from statDataAllQuiet and statDataAllMove, which contain all data points 
            # ie. filtering of >0.01 Hz has not been applied to these raw variables, so may not correspond to quietExp and moveExp
            if lineplot:

                if title == 'L5B PT':
                    linestyle=':' 
                else:
                    linestyle='-'

                utils.stats_lineplot(title, quietExp, moveExp, quietModel, moveModel, statDataAllQuiet, statDataAllMove, modelInd, 
                    figSize=(6,4), fontsize=fontsiz+8, printOutput=True, linestyle=linestyle)
                filename = root+batchLabel+'_lineplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                plt.savefig(filename, dpi=300)


            # ---------
            # line plot with std - multiple sims
            if lineplotMultipleSims:

                if title == 'L5B PT':
                    linestyle=':' 
                else:
                    linestyle='-'

                filename = root + 'v56_batch29_incK_processedModelData.pkl'
                filename = filename.replace(batchLabel, 'v56_batch29')

                with open (filename, 'rb') as f:
                    d=pickle.load(f)

                #import IPython as ipy; ipy.embed()

                dfStats1 = dfStats.replace('Model', 'Model')
                dfStats2 = d[title]['dfStats'].query('source=="Model"').replace('Model', 'Model (increased K+)')
                dfStatsCombined = pd.concat([dfStats1.query('source=="Model"'), dfStats2, dfStats1.query('source=="Experiment"')])
                
                quietModel2 = {'Model': quietModel, 'Model (increased K+)': d[title]['quietModel']}
                moveModel2 = {'Model': moveModel, 'Model (increased K+)': d[title]['moveModel']}

                my_pal = {'Model': 'royalblue', 'Model (increased K+)': 'purple', 'Experiment': 'orange'}

                utils.stats_lineplot(title, quietExp, moveExp, quietModel2, moveModel2, statDataAllQuiet, statDataAllMove, modelInd, figSize=(6,4),
                    fontsize=fontsiz+8, printOutput=True,  linestyle=linestyle, colors=my_pal)

                batchLabelCondition = 'incK'
                filename = root + batchLabel + '_' + batchLabelCondition + '_lineplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])

                fontsiz = 40
                ax=plt.gca()
                plt.ylim(-1, 25)
                plt.yticks([0, 10, 20], ['0', '10', '20'], fontsize=fontsiz)
                plt.ylabel('', fontsize=fontsiz)
                plt.tick_params(bottom = False)
                ax.axes.xaxis.set_ticklabels([])
                ax.spines['bottom'].set_visible(False)

                plt.savefig(filename, dpi=300)


            # ---------
            # bar plot with change in firing rate
            if barplot_diff:
                utils.stats_barplot_diff(quietExp, moveExp, quietModel, moveModel,figSize=(4,8), fontsize=fontsiz)
                filename = root+batchLabel+'_barplot_relative_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                plt.savefig(filename, dpi=300)
                        

            # ---------
            # 2d density plot (with gaussian KDE)
            # https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density/
            if densityplot:
                # model + exp
                utils.stats_densityplot_combined(quietExp, moveExp, quietModel, moveModel, figSize=(8,8), fontsize=fontsiz)
                filename = root + batchLabel + '_2Ddensity_%s_%d_%d' % (title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                plt.savefig(filename, dpi=300)

                # only exp
                # utils.stats_densityplot_single(quietExp, moveExp, figSize=(8,8), fontsize=fontsiz)
                # filename = root + batchLabel + '_2Ddensity_exp_%s_%d_%d' % (title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                # plt.savefig(filename, dpi=300)


            # ---------
            # scatter plot 
            if scatterplot:
                utils.stats_scatterplot(quietExp, moveExp, quietModel, moveModel, figSize = (8,8), fontsize=fontsiz)
                filename = root + batchLabel + '_scatterplot_%s_%d_%d' % (title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                plt.savefig(filename, dpi=300)


            # -------------
            #  statistical tests 
            if stattests:
                utils.stats_tests(title, quietExp, moveExp, quietModel, moveModel)


        # -------------
        #  save processed model data to file
        if saveProcessedModelData:
            filename = root + batchLabel + '_' + batchLabelCondition + '_processedModelData.pkl'
            with open (filename, 'wb') as f:
                pickle.dump(processedModelData, f)



        # ------------------------------------------------------------
        # plot combined model and experimental L5B enh vs supp distribution (for main + na_block) 
        # does not match exp
        if combined:
            # exp data
            L5BMainExpTotal = len(expData['Schi15_L5Benh5%_RS_quiet']) + len(expData['Schi15_L5Bsupp5%_RS_quiet'])
            L5BenhMainExpFrac = len(expData['Schi15_L5Benh5%_RS_quiet']) / L5BMainExpTotal * 100
            L5BsuppMainExpFrac = len(expData['Schi15_L5Bsupp5%_RS_quiet']) / L5BMainExpTotal * 100

            L5BExpNAblockTotal = len(expData['Schi15_L5Benh5%_RS_quiet_na-block']) + len(expData['Schi15_L5Bsupp5%_RS_quiet_na-block'])
            L5BenhNAblockExpFrac = len(expData['Schi15_L5Benh5%_RS_quiet_na-block']) / L5BExpNAblockTotal * 100
            L5BsuppNAblockExpFrac = len(expData['Schi15_L5Bsupp5%_RS_quiet_na-block']) / L5BExpNAblockTotal * 100

            L5BenhExpFrac = [L5BenhMainExpFrac, L5BenhNAblockExpFrac]
            L5BsuppExpFrac = [L5BsuppMainExpFrac, L5BsuppNAblockExpFrac]

            # model data
            L5BMainModelTotal = len(statDataQuietMain[4]) + len(statDataQuietMain[5])
            L5BenhMainModelFrac = len(statDataQuietMain[4]) / L5BMainModelTotal * 100
            L5BsuppMainModelFrac = len(statDataQuietMain[5]) / L5BMainModelTotal * 100

            L5BNAblockModelTotal = len(statDataQuiet[3]) + len(statDataQuiet[4])
            L5BenhNAblockModelFrac = len(statDataQuiet[3]) / L5BNAblockModelTotal * 100
            L5BsuppNAblockModelFrac = len(statDataQuiet[4]) / L5BNAblockModelTotal * 100

            L5BenhModelFrac = [L5BenhMainModelFrac, L5BenhNAblockModelFrac]
            L5BsuppModelFrac = [L5BsuppMainModelFrac, L5BsuppNAblockModelFrac]

            # cumulative bar plot exp
            ind = range(2)  
            barWidth = 0.85
            plt.figure(figsize=(8,8))
            L5BenhBar = plt.bar(ind, L5BenhExpFrac, width=barWidth, color='#f9bc86', edgecolor='white')
            L5BsuppBar = plt.bar(ind, L5BsuppExpFrac, width=barWidth, bottom=L5BenhExpFrac , color='#a3acff', edgecolor='white')

            plt.xticks(ind, ['Ctrl', 'NA-R block'], fontsize=fontsiz)
            plt.ylabel("Relative distribution (%)", fontsize=fontsiz)
            plt.legend((L5BenhBar[0], L5BsuppBar[0]), ('L5Benh', 'L5Bsupp'), fontsize=fontsiz)

            filename = root + batchLabel + '_barplotDistribExp'
            plt.savefig(filename, dpi=300)

            # cumulative bar plot exp
            ind = range(2)  
            barWidth = 0.85
            L5BenhBar = plt.bar(ind, L5BenhModelFrac, width=barWidth, color='#f9bc86', edgecolor='white')
            L5BsuppBar = plt.bar(ind, L5BsuppModelFrac, width=barWidth, bottom=L5BenhModelFrac , color='#a3acff', edgecolor='white')

            plt.xticks(ind, ['Ctrl', 'NA-R block'], fontsize=fontsiz)
            plt.ylabel("Relative distribution (%)", fontsize=fontsiz)
            plt.legend((L5BenhBar[0], L5BsuppBar[0]), ('L5Benh', 'L5Bsupp'), fontsize=fontsiz)
            #plt.title('Model')

            filename = root + batchLabel + '_barplotDistribModel'
            plt.savefig(filename, dpi=300)


            # model indiv trial data
            for i in range(25):
                L5BMainModelTotal = len(statDataAllQuietMain[i][4]) + len(statDataAllQuietMain[i][5])
                L5BenhMainModelFrac = len(statDataAllQuietMain[i][4]) / L5BMainModelTotal * 100
                L5BsuppMainModelFrac = len(statDataAllQuietMain[i][5]) / L5BMainModelTotal * 100

                L5BNAblockModelTotal = len(statDataAllQuiet[i][3]) + len(statDataAllQuiet[i][4])
                L5BenhNAblockModelFrac = len(statDataAllQuiet[i][3]) / L5BNAblockModelTotal * 100
                L5BsuppNAblockModelFrac = len(statDataAllQuiet[i][4]) / L5BNAblockModelTotal * 100

                L5BenhModelFrac = [L5BenhMainModelFrac, L5BenhNAblockModelFrac]
                L5BsuppModelFrac = [L5BsuppMainModelFrac, L5BsuppNAblockModelFrac]


                # cumulative bar plot exp
                ind = range(2)  
                barWidth = 0.85
                L5BenhBar = plt.bar(ind, L5BenhModelFrac, width=barWidth, color='#f9bc86', edgecolor='white')
                L5BsuppBar = plt.bar(ind, L5BsuppModelFrac, width=barWidth, bottom=L5BenhModelFrac , color='#a3acff', edgecolor='white')

                plt.xticks(ind, ['Ctrl', 'NA-R block'], fontsize=fontsiz)
                plt.ylabel("Relative distribution (%)", fontsize=fontsiz)
                plt.legend((L5BenhBar[0], L5BsuppBar[0]), ('L5Benh', 'L5Bsupp'), fontsize=fontsiz)
                #plt.title('Model')

                filename = root + batchLabel + '_barplotDistribModel_'+str(i)
                plt.savefig(filename, dpi=300)


            # ------------------------------------------------------------
            # SBR - L5Benh move/quiet ratio in 
            # model shows same trend; need experimental data from Fig 6

            SBR_main_exp = np.mean(expData['Schi15_L5Benh5%_RS_move']) / np.mean(expData['Schi15_L5Benh5%_RS_quiet'])
            SBR_nablock_exp = np.mean(expData['Schi15_L5Benh5%_RS_move_na-block']) / np.mean(expData['Schi15_L5Benh5%_RS_quiet_na-block'])

            SBR_main_model = np.mean(statDataMoveMain[4]) / np.mean(statDataQuietMain[4])
            SBR_nablock_model = np.mean(statDataMove[3]) / np.mean(statDataQuiet[3])


            SBR_main_model_all = [x / y for x,y in zip(statDataMoveMain[4], statDataQuietMain[4])]
            SBR_nablock_model_all = [x / y for x,y in zip(statDataMove[3], statDataQuiet[3])]

            conditions = ['Control'] * len(SBR_main_model_all) +  ['NA-block'] * len(SBR_nablock_model_all)

            dfStats = pd.DataFrame(np.array([SBR_main_model_all + SBR_nablock_model_all, conditions]).T, columns = ['rate', 'condition'])
            dfStats[['rate']] = dfStats[['rate']].apply(pd.to_numeric)  

            my_pal = {'Control': 'gray', 'NA-block': 'green'}
            fig=plt.figure(figsize=(8,8))
            ax = sb.boxplot(x='condition', y='rate', data=dfStats, palette=my_pal, showfliers=False)  #
            
            #plt.title (title, fontsize=fontsiz)
            handles, labels = ax.get_legend_handles_labels()
            #ax.legend(handles=handles, labels=labels, fontsize=fontsiz)
            plt.tight_layout()
            plt.ylabel('Signal-to-baseline ratio', fontsize=fontsiz)
            plt.xlabel('')
            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsiz) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsiz) 
            fig.subplots_adjust(left=0.13)

            filename = root + batchLabel + '_SBR_boxplot'
            plt.savefig(filename, dpi=300)



# Main code
if __name__ == '__main__': 
    fig_narblock()
