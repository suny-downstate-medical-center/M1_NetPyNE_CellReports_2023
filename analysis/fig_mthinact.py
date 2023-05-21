"""
M1 paper Figure 5

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis

# ----------------------------------------------------------------
def fig_mthinact():
    ''' Figure MTh (VL) inactivation activity: 
        - raster plot of 1-7 sec 
        - traces 1-7 sec (compare to exp?)
        - stats (boxplot,line,scatter) comparing quite+move to exp
    '''

    # ---------------------------------------------------------------------------------------------------------------
    # Config

    raster = 0
    histogram = 0
    traces = 0
    stats = 1      # boxplot of rates

    fontsiz = 26

    dataFolder = '../../data/'

    # -------------------------------------------------------------------------------------------
    # Raster plot
    if  raster:  # 2 sec N=1  
        batchLabel = 'v56_batch20'  
        simLabel = 'v56_batch20_0_1_0_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        timeRange = [4000, 7000] #[2000, 4000]
        include = allpops
        orderBy = ['pop', 'y']
        
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

        batchLabel = 'v56_batch20'  
        simLabel = 'v56_batch20_0_1_0_0_0'

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
    elif traces:
        batchLabel = 'v56_batch20'  
        simLabel = 'v56_batch20_0_1_0_0_0'

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
        plt.savefig('%s%s_traces_%d_%d_PT-%d_x1.5.png' % (root, simLabel, timeRange[0], timeRange[1], include[0]), dpi=200)


    # ---------------------------------------------------------------------------------------------------------------
    # stats plots (boxplot, line plot, scatter)

    elif stats:  # 50 sec N=1
        
        timeRangeQuiet = [1000, 5000]
        timeRangeMove = [5000, 9000]
        include = ['IT5B', 'PT5B', ['IT5B','PT5B']]
        xlim = [0, 40]  #[0,70] # ok to cut off max and flyers (Bill19)
        labelsModel = ['IT5B', 'PT5B', 'L5B', 'L5Benh', 'L5Bsupp']
        multipleTrials = True
        loadAll = 1

        # multiple conditions from same batch
        batchLabel = 'v56_batch25'
        simLabelPrefixes = ['v56_batch25_0_0_0_0_0_0_0_0',
                            'v56_batch25_1_1_1_1_1_1_1_0',
                            'v56_batch25_2_2_2_2_2_2_2_0',
                            'v56_batch25_3_3_3_3_3_3_3_0',
                           'v56_batch25_4_4_4_4_4_4_4_0',
                           'v56_batch25_5_5_5_5_5_5_5_0',
                            'v56_batch25_6_6_6_6_6_6_6_0',
                            'v56_batch25_7_7_7_7_7_7_7_0']

        batchLabelConditions = ['VL_M2_inact',
                                'VL_PO_inact',
                                'VL_M2_cM1_inact',
                                'VL_PO_M2_cM1_inact',
                                'VL_cM1_M2_OC_inact',
                               'VL_cM1_M2_OC_inact',
                               'VL_M2_OC_inact', 
                                'VL_S1h_S2h_inact',
                                'VL_cM1_M2_OC_S1h_S2h_inact'
                                 ]

        # Original
        batchLabel = 'v56_batch20'
        simLabelPrefixes = ['v56_batch20_0_1_0']
        batchLabelConditions = ['VL_inact_original']

        for simLabelPrefix, batchLabelCondition in zip(simLabelPrefixes, batchLabelConditions):
                
            # single trial
            if not multipleTrials:
                batchLabel = 'v56_batch20'
                simLabel = 'v56_batch20_0_1_0_0_0'
                sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
            
                # fix
                fig1, modelData = sim.analysis.plotSpikeStats(include=include, figSize=(8,4), timeRange=timeRangeMove, xlim=xlim,
                stats = ['rate'], includeRate0=True, legendLabels=None, fontSize = fontsiz, popColors=popColors, showFig=0, dpi=300, saveFig=False)
            
                statDataMove = modelData['statData'] 
            
            # multiple trials
            else:

                simLabel = [simLabelPrefix + '_' + str(iseed) + '_' + str(iconn) for iseed in range(5) for iconn in range(5)]

                root = dataFolder + batchLabel + '/'
                
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
                                #print(gid, end=', ')
                                _, spktcell, _ = analysis.utils.getSpktSpkid([gid], [timeRangeQuiet[0], timeRangeMove[1]], sim)
                                # slower if use pandas
                                # rateQuiet = int(dfcell.query('@timeRangeQuiet[0] < spkt < @timeRangeQuiet[1]').spkt.count()) * 1000 / (timeRangeQuiet[1] - timeRangeQuiet[0])
                                # rateMove = int(dfcell.query('@timeRangeMove[0] < spkt < @timeRangeMove[1]').spkt.count()) * 1000 / (timeRangeMove[1] - timeRangeMove[0])
                                                
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
                    filename = root+simLabel[0]+'.json'

                statDataAllMove = statData['statDataAllMove']
                statDataAllQuiet = statData['statDataAllQuiet']
                includeAll = statData['includeAll']

                # # code to combine IT5B+PT5B into L5B
                # for i in range(1):
                #     statDataAllMove[i][2] = statDataAllMove[i][0]+statDataAllMove[i][1]
                #     statDataAllQuiet[i][2] = statDataAllQuiet[i][0]+statDataAllQuiet[i][1]

                # combine data
                minValueQuiet = 0.0
                minValueMove = 0.0
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

            expData['Schi15_L5B_RS_quiet_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and code!='23'").quietFR)  # mth-inact, quiet (low VL, medium ih)
            expData['Schi15_IT2_RS_quiet_mth-inact'] = list(dfSchi.query("condition=='mth-inact' and code=='23'").quietFR)  
            expData['Schi15_IT5B_RS_quiet_mth-inact'] = list(dfSchi.query("condition=='mth-inact' and cell_class=='IT'").quietFR)  
            expData['Schi15_PT5B_RS_quiet_mth-inact']  = list(dfSchi.query("condition=='mth-inact' and cell_class=='PT'").quietFR) 
            expData['Schi15_L5Benh_RS_quiet_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='enh'").quietFR)  
            expData['Schi15_L5Bsupp_RS_quiet_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='supp'").quietFR)  
            expData['Schi15_L5Benh5%_RS_quiet_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and moveFR>1.05*quietFR").quietFR)  
            expData['Schi15_L5Bsupp5%_RS_quiet_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and moveFR<0.95*quietFR").quietFR)  

            expData['Schi15_L5B_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and code!='23'").moveFR)  # mth-inact, move (low VL, low ih)
            expData['Schi15_IT2_RS_move_mth-inact'] = list(dfSchi.query("condition=='mth-inact' and code=='23'").moveFR)  
            expData['Schi15_IT5B_RS_move_mth-inact'] = list(dfSchi.query("condition=='mth-inact' and cell_class=='IT'").moveFR)  
            expData['Schi15_PT5B_RS_move_mth-inact']  = list(dfSchi.query("condition=='mth-inact' and cell_class=='PT'").moveFR)  
            expData['Schi15_L5Benh_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='enh'").moveFR)  
            expData['Schi15_L5Bsupp_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='supp'").moveFR)  
            expData['Schi15_L5Benh_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='enh'").moveFR)  
            expData['Schi15_L5Bsupp_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and type=='supp'").moveFR)  
            expData['Schi15_L5Benh5%_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and moveFR>1.05*quietFR").moveFR)  
            expData['Schi15_L5Bsupp5%_RS_move_mth-inact'] = list(dfSchi.query("condition == 'mth-inact' and moveFR<0.95*quietFR").moveFR)  


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
            #expLabels = ['Schi15_IT5B_RS', 'Schi15_PT5B_RS', 'Schi15_L5B_RS', 'Schi15_L5Benh5%_RS', 'Schi15_L5Bsupp5%_RS']
            expLabels = ['Schi15_L5B_RS', 'Schi15_L5B_RS', 'Schi15_L5B_RS', 'Schi15_L5Benh5%_RS', 'Schi15_L5Bsupp5%_RS']

            colors = [popColors['IT5B'], popColors['PT5B'], '#800080', 'orange', 'purple']

            # what to plot
            boxplot = 0
            boxplotMultipleSims = 0
            lineplot = 0
            lineplotMultipleSims = 1
            barplot = 0
            barplot_diff = 0
            densityplot = 0
            scatterplot = 0
            stattests = 0

            processedModelData = {}
            saveProcessedModelData = 1 - int(loadAll) #True

            for title, modelInd, expLabel, color in zip(titles, modelInds, expLabels, colors):
                statDataRows = [['Quiet', 'Model', statDataQuiet[modelInd]],
                                ['Quiet', 'Experiment', expData[expLabel+'_quiet_mth-inact']],
                                ['Move', 'Model', statDataMove[modelInd]],
                                ['Move', 'Experiment', expData[expLabel+'_move_mth-inact']]]
                statDataCols = ['cond', 'source', 'rates']

                dfStats = pd.DataFrame(statDataRows, columns=statDataCols)
                
                if len(dfStats.iloc[1].rates) == 0 or len(dfStats.iloc[3].rates) == 0:
                    continue
                
                dfStats = utils.explode(dfStats, ['rates'])

                #dfStats = dfStats.convert_objects(convert_numeric=True).dropna()  # remove non-numeric 'rates' values
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
                    filename = root + batchLabel + '_' + batchLabelCondition +'_boxplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                    plt.savefig(filename, dpi=300)     


                # ---------
                # boxplot multiple sims    
                if boxplotMultipleSims:    

                    popinacts = ['VL_M2_inact',
                        'VL_PO_inact',
                        'VL_M2_cM1_inact',
                        'VL_PO_M2_cM1_inact',
                        'VL_cM1_M2_OC_inact',
                        'VL_cM1_M2_OC_inact',
                        'VL_M2_OC_inact', 
                        'VL_S1h_S2h_inact',
                        'VL_cM1_M2_OC_S1h_S2h_inact']

                    #popinacts = ['VL_M2_cM1_inact']

                    for popinact in popinacts:

                        filename = root + batchLabel + '_' +popinact + '_processedModelData.pkl'
                        filename = filename.replace(batchLabel, 'v56_batch25')

                        with open (filename, 'rb') as f:
                            d=pickle.load(f)

                        dfStats1 = dfStats.replace('Model', 'Model (inact VL)')
                        dfStats2 = d[title]['dfStats'].query('source=="Model"').replace('Model', popinact)
                        dfStatsCombined = pd.concat([dfStats1.query('source=="Model (inact VL)"'), dfStats2, dfStats1.query('source=="Experiment"')])
                        
                        quietModel2 = {'Model (inact VL)': quietModel, popinact: d[title]['quietModel']}
                        moveModel2 = {'Model (inact VL)': moveModel, popinact: d[title]['moveModel']}

                        my_pal = {'Model (inact VL)': 'royalblue', popinact: 'purple', 'Experiment': 'orange'}

                        overlayLine = True
                        utils.stats_boxplot(dfStats=dfStatsCombined, x='cond', y='rates', hue='source', color=color, 
                            overlayLine=overlayLine, quietExp=quietExp, moveExp=moveExp, quietModel=quietModel2, moveModel=moveModel2, 
                            figSize=(8,8), fontsize=fontsiz, palette=my_pal)  
                        
                        filename = root + batchLabel + '_' + batchLabelCondition + '_' + popinact + '_boxplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                        
                        plt.savefig(filename, dpi=300)     

                # ---------
                # line plot with std
                # Note: individual model line plots are taken from statDataAllQuiet and statDataAllMove, which contain all data points 
                # ie. filtering of >0.01 Hz has not been applied to these raw variables, so may not correspond to quietExp and moveExp
                if lineplot:

                    print(batchLabelCondition)
                    
                    if title == 'L5B PT' or title == 'L5B IT':
                        linestyle=':' 
                    else:
                        linestyle='-'

                    utils.stats_lineplot(title, quietExp, moveExp, quietModel, moveModel, statDataAllQuiet, statDataAllMove, modelInd, figSize=(6,4),
                        fontsize=fontsiz+8, printOutput=True,  linestyle=linestyle)
                    filename = root + batchLabel + '_' + batchLabelCondition +'_lineplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1],)
                    
                    plt.savefig(filename, dpi=300)

                # ---------
                # line plot with std - multiple sims
                if lineplotMultipleSims:

                    if title == 'L5B PT' or title == 'L5B IT':
                        linestyle=':' 
                    else:
                        linestyle='-'

                    # popinacts = ['VL_M2_inact',
                    #     'VL_PO_inact',
                    #     'VL_M2_cM1_inact',
                    #     'VL_PO_M2_cM1_inact',
                    #     'VL_cM1_M2_OC_inact',
                    #     'VL_cM1_M2_OC_inact',
                    #     'VL_M2_OC_inact', 
                    #     'VL_S1h_S2h_inact',
                    #     'VL_cM1_M2_OC_S1h_S2h_inact']

                    popinacts = ['VL_M2_cM1_inact']

                    for popinact in popinacts:

                        filename = root + batchLabel + '_' +popinact + '_processedModelData.pkl'
                        filename = filename.replace(batchLabel, 'v56_batch25')

                        with open (filename, 'rb') as f:
                            d=pickle.load(f)

                        dfStats1 = dfStats.replace('Model', 'Model (inact VL)')
                        dfStats2 = d[title]['dfStats'].query('source=="Model"').replace('Model', popinact)
                        dfStatsCombined = pd.concat([dfStats1.query('source=="Model (inact VL)"'), dfStats2, dfStats1.query('source=="Experiment"')])
                        
                        quietModel2 = {'Model (inact VL)': quietModel, popinact: d[title]['quietModel']}
                        moveModel2 = {'Model (inact VL)': moveModel, popinact: d[title]['moveModel']}

                        my_pal = {'Model (inact VL)': 'royalblue', popinact: 'purple', 'Experiment': 'orange'}

                        utils.stats_lineplot(title, quietExp, moveExp, quietModel2, moveModel2, statDataAllQuiet, statDataAllMove, modelInd, figSize=(6,4),
                            fontsize=fontsiz+8, printOutput=True,  linestyle=linestyle, colors=my_pal)

                        filename = root + batchLabel + '_' + batchLabelCondition + '_' + popinact + '_lineplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])

                        # to generate separate legend uncomment these lines:
                        # utils.stats_lineplot(title, quietExp, moveExp, quietModel2, moveModel2, [], [], modelInd, figSize=(6,4),
                        #     fontsize=fontsiz+8, printOutput=True,  linestyle=linestyle, colors=my_pal)
                        # plt.legend(['Experiment', 'Model (original)', 'Model (modified)'])                        
                        # filename = root + batchLabel + '_' + batchLabelCondition + '_' + popinact + '_boxplot_legendTemp_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])                        

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
                if barplot:
                    utils.stats_barplot(title, quietExp, moveExp, quietModel, moveModel, statDataAllQuiet, statDataAllMove, modelInd, figSize=(8,8), fontsize=fontsiz, printOutput=False)
                    filename = root+batchLabel+'_barplot_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                    plt.savefig(filename, dpi=300)


                # ---------
                # bar plot with change in firing rate
                if barplot_diff:
                    utils.stats_barplot_diff(quietExp, moveExp, quietModel, moveModel,figSize=(4,8), fontsize=fontsiz)
                    filename = root+batchLabel+'_barplot_diff_relative_%s_%d_%d'%(title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                    plt.savefig(filename, dpi=300)
                        

                # ---------
                # 2d density plot (with gaussian KDE)
                # https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density/
                if densityplot:
                    # model + exp
                    utils.stats_densityplot_combined(quietExp, moveExp, quietModel, moveModel, figSize=(8,8), fontsize=fontsiz)
                    filename = root + batchLabel + '_' + batchLabelCondition + '_2Ddensity_%s_%d_%d' % (title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
                    plt.savefig(filename, dpi=300)


                # ---------
                # scatter plot 
                if scatterplot:
                    utils.stats_scatterplot(quietExp, moveExp, quietModel, moveModel, figSize = (8,8), fontsize=fontsiz)
                    filename = root + batchLabel + '_' + batchLabelCondition + '_scatterplot_%s_%d_%d' % (title.replace(' ', '').replace('/', ''), timeRangeQuiet[0], timeRangeMove[1])
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

        

# Main code
if __name__ == '__main__': 
    fig_mthinact()
