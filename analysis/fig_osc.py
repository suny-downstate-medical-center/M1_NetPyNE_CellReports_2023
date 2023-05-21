"""
M1 paper Figure 4 

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis



# ----------------------------------------------------------------
def fig_osc_lfp(dataFolder, batchLabel, simLabels, timeRange, lfp=0, lfp_psd=0, lfp_spectrogram=0, saveToFile=1, loadFromFile=0, dataSavePSD=None, timeLabel=None):
    ''' Figure LFP filt and spectrogrma: 
        - LFP signal + filtered slow+gamma
        - LFP Spectrogram
        '''
    
    fontsiz = 16
    dataSave = []

    if lfp or lfp_psd or lfp_spectrogram:

        root = dataFolder + batchLabel + '/'

        for simLabel in simLabels:

            if not loadFromFile:
                print(dataFolder, batchLabel, simLabel)
                sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

            plt.style.use('seaborn-ticks') 

            # ---------------------------------------------------------------------------------------------------------------
            # LFP signal
            if lfp:
                
                ''' OLD electrode locations = range(200,1300,100)
                0=200 (L23), 1=300 (L23), 2=400 (L23/4 border), 3=500 (L4/5A border), 4=600 (L5A), 5=700 (L5Bu), 
                6=800 (L5Bu), 7=900 (L5Bm), 8=1000 (L5Bl), 9=1100 (L6), 10=1200 (L6), 11=1300 (L6)
                '''

                ''' NEW electrode locations = [600, 800, 1000]
                0=600 (L5A), 1=800 (L5Bu), 2=1000 (L5Bl)
                '''
                
                # time series full - 4 sec
                electrodes = list(range(11))  
                filtFreq = 200
                plots = ['timeSeries']
                filename= '%s%s_lfp_timeSeries_%d_%d_filt_%d.png' % (root, simLabel, timeRange[0], timeRange[1], filtFreq)
                fig4, figdata = sim.analysis.plotLFP(plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=(8,4), overlay=True, 
                    separation=1.5, dpi=300, colors=[[0,0,0]]*11, filtFreq=filtFreq, saveFig=0, showFig=False) 
                lfpPlot = figdata['LFP']
                t = figdata['t']
                color = 'black'
                lw = 0.5

                for elec in electrodes:
                    plt.figure(figsize=(8,4))
                    fontsiz = 16
                    plt.plot(t, -lfpPlot[:,elec], color=color, linewidth=lw)
                    ax = plt.gca()        
                    ax.invert_yaxis()
                    plt.axis('off')
                    plt.xlabel('time (ms)', fontsize=fontsiz)

                    meanSignal = np.mean(-lfpPlot[:,elec])
                    plt.ylim(meanSignal+0.4,meanSignal-0.5)
                    plt.subplots_adjust(bottom=0.0, top=0.9, left=0.1, right=0.9)

                    # calculate scalebar size and add scalebar
                    scaley = 1000.0  # values in mV but want to convert to uV
                    sizey = 100/scaley
                    labely = '%.3g $\mu$V'%(sizey*scaley)#)[1:]
                    add_scalebar(ax, hidey=True, matchy=False, hidex=True, matchx=True, sizex=None, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, 
                            unitsx='ms', loc='upper right', pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
                    
                    plt.title('LFP 0-200 Hz', fontsize=fontsiz, fontweight='bold')
                    plt.savefig(filename[:-4]+'elec_%d.png' % (elec),dpi=300)



            # ---------------------------------------------------------------------------------------------------------------
            # LFP PSD
            if lfp_psd:

                # psd
                figSiz = (8,4)
                fontsiz=12
                electrodes = [1,4,6,8] #,10]
                legendLabels = ['Electrode 300um (L2/3)', 'Electrode 600um (L5A)', 'Electrode 800um (upper L5B)', 'Electrode 1000um (lower L5B)'] #, 'Electrode 1100um (L6)']
                plots = ['PSD'] #['PSD','timeSeries', 'spectrogram']
                colors = [[0,0,0]]*11
                
                colors[0] = [x/255.0 for x in [253,116,0]]    # popColors['IT2']
                colors[1] = [x/255.0 for x in [255,225,26]]   # popColors['IT5A']
                colors[2] = [x/255.0 for x in [190, 219, 57]] #'mediumpurple'#popColors['PT5B']
                colors[3] = [x/255.0 for x in [31, 138, 112]] #'purple'

                colors[0] = 'red' #popColors['IT2']
                colors[1] = 'magenta' #popColors['IT5A']
                colors[2] = 'blue' #'mediumpurple'#popColors['PT5B']
                colors[3] = 'green' #'purple'
                #colors[4] = 'firebrick' #popColors['IT6']

                filtFreq = 200
                fig4, outData = sim.analysis.plotLFP(plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=figSiz, overlay=True,  maxFreq=80, stepFreq=1, normSignal=0, normPSD=0, separation=1.5, dpi=200, lineWidth=2.5, colors=colors, saveFig=0, showFig=False) 
                plt.legend(legendLabels, fontsize=fontsiz-3)
                ax = plt.gca()
                plt.title('')
                plt.suptitle('')
                #ax.legend_.remove()
                [plt.setp(l,linewidth=2.5) for l in ax.lines] 
                #plt.xscale('log')
                plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
                plt.ylabel('Normalized LFP power',fontsize=fontsiz)
                plt.setp(ax.get_xticklabels(),fontsize=fontsiz-2)
                plt.setp(ax.get_yticklabels(),fontsize=fontsiz-2)
                
                #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
                plt.subplots_adjust(bottom=0.17, top=0.95, right=0.9, left=0.1)
                plt.savefig('%s%s_lfp_psd_morlet_notnormx2_%d_%d_%d.png' % (root, simLabel, timeRange[0], timeRange[1],filtFreq), dpi=300)


                signal = {}
                freq = {}
                for ilabel,label in enumerate(legendLabels):
                    signal[label] = outData['allSignal'][ilabel]
                    freq[label] = outData['allFreqs'][ilabel]
                
                if dataSavePSD:
                    dataSavePSD[timeLabel].append(signal)
                    dataSavePSD['freq'] = freq

                    if saveToFile:
                        with open('%s/combined_psd_lfp_norm_stats_data.pkl' % (root), 'wb') as f:
                            pickle.dump(dataSavePSD, f)



            # ---------------------------------------------------------------------------------------------------------------
            # LFP Spectrogram
            if lfp_spectrogram:
                electrodesList = [[4],[5],[6],[7],[8]]

                for electrodes in electrodesList:
                    plots = ['spectrogram'] #['PSD','timeSeries', 'spectrogram']
                    filtFreq = 200

                    fig4, outdata = sim.analysis.plotLFP(plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=(8,4), overlay=True, maxFreq=80,
                        normSpec=True,  dpi=300, saveFig=0, showFig=False, logy=1) 

                    ax = plt.gca()
                    plt.title('')
                    plt.suptitle('')
                    
                    f = outdata['freqs']

                    plt.ylim(np.min(f), np.max(f))
                    ystep = 20
                    yrange = np.arange(ystep, np.max(f)+1, ystep)
                    ax.set_yticks(yrange)
                    ax.set_yticklabels(['%.0f' % (f[y]) for y in [-3*ystep, -2*ystep, -ystep, -1]])

                    ax.set_xticklabels(range(int(timeRange[0]/1000.), int(timeRange[1]/1000.)))
                    ax.set_xticks(range(int(timeRange[0]), int(timeRange[1]), 1000))
                    plt.setp(ax.get_xticklabels(),fontsize=fontsiz-2)
                    plt.setp(ax.get_yticklabels(),fontsize=fontsiz-2)
                    plt.ylabel('Frequency (Hz)',fontsize=fontsiz)
                    plt.xlabel('Time (s)',fontsize=fontsiz)
                    plt.subplots_adjust(bottom=0.25, top=0.9, right=1.0, left=0.1)
                    plt.savefig('%s%s_lfp_spec_ylog_%d_%d_filt_%d_elec_%d.png' % (root, simLabel, timeRange[0], timeRange[1], filtFreq, electrodes[0]), dpi=300)



def plot_combined_psd_scatter(dataFolder, dataSavePSD, filenameModifier='', norm=True):
    # ---------------------------------------------------------------------------------------------------------------
    # LFP PSD scatter peaks by freq

    for k in dataSavePSD['quiet'][0].keys():

        nsims = 25*3 if k == 'L5' else 25

        quietData = [dataSavePSD['quiet'][i][k] for i in range(nsims)] 
        moveData = [dataSavePSD['move'][i][k] for i in range(nsims)] 
        
        # 1–4 Hz, Theta: 4–8 Hz, Alpha: 8–12 Hz, Beta: 13–30 Hz, Gamma: 30–80 Hz
        freqRanges = [[0,4], [4,8], [8,12], [12,30], [30,80]]

        freqLabels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']

        quietPeaks = []
        movePeaks = []
        changePeaks = []

        for freq in freqRanges:
            quietPeaks.append([np.sum(x[freq[0]:freq[1]])/np.sum(x) for x in quietData])
            movePeaks.append([np.sum(x[freq[0]:freq[1]]/np.sum(x)) for x in moveData])
            changePeaks.append([(m - q) / q * 100 for q, m in zip(quietPeaks[-1], movePeaks[-1])])


        quietPeaksT = np.array(quietPeaks).T
        pdQuiet = pd.DataFrame(quietPeaksT)
        pdQuiet.columns = freqLabels

        movePeaksT = np.array(movePeaks).T
        pdMove = pd.DataFrame(movePeaksT)
        pdMove.columns = freqLabels

        pdChange = pd.DataFrame(np.array(changePeaks).T)
        pdChange.columns = freqLabels

        freqs = dataSavePSD['freq']
        legendLabels = ['quiet', 'move']
        fontsiz = 18


        plotQuiet = 0
        plotMove = 0
        plotCombined = 1
        plotChange = 0


        # ---------------
        # plot pdQuiet
        if plotQuiet:
            plt.figure(figsize=(8*2,4*2))
            sb.boxplot(data=pdQuiet)
            sb.stripplot(data=pdQuiet, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

            #plt.legend(legendLabels, fontsize=fontsiz)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(freqLabels)
            ax.set_xticks(range(0, 5))
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_%s_scatter_quiet%s.png' % (dataFolder, k.replace(' ', '').replace('/',''), filenameModifier), dpi=300)

        # ---------------
        # plot pdMove
        if plotMove:
            plt.figure(figsize=(8*2,4*2))
            sb.boxplot(data=pdMove)
            sb.stripplot(data=pdMove, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(freqLabels)
            ax.set_xticks(range(0, 5))
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_%s_scatter_move%s.png' % (dataFolder, k.replace(' ', '').replace('/',''), filenameModifier), dpi=300)


        # ---------------
        # plot combined quiet+move
        if plotCombined:
            plt.figure(figsize=(8*3,4*2))

            pdQuiet['state'] = ['Quiet'] * len(pdQuiet)
            pdMove['state'] = ['Movement'] * len(pdMove)
            pdCombined = pd.melt(pd.concat([pdQuiet,pdMove]), id_vars=['state'])

            ax = sb.boxplot(data=pdCombined, x='variable', y='value', hue='state', width=0.8)
            sb.stripplot(data=pdCombined, x='variable', y='value', hue='state', dodge=True, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

            handles, _ = ax.get_legend_handles_labels() 
            plt.legend(handles, ['Quiet', 'Movement'], loc='best', fontsize=fontsiz)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_%s_scatter_quietmove%s.png' % (dataFolder, k.replace(' ', '').replace('/',''), filenameModifier), dpi=300)

            filename = '%s/combined_lfp_psd_morlet_notnorm_%s_scatter_quietmove%s.pkl' % (dataFolder, k.replace(' ', '').replace('/',''), filenameModifier)
            with open(filename, 'wb') as f:
                pickle.dump(pdCombined, f)                


        # ---------------
        # plot pdChange
        if plotChange:
            plt.figure(figsize=(8*2,4*2))
            sb.boxplot(data=pdChange)
            sb.stripplot(data=pdChange, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Change in LFP power (%)', fontsize=fontsiz)
            ax.set_xticklabels(freqLabels)
            ax.set_xticks(range(0, 5))
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_%s_scatter_change%s.png' % (dataFolder, k.replace(' ', '').replace('/',''), filenameModifier), dpi=300)



def plot_combined_psd_scatter_expVsModel(dataFolder):
    # ---------------------------------------------------------------------------------------------------------------
    # Plot combined scatter experiment vs model

        # load model data
        with open ('%s/combined_lfp_psd_morlet_notnorm_L5_scatter_quietmove.pkl' % (dataFolder), 'rb') as f:
            pdModel = pickle.load(f)

        # load exp data
        with open ('../../data/Schi15/LFP_analysis_4sec_fixedDur_step1_notCombined_freqsFixed/combined_lfp_psd_morlet_notnorm_scatter_quietmove_clusters.pkl', 'rb') as f:
            pdExp = pickle.load(f)

        # load exp data info
        with open('../../data/Schi15/LFP_analysis_4sec_fixedDur_step1_notCombined_freqsFixed/combined_lfp_psd_dinfo_scatter_quietmove_clusters.pkl', 'rb') as f:
            pdExpInfo = pickle.load(f)

        pdExpInfoQuiet = pdExpInfo['quiet'].drop(columns = ['counter', 'kmeans', 'psd'])
        pdExpInfoMove = pdExpInfo['move'].drop(columns = ['counter', 'kmeans', 'psd'])

        # add source column
        pdModel['source']=['Model']*len(pdModel)
        pdExp['source']=['Experiment']*len(pdExp)

        # generate combined quiet vs move
        pdQuiet = pd.concat([pdExp.query('state=="Quiet"'), pdModel.query('state=="Quiet"')])
        pdMove = pd.concat([pdExp.query('state=="Movement"'), pdModel.query('state=="Movement"')])


        # ---------------------------------------------
        # plot quiet
        plt.figure(figsize=(8*3,4*2))
        fontsiz = 18

        ax = sb.boxplot(data=pdQuiet, x='variable', y='value', hue='source', width=0.8)
        sb.stripplot(data=pdQuiet, x='variable', y='value', hue='source', dodge=True, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        plt.legend(handles, ['Experiment', 'Model'], loc='best', fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_scatter_quiet_modelVsExp.png' % (dataFolder ), dpi=300)


        # ---------------------------------------------
        # plot move
        plt.figure(figsize=(8*3,4*2))
        fontsiz = 18

        ax = sb.boxplot(data=pdMove, x='variable', y='value', hue='source', width=0.8)
        sb.stripplot(data=pdMove, x='variable', y='value', hue='source', dodge=True, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        plt.legend(handles, ['Experiment', 'Model'], loc='best', fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_scatter_move_modelVsExp.png' % (dataFolder ), dpi=300)


        # ---------------------------------------------
        # plot quiet+move combined
        plt.figure(figsize=(2*8*3,4*2))
        fontsiz = 18

        pdQuietMove = pd.concat([pdQuiet, pdMove])
        pdQuietMove['Source_State'] = pdQuietMove['source'] + ', ' + pdQuietMove['state']

        ax = sb.boxplot(data=pdQuietMove, x='variable', y='value', hue='Source_State', width=0.8)
        plt.legend(loc='best', fontsize=fontsiz)
        sb.stripplot(data=pdQuietMove, x='variable', y='value', hue='Source_State', dodge=True, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_scatter_quietmove_modelVsExp.png' % (dataFolder ), dpi=300)


        # ---------------------------------------------
        # plot exp+model quiet+move combined (4 bars)

        plt.figure(figsize=(8*3,4*2))
        fontsiz = 24

        pdQuietMove = pd.concat([pdExp, pdModel])
        pdQuietMove['Source_State'] = pdQuietMove['source'] + ', ' + pdQuietMove['state']

        colors = ['#ffd1a8', '#ff7f0e', '#81beea', '#1f77b4']

        ax = sb.boxplot(data=pdQuietMove, x='variable', y='value', hue='Source_State', width=0.8, palette=colors)
        plt.legend(loc='best', fontsize=fontsiz)
        sb.stripplot(data=pdQuietMove, x='variable', y='value', hue='Source_State', dodge=True, linewidth=0.5, edgecolor='gray', palette=colors)#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        plt.legend(handles, ['Experiment Quiet', 'Experiment Move', 'Model Quiet', 'Model Move'], loc='best', fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_scatter_modelVsExp_quietmove_narrow.png' % (dataFolder ), dpi=300)

        # print results
        for source in ['Experiment', 'Model']:
            for state in ['Quiet', 'Movement']:
                for band in ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']:
                    print('%s %s %s: %.2f +- %.2f ' % (source, state, band,
                    np.median(pdQuietMove.query('source==@source and state==@state and variable==@band')['value']),
                    scipy.stats.iqr(pdQuietMove.query('source==@source and state==@state and variable==@band')['value'])))


        # ---------------------------------------------
        # plot change in exp+model

        # change in model
        pdChangeModel = pd.DataFrame(pdModel.query('state=="Quiet"'))
        pdChangeModel['value'] = np.array(pdMove.query('source=="Model"')['value'])-np.array(pdQuiet.query('source=="Model"')['value'])
        pdChangeModel['state'] = ['Change'] * len(pdChangeModel)

        # change in experiment
        # note: for exp there is no one-to-one mapping of quiet vs move
        # so for each file, found quiet samples before move samples, and averaged each set 

        uniqueFiles = uniqueFiles=sorted(list(set(pdExpInfoQuiet['file'])))

        sampleStep = 1000.

        pdChangeExp = pd.DataFrame()

        for file in uniqueFiles:
            pdFileQuiet = pdExpInfoQuiet.query('file==@file')
            pdFileMove = pdExpInfoMove.query('file==@file')
            
            groupedQuiet = []
            groupedMove = []
            
            prevTimeQuiet = pdFileQuiet.iloc[0]['timeRange'][0]
            prevTimeMove = pdFileMove.iloc[0]['timeRange'][0]
            
            for irowQuiet, (_, rowQuiet) in enumerate(pdFileQuiet.iterrows()):
                print('\nReading quiet row: ', prevTimeQuiet, rowQuiet['timeRange'][0])

                if (rowQuiet['timeRange'][0] <= prevTimeQuiet + (sampleStep+1) and (irowQuiet < len(pdFileQuiet)-1)):
                    prevTimeQuiet = rowQuiet['timeRange'][0]
                    groupedQuiet.append(rowQuiet)
                    print('Added to quiet group')
                else:
                    if irowQuiet == len(pdFileQuiet)-1:
                        maxTime = pdFileMove.iloc[-1]['timeRange'][0]
                    else: 
                        maxTime = rowQuiet['timeRange'][0]
                        
                    print('\nEND OF QUIET GROUP!\n')

                    # check if move after
                    for irowMove, rowMove in pdFileMove.iterrows():

                        print('\nReading move row: ', groupedQuiet[-1]['timeRange'][1], prevTimeMove, rowMove['timeRange'][0], maxTime)
                        
                        if groupedQuiet[-1]['timeRange'][1] < rowMove['timeRange'][0] <= maxTime:   
                            groupedMove.append(rowMove)
                            print('Added to move group')
                        
                        prevTimeMove = rowMove['timeRange'][0]
                    
                    print('\nEND OF MOVE GROUP!\n')

                    if len(groupedMove) > 0:
                        # average value of grouped +  add exp data point
                        pdChangeExpRow = pd.DataFrame(groupedMove).mean() - pd.DataFrame(groupedQuiet).mean()
                        pdChangeExp = pdChangeExp.append(pdChangeExpRow, ignore_index=True)
                        # calculate average of model + add exp data point

                    # keep track of last rowQuiet -- reset grouped
                    prevTimeQuiet = rowQuiet['timeRange'][0]
                    groupedQuiet = [rowQuiet]
                    groupedMove = [] 


        freqLabels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']
        pdChangeExp=pdChangeExp.drop(columns=['index'])
        pdChangeExp = pdChangeExp[freqLabels]
        pdChangeExp['state'] = ['Change'] * len(pdChangeExp)
        pdChangeExp['source'] = ['Experiment'] * len(pdChangeExp)
        pdChangeExp = pd.melt(pdChangeExp, id_vars=['state', 'source'])
                                

        pdChange = pd.concat([pdChangeExp, pdChangeModel])

        # plot difference
        plt.figure(figsize=(8*3, 4*2))
        fontsiz = 24

        colors = [ '#ff7f0e', '#1f77b4']

        ax = sb.boxplot(data=pdChange, x='variable', y='value', hue='source', order=freqLabels, width=0.6, palette=colors)
        sb.stripplot(data=pdChange, x='variable', y='value', hue='source', order=freqLabels, dodge=True, linewidth=0.5, edgecolor='gray', palette=colors)#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        plt.legend(handles, ['Experiment', 'Model'], loc='best', fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Change in normalized LFP power', fontsize=fontsiz)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_scatter_change_modelVsExp_v2_narrow.png' % (dataFolder ), dpi=300)

        # print stats
        for source in ['Experiment', 'Model']:
            for band in ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']:
                print('%s %s: %.2f +- %.2f ' % (source, band,
                np.median(pdChange.query('source==@source and variable==@band')['value']),
                scipy.stats.iqr(pdChange.query('source==@source and variable==@band')['value'])))



def addCombinedL5(dataSavePSD):
    # add L5 combining L5A+L5Bu+L5Bl

    nsims = len(dataSavePSD['quiet'])
    L5elecs = ['Electrode 600um (L5A)', 'Electrode 800um (upper L5B)', 'Electrode 1000um (lower L5B)']

    method = 'concat'

    if method == 'mean':
        for i in range(nsims):
            dataSavePSD['quiet'][i]['L5'] = np.mean([dataSavePSD['quiet'][i][L5elec] for L5elec in L5elecs],0)
            dataSavePSD['move'][i]['L5'] = np.mean([dataSavePSD['move'][i][L5elec] for L5elec in L5elecs],0)

    elif method == 'concat':
        for ielec, elec in enumerate(L5elecs):
            for i in range(nsims):    
                index = (ielec*25)+i
                if index >= 25:
                    dataSavePSD['quiet'].append({})
                    dataSavePSD['move'].append({})
                
                dataSavePSD['quiet'][index]['L5'] = dataSavePSD['quiet'][i][elec]
                dataSavePSD['move'][index]['L5'] = dataSavePSD['move'][i][elec]

    dataSavePSD['freq']['L5'] = dataSavePSD['freq'][L5elecs[0]]

    return dataSavePSD


# Main code
if __name__ == '__main__':
    #     ''' Figure osc activity for 2 states (quiet, movement) x 3 conditions (main, mth-inactive, na-block): 
    #         - LFP signal, PSD, spectrogram (by elec depth)
    #         - Spike time histogram, PSD (by population)

    dataFolder = '../../data/'

    batchLabelList = ['v56_batch23'] * 2 #, 'v56_batch23', 'v56_batch23', 'v56_batch23']  
    simLabelsList = [['v56_batch23_0_0_%d_%d' % (i, j) for i in range(5) for j in range(5)]] * 2
    timeRangeList = [[1000, 5000], [5000, 9000]] #* 25   #[1000, 13000], [1000, 9000], ] 
    timeLabelList = ['quiet', 'move'] #* 25
    
    loadFromFile = 1 #1
    saveToFile = 1 - loadFromFile


    if loadFromFile:
        with open('%s/%s/%s_combined_psd_lfp_norm_stats_data.pkl' % (dataFolder, batchLabelList[0], batchLabelList[0]), 'rb') as f:
            dataSavePSD = pickle.load(f)

    else:
        dataSavePSD = {'quiet': [], 'move': [], 'freq': None}
        for batchLabel, simLabels, timeRange, timeLabel in zip(batchLabelList, simLabelsList, timeRangeList, timeLabelList):
            fig_osc_lfp(dataFolder, batchLabel, simLabels, timeRange, lfp=0, lfp_psd=0, lfp_spectrogram=1, saveToFile=saveToFile, 
                        loadFromFile=loadFromFile, dataSavePSD=dataSavePSD, timeLabel=timeLabel)  # lfp signal
   
    # combined analysis
    # dataSavePSD = addCombinedL5(dataSavePSD)
    # plot_combined_psd_scatter(dataFolder=dataFolder+'/v56_batch23', dataSavePSD=dataSavePSD, norm=True)
    plot_combined_psd_scatter_expVsModel(dataFolder=dataFolder+'/v56_batch23')

