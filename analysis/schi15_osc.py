"""
paper.py 

Paper figures

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis
from scipy.io import loadmat
import pickle as pkl


# ----------------------------------------------------------------
def fig_osc_lfp(inputDataFolder, outputDataFolder, dataFile, lfp=0, lfp_psd=0, lfp_spectrogram=0, durMin=2, durMax=None, 
            fixedDurSegmentsStep=False, dataSavePSD=None, counter={}, saveToFile=1, loadFromFile=0):
    
    ''' Figure LFP filt and spectrogrma: 
        - LFP signal + filtered 
        - LFP Spectrogram
        '''

    fontsiz = 16
    dataSave = []

    from netpyne import sim, specs

    if lfp or lfp_psd or lfp_spectrogram:

        if not loadFromFile:
            try:
                data = loadmat(inputDataFolder + dataFile)  # use simplify_cells = True to get proper cells structure
            except:
                print('Error loading file %s' % inputDataFolder + dataFile)
                return -1
        
        print('\n %s' % (inputDataFolder+dataFile))

        onlimits = [float(x) for x in data['all_var']['movtrig'][0][0][0]['onlimitbigbin'][0]]
        offlimits = [float(x) for x in data['all_var']['movtrig'][0][0][0]['offlimitbigbin'][0]]
        expLFP = data['all_var']['basic'][0][0]['trace'][0][0]
        dt = float(data['all_var']['basic'][0][0]['dt'][0][0])
        
        #import IPython; IPython.embed()

        plt.style.use('seaborn-ticks')

        if onlimits[0] < offlimits[0]:
            timeRanges =  [[0, onlimits[0]]] \
                        + [[on, off] for on, off in zip(onlimits, offlimits)] \
                        + [[off, on] for on, off in zip(onlimits[1:], offlimits)] \
                        + [[offlimits[-1], len(expLFP)]]
                        
        elif offlimits[0] < onlimits[0]:
            timeRanges =  [[0, offlimits[0]]] \
                        + [[off, on] for on, off in zip(onlimits, offlimits)] \
                        + [[on, off] for on, off in zip(onlimits, offlimits[1:])] \
                        + [[onlimits[-1], len(expLFP)]]


        # remove short ones (< threshold)
        print('Keeping only segments with duration > %.1f secs...' % (durMin))
        durThreshold = durMin * 1./dt      # 2 secs ()                      
        timeRanges = [t for t in timeRanges if (t[1]-t[0]) > durThreshold ]

        durThresholdMax = durMax * 1./dt if durMax is not None else None
        if durThresholdMax:
            print('Keeping only segments with duration < %.1f secs...' % (durMax))
            timeRanges = [t for t in timeRanges if  (t[1]-t[0]) < durThresholdMax ]

        durs = [(t[1]-t[0]) for t in timeRanges]

        # add 'move' or 'quiet' labels
        timeLabels = ['move' if t[0] in onlimits else 'quiet' for t in timeRanges]  # list of boolean indicating if move period

        
                
        # get all segments of fixed duration x secs  
        if fixedDurSegmentsStep:
            print('Adding all segments with fixed duration %.1f secs at interval %.1f...' % (durMin, fixedDurSegmentsStep))
            step = fixedDurSegmentsStep * 1./dt
            timeRangesNew = []
            timeLabelsNew = []
            
            for i,timeRange in enumerate(timeRanges):
                for startTime in range(int(timeRange[0]), int(timeRange[1])-int(durThreshold), int(step)):
                    timeRangesNew.append([startTime, startTime+durThreshold])
                    timeLabelsNew.append(timeLabels[i])
            
            timeRanges = timeRangesNew
            timeLabels = timeLabelsNew

        # convert to ms
        timeRanges = [[t[0] * 1000 * dt, t[1] * 1000 * dt] for t in timeRanges]

        print(timeRanges)
        print(timeLabels)

        dataSavePSDEachFile = {'quiet': [], 'move': [], 'freq': []}

        counter = {'quiet': len(dataSavePSD['quiet']) - 1, 'move': len(dataSavePSD['move']) - 1}
        
        print(len(dataSavePSD['quiet']))

        for timeRange, timeLabel in zip(timeRanges, timeLabels):
            
            counter[timeLabel] += 1
            print(counter)

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

                electrodes = [0] 
                filtFreq = 200
                plots = ['timeSeries']
                filename = '%s%s_lfp_timeSeries_%s-%d_%d_%d_filt_%d.png' % (outputDataFolder, dataFile, timeLabel, counter[timeLabel], timeRange[0], timeRange[1], filtFreq)
                sim.cfg = specs.SimConfig()

                sim.cfg.recordStep = 0.1
                fig4, figdata = sim.analysis.plotLFP(inputLFP=expLFP, plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=(8,4), overlay=True, 
                    separation=1.5, dpi=300, colors=[[0,0,0]]*1, filtFreq=filtFreq, saveFig=0, showFig=False) 
                lfpPlot = figdata['LFP']
                t = figdata['t']
                color = 'black'
                lw = 0.5

                plt.figure(figsize=(8,4))
                fontsiz = 16
                
                if len(t) < len(lfpPlot):
                    lfpPlot = lfpPlot[:len(t)]

                plt.plot(t[0:len(lfpPlot)], -lfpPlot, color=color, linewidth=lw)

                ax = plt.gca()        
                ax.invert_yaxis()
                plt.axis('off')
                plt.xlabel('time (ms)', fontsize=fontsiz)
                plt.ylim(0.4,-0.5)
                plt.subplots_adjust(bottom=0.0, top=0.9, left=0.1, right=0.9)

                # calculate scalebar size and add scalebar
                scaley = 1000.0  # values in mV but want to convert to uV
                sizey = 100/scaley
                labely = '%.3g $\mu$V'%(sizey*scaley)#)[1:]
                add_scalebar(ax, hidey=True, matchy=False, hidex=True, matchx=True, sizex=None, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, 
                        unitsx='ms', loc='upper right', pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
                
                plt.title('LFP 0-200 Hz', fontsize=fontsiz, fontweight='bold')
                plt.savefig(filename,dpi=300)

            # ---------------------------------------------------------------------------------------------------------------
            # LFP Spectrogram
            if lfp_spectrogram:
                # spectrogram
                electrodes = [0]
                plots = ['spectrogram'] #['PSD','timeSeries', 'spectrogram']
                filtFreq = 200

                # log scale
                fig4, outdata = sim.analysis.plotLFP(inputLFP=expLFP, plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=(8,4), overlay=True, maxFreq=80,
                    normSpec=True,  dpi=300, saveFig=0, showFig=False, logy=1) 

                ax = plt.gca()
                plt.title('')
                plt.suptitle('')
                #plt.suptitle('LFP spectrogram', fontsize=fontsiz,  fontweight='bold')
            
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
                plt.savefig('%s%s_lfp_spec_log_%s-%d_%d_%d_filt_%d.png' % (outputDataFolder, dataFile, timeLabel, counter[timeLabel], timeRange[0], timeRange[1], filtFreq), dpi=300)


                # linear scale
                fig4, outdata = sim.analysis.plotLFP(inputLFP=expLFP, plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=(8,4), overlay=True, maxFreq=80,
                    normSpec=True,  dpi=300, saveFig=0, showFig=False, logy=0) 

                ax = plt.gca()
                plt.title('')
                plt.suptitle('')
                #plt.suptitle('LFP spectrogram', fontsize=fontsiz,  fontweight='bold')
            
                # f = outdata['freqs']
                # ystep = 20
                # ax.set_yticks(range(1, len(f), ystep))
                # ax.set_yticklabels(['%d' % (x) for x in f[::ystep]])

                ax.set_xticklabels(range(int(timeRange[0]/1000.), int(timeRange[1]/1000.)))
                ax.set_xticks(range(int(timeRange[0]), int(timeRange[1]), 1000))
                plt.setp(ax.get_xticklabels(),fontsize=fontsiz-2)
                plt.setp(ax.get_yticklabels(),fontsize=fontsiz-2)
                plt.ylabel('Frequency (Hz)',fontsize=fontsiz)
                plt.xlabel('Time (s)',fontsize=fontsiz)
                plt.subplots_adjust(bottom=0.25, top=0.9, right=1.0, left=0.1)
                plt.savefig('%s%s_lfp_spec_lin_%s-%d_%d_%d_filt_%d.png' % (outputDataFolder, dataFile, timeLabel, counter[timeLabel], timeRange[0], timeRange[1], filtFreq), dpi=300)


            # ---------------------------------------------------------------------------------------------------------------
            # LFP PSD
            if lfp_psd:

                # psd
                figSiz = (8,4)
                fontsiz=12
                electrodes = [0] # [1,4,6,8] #,10]
                legendLabels = ['Experiment LFP'] # ['Electrode 300um (L2/3)', 'Electrode 600um (L5A)', 'Electrode 800um (upper L5B)', 'Electrode 1000um (lower L5B)'] #, 'Electrode 1100um (L6)']
                plots = ['PSD'] #['PSD','timeSeries', 'spectrogram']
                colors = [[0,0,0]]*11
                
                colors[0] = [x/255.0 for x in [253,116,0]]    # popColors['IT2']
                colors[1] = [x/255.0 for x in [255,225,26]]   # popColors['IT5A']
                colors[0] = [x/255.0 for x in [190, 219, 57]] #'mediumpurple'#popColors['PT5B']
                colors[3] = [x/255.0 for x in [31, 138, 112]] #'purple'

                colors[0] = 'red' #popColors['IT2']
                colors[1] = 'magenta' #popColors['IT5A']
                colors[2] = 'blue' #'mediumpurple'#popColors['PT5B']
                colors[3] = 'green' #'purple'
                #colors[4] = 'firebrick' #popColors['IT6']

                filtFreq = 200
                fig4, outData = sim.analysis.plotLFP(inputLFP=expLFP, plots=plots, electrodes=electrodes, timeRange=timeRange, figSize=figSiz, overlay=True,  maxFreq=80, stepFreq=1, normSignal=0, normPSD=0, separation=1.5, dpi=200, lineWidth=2.5, colors=colors, saveFig=0, showFig=False) 
                
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
                plt.savefig('%s%s_lfp_psd_morlet_notnormx2_%s-%d_%d_%d_%d.png' % (outputDataFolder, dataFile, timeLabel, counter[timeLabel], timeRange[0], timeRange[1],filtFreq), dpi=300)

                signal = {}
                freq = {}
                for ilabel,label in enumerate(legendLabels):
                    signal[label] = outData['allSignal'][ilabel]
                    freq[label] = outData['allFreqs'][ilabel]

                if dataSavePSD:
                    psdInfo = {'psd': signal['Experiment LFP'], 'file': dataFile, 'timeRange': timeRange, 'counter': counter[timeLabel], 'timeLabel': timeLabel}
                    
                    dataSavePSDEachFile[timeLabel].append(signal['Experiment LFP'])
                    dataSavePSDEachFile['freq'] = freq['Experiment LFP']

                    dataSavePSD[timeLabel].append(signal['Experiment LFP'])
                    dataSavePSD['info'].append(psdInfo)
                    dataSavePSD['freq'] = freq['Experiment LFP']

                    #print(dataSavePSD)
                    if saveToFile:
                        with open('%s/%s_psd_lfp_norm_stats_data.pkl' % (outputDataFolder, 'combined'), 'wb') as f:
                            pickle.dump(dataSavePSD, f)
                    


        # ---------------------------------------------------------------------------------------------------------------
        # plot combined PSD for this file
        if lfp_psd:

            CVQuiet = np.max(scipy.stats.variation(dataSavePSDEachFile['quiet']))
            CVMove = np.max(scipy.stats.variation(dataSavePSDEachFile['move']))

            print('Max CV quiet = %.2f, move = %.2f' % (CVQuiet, CVMove))                

            try:
                plot_combined_psd(inputDataFolder, outputDataFolder, dataSavePSDEachFile, filenameModifier='_'+dataFile)
            except:
                pass



def plot_combined_psd(inputDataFolder, outputDataFolder, dataSavePSD, filenameModifier=''):
    # ---------------------------------------------------------------------------------------------------------------
    # LFP PSD

    quietPSD = np.mean(dataSavePSD['quiet'], 0)
    movePSD = np.mean(dataSavePSD['move'], 0)
    quietPSDstd = np.std(dataSavePSD['quiet'], 0)
    movePSDstd = np.std(dataSavePSD['move'], 0)

    freqs = dataSavePSD['freq']
    legendLabels = ['quiet', 'move']
    fontsiz = 14

    plt.figure(figsize=(8,4))

    plt.plot(freqs, quietPSD) #, color='blue')
    plt.plot(freqs, movePSD)
    plt.fill_between(freqs, quietPSD-quietPSDstd, quietPSD+quietPSDstd, alpha=0.4, edgecolor=None, facecolor='#1f77b4')
    plt.fill_between(freqs, movePSD-movePSDstd, movePSD+movePSDstd, alpha=0.4, edgecolor=None, facecolor='#ff7f0e')
                                                                                                                                  
    plt.legend(legendLabels, fontsize=fontsiz)
    ax = plt.gca()
    plt.title('')
    plt.suptitle('')
    #ax.legend_.remove()
    [plt.setp(l,linewidth=2.5) for l in ax.lines] 
    #plt.xscale('log')
    plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
    plt.ylabel('LFP power (mV^2/Hz)', fontsize=fontsiz)
    ax.set_xticklabels(range(0,90,10))
    ax.set_xticks(range(0, 90, 10))
    #plt.ylim(-0.003, 0.005)
    plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
    plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
    
    plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
    plt.savefig('%s/combined_lfp_psd_morlet_notnorm%s.png' % (outputDataFolder, filenameModifier), dpi=300)


def plot_combined_psd_cluster(inputDataFolder, outputDataFolder, dataSavePSD, filenameModifier='', n_clusters=4):
    # ---------------------------------------------------------------------------------------------------------------
    # LFP PSD

    from sklearn import cluster

    quietData = np.array(dataSavePSD['quiet'])
    quietDataNorm = np.array([x/np.sum(x) for x in quietData])
    moveData = np.array(dataSavePSD['move'])
    moveDataNorm = np.array([x/np.sum(x) for x in moveData])
    
    # 1–4 Hz, Theta: 4–7 Hz, Alpha: 8–12 Hz, Beta: 13–30 Hz, Gamma: 30–80 Hz
    freqRanges = [[0,4], [4,8], [8,13], [13,30], [30,80]]
    freqLabels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']

    
    quietPeaks = []
    movePeaks = []
    
    for freq in freqRanges:
        quietPeaks.append([np.sum(x[freq[0]:freq[1]]) / np.sum(x) for x in quietData])
        movePeaks.append([np.sum(x[freq[0]:freq[1]] / np.sum(x)) for x in moveData])

    quietPeaksT = np.array(quietPeaks).T
    pdQuiet = pd.DataFrame(quietPeaksT)
    pdQuiet.columns = freqLabels

    movePeaksT = np.array(movePeaks).T
    pdMove = pd.DataFrame(np.array(movePeaks).T)
    pdMove.columns = freqLabels


    # k-means clustering
    kmeansQuiet = cluster.KMeans(n_clusters=n_clusters).fit(quietPeaksT)
    pdQuiet['kmeans'] = kmeansQuiet.labels_
    kmeansMove = cluster.KMeans(n_clusters=n_clusters).fit(movePeaksT)
    pdMove['kmeans'] = kmeansMove.labels_

    # UMAP
    import umap

    n_neighbors = 30

    reducerQuiet = umap.UMAP(random_state=1, n_neighbors=n_neighbors, min_dist=0.0, n_components=2)  # 15, 0.1
    umapQuiet = reducerQuiet.fit_transform(quietPeaksT)

    reducerMove = umap.UMAP(random_state=1, n_neighbors=n_neighbors, min_dist=0.0, n_components=2)
    umapMove = reducerMove.fit_transform(movePeaksT)

    yQuiet = list(pdQuiet.Gamma)
    yMove = list(pdMove.Gamma)

    plt.figure()
    plt.scatter(umapQuiet[:, 0], umapQuiet[:, 1], s= 5, c=yQuiet, cmap='Spectral')
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.title('Quiet UMAP')
    plt.savefig('%s/UMAP_%d_quiet.png' % (outputDataFolder, n_neighbors))

    plt.figure()
    for i, c in zip(range(n_clusters), 'grbcmykw'):
        plt.scatter(umapQuiet[kmeansQuiet.labels_==i, 0], umapQuiet[kmeansQuiet.labels_==i, 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')
    plt.savefig('%s/UMAP_vs_kmeans_%d_quiet.png' % (outputDataFolder, n_neighbors))


    plt.figure()
    plt.scatter(umapMove[:, 0], umapMove[:, 1], s=5, c=yMove, cmap='Spectral')
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.title('Move UMAP')
    plt.savefig('%s/UMAP_%d_move.png' % (outputDataFolder, n_neighbors))

    plt.figure()
    for i, c in zip(range(n_clusters), 'grbcmykw'):
        plt.scatter(umapMove[kmeansMove.labels_==i, 0], umapMove[kmeansMove.labels_==i, 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')
    plt.savefig('%s/UMAP_vs_kmeans_%d_move.png' % (outputDataFolder, n_neighbors))

    # apply k-means to UMAP
    n_clusters = 2
    kmeansUMAPQuiet = cluster.KMeans(n_clusters=n_clusters).fit(umapQuiet)
    kmeansUMAPMove = cluster.KMeans(n_clusters=n_clusters).fit(umapMove)

    plt.figure()
    for i, c in zip(range(n_clusters), 'grbcmykw'):
        plt.scatter(umapQuiet[kmeansUMAPQuiet.labels_==i, 0], umapQuiet[kmeansUMAPQuiet.labels_==i, 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')
    plt.savefig('%s/kmeans_UMAP_%d_quiet.png' % (outputDataFolder, n_neighbors))

    plt.figure()
    for i, c in zip(range(n_clusters), 'grbcmykw'):
        plt.scatter(umapMove[kmeansUMAPMove.labels_==i, 0], umapMove[kmeansUMAPMove.labels_==i, 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')
    plt.savefig('%s/kmeans_UMAP_%d_move.png' % (outputDataFolder, n_neighbors))


    # plot quiet and move on same space
    allPeaksT=np.concatenate([quietPeaksT,movePeaksT])
    reducerAll = umap.UMAP(random_state=1, n_neighbors=n_neighbors, min_dist=0.0, n_components=2)  # 15, 0.1
    umapAll = reducerQuiet.fit_transform(allPeaksT)
    
    plt.figure()
    for i, c in zip(range(n_clusters), 'grbcmykw'):
        plt.scatter(umapAll[np.array(list(kmeansQuiet.labels_==i) + [False]*len(movePeaksT)), 0], 
                    umapAll[np.array(list(kmeansQuiet.labels_==i) + [False]*len(movePeaksT)), 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')
    
    for i, c in zip(range(n_clusters), 'bmykw'):
        plt.scatter(umapAll[np.array([False]*len(quietPeaksT) + list(kmeansMove.labels_==i)), 0], 
                    umapAll[np.array([False]*len(quietPeaksT) + list(kmeansMove.labels_==i)), 1], c=c, s=5)
        plt.gca().set_aspect('equal', 'datalim')

    plt.savefig('%s/ALL_UMAP_vs_kmeans_%d_quiet.png' % (outputDataFolder, n_neighbors))

    # calculate distance between clusters 
    meanQuietCluster1 = np.mean(umapAll[np.array(list(kmeansQuiet.labels_==0) + [False]*len(movePeaksT))])
    meanQuietCluster2 = np.mean(umapAll[np.array(list(kmeansQuiet.labels_==1) + [False]*len(movePeaksT))])
    meanMoveCluster1 = np.mean(umapAll[np.array([False]*len(quietPeaksT) + list(kmeansMove.labels_==0))])
    meanMoveCluster2 = np.mean(umapAll[np.array([False]*len(quietPeaksT) + list(kmeansMove.labels_==1))])

    distQuiet1Quiet2 = np.abs(meanQuietCluster1-meanQuietCluster2)
    distMove1Move2 = np.abs(meanMoveCluster1-meanMoveCluster2)
    distQuiet1Move1 = np.abs(meanQuietCluster1-meanMoveCluster1)
    distQuiet1Move2 = np.abs(meanQuietCluster1-meanMoveCluster2)
    
    distQuiet2Move1 = np.abs(meanQuietCluster2-meanMoveCluster1)
    distQuiet2Move2 = np.abs(meanQuietCluster2-meanMoveCluster2)

    print("QS-QL: %.1f, MS-ML: %.1f, QS-ML: %.1f, QS-MS: %.1f, QL-ML: %.1f, QL-MS: %.1f" 
        % (distQuiet1Quiet2, distMove1Move2, distQuiet1Move1, distQuiet1Move2, distQuiet2Move1, distQuiet2Move2))


    # calculate CV before clustering
    CVQuiet = scipy.stats.variation(pdQuiet)
    CVMove = scipy.stats.variation(pdMove)

    print('\n Unclustered dataset: Mean CV quiet = %.2f, move = %.2f' % (np.mean(CVQuiet), np.mean(CVMove)))
    print('\n Unclustered dataset: Max CV quiet = %.2f, move = %.2f' % (np.max(CVQuiet), np.max(CVMove)))


    freqs = dataSavePSD['freq']
    legendLabels = ['quiet', 'move']
    fontsiz = 18

    # PCA
    from sklearn.decomposition import PCA
    from mpl_toolkits.mplot3d import Axes3D

    pcaQuiet = PCA(n_components=2)
    pcaQuiet.fit(quietPeaksT)
    X_pcaQuiet = pcaQuiet.transform(quietPeaksT)

    pcaMove = PCA(n_components=2)
    pcaMove.fit(movePeaksT)
    X_pcaMove = pcaMove.transform(movePeaksT)

    print(pcaQuiet.explained_variance_ratio_)
    print(pcaMove.explained_variance_ratio_)

    fig = plt.figure(figsize=(6, 5))
    #ax = Axes3D(fig, rect=[0, 0, .95, 1])#, elev=48, azim=134)
    for i, c in zip(range(n_clusters), 'rgbcmykw'):
        plt.scatter(X_pcaQuiet[kmeansQuiet.labels_==i, 0], X_pcaQuiet[kmeansQuiet.labels_==i, 1], c=c)
    plt.legend()
    plt.savefig('%s/PCA2_quiet.png' % (outputDataFolder))

    fig = plt.figure(figsize=(6, 5))
    #ax = Axes3D(fig, rect=[0, 0, .95, 1])#, elev=48, azim=134)
    for i, c in zip(range(n_clusters), 'rgbcmykw'):
        plt.scatter(X_pcaMove[kmeansMove.labels_==i, 0], X_pcaMove[kmeansMove.labels_==i, 1], c=c)
    plt.legend()
    plt.savefig('%s/PCA2_move.png' % (outputDataFolder))


    # add info of each sample (file, time, etc)
    dinfo=pd.DataFrame(dataSavePSD['info'])
    dinfoQuiet=dinfo.query("timeLabel=='quiet'")
    dinfoMove=dinfo.query("timeLabel=='move'")

    # add missing row in quiet (544) -- due to memory error when processing all files, and continued on wrong index
    # row = pd.DataFrame([{'psd': quietData[544],
    #                     'file': '25_LFP.mat',
    #                     'timeRange': [55920.5, 59920.5],
    #                     'counter': 544,
    #                     'timeLabel': 'quiet'}],
    #                      index=[544])
    # dinfoQuiet = pd.concat([dinfoQuiet.iloc[:544], row, dinfoQuiet.iloc[544:]]).reset_index(drop=True)

    # combine with psd power per freq band
    dinfoQuiet=pd.concat([dinfoQuiet, pdQuiet], axis=1)
    dinfoMove=pd.concat([dinfoMove.reset_index(), pdMove], axis=1)


    plotEachCluster = 1
    if plotEachCluster:

        for k in range(n_clusters):
            # ------------------
            # quiet cluster PSD

            plt.figure(figsize=(8,4))

            quietPSD = np.mean(quietDataNorm[kmeansQuiet.labels_==k], 0)
            quietPSDstd = np.std(quietDataNorm[kmeansQuiet.labels_==k], 0)
            plt.plot(freqs, quietPSD, color='#1f77b4') #, color='blue')
            plt.fill_between(freqs, quietPSD-quietPSDstd, quietPSD+quietPSDstd, alpha=0.4, edgecolor=None, facecolor='#1f77b4')

            plt.legend(legendLabels, fontsize=fontsiz)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            #ax.legend_.remove()
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            #plt.xscale('log')
            plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(range(0,90,10))
            ax.set_xticks(range(0, 90, 10))
            #plt.ylim(-0.003, 0.005)
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_quiet_cluster_%d_%d%s.png' % (outputDataFolder, k, n_clusters, filenameModifier), dpi=300)

            # ------------------
            # quiet cluster boxplot
            plt.figure(figsize=(8*2,4*2))

            pdQuietCluster = pdQuiet.query('kmeans == @k').drop(columns='kmeans')
            sb.boxplot(data=pdQuietCluster)
            sb.stripplot(data=pdQuietCluster, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(freqLabels)
            ax.set_xticks(range(0, 5))
            #plt.ylim(-150, 1100)
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_scatter_quiet_cluster%d_%d%s.png' % (outputDataFolder, k, n_clusters, filenameModifier), dpi=300)

            print('\n Quiet cluster %d of %d includes %d / %d (%.1f%%) data points' % 
                (k+1, n_clusters, len(pdQuietCluster), len(pdQuiet), len(pdQuietCluster)/len(pdQuiet) *100 ))

            # ------------------
            # move clusters PSD
            plt.figure(figsize=(8,4))
            movePSD = np.mean(moveDataNorm[kmeansMove.labels_==k], 0)
            movePSDstd = np.std(moveDataNorm[kmeansMove.labels_==k], 0)
            plt.plot(freqs, movePSD, color ='#ff7f0e')
            plt.fill_between(freqs, movePSD-movePSDstd, movePSD+movePSDstd, alpha=0.4, edgecolor=None, facecolor='#ff7f0e')
                                                                                                                                    
            plt.legend(legendLabels, fontsize=fontsiz)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            #ax.legend_.remove()
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            #plt.xscale('log')
            plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(range(0,90,10))
            ax.set_xticks(range(0, 90, 10))
            #plt.ylim(-0.003, 0.005)
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_move_cluster_%d_%d%s.png' % (outputDataFolder, k, n_clusters, filenameModifier), dpi=300)


            # ------------------
            # move cluster boxplot
            plt.figure(figsize=(8*2,4*2))

            pdMoveCluster = pdMove.query('kmeans == @k').drop(columns='kmeans')
            sb.boxplot(data=pdMoveCluster)
            sb.stripplot(data=pdMoveCluster, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)
            ax = plt.gca()
            plt.title('')
            plt.suptitle('')
            [plt.setp(l,linewidth=2.5) for l in ax.lines] 
            plt.xlabel('Frequency band', fontsize=fontsiz) 
            plt.ylabel('Normalized LFP power', fontsize=fontsiz)
            ax.set_xticklabels(freqLabels)
            ax.set_xticks(range(0, 5))
            #plt.ylim(0,0.6)
            plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
            plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
            
            #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
            plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
            plt.savefig('%s/combined_lfp_psd_morlet_notnorm_scatter_move_cluster%d_%d%s.png' % (outputDataFolder, k, n_clusters, filenameModifier), dpi=300)

            print('\n Move cluster %d of %d includes %d / %d (%.1f%%) data points' % 
                (k+1, n_clusters, len(pdMoveCluster), len(pdMove), len(pdMoveCluster)/len(pdMove) * 100))

            # -----------------------------
            # CVs
            CVQuiet = scipy.stats.variation(pdQuietCluster)
            CVMove = scipy.stats.variation(pdMoveCluster)

            print('\n Mean CV quiet = %.2f, move = %.2f' % (np.mean(CVQuiet), np.mean(CVMove)))
            print('\n Max CV quiet = %.2f, move = %.2f' % (np.max(CVQuiet), np.max(CVMove)))


    # ---------------
    # plot combined quiet+move scatter + PSD
    plotCombined = 1
    if plotCombined:
        
        # scatter
        quietK = np.argmax(np.bincount(kmeansQuiet.labels_))
        moveK = np.argmax(np.bincount(kmeansMove.labels_))

        plt.figure(figsize=(8*3,4*2))

        dinfoQuietCluster = dinfoQuiet[kmeansQuiet.labels_==quietK]
        dinfoMoveCluster = dinfoMove[kmeansMove.labels_==moveK]
        #dinfoCombined = pd.melt(pd.concat([dinfoQuietCluster, dinfoMoveCluster]), id_vars=['file', 'timeRange', 'counter', 'timeLabel'])

        pdQuietCluster = pdQuiet.query('kmeans == @quietK').drop(columns='kmeans')
        pdMoveCluster = pdMove.query('kmeans == @moveK').drop(columns='kmeans')

        pdQuietCluster['state'] = ['Quiet'] * len(pdQuietCluster)
        pdMoveCluster['state'] = ['Movement'] * len(pdMoveCluster)
        pdCombined = pd.melt(pd.concat([pdQuietCluster, pdMoveCluster]), id_vars=['state'])

        ax = sb.boxplot(data=pdCombined, x='variable', y='value', hue='state', width=0.8)
        sb.stripplot(data=pdCombined, x='variable', y='value', hue='state', dodge=True, linewidth=0.5, edgecolor='gray')#, size=4, color=".3", linewidth=0)

        handles, _ = ax.get_legend_handles_labels() 
        plt.legend(handles, ['Quiet', 'Movement'], loc='best', fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        #ax.legend_.remove()
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        #plt.xscale('log')
        plt.xlabel('Frequency band', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        #ax.set_xticklabels(freqLabels)
        #ax.set_xticks(range(0, 5))
        #plt.ylim(-150, 1100)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_notnorm_scatter_quietmove_clusters_%d_%d_of_%d%s.png' % (outputDataFolder, quietK, moveK, n_clusters, filenameModifier), dpi=300)

        filename = '%s/combined_lfp_psd_morlet_notnorm_scatter_quietmove_clusters%s.pkl' % (outputDataFolder, filenameModifier)
        with open(filename, 'wb') as f:
            pkl.dump(pdCombined, f)

        filename = '%s/combined_lfp_psd_dinfo_scatter_quietmove_clusters%s.pkl' % (outputDataFolder,filenameModifier)
        with open(filename, 'wb') as f:
            pkl.dump({'quiet': dinfoQuietCluster, 'move':  dinfoMoveCluster}, f)

        # PSD norm
        plt.figure(figsize=(8,4))

        quietPSD = np.mean(quietDataNorm[kmeansQuiet.labels_==quietK], 0)
        quietPSDstd = np.std(quietDataNorm[kmeansQuiet.labels_==quietK], 0)
        plt.plot(freqs, quietPSD, color='#1f77b4') #, color='blue')
        plt.fill_between(freqs, quietPSD-quietPSDstd, quietPSD+quietPSDstd, alpha=0.4, edgecolor=None, facecolor='#1f77b4')

        movePSD = np.mean(moveDataNorm[kmeansMove.labels_==moveK], 0)
        movePSDstd = np.std(moveDataNorm[kmeansMove.labels_==moveK], 0)
        plt.plot(freqs, movePSD, color='#ff7f0e') #, color='blue')
        plt.fill_between(freqs, movePSD-movePSDstd, movePSD+movePSDstd, alpha=0.4, edgecolor=None, facecolor='#ff7f0e')

        plt.legend(legendLabels, fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        #ax.legend_.remove()
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        #plt.xscale('log')
        plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
        plt.ylabel('Normalized LFP power', fontsize=fontsiz)
        ax.set_xticklabels(range(0,90,10))
        ax.set_xticks(range(0, 90, 10))
        #plt.ylim(-0.003, 0.005)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_norm_quietmove_clusters_%d_%d_of_%d%s.png' % (outputDataFolder, quietK, moveK, n_clusters, filenameModifier), dpi=300)


        # PSD absolute
        plt.figure(figsize=(8,4))

        quietPSD = np.mean(quietData[kmeansQuiet.labels_==quietK], 0)
        quietPSDstd = np.std(quietData[kmeansQuiet.labels_==quietK], 0)
        plt.plot(freqs, quietPSD, color='#1f77b4') #, color='blue')
        plt.fill_between(freqs, quietPSD-quietPSDstd, quietPSD+quietPSDstd, alpha=0.4, edgecolor=None, facecolor='#1f77b4')

        movePSD = np.mean(moveData[kmeansMove.labels_==moveK], 0)
        movePSDstd = np.std(moveData[kmeansMove.labels_==moveK], 0)
        plt.plot(freqs, movePSD, color='#ff7f0e') #, color='blue')
        plt.fill_between(freqs, movePSD-movePSDstd, movePSD+movePSDstd, alpha=0.4, edgecolor=None, facecolor='#ff7f0e')

        plt.legend(legendLabels, fontsize=fontsiz)
        ax = plt.gca()
        plt.title('')
        plt.suptitle('')
        #ax.legend_.remove()
        [plt.setp(l,linewidth=2.5) for l in ax.lines] 
        #plt.xscale('log')
        plt.xlabel('Frequency (Hz)', fontsize=fontsiz) 
        plt.ylabel('LFP power (mV^2/Hz)', fontsize=fontsiz)
        ax.set_xticklabels(range(0,90,10))
        ax.set_xticks(range(0, 90, 10))
        #plt.ylim(-0.003, 0.005)
        plt.setp(ax.get_xticklabels(),fontsize=fontsiz)
        plt.setp(ax.get_yticklabels(), fontsize=fontsiz)
        
        #fplt.suptitle('LFP power spectral density', fontsize=fontsiz, fontweight='bold')
        plt.subplots_adjust(bottom=0.15, top=0.95, right=0.95, left=0.15)
        plt.savefig('%s/combined_lfp_psd_morlet_notnorm_quietmove_clusters_%d_%d_of_%d%s.png' % (outputDataFolder, quietK, moveK, n_clusters, filenameModifier), dpi=300)






# ----------------------------------------------------------------------------------------------------
# Main code
if __name__ == '__main__':
    
    inputDataFolder = '../../data/Schi15/LFP_raw/'
    outputDataFolder = '../../data/Schi15/LFP_analysis/'
    try:
        os.mkdir(outputDataFolder)
    except:
        pass

    dataFileList = ['%d_LFP.mat' % (i) for i in range(10,41)] # range(26,41)] #range(35,41)] # range(10,41)] 
    
    loadFromFile = True
    saveToFile = 0 #1 - loadFromFile

    if loadFromFile:
        with open('%s/%s_psd_lfp_norm_stats_data.pkl' % (outputDataFolder, 'combined'), 'rb') as f:
            dataSavePSD = pickle.load(f)
    else:
        dataSavePSD = {'quiet': [], 'move': [], 'info': [], 'freq': None}
        

        for dataFile in dataFileList:
            fig_osc_lfp(inputDataFolder, outputDataFolder, dataFile, 
                        lfp=1, 
                        lfp_psd=1, 
                        lfp_spectrogram=1, 
                        lfp_cfc=0, 
                        durMin=4, 
                        durMax=None, 
                        fixedDurSegmentsStep=1,
                        dataSavePSD=dataSavePSD, 
                        saveToFile=saveToFile, 
                        loadFromFile=loadFromFile)  # lfp 
    
    plot_combined_psd_cluster(inputDataFolder, outputDataFolder, dataSavePSD, n_clusters=2)
