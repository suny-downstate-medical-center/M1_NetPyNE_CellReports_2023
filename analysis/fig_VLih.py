"""
paper.py 

Paper figures

Contributors: salvadordura@gmail.com
"""

from shared import *
from netpyne import analysis



# ----------------------------------------------------------------
def fig_VLvsIh_test(vmin=None, vmax=None, method = 'mean_>=0.01', load=False, interpMissing=False):

    dataFolder = '../../data/'


    '''
    1) quiet state, normal condition -> medium VL, medium ih; 2) movement state, normal condition -> high VL, low ih, low ih; 
    3) quiet state, motor thalamus inactive -> low VL, medium ih; 4) movement state, motor thalamus inactive -> low VL, low ih;
    5) quiet state, NA blocked -> medium VL, high ih; 6) movement state, NA blocked -> high VL, high ih.

    ih = [0.0, 0.25, 0.5, 0.75, 1.0]
    VL = [0.01, 2.5, 5, 10, 15] # /2 so equivalent to ratesLong uniform distrib (seems more similar to max)
    
    main, quiet: ih = 0.75, VL = 2.5  -> 3_1
    main, move : ih = 0.25, VL = 10   -> 1_3
    mth, quiet : ih = 0.75, VL = 0.01 -> 3_0
    mth, move  : ih = 0.25, VL = 0.01 -> 1_0
    na, quiet  : ih = 1.0, VL = 2.5   -> 4_1
    na, move   : ih = 1.0, VL = 10    -> 4_3
    '''

    # v56_batch2
    # batchLabel = 'v56_batch2'  #'v53_batch12' 
    
    # ihValues = [0, 0.25, 0.5, 1.0, 1.25]#, 1.5]   
    # VLValues = [0.01, 2.5, 5, 10, 15]
    # timeRange = [1000, 5000]

    # v56_batch30
    batchLabel = 'v56_batch30'  #'v53_batch12' 
    timeRange = [1000, 5000]  # 1-5  
    
    ihValues = [0, 0.25, 0.5, 0.75, 1.0]#, 1.25, 1.5]  # [0.25, 0.75, 1.0] # 
    VLValues = [0.01, 2.5, 5, 7.5, 10]#, 12.5, 15]  # [0.01, 2.5, 15] # 

    ihInds = list(range(len(ihValues)))
    VLInds = list(range(len(VLValues)))

    #KgbarFactor = [1.0, 1.5]
    #ratesLongM2cM1 = [0.01,  2.5]

    ihSubset = [ihValues[i] for i in ihInds] #[0:5]]
    VLSubset = [VLValues[i] for i in VLInds] #[0:5]]


    popLabels = ['PT5B', 'IT5A', 'IT5B', 'upperPT5B', 'lowerPT5B', 'L5B']

    plt.style.use('seaborn-ticks') 

    rates = {l: np.zeros((len(ihValues), len(VLValues))) for l in popLabels}

    #load = True

    if load:
        simLabel = batchLabel
        with open('%s%s/%s_VL_ih_data_%s.pkl' % (dataFolder, batchLabel, batchLabel, method), 'rb') as f:
            rates=pickle.load(f)
        for iih, ih in enumerate(ihValues):
            for iVL, VL in enumerate(VLValues):
                print('')
                for popLabel in popLabels: 
                    print('ih=%.2f, VL=%.2f: pop %s=%.2f Hz' % (ih, VL, popLabel, rates[popLabel][iih, iVL]))
    else:

        for iih, ih in enumerate(ihValues):
            for iVL, VL in enumerate(VLValues):


                iratesLongM2cM1 = 1                 
                iKgbarFactor = 0 

                #iratesLongM2cM1 = 0 if iVL==0 else 1
                #iKgbarFactor = 0 if iih <=3 else 1

                # iratesLongM2cM1 = iVL
                # iKgbarFactor = iih

                #try:
                simLabel = '%s_%d_%d_%d_%d_%d' % (batchLabel, iratesLongM2cM1, iratesLongM2cM1, iih, iVL, iKgbarFactor) # v56_batch30
                #simLabel = '%s_%d_%d_%d_%d' % (batchLabel, iratesLongM2cM1, iratesLongM2cM1, iVL, iih) #, iKgbar) # v45_batch31,33
                #simLabel = '%s_%d_%d_%d' % (batchLabel, iih, iKgbarFactor, iVL) #, iKgbar)  # # v45_batch32

                # only do this preprocessing once
                if iih == 0 and iVL == 0:
                    from netpyne import sim
                    with open(dataFolder + 'v56_batch2_net/v56_batch2_net_0_0.pkl', 'rb') as f:
                        netLoad = pickle.load(f)

                    allCells = netLoad['net']['cells']
                    popGids, popNumCells = loadPopGids(dataFolder, popLabels)
                    
                # load just simData   
                try: 
                    with open(dataFolder+batchLabel+'/'+simLabel+'.json', 'r') as f:
                        data = json.load(f)
                    sim.allSimData = data['simData']
                    spkts = sim.allSimData['spkt']
                    spkids = sim.allSimData['spkid']
                    spkids,spkts = zip(*[(spkid,spkt-timeRange[0]) for spkid,spkt in zip(spkids,spkts) 
                        if timeRange[0] <= spkt <= timeRange[1]])

                    print('ih=%.2f, VL=%.2f' % (ih, VL))

                    # convert to dataframe
                    df=pd.DataFrame(np.array([spkts, spkids]).T, columns=['spkt', 'spkid'])
                    dfcount=df.groupby('spkid').count()
                    dfrates=dfcount / (timeRange[1] - timeRange[0]) * 1000

                    # calculate avg firing rates from spk times
                    if method == 'mean_>=0.01':
                        methodLabel = 'mean_>=0.01'

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            popRates = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(popRates.mean())
                            print(popLabel, rates[popLabel][iih, iVL])


                    ## mean, >=0Hz
                    elif method == 'mean_>=0.0':    
                        methodLabel = 'mean_>=0.0'
                        dfrates=dfrates.reindex(list(range(0, len(allCells))), fill_value=0)

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            popRates = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(popRates.mean())
                            print(popLabel, rates[popLabel][iih, iVL])
                        
                        #if iih == 3 and iVL == 1:
                        #    ipy.embed()

                    ## median, >=0.1Hz
                    elif method == 'median_>=0.01':
                        methodLabel = 'median_>=0.01'

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            spks = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(spks.median())
                            print(popLabel, rates[popLabel][iih, iVL])


                    ## median, >=0Hz
                    elif method == 'median_>=0.0':
                        methodLabel = 'median_>=0.0'
                        dfrates=dfrates.reindex(list(range(0, len(allCells))), fill_value=0)

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            spks = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(spks.median())
                            print(popLabel, rates[popLabel][iih, iVL])
                        
                except:
                    print('Error in file: %s' % (simLabel))
                    for popLabel in popLabels:
                        rates[popLabel][iih, iVL] = None

        with open('%s%s/%s_VL_ih_data_%s.pkl' % (dataFolder, batchLabel, batchLabel, methodLabel), 'wb') as f:
           pickle.dump(rates, f)

    # set subset
    ratesSub = {l: np.zeros((len(ihInds), len(VLInds))) for l in popLabels}

    for l in popLabels:
        for iihSub,ihSub in enumerate(ihInds):
            ratesSub[l][iihSub,:] = rates[l][ihSub, VLInds]



    # plot
    fontsiz = 22
    for label in popLabels: 

        if interpMissing:
            ihValues = [0, 0.25, 0.5, 0.75, 1.0]#, 1.5]   
            VLValues = [0.0, 2.5, 5, 7.5, 10, 12.5, 15]

            # add VL 7.5
            rates[label] = np.hstack((rates[label][:, 0:3], np.reshape(np.mean(rates[label][:, 2:4],1), (5,1)), rates[label][:, 3:5]))
            
            # add VL 12.5
            rates[label] = np.hstack((rates[label][:, 0:5], np.reshape(np.mean(rates[label][:, 4:6],1), (5,1)), rates[label][:, 5:7]))

        # ----------------------------------
        # full 
        plt.figure(figsize=(10, 8))
        plt.plot(rates[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLValues)), VLValues)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLValues, title='VL firing rate (Hz)')
        plt.legend(ihValues, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot.png' % (dataFolder, batchLabel, batchLabel, label))

        plt.figure(figsize=(14*0.65, 12*0.65))

        plt.imshow(rates[label].T, cmap = 'viridis', origin='lower', interpolation = 'spline16',vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([x for x in range(len(ihValues))], ihValues)
        plt.yticks([0, 2, 4, 6], [0, 5, 10, 15])  #([x for x in range(len(VLValues))], VLValues)
        #plt.xlim(0 - 0.5, len(ihValues) - 0.5)
        plt.xlim(0, max(ihValues))
        plt.ylim(0, max(VLValues))
        
        # plt.pcolor(rates[label].T, cmap='viridis')
        # plt.xticks([x+0.5 for x in range(len(ihValues))], ihValues)
        # plt.yticks([x+0.5 for x in range(len(VLValues))], VLValues)

        plt.ylabel('VL firing rate (Hz)', fontsize=fontsiz)
        plt.xlabel('PT5B Ih level (NA-R block)', fontsize=fontsiz)
        axisFontSize(plt.gca(), fontsiz)
        if label == 'L5B':
            cbar = plt.colorbar(ticks=[4, 5, 6, 7, 8])
        else:
            cbar = plt.colorbar()
        cbar.set_label('%s firing rate (Hz)' % label, fontsize=fontsiz)
        cbar.ax.tick_params(labelsize=fontsiz)

        #plt.scatter([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
        #            [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4], marker = 'o', s = 20, color="none", edgecolor = 'gray')

        plt.savefig('%s%s/%s_VLvsIh_%s_matrix_interpolated.png' % (dataFolder, batchLabel, batchLabel, label))

        plt.figure(figsize=(10, 8))
        plt.plot(rates[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLSubset)), VLSubset)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLSubset, title='VL firing rate (Hz)')
        plt.legend(ihSubset, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot.png' % (dataFolder, batchLabel, batchLabel, label))


        # ----------------------------------
        # subset (3x3)
        plt.figure(figsize=(14, 8))

        ratesSub_reverse = np.flipud(ratesSub[label])
        plt.imshow(ratesSub_reverse.T, cmap = 'viridis', origin='lower', interpolation = 'spline16', vmin=vmin, vmax=vmax)
        #plt.imshow(ratesSub[label].T, cmap = 'viridis', origin='lower', interpolation = 'bilinear', vmin=vmin, vmax=vmax)
        plt.xticks([x for x in range(len(ihSubset))], [int(x*100) for x in ihSubset])
        plt.yticks([x for x in range(len(VLSubset))], VLSubset)
        plt.xlim(0 - 0.5, len(ihSubset) - 0.5)
        #plt.xlim(0, len(ihSubset)-1)
        #plt.ylim(0, len(VLSubset)-1)
        
        # plt.pcolor(ratesSub[label].T, cmap='viridis')
        # plt.xticks([x+0.5 for x in range(len(ihValues))], ihValues)
        # plt.yticks([x+0.5 for x in range(len(VLSubset))], VLSubset)

        plt.ylabel('MTh input (Hz)', fontsize=fontsiz, fontweight='bold')
        #plt.xlabel('PT5B Ih level (NA-R block)', fontsize=fontsiz)
        plt.xlabel('NA input (%)', fontsize=fontsiz, fontweight='bold')
        axisFontSize(plt.gca(), fontsiz)
        cbar = plt.colorbar()
        cbar.set_label('%s firing rate (Hz)' % label, fontsize=fontsiz, fontweight='bold', labelpad=20)
        cbar.ax.tick_params(labelsize=fontsiz)

        #plt.scatter([0,0,0,1,1,1,2,2,2], [0,1,2,0,1,2,0,1,2], marker = 'o', s = 20, color="none", edgecolor = 'gray')
        print(label, ratesSub[label].min(), ratesSub[label].max())

        plt.savefig('%s%s/%s_VLvsIh_%s_matrix_interpolated_subset.png' % (dataFolder, batchLabel, batchLabel, label))

        #ipy.embed()

        plt.figure(figsize=(10, 8))
        plt.plot(ratesSub[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLSubset)), VLSubset)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('MTh firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLSubset, title='VL firing rate (Hz)')
        plt.legend(ihSubset, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot_subset.png' % (dataFolder, batchLabel, batchLabel, label))


    #ipy.embed()


# ----------------------------------------------------------------
def fig_VLvsIh(vmin=None, vmax=None, method = 'mean_>=0.01', interpMissing=False):

    dataFolder = '../../data/'


    '''
    1) quiet state, normal condition -> medium VL, medium ih; 2) movement state, normal condition -> high VL, low ih, low ih; 
    3) quiet state, motor thalamus inactive -> low VL, medium ih; 4) movement state, motor thalamus inactive -> low VL, low ih;
    5) quiet state, NA blocked -> medium VL, high ih; 6) movement state, NA blocked -> high VL, high ih.

    ih = [0.0, 0.25, 0.5, 0.75, 1.0]
    VL = [0.01, 2.5, 5, 10, 15] # /2 so equivalent to ratesLong uniform distrib (seems more similar to max)
    
    main, quiet: ih = 0.75, VL = 2.5  -> 3_1
    main, move : ih = 0.25, VL = 10   -> 1_3
    mth, quiet : ih = 0.75, VL = 0.01 -> 3_0
    mth, move  : ih = 0.25, VL = 0.01 -> 1_0
    na, quiet  : ih = 1.0, VL = 2.5   -> 4_1
    na, move   : ih = 1.0, VL = 10    -> 4_3
    '''

    # v56_batch2
    # batchLabel = 'v56_batch2'  #'v53_batch12' 
    
    # ihValues = [0, 0.25, 0.5, 1.0, 1.25]#, 1.5]   
    # VLValues = [0.01, 2.5, 5, 10, 15]
    # timeRange = [1000, 5000]

    # v56_batch5
    batchLabel = 'v56_batch5b'  #'v53_batch12' 
    timeRange = [2000, 3000]  # 1-5  
    
    ihValues = [0, 0.25, 0.5, 0.75, 1.0]  # [0.25, 0.75, 1.0] # 
    VLValues = [0.01, 2.5, 5, 10, 15]  # [0.01, 2.5, 15] # 

    # ihValues = [0.25, 0.75, 1.0]  # [0.25, 0.75, 1.0] # 
    # VLValues = [0.01, 2.5, 10]  # [0.01, 2.5, 15] # 

    ihInds = [1, 3, 4]
    VLInds = [0, 1, 3]

    ihSubset = [ihValues[i] for i in ihInds]
    VLSubset = [VLValues[i] for i in VLInds]


    popLabels = ['PT5B', 'IT5A', 'IT5B', 'upperPT5B', 'lowerPT5B', 'L5B']

    plt.style.use('seaborn-ticks') 

    rates = {l: np.zeros((len(ihValues), len(VLValues))) for l in popLabels}

    load = False

    if load:
        simLabel = batchLabel
        with open('%s%s/%s_VL_ih_data_%s.pkl' % (dataFolder, batchLabel, batchLabel, method), 'rb') as f:
            rates=pickle.load(f)
        for iih, ih in enumerate(ihValues):
            for iVL, VL in enumerate(VLValues):
                print('')
                for popLabel in popLabels: 
                    print('ih=%.2f, VL=%.2f: pop %s=%.2f Hz' % (ih, VL, popLabel, rates[popLabel][iih, iVL]))
    else:

        for iih, ih in enumerate(ihValues):
            for iVL, VL in enumerate(VLValues):

                #try:
                simLabel = '%s_%d_%d' % (batchLabel, iih, iVL)

                # only do this preprocessing once
                if iih == 0 and iVL == 0:
                    from netpyne import sim
                    with open(dataFolder + 'v56_batch2_net/v56_batch2_net_0_0.pkl', 'rb') as f:
                        netLoad = pickle.load(f)

                    allCells = netLoad['net']['cells']
                    simLabel = '%s_%d_%d' % (batchLabel, iih, iVL)

                    popGids, popNumCells = loadPopGids(dataFolder, popLabels)
                    
                # load just simData   
                try: 
                    with open(dataFolder+batchLabel+'/'+simLabel+'.json', 'r') as f:
                        data = json.load(f)
                    sim.allSimData = data['simData']
                    spkts = sim.allSimData['spkt']
                    spkids = sim.allSimData['spkid']
                    spkids,spkts = zip(*[(spkid,spkt-timeRange[0]) for spkid,spkt in zip(spkids,spkts) 
                        if timeRange[0] <= spkt <= timeRange[1]])

                    print('ih=%.2f, VL=%.2f' % (ih, VL))

                    # convert to dataframe
                    df=pd.DataFrame(np.array([spkts, spkids]).T, columns=['spkt', 'spkid'])
                    dfcount=df.groupby('spkid').count()
                    dfrates=dfcount / (timeRange[1] - timeRange[0]) * 1000

                    # calculate avg firing rates from spk times
                    if method == 'mean_>=0.01':
                        methodLabel = 'mean_>=0.01'

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            popRates = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(popRates.mean())
                            print(popLabel, rates[popLabel][iih, iVL])


                    ## mean, >=0Hz
                    elif method == 'mean_>=0.0':    
                        methodLabel = 'mean_>=0.0'
                        dfrates=dfrates.reindex(list(range(0, len(allCells))), fill_value=0)

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            popRates = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(popRates.mean())
                            print(popLabel, rates[popLabel][iih, iVL])
                        
                        #if iih == 3 and iVL == 1:
                        #    ipy.embed()

                    ## median, >=0.1Hz
                    elif method == 'median_>=0.01':
                        methodLabel = 'median_>=0.01'

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            spks = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(spks.median())
                            print(popLabel, rates[popLabel][iih, iVL])


                    ## median, >=0Hz
                    elif method == 'median_>=0.0':
                        methodLabel = 'median_>=0.0'
                        dfrates=dfrates.reindex(list(range(0, len(allCells))), fill_value=0)

                        for popLabel in popLabels:
                            popGidsList = popGids[popLabel]
                            spks = dfrates.query('spkid in @popGidsList')
                            rates[popLabel][iih, iVL] = float(spks.median())
                            print(popLabel, rates[popLabel][iih, iVL])
                        
                except:
                    print('Error in file: %s' % (simLabel))
                    for popLabel in popLabels:
                        rates[popLabel][iih, iVL] = None

        with open('%s%s/%s_VL_ih_data_%s.pkl' % (dataFolder, batchLabel, batchLabel, methodLabel), 'wb') as f:
           pickle.dump(rates, f)

    # set subset
    ratesSub = {l: np.zeros((len(ihInds), len(VLInds))) for l in popLabels}

    for l in popLabels:
        for iihSub,ihSub in enumerate(ihInds):
            ratesSub[l][iihSub,:] = rates[l][ihSub, VLInds]



    # plot
    fontsiz = 20
    for label in popLabels: 

        if interpMissing:
            ihValues = [0, 0.25, 0.5, 0.75, 1.0]#, 1.5]   
            VLValues = [0.0, 2.5, 5, 7.5, 10, 12.5, 15]

            # add VL 7.5
            rates[label] = np.hstack((rates[label][:, 0:3], np.reshape(np.mean(rates[label][:, 2:4],1), (5,1)), rates[label][:, 3:5]))
            
            # add VL 12.5
            rates[label] = np.hstack((rates[label][:, 0:5], np.reshape(np.mean(rates[label][:, 4:6],1), (5,1)), rates[label][:, 5:7]))

        # ----------------------------------
        # full 
        plt.figure(figsize=(10, 8))
        plt.plot(rates[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLValues)), VLValues)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLValues, title='VL firing rate (Hz)')
        plt.legend(ihValues, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot.png' % (dataFolder, batchLabel, batchLabel, label))

        plt.figure(figsize=(14*0.65, 12*0.65))

        plt.imshow(rates[label].T, cmap = 'viridis', origin='lower', interpolation = 'spline16',vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([x for x in range(len(ihValues))], ihValues)
        plt.yticks([0, 2, 4, 6], [0, 5, 10, 15])  #([x for x in range(len(VLValues))], VLValues)
        plt.xlim(0 - 0.5, len(ihValues) - 0.5)
        
        # plt.pcolor(rates[label].T, cmap='viridis')
        # plt.xticks([x+0.5 for x in range(len(ihValues))], ihValues)
        # plt.yticks([x+0.5 for x in range(len(VLValues))], VLValues)

        plt.ylabel('VL firing rate (Hz)', fontsize=fontsiz)
        plt.xlabel('PT5B Ih level (NA-R block)', fontsize=fontsiz)
        axisFontSize(plt.gca(), fontsiz)
        if label == 'L5B':
            cbar = plt.colorbar(ticks=[4, 5, 6, 7, 8])
        else:
            cbar = plt.colorbar()
        cbar.set_label('%s firing rate (Hz)' % label, fontsize=fontsiz)
        cbar.ax.tick_params(labelsize=fontsiz)

        #plt.scatter([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
        #            [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4], marker = 'o', s = 20, color="none", edgecolor = 'gray')

        plt.savefig('%s%s/%s_VLvsIh_%s_matrix_interpolated.png' % (dataFolder, batchLabel, batchLabel, label))

        plt.figure(figsize=(10, 8))
        plt.plot(rates[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLSubset)), VLSubset)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLSubset, title='VL firing rate (Hz)')
        plt.legend(ihSubset, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot.png' % (dataFolder, batchLabel, batchLabel, label))


        # ----------------------------------
        # subset (3x3)
        plt.figure(figsize=(14, 8))

        plt.imshow(ratesSub[label].T, cmap = 'viridis', origin='lower', interpolation = 'bilinear', vmin=vmin, vmax=vmax)
        plt.xticks([x for x in range(len(ihSubset))], ihSubset)
        plt.yticks([x for x in range(len(VLSubset))], VLSubset)
        plt.xlim(0 - 0.5, len(ihSubset) - 0.5)
        
        # plt.pcolor(ratesSub[label].T, cmap='viridis')
        # plt.xticks([x+0.5 for x in range(len(ihValues))], ihValues)
        # plt.yticks([x+0.5 for x in range(len(VLSubset))], VLSubset)

        plt.ylabel('VL firing rate (Hz)', fontsize=fontsiz)
        plt.xlabel('PT5B Ih level (NA-R block)', fontsize=fontsiz)
        axisFontSize(plt.gca(), fontsiz)
        cbar = plt.colorbar()
        cbar.set_label('%s firing rate (Hz)' % label, fontsize=fontsiz)
        cbar.ax.tick_params(labelsize=fontsiz)

        #plt.scatter([0,0,0,1,1,1,2,2,2], [0,1,2,0,1,2,0,1,2], marker = 'o', s = 20, color="none", edgecolor = 'gray')

        plt.savefig('%s%s/%s_VLvsIh_%s_matrix_interpolated_subset.png' % (dataFolder, batchLabel, batchLabel, label))

        plt.figure(figsize=(10, 8))
        plt.plot(ratesSub[label].T, '-o')
        #plt.xticks(range(len(ihValues)), ihValues)
        plt.xticks(range(len(VLSubset)), VLSubset)
        #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
        plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
        plt.ylabel('%s firing rate (Hz)' % label, fontsize=fontsiz)
        #plt.legend(VLSubset, title='VL firing rate (Hz)')
        plt.legend(ihSubset, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
        plt.savefig('%s%s/%s_VLvsIh_%s_lineplot_subset.png' % (dataFolder, batchLabel, batchLabel, label))


    #ipy.embed()



# ----------------------------------------------------------------
def fig_VLvsIh_exp(vmin=None, vmax=None, avg='mean', minRate=0.0):
    # plot experimental data
    '''
    - States to reproduce:
    -- quiet (wakefulness) --> ih = medium (1.0, 0.75, 0.25?); only bkg inputs (~5Hz)
    -- movement (self-paced, voluntary) --> ih = low (0.25) + increased VL or other inputs (bimodal: ~2Hz vs ~13Hz)

    -- inactive VL; quiet --> ih medium (?); only bkg inputs; VL = 0 (~1Hz)
    -- inactive VL; movement -->  ih low (0.25); only bkg inputs; VL = 0 (~2Hz)

    -- NA-R antagonist; quiet --> ih high (?); only bkg inputs0 (~2Hz)
    - -NA - R antagonist; movement - ->ih high(?); bkg inputs + high VL(bimodal: ~ 1 HZ vs ~ 7 Hz)
    '''

    dataFolder = '../../data/'
    batchLabel = 'v56_batch5b'  #'v53_batch12' 
    
    plt.style.use('seaborn-ticks') 

    ihValues = ['low', 'medium', 'high']
    VLValues = ['low', 'medium', 'high']

    import schi15
    # df = schi15.readData()
    # expRates = schi15.dataFromPaperFigs()  # getIhVsVLRates(df)
    df = schi15.readExcelFiringRatesAndMetada()
    df.dropna()

    #import IPython; IPython.embed()

    df = df.query('quietFR >= @minRate and moveFR >= @minRate')
    expRates = schi15.getIhVsVLRates(df, avg=avg)

    # extrapolate missing values
    # from scipy import interpolate
    # x = np.arange(0, 2.01, 1)
    # y = np.arange(0, 2.01, 1)
    # xx, yy = np.meshgrid(x, y)
    # z = np.sin(xx**2+yy**2)
    # f = interpolate.interp2d(x, y, z, kind='cubic')

    # xnew = np.arange(-5.01, 5.01, 1e-2)
    # ynew = np.arange(-5.01, 5.01, 1e-2)
    # znew = f(xnew, ynew)
    # plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
   
    expRates[0, 1] = (expRates[0, 0] + expRates[0, 2]) / 2.0
    expRates[1, 2] = (expRates[0, 2] + expRates[2, 2]) / 2.0

    from sklearn.linear_model import LinearRegression  
    x = np.array([2,3]).reshape((-1, 1))

    y = [expRates[1, 0], expRates[0, 0]]
    model = LinearRegression().fit(x,y)   
    x_new1 = model.predict(np.array([1]).reshape((-1,1)))  

    y = [expRates[2, 1], expRates[2, 2]]
    model = LinearRegression().fit(x,y)   
    x_new2 = model.predict(np.array([1]).reshape((-1,1)))  

    expRates[2, 0] = (x_new1 + x_new2) / 2.0  # based on trend

    #import IPython; IPython.embed()
    
    fontsiz = 22

    ## line plot
    plt.figure(figsize=(10, 8))
    plt.plot(expRates.T, '-o')
    #plt.xticks(range(len(ihValues)), ihValues)
    plt.xticks(range(len(VLValues)), VLValues)
    #plt.xlabel('PT5B Ih level (1.0 = in vitro)')
    plt.xlabel('VL firing rate (Hz)', fontsize=fontsiz) 
    plt.ylabel('L5B firing rate (Hz)', fontsize=fontsiz)
    #plt.legend(VLValues, title='VL firing rate (Hz)')
    plt.legend(ihValues, title='PT5B Ih level (1.0 = in vitro)', fontsize=fontsiz)
    plt.savefig('%s%s/%s_VLvsIh_experiment_lineplot_%s_%.2f.png' % (dataFolder, batchLabel, batchLabel, avg, minRate))


    ## 2D plot
    # plot first real points 
    plt.figure(figsize = (14, 8))

    # plt.pcolor(expRates.T, cmap = 'viridis')
    # plt.xticks([x+0.5 for x in range(len(ihValues))], ihValues)
    # plt.yticks([x+0.5 for x in range(len(VLValues))], VLValues)

    expRates_reverse = np.flipud(expRates)

    plt.imshow(expRates_reverse.T, cmap = 'viridis', origin='lower', interpolation = 'spline16', vmin=vmin, vmax=vmax)
    plt.xticks([x for x in range(len(ihValues))], ihValues)
    plt.yticks([x for x in range(len(VLValues))], VLValues)   
    #plt.xlim(0,len(ihValues)-1) 
    #plt.ylim(0,len(VLValues)-1) 
    
    plt.ylabel('MTh input', fontsize=fontsiz, fontweight='bold')
    #plt.xlabel('PT5B Ih level (NA-R block)', fontsize = fontsiz)
    plt.xlabel('NA input', fontsize = fontsiz, fontweight='bold')
  
    axisFontSize(plt.gca(), fontsiz)
    cbar = plt.colorbar()
    cbar.set_label('L5B firing rate (Hz)', fontsize = fontsiz, fontweight='bold', labelpad=20)
    cbar.ax.tick_params(labelsize=fontsiz)
    #plt.scatter([1.5, 0.5, 1.5, 0.5, 2.5, 2.5], [1.5, 2.5, 0.5, 0.5, 1.5, 2.5], marker='x', color='black')
    plt.scatter([1, 0, 1, 0, 2, 2], [1, 2, 0, 0, 1, 2], marker = 'o', s = 70, color="none", edgecolor = 'red', lw=3)
    
    #plt.contour(ihValues, VLValues, expRates.T, corner_mask=True,  extend='both') #cmap='viridis',

    print('Experiment ', expRates.min(), expRates.max())

    plt.savefig('%s%s/%s_VLvsIh_experiment_matrix_interpolated_reverse_%s_%.2f.png' % (dataFolder, batchLabel, batchLabel, avg, minRate))


    #ipy.embed()

    return expRates



def plot_corner_rasters():

    # -------------------------------------------------------------------------------------------
    # Raster plot
    dataFolder = '../../data/'
    batchLabel = 'v56_batch30'  
    # iratesLongM2cM1, iratesLongM2cM1, iih, iVL, iKgbarFactor
    simLabels = ['v56_batch30_1_1_0_0_0', 'v56_batch30_1_1_0_4_0', 'v56_batch30_1_1_4_0_0', 'v56_batch30_1_1_4_4_0']

    # if includeiratesLongM2cM1 and iKgbarFactor effects
    #iratesLongM2cM1 = 0 if iVL==0 else 1
    #iKgbarFactor = 0 if iih <=3 else 1
    simLabels = ['v56_batch30_0_0_0_0_0', 'v56_batch30_1_1_0_4_0', 'v56_batch30_0_0_4_0_1', 'v56_batch30_1_1_4_4_1']


    for simLabel in simLabels:
        #simLabel = 'v56_batch19_0_0_0_0'

        sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

        timeRange = [1000, 5000] #[2000, 4000]
        include = allpops
        orderBy = ['pop', 'y']
        #filename = '%s%s_raster_%d_%d_%s.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
        fig1 = sim.analysis.plotRaster(include=include, timeRange=timeRange, labels='overlay', 
            popRates=0, orderInverse=True, lw=0, markerSize=3.5, marker='.', popColors=popColors, 
            showFig=0, saveFig=0, figSize=(8.5*1,10), orderBy=orderBy)# 
        ax = plt.gca()

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner
        #plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')  #remove ticks
        plt.xticks([1000, 5000], ['1', '5'])
        plt.yticks([0, 5000, 10000], [0, 5000, 10000])
        
        fontsiz=16
        plt.ylabel('Neuron ID', fontsize=fontsiz) #Neurons (ordered by NCD within each pop)')
        plt.xlabel('Time (s)', fontsize=fontsiz)
        axisFontSize(plt.gca(), fontsiz)

        plt.title('')
        filename='%s%s_raster_%d_%d_%s_v2.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
        plt.savefig(filename, dpi=600)


# Main code
if __name__ == '__main__': 
    #for method in ['mean_>=0.01', 'mean_>=0.0', 'median_>=0.01', 'median_>=0.0']:
    #    fig_VLvsIh(method=method) #vmin=3.0, vmax=7.5)#vmin=3.0, vmax=6.5)
    #    input("Press Enter to continue...")

    #fig_VLvsIh_test(method='mean_>=0.0', load=False , interpMissing=False)

    fig_VLvsIh_exp(avg='mean', minRate=0.0) #vmin=3.0, vmax=7.5)# vmin=3.0, vmax=6.5)

    #plot_corner_rasters()
