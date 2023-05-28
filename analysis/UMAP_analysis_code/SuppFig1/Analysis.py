import pickle, json
import matplotlib.pyplot as plt
from netpyne import specs,sim
import numpy as np
import scipy
import os
import umap
import random
from matplotlib import cm

Trials = ['0_0','0_1','0_2','0_3','0_4',
          '1_0','1_1','1_2','1_3','1_4',
          '2_0','2_1','2_2','2_3','2_4',
          '3_0','3_1','3_2','3_3','3_4',
          '4_0','4_1','4_2','4_3','4_4']

Ndata = 1  # One condition
reading_folder  = ['Data_Control']
pick_label_pre = ['v56_batch23_0_0_']

pick_label = []
for label in pick_label_pre:
    for trial in Trials: pick_label.append(label+trial)

Twindow = 25
WindowOverlap = 0
Tstart = 3000
Ton = 5000
Tend = 7000

## DATA PREPROCESSING
data = {}
spkt = {}
spkid = {}

NcellsSample = 1000

## UMAP
n_neighbors = 100
min_dist = 0.95
metric = 'euclidean'


###############################################################################
# Structuring analysis
if isinstance(Twindow,int): Twindow = [Twindow]

for label in pick_label:
    if label.startswith('v56_batch23'):
        folder = os.path.join(os.getcwd(),reading_folder[0])
        filename = label+'.pkl'
        filename = os.path.join(folder,filename)
        print(('Loading file %s ... ' % (filename)))
        with open(filename, 'rb') as fileObj:
            data[label] = pickle.load(fileObj, encoding='latin1')

    spkt[label] = np.array(data[label]['simData']['spkt'])
    spkid[label] = np.array(data[label]['simData']['spkid'])


###############################################################################
# Creating the populations
L5Bmin=0.47
L5Bmax=0.8
L5Bmid = L5Bmin + (L5Bmax-L5Bmin)/2
pops_labels = ['IT2','SOM2','PV2','IT4','IT5A','SOM5A','PV5A','IT5B','upperPT5B','lowerPT5B','SOM5B','PV5B','IT6','CT6','SOM6','PV6']
Npops = len(pops_labels)

# Create the network according to batch23
if  pick_label[0].startswith('v56_batch23'):
    cfg = specs.SimConfig(data[pick_label[0]]['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data[pick_label[0]], instantiate=False)
    sim.setSimCfg(cfg)

    sim.net.createPops()
    sim.net.createCells()

    pops = specs.ODict({'IT2'      :[c.gid for c in sim.net.cells if c.tags['pop']=='IT2'],
                        'SOM2'     :[c.gid for c in sim.net.cells if c.tags['pop']=='SOM2'],
                        'PV2'      :[c.gid for c in sim.net.cells if c.tags['pop']=='PV2'],
                        'IT4'      :[c.gid for c in sim.net.cells if c.tags['pop']=='IT4'],
                        'IT5A'     :[c.gid for c in sim.net.cells if c.tags['pop']=='IT5A'],
                        'SOM5A'    :[c.gid for c in sim.net.cells if c.tags['pop']=='SOM5A'],
                        'PV5A'     :[c.gid for c in sim.net.cells if c.tags['pop']=='PV5A'],
                        'IT5B'     :[c.gid for c in sim.net.cells if c.tags['pop']=='IT5B'],
                        'upperPT5B':[c.gid for c in sim.net.cells if L5Bmin <= c.tags['ynorm'] <= L5Bmid and c.tags['pop']=='PT5B'],
                        'lowerPT5B':[c.gid for c in sim.net.cells if L5Bmid <= c.tags['ynorm'] <= L5Bmax and c.tags['pop']=='PT5B'],
                        'SOM5B'    :[c.gid for c in sim.net.cells if c.tags['pop']=='SOM5B'],
                        'PV5B'     :[c.gid for c in sim.net.cells if c.tags['pop']=='PV5B'],
                        'IT6'      :[c.gid for c in sim.net.cells if c.tags['pop']=='IT6'],
                        'CT6'      :[c.gid for c in sim.net.cells if c.tags['pop']=='CT6'],
                        'SOM6'     :[c.gid for c in sim.net.cells if c.tags['pop']=='SOM6'],
                        'PV6'      :[c.gid for c in sim.net.cells if c.tags['pop']== 'PV6'],
                        # External inputs
                        'TPO'      :[c.gid for c in sim.net.cells if c.tags['pop']== 'TPO'],
                        'TVL'      :[c.gid for c in sim.net.cells if c.tags['pop']== 'TVL'],
                        'S1'       :[c.gid for c in sim.net.cells if c.tags['pop']== 'S1'],
                        'S2'       :[c.gid for c in sim.net.cells if c.tags['pop']== 'S2'],
                        'cM1'      :[c.gid for c in sim.net.cells if c.tags['pop']== 'cM1'],
                        'M2'       :[c.gid for c in sim.net.cells if c.tags['pop']== 'M2'],
                        'OC'       :[c.gid for c in sim.net.cells if c.tags['pop']== 'OC']
                        })
    NcellsTotal = sum([len(pops[pop]) for pop in pops])
    cells_pop = [0]*NcellsTotal
    for pop in pops:
        for nn in range(len(pops[pop])):
            cells_pop[pops[pop][nn]] = pop

Tduration = data[pick_label[0]]['simConfig']['duration']
Ncells = sum(len(pops[pop]) for pop in pops_labels)       # excluding stims

###############################################################################

# Selected neurons
validPops = pops_labels
Id_validCells = [id for id in range(NcellsTotal) if cells_pop[id] in validPops]

activityVector = {}

for Tw in Twindow:
    print(Tw)

    t_stamps = np.linspace(Tw/2.0,Tduration-Tw/2.0,int((Tduration-Tw)/((1.0-WindowOverlap)*Tw))+1,endpoint=True)
    t_limits = [(t_stamps[nt]-Tw/2.0,t_stamps[nt]+Tw/2.0) for nt in range(len(t_stamps))]
    index_start = [Tstart<=t_stamps[nt] for nt in range(len(t_stamps))].index(True)
    index_end =   [Tend<=t_stamps[nt] for nt in range(len(t_stamps))].index(True)

    activityVector[str(Tw)] = {}
    for label in pick_label:

        activityVector[str(Tw)][label] = np.zeros((Ncells,len(t_stamps)))

        for n in range(len(spkt[label])):
            neuron = int(spkid[label][n])
            if cells_pop[neuron] in pops_labels:   # cells and not stimuli
                nt = [t_limits[nt][0] <= spkt[label][n] < t_limits[nt][1] for nt in range(len(t_limits))].index(True)
                activityVector[str(Tw)][label][neuron][nt] = activityVector[str(Tw)][label][neuron][nt] + 1000.0/Tw

        if label == pick_label[0]:
            Matrix = activityVector[str(Tw)][label][Id_validCells,index_start:index_end]
        else:
            Matrix = np.concatenate((Matrix,activityVector[str(Tw)][label][Id_validCells,index_start:index_end]),axis=1)

    mean_activity = np.sum(Matrix,1)*(Tw/1000.0) / ((Tend-Tstart)/1000.0)   # number of spikes in [Tstart,Tend] / T
    Valid_Ids = [n for n in range(len(mean_activity)) if mean_activity[n]!=0]
    Ids = random.sample(Valid_Ids, NcellsSample)
    Ids.sort()

    x = np.transpose(Matrix[Ids,:])

    for n_components in [3]: #[2,3]:    
        umap_reduction = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,n_components=n_components,metric=metric).fit(x)
        umap_representation = umap_reduction.transform(x)
        umap_representation_back = umap_reduction.inverse_transform(umap_representation)

        # Plot as in the paper (not exactly the same figure because sampled cells are different)
        if n_components==3 and Tw==25:

            time = t_stamps[index_start:index_end]
            index_on = [Ton<=time[nt] for nt in range(len(time))].index(True)
            index_end = len(time)

            orig_data = {}
            umap_data = {}
            umap_back_data = {}
            for index,label in enumerate(pick_label):
                orig_data[label] = x[index*len(time):(index+1)*len(time),:]
                umap_data[label] = umap_representation[index*len(time):(index+1)*len(time),:]
                umap_back_data[label] = umap_representation_back[index*len(time):(index+1)*len(time),:]

            # Plotting
            tab20 = cm.get_cmap('tab20', 20)
            cmap1 = tab20(range(20))[[1,1,1,1,1, 3,3,3,3,3, 5,5,5,5,5, 7,7,7,7,7, 13,13,13,13,13]]
            cmap2 = tab20(range(20))[[0,0,0,0,0, 2,2,2,2,2, 4,4,4,4,4, 6,6,6,6,6, 12,12,12,12,12]]

            # All 25 (5 animals x 5 trials)
            trials_transparent = ['0_1','0_2','0_3','0_4',
                                  '1_1','1_2','1_3','1_4',
                                  '2_1','2_2','2_3','2_4',
                                  '3_1','3_2','3_3','3_4',
                                  '4_1','4_2','4_3','4_4']

            for az in np.linspace(0,360,37):
                fig = plt.figure(figsize=(20, 16), dpi=300)
                ax = fig.gca(projection='3d')
                ax.azim = az
                ax.dist = 10
                ax.elev = 25

                lab = pick_label_pre[0]
                for j,trial in enumerate(Trials):
                    label = str(lab) + str(trial)

                    transparency = 0.75
                    if trial in trials_transparent: 
                        transparency = 0.1

                    ax.plot3D(umap_data[label][0:index_on+1,0], umap_data[label][0:index_on+1,1], umap_data[label][0:index_on+1,2], color=cmap1[j], linewidth=1, alpha = transparency)
                    ax.plot3D(umap_data[label][index_on:index_end,0], umap_data[label][index_on:index_end,1], umap_data[label][index_on:index_end,2], color=cmap2[j], linewidth=1, alpha = transparency)
                    ax.scatter3D(umap_data[label][0:index_on,0], umap_data[label][0:index_on,1], umap_data[label][0:index_on,2], color=cmap1[j], s=50, alpha = transparency)
                    ax.scatter3D(umap_data[label][index_on:index_end,0], umap_data[label][index_on:index_end,1], umap_data[label][index_on:index_end,2], color=cmap2[j], s=50, edgecolors = (0,0,0,1), alpha = transparency)

                ax.set_top_view()
                    
                ax.set_xlabel('$C_1$', fontsize=40, labelpad=40)
                ax.set_ylabel('$C_2$', fontsize=40, labelpad=40)
                ax.set_zlabel('$C_3$', fontsize=40, labelpad=40)

                ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

                plt.savefig('ReducedDim-'+str(az)+'.png')
                plt.close(fig)