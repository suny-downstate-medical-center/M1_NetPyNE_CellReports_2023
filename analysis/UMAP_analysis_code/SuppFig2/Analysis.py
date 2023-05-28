import pickle, json
import matplotlib.pyplot as plt
from netpyne import specs,sim
import numpy as np
import scipy
import os
import umap
import random

Trials = ['0_0']

Ndata = 3   # Three conditions
reading_folder  = ['Data_Control', 'Data_MotorThalamusInactivation', 'Data_NoradrenalineBlock']
pick_label_pre = ['v56_batch23_0_0_', 'v56_batch25_2_2_2_2_2_2_2_0_', 'v56_batch29_7_']

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
            
    if label.startswith('v56_batch25'):
        folder = os.path.join(os.getcwd(),reading_folder[1])
        filename = label+'.json'
        filename = os.path.join(folder,filename)
        print(('Loading file %s ... ' % (filename)))
        with open(filename, 'rb') as fileObj:
            data[label] = json.load(fileObj, encoding='latin1')

    if label.startswith('v56_batch29'):
        folder = os.path.join(os.getcwd(),reading_folder[2])
        filename = label+'.json'
        filename = os.path.join(folder,filename)
        print(('Loading file %s ... ' % (filename)))
        with open(filename, 'rb') as fileObj:
            data[label] = json.load(fileObj, encoding='latin1')

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


            # Plot data
            color1 = [ (0.175, 0.35, 0.66, 1), (0.47, 0.81, 0.92, 1), (0.55, 0.76, 0.36, 1)]
            color2 = [ (0.90, 0.21, 0.16, 1), (0.61, 0.34, 0.67, 1), (0.95, 0.70, 0.36, 1)]
            
            for az in np.linspace(0,360,37):
                fig = plt.figure(figsize=(20, 16), dpi=300)
                ax = fig.gca(projection='3d')
                ax.azim = az
                ax.dist = 10
                ax.elev = 25

                for i,lab in enumerate(pick_label_pre):
                    for trial in Trials:
                        label = str(lab) + str(trial)

                        import csaps
                        theta = np.arange(len(umap_data[label]))
                        data = np.transpose(umap_data[label])
                        sp_theta = csaps.MultivariateCubicSmoothingSpline(data, theta, smooth=1.0)
                        theta_i = np.linspace(0, len(umap_data[label])-1, 100*len(umap_data[label]))
                        data_i = np.transpose(sp_theta(theta_i))

                        ax.scatter3D(umap_data[label][0:index_on,0], umap_data[label][0:index_on,1], umap_data[label][0:index_on,2], color=color1[i], s=40, alpha=1.0)
                        ax.scatter3D(umap_data[label][index_on:index_end,0], umap_data[label][index_on:index_end,1], umap_data[label][index_on:index_end,2], color=color2[i], s=40, alpha=1.0)
                    
                        ax.plot3D(data_i[0:int(len(data_i)/2),0], data_i[0:int(len(data_i)/2),1], data_i[0:int(len(data_i)/2),2], color=color1[i], linewidth=1.0)
                        ax.plot3D(data_i[int(len(data_i)/2):len(data_i),0], data_i[int(len(data_i)/2):len(data_i),1], data_i[int(len(data_i)/2):len(data_i),2], color=color2[i], linewidth=1.0)

                        startpoint = ax.scatter3D(umap_data[label][0,0], umap_data[label][0,1], umap_data[label][0,2], color=color1[i], s=150, alpha = 1.0, label='Start (Control)')
                        transitionpoint = ax.scatter3D(umap_data[label][index_on,0], umap_data[label][index_on,1], umap_data[label][index_on,2], color=color2[i], edgecolor = 'k', linewidth=2, s=150, alpha = 1.0, label='Transition (Control)')
                        endpoint = ax.scatter3D(umap_data[label][index_end-1,0], umap_data[label][index_end-1,1], umap_data[label][index_end-1,2], color=color2[i], s=150, alpha = 1.0, label='End (Control)')

                #ax.set_title('$\Delta T = $'+str(Tw)+'$~ ms$', fontsize=40)
                ax.set_top_view()

                ax.set_xlabel('$C_1$', fontsize=40, labelpad=25)
                ax.set_ylabel('$C_2$', fontsize=40, labelpad=25)
                ax.set_zlabel('$C_3$', fontsize=40, labelpad=25)

                ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

                plt.savefig('ReducedDim-'+str(az)+'.png')
                plt.close(fig)
