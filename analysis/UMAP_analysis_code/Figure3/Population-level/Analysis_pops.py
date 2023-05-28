import pickle, json
import matplotlib.pyplot as plt
from netpyne import specs,sim
import numpy as np
import scipy
import os
import umap

Nseeds_conn = 1
Nseeds_stim = 1

Ndata = 3   # Three conditions
reading_folder  = ['Data_Control', 'Data_MotorThalamusInactivation', 'Data_NoradrenalineBlock']
pick_label = ['v56_batch23_0_0_0_0', 'v56_batch20_0_1_0_0_0', 'v56_batch22_0_0']

Twindow = 10
WindowOverlap = 0
Tstart = 3000
Ton = 5000
T0 = 5660
T1 = 6000
T2 = 6040
Tend = 7000

## DATA PREPROCESSING
data = {}
spkt = {}
spkid = {}
spkt_pops = {}
spkid_pops = {}

## UMAP
n_neighbors = 15
min_dist = 0.30
metric = 'euclidean'


###############################################################################
# Structuring analysis
if isinstance(Twindow,int): Twindow = [Twindow]

for i,label in enumerate(pick_label):
    if label.startswith('v56_batch23'): # pkl file
        folder = os.path.join(os.getcwd(),reading_folder[0])
        filename = label+'.pkl'
        filename = os.path.join(folder,filename)
        print(('Loading file %s ... ' % (filename)))
        with open(filename, 'rb') as fileObj:
            data[label] = pickle.load(fileObj, encoding='latin1')
            
    if label.startswith('v56_batch20'): # json file
        folder = os.path.join(os.getcwd(),reading_folder[1])
        filename = label+'.json'
        filename = os.path.join(folder,filename)
        print(('Loading file %s ... ' % (filename)))
        with open(filename, 'rb') as fileObj:
            data[label] = json.load(fileObj, encoding='latin1')
            
    if label.startswith('v56_batch22'): # json file
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

Tduration = data[label]['simConfig']['duration']
Ncells = sum(len(pops[pop]) for pop in pops_labels)       # excluding stims

###########################################################################
## Structuring spiking data
for label in pick_label:
    spkt_pops[label] = {}
    spkid_pops[label] = {}
    for pop in pops:
        spkt_pops[label][pop] = []
        spkid_pops[label][pop] = []

    for n in range(len(spkt[label])):
        pop = cells_pop[int(spkid[label][n])]
        spkt_pops[label][pop].append(spkt[label][n])
        spkid_pops[label][pop].append(int(spkid[label][n]))


## Generating activity vectors along time
actVector_pops = {}
for Tw in Twindow:
    print(Tw)

    t_stamps = np.linspace(Tw/2.0,Tduration-Tw/2.0,int((Tduration-Tw)/((1.0-WindowOverlap)*Tw))+1,endpoint=True)
    t_limits = [(t_stamps[nt]-Tw/2.0,t_stamps[nt]+Tw/2.0) for nt in range(len(t_stamps))]
    index_start = [Tstart<=t_stamps[nt] for nt in range(len(t_stamps))].index(True)
    index_end =   [Tend<=t_stamps[nt] for nt in range(len(t_stamps))].index(True)
    
    actVector_pops[str(Tw)] = {}
    for label in pick_label:

        actVector_pops[str(Tw)][label] = np.zeros((Npops,len(t_stamps)))

        for n in range(len(spkt[label])):
            neuron = int(spkid[label][n])
            if cells_pop[neuron] in pops_labels:   # cells and not stimuli
                nt = [t_limits[nt][0] <= spkt[label][n] < t_limits[nt][1] for nt in range(len(t_limits))].index(True)
                pop = list(pops.keys()).index(cells_pop[neuron])
                actVector_pops[str(Tw)][label][pop][nt] = actVector_pops[str(Tw)][label][pop][nt] + 1000.0/(Tw*len(pops[cells_pop[neuron]]))

        if label == pick_label[0]:
            Matrix = actVector_pops[str(Tw)][label][:,index_start:index_end]
        else:
            Matrix = np.concatenate((Matrix,actVector_pops[str(Tw)][label][:,index_start:index_end]),axis=1)

    x = np.transpose(Matrix)
    
    for n_components in [3]: #[2,3]:  
        umap_reduction = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,n_components=n_components,metric=metric).fit(x)
        umap_representation = umap_reduction.transform(x)
        umap_representation_back = umap_reduction.inverse_transform(umap_representation)

        # Plot as in the paper
        if n_components==3 and Tw==10:

            time = t_stamps[index_start:index_end]
            index_on = [Ton<=time[nt] for nt in range(len(time))].index(True)
            index0 = [T0<=time[nt] for nt in range(len(time))].index(True)
            index1 = [T1<=time[nt] for nt in range(len(time))].index(True)
            index2 = [T2<=time[nt] for nt in range(len(time))].index(True)
            index_end = len(time)

            orig_data = {}
            umap_data = {}
            umap_back_data = {}
            for index,label in enumerate(pick_label):
                orig_data[label] = x[index*len(time):(index+1)*len(time),:]
                umap_data[label] = umap_representation[index*len(time):(index+1)*len(time),:]
                umap_back_data[label] = umap_representation_back[index*len(time):(index+1)*len(time),:]

  
            # Plot data
            color1 = (0.175, 0.35, 0.66, 1)
            color2 = (0.90, 0.21, 0.16, 1)

            # only first trial - control condition
            label = pick_label[0]
            
            import csaps
            theta = np.arange(len(umap_data[label]))
            data_traj = np.transpose(umap_data[label])
            sp_theta = csaps.MultivariateCubicSmoothingSpline(data_traj, theta, smooth=1.0)                        
            theta_i = np.linspace(0, len(umap_data[label])-1, 100*len(umap_data[label]))
            data_i = np.transpose(sp_theta(theta_i))

            for az in np.linspace(0,360,37):
                fig = plt.figure(figsize=(20, 16), dpi=300)
                ax = fig.gca(projection='3d')
                ax.azim = az
                ax.dist = 10
                ax.elev = 20

                ax.scatter3D(umap_data[label][0:index_on,0], umap_data[label][0:index_on,1], umap_data[label][0:index_on,2], color=color1, s=25, alpha = 1.0)
                ax.scatter3D(umap_data[label][index_on:index0,0], umap_data[label][index_on:index0,1], umap_data[label][index_on:index0,2], color=color2, s=25, edgecolors = color1, alpha = 1.0)
                ax.scatter3D(umap_data[label][index0:index_end,0], umap_data[label][index0:index_end,1], umap_data[label][index0:index_end,2], color=color2, s=25, alpha = 1.0)

                ax.plot3D(data_i[0:int(len(data_i)/2),0], data_i[0:int(len(data_i)/2),1], data_i[0:int(len(data_i)/2),2], color=color1, linewidth=1.0)
                ax.plot3D(data_i[int(len(data_i)/2):len(data_i),0], data_i[int(len(data_i)/2):len(data_i),1], data_i[int(len(data_i)/2):len(data_i),2], color=color2, linewidth=1.0)

                startpoint = ax.scatter3D(umap_data[label][0,0], umap_data[label][0,1], umap_data[label][0,2], color=color1, s=150, alpha = 1.0, label='Start')
                transitionpoint = ax.scatter3D(umap_data[label][index_on,0], umap_data[label][index_on,1], umap_data[label][index_on,2], color=color2, edgecolor = 'k', linewidth=2, s=150, alpha = 1.0, label='Transition')
                endpoint = ax.scatter3D(umap_data[label][index_end-1,0], umap_data[label][index_end-1,1], umap_data[label][index_end-1,2], color=color2, s=150, alpha = 1.0, label='End')
                #ax.legend(handles=[startpoint,endpoint],fontsize = 25)
                        
                #ax.plot3D(umap_data[str(Tw)][label][index1:index2,0], umap_data[str(Tw)][label][index1:index2,1], umap_data[str(Tw)][label][index1:index2,2], color='k', linestyle='--', linewidth=4)
                ax.scatter3D(umap_data[label][index1:index2,0], umap_data[label][index1:index2,1], umap_data[label][index1:index2,2], color='k', s=100)
                index1 = 30074
                index2 = 30375
                ax.plot3D(data_i[index1:index2,0], data_i[index1:index2,1], data_i[index1:index2,2], color='k', linestyle='--', linewidth=4)

                ax.set_top_view()
                    
                ax.set_xlabel('$C_1$', fontsize=40, labelpad=35)
                ax.set_xticks([4, 6, 8, 10, 12])
                ax.set_xticklabels([])
                    
                ax.set_ylabel('$C_2$', fontsize=40, labelpad=35)
                ax.set_yticks([0, 2, 4, 6, 8, 10, 12])
                ax.set_yticklabels([])

                ax.set_zlabel('$C_3$', fontsize=40, labelpad=35)
                ax.set_zticks([3, 5, 7, 9])
                ax.set_zticklabels([])

                ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                    
                ax.set_xlim([3, 13])
                ax.set_ylim([-0.5, 12.5])
                ax.set_zlim([2, 10])
                    
                plt.savefig('ReducedDim-'+str(az)+'.png')
                plt.close(fig)

            # Reconstruction
            fig, (ax1, ax2) = plt.subplots(1,2, dpi=300)
            
            minSpikeCount = -5
            maxSpikeCount = 60

            (pearsonCorr,pvalue) = scipy.stats.pearsonr(orig_data[label].flatten(), umap_back_data[label].flatten())

            fig.set_figheight(7.5)
            fig.set_figwidth(22.5)
        
            ax1.set_title('Original activity data', fontsize=40, pad=25)
            im1 = ax1.imshow(np.transpose(orig_data[label][:,:]),cmap='viridis', aspect='auto')
            cbar1 = plt.colorbar(im1,ax=ax1)
            cbar1.ax.tick_params(labelsize=35)
            ax1.set_xticks([-0.5,len(time)/2 - 0.5,len(time) - 0.5])
            ax1.set_xticklabels([0 , 2, 4])
            ax1.set_yticks([0,3,6,9,12,15])
            ax1.tick_params(labelsize=35)
            ax1.set_xlabel('Time (s)', fontsize=40)
            ax1.set_ylabel('Population ID', fontsize=40)
            im1.set_clim(minSpikeCount, maxSpikeCount)
            ax1.axvline(len(time)/2,color='white',dashes=[8, 4],linewidth=4)
            ax1.axvline(index0,color='silver',dashes=[4, 4],linewidth=4)
            
            
            ax2.set_title(r'Reconstruction ($\rho$ = %.3f)' % (pearsonCorr), fontsize=40, pad=25)
            im2 = ax2.imshow(np.transpose(umap_back_data[label][:,:]),cmap='viridis', aspect='auto')
            cbar2 = plt.colorbar(im2,ax=ax2)
            cbar2.ax.tick_params(labelsize=35)
            ax2.set_xticks([-0.5,len(time)/2 - 0.5,len(time) - 0.5])
            ax2.set_xticklabels([0 , 2, 4])
            ax2.set_yticks([0,3,6,9,12,15])
            ax2.tick_params(labelsize=35)
            ax2.set_xlabel('Time (s)', fontsize=40)
            cbar2.set_label('Rate (Hz)', fontsize=40)
            im2.set_clim(minSpikeCount, maxSpikeCount)
            ax2.axvline(len(time)/2,color='white',dashes=[8, 4],linewidth=4)
            ax2.axvline(index0,color='silver',dashes=[4, 4],linewidth=4)
        
            plt.tight_layout()
        
            plt.savefig('ReconstructionMap.png')
            plt.close(fig)

