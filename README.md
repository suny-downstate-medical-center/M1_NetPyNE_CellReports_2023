# Primary motor cortex (M1) circuits model
## Description
Multiscale model of mouse primary motor cortex (M1) developed using NetPyNE (www.netpyne.org).

The model is described in the following paper: 

- Dura-Bernal S, Neymotin SA, Suter BA, Dacre J, Joao VS Moreira, Eugenio Urdapilleta, Schiemann J, Duguid I, Shepherd GMG, Lytton WW. **Multiscale model of primary motor cortex circuits predicts in vivo cell type-specific, behavioral state-dependent dynamics.** BioRxiv 2022.02.03.479040; doi: https://doi.org/10.1101/2022.02.03.479040. (Under review in Cell Reports) 

A previous version of the model was also described in the following papers:

- Dura-Bernal S, Suter B, Gleeson P, Cantarelli M, Quintana A, Rodriguez F, Kedziora DJ, Chadderdon GL, Kerr CC, Neymotin SA, McDougal R, Hines M, Shepherd GMG, Lytton WW. **NetPyNE: a tool for data-driven multiscale modeling of brain circuits.** eLife 2019;8:e44494 

- Sivagnanam S, Gorman W, Doherty D, Neymotin SA, Fang S, Hovhannisyan H, Lytton WW, Dura-Bernal S. **Simulating large-scale models of brain neuronal circuits using Google Cloud Platform.** Practice and Experience in Advanced Research Computing, PEARC2020. 10.1145/3311790.3399621

A modified version of the model was used to benchmark CoreNEURON in the following paper:

- Awile O, Kumbhar P, Cornu N, Dura-Bernal S, Gonzalo JK, Lupton O, Magkanaris I, McDougal R, Newton AJH, Pereira A, Savulescu A, Carnevale NT, Hines M, Lytton WW, Schurmann F. **Modernizing the NEURON Simulator for Sustainability, Portability, and Performance**. Frontiers in Neuroinformatics (Under Revision). Research Topic: "Neuroscience, Computing, Performance, and Benchmarks: Why It Matters to Neuroscience How Fast We Can Compute." 


## Setup and execution

Requires NEURON with Python and MPI support. 

1. From /sim run `nrnivmodl ../mod`. This should create a directory called x86_64. 
2. To run type: `./runsim [num_proc]' or the equivalent `mpiexec -np [num_proc] nrniv -python -mpi init.py`

## Overview of file structure:

* /sim/init.py: Main executable; calls functions from other modules. Sets what parameter file to use.

* /sim/netParams.py: Network parameters

* /sim/cfg.py: Simulation configuration

* /cells: source .py, .hoc, .json or .pkl files for the different cell types used in the model; these will be imported into netpyne

* /conns: source .py, .json or .pkl files for the network connectivity; these will be imported into netpyne

* /mod: NMODL files containing the ionic channel and synaptic mechanisms used in the model 

* /data: where the model and simulation data is stored 


For further information please contact: salvadordura@gmail.com 

