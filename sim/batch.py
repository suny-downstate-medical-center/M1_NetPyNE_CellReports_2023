"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Weight Normalization Exc
# ----------------------------------------------------------------------------------------------
def weightNormE(pops=['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2'], 
    segs = None, allSegs = True, rule = 'IT2_reduced', weights=list(np.arange(0.01, 0.2, 0.01)/100.0)):

    # Add params
    from cfg_cell import cfg
    from netParams_cell import netParams

    excludeSegs = ['axon']
    if not segs:
        secs = []
        locs = []
        for secName,sec in netParams.cellParams[rule]['secs'].iteritems():
            if secName not in excludeSegs:
                if allSegs:
                    nseg = sec['geom']['nseg']
                    for iseg in range(nseg):
                        secs.append(secName) 
                        locs.append((iseg+1)*(1.0/(nseg+1)))
                else:
                    secs.append(secName) 
                    locs.append(0.5)

    params = specs.ODict()
    params[('NetStim1', 'pop')] = pops
    params[('NetStim1', 'sec')] = secs
    params[('NetStim1', 'loc')] = locs
    params[('NetStim1', 'weight')] = weights

    groupedParams = [('NetStim1', 'sec'), ('NetStim1', 'loc')] 

    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    initCfg['weightNorm'] = False
    initCfg['stimSubConn'] = False
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'synMech')] = ['AMPA','NMDA']
    initCfg[('NetStim1','synMechWeightFactor')] = [0.5,0.5]
    initCfg[('NetStim1', 'start')] = 700
    initCfg[('NetStim1', 'interval')] = 1000
    initCfg[('NetStim1','ynorm')] = [0.0, 1.0]

    initCfg[('NetStim1', 'noise')] = 0
    initCfg[('NetStim1', 'number')] = 1
    initCfg[('NetStim1', 'delay')] = 1
    #initCfg[('GroupNetStimW1', 'pop')] = 'None'
    initCfg[('NetStim1', 'delay')] = 1
    initCfg['addIClamp'] = 0
    
    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# EPSPs via NetStim
# ----------------------------------------------------------------------------------------------
def EPSPs():
    params = specs.ODict()

    params['groupWeight'] = [x*0.05 for x in np.arange(1, 8, 1)]
    params['ihGbar'] = [0.0, 1.0]
 
    
    # initial config
    initCfg = {}
    initCfg['duration'] = 0.5*1e3
    initCfg['addIClamp'] = False
    initCfg['addNetStim'] = True
    initCfg[('GroupNetStimW1', 'pop')] = 'PT5B'
    initCfg[('analysis','plotTraces','timeRange')] = [0, 500]
    initCfg['excTau2Factor'] = 2.0
    initCfg['weightNorm'] = True
    initCfg['stimSubConn'] = False
    initCfg['ihGbarZD'] = None

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve():
    params = specs.ODict()

    params[('IClamp1', 'pop')] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2']
    params[('IClamp1', 'amp')] = list(np.arange(0.0, 6.5, 0.5)/10.0) 
    #params['ihGbar'] = [0.0, 1.0, 2.0]
    # params['axonNa'] = [5, 6, 7, 8] 
    # params['gpas'] = [0.6, 0.65, 0.70, 0.75] 
    # params['epas'] = [1.0, 1.05] 
    # params['ihLkcBasal'] = [0.0, 0.01, 0.1, 0.5, 1.0] 

    # initial config
    initCfg = {}
    initCfg['duration'] = 1.5*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 500
    initCfg[('IClamp1','dur')] = 1000
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1500]

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b



# ----------------------------------------------------------------------------------------------
# v56_batch7 - Figure 2 (Control Quiet) (5*5 seeds)
#
# Note: v56_batch18 is the same but recording only 1 seed and more voltage traces
# ----------------------------------------------------------------------------------------------
def v56_batch7():

    params = specs.ODict()

    params['ihGbar'] = [1e-6, 0.25, 0.5, 0.75, 1.0] #, 0.25, 0.5, 0.75, 1.0]   #[0.25] # [0, 0.25] 
    
    params[('seeds', 'conn')] = [4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234+(17*i) for i in range(5)]
    
    groupedParams = []

    # initial config
    initCfg = {}
    duration = 5*1e3
    initCfg['duration'] = duration
    initCfg[('analysis', 'plotRaster', 'timeRange',0)] = 1000
    initCfg[('analysis', 'plotLFP', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotTraces', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotRaster', 'timeRange',1)] = duration
    initCfg[('analysis', 'plotLFP', 'timeRange', 1)] = duration
    initCfg[('analysis', 'plotTraces', 'timeRange', 1)] = duration

    initCfg['ihGbar'] = 0.75
    # initCfg[('modifyMechs', 'startTime')] = 1500
    # initCfg[('modifyMechs', 'endTime')] = 3000
    # initCfg[('modifyMechs', 'origFactor')] = initCfg['ihGbar']

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5
    
    # initCfg[('pulse', 'pop')] = 'TVL'
    # initCfg[('pulse', 'rate')] = [0, 10.0]
    # initCfg[('pulse', 'start')] = 1500.0
    # initCfg[('pulse', 'end')] = 3000.0
    # initCfg[('pulse', 'noise')] = 1.0

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [] #('IEweights',0), ('IIweights',0), ('IEweights',1), ('IIweights',1), ('IEweights',2), ('IIweights',2)]

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)
    b.batchLabel = 'v56_batch7'

    return b


# ----------------------------------------------------------------------------------------------
# v56_batch19 - Figure 3 (Control Quiet+Move)
#
# Note: v56_batch23 is the same but recording LFPs
# Note: v56_batch39 is the same but recording population LFPs and with recordStep=0.025
#
# ----------------------------------------------------------------------------------------------
def v56_batch19():
    params = specs.ODict()
    
    params[('modifyMechs', 'newFactor')] = [0.25]
    params[('pulse', 'rate', 1)] = [10]

    params[('seeds', 'conn')] = [4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234+(17*i) for i in range(5)]

    groupedParams = []

    # initial config
    initCfg = {}
    duration = 13*1e3
    initCfg['duration'] = duration
    initCfg[('analysis', 'plotRaster', 'timeRange',0)] = 1000
    initCfg[('analysis', 'plotLFP', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotTraces', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotRaster', 'timeRange',1)] = duration
    initCfg[('analysis', 'plotLFP', 'timeRange', 1)] = duration
    initCfg[('analysis', 'plotTraces', 'timeRange', 1)] = duration

    initCfg['ihGbar'] = 0.75
    initCfg[('modifyMechs', 'startTime')] = 5000
    initCfg[('modifyMechs', 'endTime')] = 9000
    initCfg[('modifyMechs', 'origFactor')] = initCfg['ihGbar']

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5
    
    initCfg[('pulse', 'pop')] = 'TVL'
    initCfg[('pulse', 'rate')] = [0, 10.0]
    initCfg[('pulse', 'start')] = 5000.0
    initCfg[('pulse', 'end')] = 9000.0
    initCfg[('pulse', 'noise')] = 1.0

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [] #('IEweights',0), ('IIweights',0), ('IEweights',1), ('IIweights',1), ('IEweights',2), ('IIweights',2)]

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)
    b.batchLabel = 'v56_batch19'

    return b


# ----------------------------------------------------------------------------------------------
# v56_batch20 - Figure 5 (MTh Inact Quiet+Move)
# ----------------------------------------------------------------------------------------------

def v56_batch20():
    params = specs.ODict()

    params[('ratesLong', 'TVL', 1)] = [0.01]  #[2.5] # [0.01, 2.5] #[2.5]  #
    params[('ratesLong', 'M2', 1)] = [0.01, 2.5]
    params[('modifyMechs', 'newFactor')] = [0.25] # , 0.5, 0.75, 1.0]  #0.5, 0.75, 1.0]

    params[('seeds', 'conn')] = [4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234+(17*i) for i in range(5)]

    groupedParams = []

    # initial config
    initCfg = {}
    duration = 13*1e3
    initCfg['duration'] = duration
    initCfg[('analysis', 'plotRaster', 'timeRange',0)] = 1000
    initCfg[('analysis', 'plotLFP', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotTraces', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotRaster', 'timeRange',1)] = duration
    initCfg[('analysis', 'plotLFP', 'timeRange', 1)] = duration
    initCfg[('analysis', 'plotTraces', 'timeRange', 1)] = duration

    initCfg['ihGbar'] = 0.75
    initCfg[('modifyMechs', 'startTime')] = 5000
    initCfg[('modifyMechs', 'endTime')] = 9000
    initCfg[('modifyMechs', 'origFactor')] = initCfg['ihGbar']

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 0.01
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5
    
    initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'rate')] = [0, 10.0]
    initCfg[('pulse', 'start')] = 5000.0
    initCfg[('pulse', 'end')] = 9000.0
    initCfg[('pulse', 'noise')] = 1.0

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)
    b.batchLabel = 'v56_batch20'
    
    return b


# ----------------------------------------------------------------------------------------------
# v56_batch22 - Figure 5 (MTh Inact Quiet+Move)
# ----------------------------------------------------------------------------------------------

def v56_batch22():
    params = specs.ODict()

    params[('seeds', 'conn')] = [4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234+(17*i) for i in range(5)]

    groupedParams = []

    # initial config
    initCfg = {}
    duration = 13*1e3
    initCfg['duration'] = duration
    initCfg[('analysis', 'plotRaster', 'timeRange',0)] = 1000
    initCfg[('analysis', 'plotLFP', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotTraces', 'timeRange', 0)] = 1000
    initCfg[('analysis', 'plotRaster', 'timeRange',1)] = duration
    initCfg[('analysis', 'plotLFP', 'timeRange', 1)] = duration
    initCfg[('analysis', 'plotTraces', 'timeRange', 1)] = duration

    initCfg['ihGbar'] = 1.0

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5
    
    initCfg[('pulse', 'pop')] = 'TVL'
    initCfg[('pulse', 'rate')] = [0, 10.0]
    initCfg[('pulse', 'start')] = 5000.0
    initCfg[('pulse', 'end')] = 9000.0
    initCfg[('pulse', 'noise')] = 1.0

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [] #('IEweights',0), ('IIweights',0), ('IEweights',1), ('IIweights',1), ('IEweights',2), ('IIweights',2)]

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)
    b.batchLabel = 'v56_batch22'

    return b



# ----------------------------------------------------------------------------------------------
# v56_batch5b - Figure 6 (VL vs Ih Quiet+Move)
# ----------------------------------------------------------------------------------------------
def v56_batch5b():
    params = specs.ODict()

    params[('modifyMechs', 'newFactor')] = [1e-6, 0.25, 0.5, 0.75, 1.0]  
    params[('pulse', 'rate', 1)] = [0.01, 2.5, 5, 10, 15, 20, 30] 

    groupedParams = []

    # initial config
    initCfg = {}
    duration = 5*1e3
    initCfg['duration'] = duration
    initCfg[('analysis', 'plotRaster', 'timeRange',0)] = 0
    #initCfg[('analysis', 'plotLFP', 'timeRange', 0)] = 0
    initCfg[('analysis', 'plotTraces', 'timeRange', 0)] = 0
    initCfg[('analysis', 'plotRaster', 'timeRange',1)] = duration
    #initCfg[('analysis', 'plotLFP', 'timeRange', 1)] = duration
    initCfg[('analysis', 'plotTraces', 'timeRange', 1)] = duration

    initCfg['ihGbar'] = 0.75
    initCfg[('modifyMechs', 'startTime')] = 1500
    initCfg[('modifyMechs', 'endTime')] = 3000
    initCfg[('modifyMechs', 'origFactor')] = initCfg['ihGbar']

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5
    
    initCfg[('pulse', 'pop')] = 'TVL'
    initCfg[('pulse', 'rate')] = [0, 10.0]
    initCfg[('pulse', 'start')] = 1500.0
    initCfg[('pulse', 'end')] = 3000.0
    initCfg[('pulse', 'noise')] = 1.0

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [] #('IEweights',0), ('IIweights',0), ('IEweights',1), ('IIweights',1), ('IEweights',2), ('IIweights',2)]

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin', nodes=1, coresPerNode=8):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 4,
            'script': 'init_cell.py',
            'mpiCommand': 'mpirun',
            'skip': True}

    elif type=='hpc_torque':
        b.runCfg = {'type': 'hpc_torque',
             'script': 'init.py',
             'nodes': 3,
             'ppn': 8,
             'walltime': "12:00:00",
             'queueName': 'longerq',
             'sleepInterval': 5,
             'skip': True}

    elif type=='hpc_slurm_comet':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'shs100', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            #'reservation': 'salva1',
            'walltime': '6:00:00',
            'nodes': 4,
            'coresPerNode': 24,  # comet=24, bridges=28
            'email': 'salvadordura@gmail.com',
            'folder': '/home/salvadord/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'ibrun', # comet='ibrun', bridges='mpirun'
            'skipCustom': '_raster.png'}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403', gcp='default'
            'walltime': '240:00:00', #'48:00:00',
            'nodes': nodes,
            'coresPerNode': coresPerNode,  # comet=24, bridges=28, gcp=32
            'email': 'salvadordura@gmail.com',
            'folder': '/home/ext_salvadordura_gmail_com/m1_paper2019_py3/sim/', # '/home/salvadord/m1/sim/',  # comet,gcp='/salvadord', bridges='/salvi82' 
            'script': 'init.py', 
            'mpiCommand': 'mpirun', # comet='ibrun', bridges,gcp='mpirun' 
            'skipCustom': '_raster.png',
            'sleepInterval': 4}
            #'custom': 'cd /home/ext_salvadordura_gmail_com/m1_paper2019_py3/sim/; cp ../data/tech/help_data.dat /usr/local/lib64/python3.6/site-packages/neuron/help_data.dat'}#,
            #'custom': '#SBATCH --partition="debug"'} # only use first 16 nodes (non-preemptible for long runs )
            # --nodelist=compute1


    elif type=='hpc_slurm_bridges':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'ib4iflp', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            'walltime': '06:00:00',
            'nodes': 2,
            'coresPerNode': 28,  # comet=24, bridges=28
            'email': 'salvadordura@gmail.com',
            'folder': '/home/salvi82/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'mpirun', # comet='ibrun', bridges='mpirun'
            'skip': True}


# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__': 

    # Figure 2 (Control Quiet) 
    # b = v56_batch7()
    
    # Figures 3, 4 (Control Quiet+Move)
    b = v56_batch19()

    # Figure 5 (MTh Inact Quiet+Move)
    # b = v56_batch20()

    # Figure 5 (NA-R block Quiet+Move)
    # b = v56_batch22()

    # Figure 6 (VL vs Ih Quiet+Move)
    # b = v56_batch5b()

    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'  # evol
    setRunCfg(b, 'hpc_slurm_gcp', nodes=1, coresPerNode=96)  # cores = nodes * 8 
    b.run() # run batch 