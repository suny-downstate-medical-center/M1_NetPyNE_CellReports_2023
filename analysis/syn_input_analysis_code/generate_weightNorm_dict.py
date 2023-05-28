#! /Users/joao/opt/anaconda3/bin/python

from fileinput import filename
import pickle

from netpyne import specs
netParams = specs.NetParams()

# --- Cell params to be loaded
loadCellParamLabels = ['PT5B_full', 'CT6_reduced', 'IT5B_reduced', 'IT6_reduced'] #  # list of cell rules to load from file

for cell_label in loadCellParamLabels:
    cell_params_fileName = cell_label+'_cellParams'
    cell_params_filePath='../cells/'+cell_params_fileName+'.pkl'
    netParams.loadCellParams(label=cell_label, fileName=cell_params_filePath)

# --- Generate a dictionary with values stored as dict[cellLabel][sec][weightNorm]
weightNorm_dict={}
for cellLabel in loadCellParamLabels:
    weightNorm_dict.update({cellLabel:{}})
    for sec in netParams.cellParams[cellLabel]['secs'].keys():
        wN = netParams.cellParams[cellLabel]['secs'][sec]['weightNorm'][0]
        print(cellLabel,sec,wN)
        weightNorm_dict[cellLabel].update({sec:wN})


# --- Filename
weightNorm_dict_filename = '../weightNorm/weightNorm_dict.pkl'

# --- Saving the weightNorm in a pkl file
with open(weightNorm_dict_filename, 'wb') as f:
    pickle.dump(weightNorm_dict, f)


