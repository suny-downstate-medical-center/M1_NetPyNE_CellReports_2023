
import os
import json
import scipy
from matplotlib import pyplot as plt

# --- Snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
import sys
from pydoc import source_synopsis
from matplotlib  import  pyplot  as plt
from numpy import size

print("Matplotlib backend (default): %s" %plt.get_backend())
modules = []
for module in sys.modules:
    if module.startswith('matplotlib'):
        modules.append(module)
for module in modules:
    sys.modules.pop(module)
import matplotlib
matplotlib.use("MacOSX")
from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (dynamic): %s" %plt.get_backend())

folderPath = '/Users/joao/Research/Models/M1Salva_paper_analysis/figs/2023_02_20/windowAnalysis'
# allFilesList = os.listdir(folderPath)

allFilesList = ['PT5B_spth_2_bar_cellType_H_6000_9000ms_allPops_enhanced_cells.json']

filesList_JSON=[ fileName for fileName in allFilesList if fileName.endswith('.json')]
for fileName in filesList_JSON:
    saveData_dict={}
    fObj = open(folderPath+'/'+fileName)
    fileData =json.load(fObj)
    saveData_dict.update({fileName:{}})
    try:
        State_1='Quiet'
        State_2='Movement'
        State_1_data = fileData[State_1]
        State_2_data = fileData[State_2]
    except:
        State_1='Cluster1'
        State_2='Cluster2'
        State_1_data = fileData[State_1]
        State_2_data = fileData[State_2]
    for mech_type in fileData[State_1].keys():
        for pre_pop in fileData[State_1][mech_type].keys():
            # print('=',mech_type,pre_pop)
            if len(fileData[State_1][mech_type][pre_pop]['vals'])>0:
                plt.figure()
                print('>',mech_type,pre_pop,'\n')
                Data1=State_1_data[mech_type][pre_pop]['vals']
                Data2=State_2_data[mech_type][pre_pop]['vals']
                # print(State_1,len(Data1))
                # print(State_2,len(Data2),'\n')
                # stats_test = scipy.stats.mannwhitneyu(State_1_data,State_2_data) 
                stats_test = scipy.stats.mannwhitneyu(Data1, Data2) 
                plt.hist(Data1,histtype=u'step',bins=100)
                plt.hist(Data2,histtype=u'step',bins=100)
                plt.title(pre_pop+'_'+mech_type)
                # print(stats_test,'\n\n')
                print(pre_pop,'_',mech_type,'\tpvalue: ',stats_test.pvalue,' | statistic: ',stats_test.statistic)
                saveData_dict[fileName].update({(pre_pop+'_'+mech_type):{'pvalue':stats_test.pvalue,'statistic':stats_test.statistic}})
            # else:
            #     print('skipping ',mech_type,pre_pop,'\n')
    json_object = json.dumps(saveData_dict, indent=4)
    with open(folderPath+'/statistics_'+fileName, "w") as outfile:
        outfile.write(json_object)

plt.show()