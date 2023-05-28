<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>

<!-- ABOUT THE PROJECT -->
# About The Project

This code is used to generate the images that comprise the Figure 7 of the paper Dura-Bernal et al 2023 "Multiscale model of primary motor cortex circuits predicts in vivo cell type-specific, behavioral state-dependent dynamics" (Cell Reports).
It is developed in Python using the [NetPyNE](http://www.netpyne.org/) package and the [NEURON](https://www.neuron.yale.edu/neuron/) simulator
<br/><br/>

# Getting Started

The code provided includes all the scripts and the processed Spike Histogram dataset used to make the figures. 

The processed dataset is provided for file size reasons, but if you wish to download the original files and re-build the analysis, you should add the required files in the `./data` folder and run the code using these flags set to `True` in the `runAnalysis.py` file: 
* `updateConnectivity = True`
* `updateHistDict = True`

The `updateConnectivity` flag will build connectivity `.pkl` files for each cell in the folder `./data/conn_info`, named as `PopName_CellGID` (e.g. `PT5B_5132.pkl`). It contains the individual connections for each cell, which are loaded to buid the Spike Histograms.

The `updateHistDict` will generate the Spike Histogram for each population, which is already provided in the repository as the processed dataset.

The data will be saved in the `./figs/histogram_figures/windowAnalysis/<POP_NAME>/hist_window_<WINDOW_SIZE>_ms/9_post_analysis` folder.
We also provide the code for other analysis and figures developed during the study, which were omitted in the publication for simplicity.
<br/><br/>

# Prerequisites

The code requires the [NetPyNE](http://www.netpyne.org/) package and the [NEURON](https://www.neuron.yale.edu/neuron/) simulator, along with common data processing and plotting packages, such as [matplotlib](https://matplotlib.org/) and [numpy](https://numpy.org/).

A full list of required packages can be found in the headed or the `runAnalysis.py` and `AnalyzeData.py` files.

It also requires [scikit-learn](https://scikit-learn.org/stable/) and [umap](https://umap-learn.readthedocs.io/en/latest/) for multivariate analysis.
<br/><br/>

# Installation

The project requires no installation other than the necessary packages.
<br/><br/>

<!-- CONTACT -->
# Contact

Joao Moreira - [ResearchGate](https://www.researchgate.net/profile/Joao-Moreira) - [LinkedIn](https://www.linkedin.com/in/joaovvitor/) - joao.moreira@downstate.edu

Salvador Dura-Bernal - salvador.dura-bernal@downstate.edu

Project Link: [M1_NetPyNE_CellReports_2023](https://github.com/suny-downstate-medical-center/M1_NetPyNE_CellReports_2023)
<br/><br/>

<p align="right">(<a href="#readme-top">back to top</a>)</p>
