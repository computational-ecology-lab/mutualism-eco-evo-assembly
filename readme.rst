This repository contains the code developed to implement the model and run numerical simulations to produce the results presented in the paper: 

**Mutualism provides a basis for biodiversity in eco-evolutionary community assembly**
by Gui Araujo & Miguel Lurgi
Published in [**PLoS Computational Biology**](https://doi.org/10.1371/journal.pcbi.1013402)
(2025), 21(9): e1013402

**#1. fevomodelM:**
Run this code to simulate the model and save the generated community.

**#2. Mut_Dataframe:**
Run this code to load communities and generate a dataframe with the histories of variables

**#3. Mut_Analysis:**
Run this code to load dataframes and produce figures

**#4. Mut_Eigenval:**
Run this code to load dataframes to calculate eigenvalues and produce the figures for distributions of eigenvalues


The idea is to generate a series of samples for each scenario using #1, separating the files by folders. Then using #2, read all samples of the same scenario into a single dataframe, for all types. Then using #3, read multiple dataframes to generate the figures comparing scenarios. Using #4, read dataframes to generate the figures of eigenvalue distributions.

The main 2 options to choose are the variables isTest and community_type. Choose isTest=1 to test custom parameters and run your own simulations.

If you choose isTest=0, then you will run pre-made scenarios to reproduce the results of the paper. There are 3 types of scenarios:

**community_type = 'main' :** runs the main scenarios in the paper, evolution and invasion models, with or without mutualisms, and also invasion with high proportions of mutualisms.

**community_type = 'onoff' :** runs the scenarios where mutualism is turned on or off during the simulation.

**community_type = 'mixed' :** runs the scenarios going from pure evolution to pure invasion, with mixed models ranging from 0.2, 0.5, and 0.8 proportions of invasion or evolution (as in the paper).

All the folders are ready to store the generated data, inside the folders 'data' and 'df'. For example, to run the main scenarios in the paper, just choose isTest=0 and community_type='main' in all scripts and run them sequentially, from #1 to #4. In this case, after running #1, the files inside data/mainData_1 will have the simulation data. After running #2, df/main_1 will have the dataframes. Then, running #3 and #4 will produce the figures.

To change any variable, follow the comments explaining what they are. Be sure to change the affected variables in other scripts. For example, to change the number of assembly events, change the variable n_evolutions in fevomodelM. Then, you will have to change accordingly the tvec variable in the other scripts.

This material includes the df's used in the paper for all scenarios: main, onoff,and mixed. Therefore, it is possible to run only scripts #3 and #4 to reproduce the same analyses used in the paper. If you run #1 + #2 for new simulations of these scenarios, the original df's will be overwritten. (these df's are larger files not included in the github repository, only in zenodo:  https://doi.org/10.5281/zenodo.15439137)

When downloading from github, first include the contents from **data_folders.zip** into the same folder containing the scripts.

The file **spec-file.txt** contains all the information about the conda environment that hosted this project.


This project was supported by the Leverhulme Trust through Research Project Grant **\# RPG-2022-114**.

