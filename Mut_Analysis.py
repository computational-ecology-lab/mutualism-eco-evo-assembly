## ---------------------------
##
## Script name: Mut_Analysis.py
##
## Purpose of script: This code loads the dataframes of simulations and creates the figures from the paper (apart from the eigenvalue analysis).
## after dataframes are produced by Mut_Dataframe, they can be loaded into this code to be processed for figures.
## Each dataframe contains a collection of simulations (samples).
##
## Author: gui araujo
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
##
## python 3.11
##
## Date Created: November-2024
##
## Email: gui.dav.araujo@gmail.com
##
## ---------------------------


# Packages used
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import math
import seaborn as sns
import community
import matplotlib as mpl

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

# Set parameters for figures
fsize=14
params = { 'figure.dpi': 600,
        'legend.fontsize': fsize-3,
          #'figure.figsize': (3.54, 3.54),
         'axes.labelsize': fsize,
         'axes.titlesize':fsize,
         'xtick.labelsize':fsize-2,
         'ytick.labelsize':fsize-2}
mpl.rcParams.update(params)




def Modularity(Gr):
# Goal: calculate unweighted modularity of network
# Inputs: network (object from networkx)
# Outputs: modularity value

    G = Gr.to_undirected()
    part = community.best_partition(G)
    mod = community.modularity(part, G)
    return mod


def Entropy(G,degrees):
# Goal: calculate degree entropy of network
# Inputs: network (object from networkx), list of degrees
# Outputs: entropy value
    
    vk = degrees
    try:
        maxk = np.max(vk)
    
        #mink = np.min(vk)
        kvalues= np.arange(0,maxk+1) # possible values of k
        Pk = np.zeros(maxk+1) # P(k)
        for k in vk:
            Pk[k] = Pk[k] + 1
        Pk = Pk/sum(Pk) # the sum of the elements of P(k) must to be equal to one
        
        H = 0
        for p in Pk:
            if(p > 0):
                H = H - p*math.log(p, 2)
    except:
        H=0
    
    return H


def NetworkMetrics(G,degrees):
# Goal: calculate connectance, number of nodes (richness), and entropy (using the Entropy function)
# Inputs: network (object from networkx), list of degrees
# Outputs: values of connectance, number of nodes, and entropy

    u_G = G.to_undirected()
    
    if(u_G.number_of_nodes()!=0):
        Cr = 2*u_G.number_of_edges()/u_G.number_of_nodes()**2 # C = L / (S^2/2)
    else:
        Cr=0
    n_nodes = G.number_of_nodes()
    

    Hr = Entropy(u_G,degrees)

    
    return Cr, n_nodes, Hr


def GenerateNetworks(aa):
# Goal: generate the newtorkx network
# Inputs: matrix with all interactions (both positive and negative)
# Outputs: the network object from networkx    

    aa = np.absolute(aa) # it doesn't care about interaction type

    RN = nx.from_numpy_array(aa)
    
    RN.remove_nodes_from(list(nx.isolates(RN))) # makes sure there are no isolated nodes left

    return RN









tvec = [1,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,999] # times processed (every 50 assembly events), same as in Mut_Analysis


################################################## CHANGE THESE VARIABLES TO RUN THE PRE-MADE SCENARIOS
isTest=0 # change this to 1 to run a custom test simulation, keep at 0 to run the pre-made simulations
community_type = 'main' # this control the type of simulation... the paper features 3 types: main, onoff, mixed
num_samples=10 # number of samples, to calculate the standard error


# TO USE WITH TEST DATAFRAMES
if(isTest==1):
    filecodes = [[0,0]]
    #dicKind = {0:'evo_$\sigma$=0.4',1:'evo_$\sigma$=0.6',2:'evo_$\sigma$=0.8',3:'inv_$\sigma$=0.4',4:'inv_$\sigma$=0.6',5:'inv_$\sigma$=0.8'}
    #dicKind = {0:'evo_1/25',1:'evo_1/100',2:'evo_1/1000',3:'inv_1/25',4:'inv_1/100',5:'inv_1/1000'}
    #dicKind = {0:'evo_$\Delta$=15',1:'evo_$\Delta$=25',2:'evo_$\Delta$=35'}
    dicKind = {0:'evo_abund.'}
    #colors=['#377eb8','#ff7f00','#4daf4a','#f781bf','#a65628','#984ea3']
    #colors=['#377eb8','#ff7f00','#4daf4a']
    colors=['#377eb8']
    folder='abundance_df/'

# THESE ARE FOR DATAFRAMES OF PRE-MADE SIMULATIONS, EXACTLY AS THE PAPER
if(isTest==0):
    
    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE MAIN SCENARIOS 
    if(community_type == 'main'):
        filecodes = [[0,0],[1,0],[2,0],[3,0],[4,0]] # this is to load the dataframes based on their kind numbers
        # Set the names of networks based on kinds and their respective colors
        dicKind = {0:'Evo.', 1:'Evo. No-M', 2:'Inv.', 3:'Inv. No-M.', 4:'Inv. High-M.'} # to determine the name of scenarios, based on their kind number
        colors=['#377eb8','#ff7f00','#4daf4a','#f781bf','#a65628'] # to detemrine the colours of scenarios
        folder='main_1/'
    
    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE MIXED SCENARIOS 
    if(community_type == 'mixed'):
        filecodes = [[0,0],[1,0],[2,0],[3,0],[4,0]]
        dicKind = {0:'Evo.', 1:'Mixed 0.2', 2:'Mixed 0.5', 3:'Mixed 0.8', 4:'Inv.'}
        colors=['#377eb8','#17a6e3','#2cc1d1','#44c2ad','#53ad8b']
        folder='mixed_1/'
        
    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE MIXED SCENARIOS 
    if(community_type == 'onoff'):
        filecodes = [[0,0],[1,0]]
        dicKind = {0:'Inv. - Mut. Introduced', 1:'Inv. - Mut. Stopped' }
        colors=['black','gray']
        folder='onoff_1/'
    

# Read the dataframes listed in filecodes
df = pd.DataFrame()
evo_history = []
for i in filecodes:
    path=dir_path+'/df/'+folder+str(i[0])+'-df-'+str(i[1])+'.pkl'
    if(isTest==1):
        path=dir_path+'/SI_tests/df/'+folder+str(i[0])+'-df-'+str(i[1])+'.pkl'
    with open(path, 'rb') as file:
          
        sollist=pickle.load(file)
    
        df = pd.concat([df,sollist[0]], ignore_index = True)
        evo_history.extend(sollist[1])


# This generates the graphics of proportions of interaction types
# The codes don't work properly if there's just one sample in a dataframe,
# because they need to calculate variances
############################# PROPORTIONS ##################################
#########################################################################
df_types = []
for k in df['kind'].unique():
    df_types.append(df[df['kind']==k].copy())


for name in ['propM','propP','propC']:
    parl = []
    parls=[]
    times = tvec
    for k in df['kind'].unique():
        parl.append([])
        parls.append([])
    for i in times:
        par = name+str(i)

        for k in range(len(df_types)):
            parl[k].append(df_types[k][par].mean())
            parls[k].append(df_types[k][par].std())



    for k in range(len(df_types)):
        lab = dicKind
        plt.plot(times,parl[k],'-o',c=colors[k],label=lab[k])
        plt.errorbar(times[1:],np.array(parl[k])[1:],yerr=np.array(parls[k])[1:]/np.sqrt(num_samples),fmt='o',c=colors[k],alpha=0.3)


    
    plt.xlabel('Assembly event')
    if(name=='propM'):
        if(community_type=='main' or community_type=='onoff'):
            plt.legend(loc=(0.65,0.1))
        plt.ylabel('Proportion of mutualism')
    if(name=='propP'):
        plt.ylabel('Proportion of consumer-resource')
    if(name=='propC'):
        plt.ylabel('Proportion of competition')
    plt.show()
########################################################################


# This is the same as above, but for the other variables
############################# OTHER VARIABLES ##################################
#########################################################################
df_types = []
for k in df['kind'].unique():
    df_types.append(df[df['kind']==k].copy())

for name in ['n_species','Con','r_dens','CS','r_P']: # Richness, connectance, average density, complexity, consumer-resource per species
    parl = []
    parls=[]
    times = tvec
    for k in df['kind'].unique():
        parl.append([])
        parls.append([])
    for i in times:
        if(name=='r_dens'):
            par2= 'n_species'+str(i)
            par1= 'tot_density'+str(i)

            for k in range(len(df_types)):
                parl[k].append((df_types[k][par1]/df_types[k][par2]).mean())
                parls[k].append((df_types[k][par1]/df_types[k][par2]).std())
        
        elif(name=='CS'):
            par2= 'n_species'+str(i)
            par1= 'Con'+str(i)

            for k in range(len(df_types)):
                parl[k].append((df_types[k][par1]*df_types[k][par2]).mean())
                parls[k].append((df_types[k][par1]*df_types[k][par2]).std())
        
        elif(name=='r_P'):
            par2= 'n_species'+str(i)
            par1= 'P'+str(i)

            for k in range(len(df_types)):
                parl[k].append((df_types[k][par1]/df_types[k][par2]).mean())
                parls[k].append((df_types[k][par1]/df_types[k][par2]).std())
                
                
        else:    
            par = name+str(i)
        

            for k in range(len(df_types)):
                parl[k].append(df_types[k][par].mean())
                parls[k].append(df_types[k][par].std())



    for k in range(len(df_types)):
        lab = dicKind
        plt.plot(times,parl[k],'-o',c=colors[k],label=lab[k])
        plt.errorbar(times[1:],np.array(parl[k])[1:],yerr=np.array(parls[k])[1:]/np.sqrt(num_samples),fmt='o',c=colors[k],alpha=0.3)


    if(name=='n_species'):
        plt.ylabel('Species richness (S)')
        if(community_type=='mixed' or community_type=='main'):
            plt.legend()
    if(name=='Con'):
        plt.ylabel('Connectance (C)')
    if(name=='r_dens'):
        plt.ylabel('Average abundance')
        plt.yscale('log')
    if(name=='CS'):
        plt.ylabel('Complexity (S*C)')
    if(name=='r_P'):
        plt.ylabel('Consumer-resource per species')
        if(community_type=='main'):
            plt.legend()
    plt.xlabel('Assembly event')
    #plt.legend()
    plt.show()

    
for kind_idx in df['kind'].unique():
    for k in df[df['kind']==kind_idx].index: 
        if(k==df[df['kind']==kind_idx].index[-1]):
            plt.scatter(df.loc[k,'n_species'+str(tvec[-1])],df.loc[k,'Con'+str(tvec[-1])],color=colors[df.loc[k,'kind']],label=dicKind[df.loc[k,'kind']])
        else:
            plt.scatter(df.loc[k,'n_species'+str(tvec[-1])],df.loc[k,'Con'+str(tvec[-1])],color=colors[df.loc[k,'kind']])
plt.xlabel('Species richness')
plt.ylabel('Connectance')
if(community_type=='mixed'):
    plt.legend()
plt.show()







# This is to calculate values of entropy and modularity, then generate the graphics
# Used only for the main scenarios
######### DEGREE ENTROPY AND MODULARITY
######################################################################################



# Create networks for each sample
nets=[]
for k in (df.index):


    RNm = GenerateNetworks(evo_history[k][0][-1])

    nets.append(RNm)


# Define a sublist for each kind of sample
SEntropy=[]
Modul=[]
for i in df['kind'].unique():
    SEntropy.append([])
    Modul.append([])
    
    

for i in df['kind'].unique():
    for k in df.index:
        if(i==df.loc[k,'kind']):
            print('Start:',i,k)
            
            net = nets[k]
            
            degrees = np.array([val for (node, val) in net.degree()])
            Cr, n_nodes, H = NetworkMetrics(net,degrees)
            
            
            mod = Modularity(net)

            
            # This part calculates the values of the random network with same C and S, then their average, then subtracts it from the original value
            randomList = []
            for z in range(50):
                randomList.append(nx.erdos_renyi_graph(n_nodes,Cr))

            r_sentropy=[]
            r_modul=[]
            for z in range(len(randomList)):
                rdegrees = np.array([val for (node, val) in randomList[z].degree()])
                rCr, rn_nodes, rH = NetworkMetrics(randomList[z],rdegrees)
                r_sentropy.append(rH)
                r_modul.append(Modularity(randomList[z]))
            
            SEntropy[i].append(H-np.array(r_sentropy).mean())
            Modul[i].append(mod-np.array(r_modul).mean())
        




# Generate entropy graphic
metric = 'Entropy increase' 
mdata=SEntropy.copy()
      
rows_list=[]
for kind in range(len(mdata)):
    for i in range(len(mdata[kind])):
        rows_list.append({metric:mdata[kind][i],'Scenario':dicKind[kind]})
mdatadf = pd.DataFrame(rows_list)
palette = sns.color_palette(colors)
if(isTest==0):
    f = plt.figure(figsize=[4,4])
if(isTest==1):
    f = plt.figure(figsize=[4,4])
ax = f.add_subplot(111)
sns.boxplot(x='Scenario',y=metric,data=mdatadf,palette=palette,ax=ax, width=0.3,linewidth=1)
if(community_type=='main'):
    if(isTest==0):
        plt.xticks([0,1,2,3,4],['Evo.','Evo.\nNo-M','Inv.','Inv.\nNo-M','Inv.\nHigh-M'])
if(isTest==1):
    #plt.xticks([0,1,2,3,4,5],['evo\n$\sigma$=0.4','evo\n$\sigma$=0.6','evo\n$\sigma$=0.8','inv\n$\sigma$=0.4','inv\n$\sigma$=0.6','inv\n$\sigma$=0.8'])
    #plt.xticks([0,1,2,3,4,5],['evo\n1/25','evo\n1/100','evo\n1/1000','inv\n1/25','inv\n1/100','inv\n1/1000'])
    plt.xticks([0],['evo_abund.'])
    #plt.xticks([0,1,2],['evo\n$\Delta$=15','evo\n$\Delta$=25','evo\n$\Delta$=35'])
f.tight_layout()




# Generate modularity graphic
metric = 'Modularity increase'
mdata=Modul.copy() 

rows_list=[]
for kind in range(len(mdata)):
    for i in range(len(mdata[kind])):
        rows_list.append({metric:mdata[kind][i],'Scenario':dicKind[kind]})
mdatadf = pd.DataFrame(rows_list)
palette = sns.color_palette(colors)
if(isTest==0):
    f = plt.figure(figsize=[4,4])
if(isTest==1):
    f = plt.figure(figsize=[4,4])
ax = f.add_subplot(111)
sns.boxplot(x='Scenario',y=metric,data=mdatadf,palette=palette,ax=ax, width=0.3,linewidth=1)
if(community_type=='main'):
    if(isTest==0):
        plt.xticks([0,1,2,3,4],['Evo.','Evo.\nNo-M','Inv.','Inv.\nNo-M','Inv.\nHigh-M'])
if(isTest==1):
    #plt.xticks([0,1,2,3,4,5],['evo\n$\sigma$=0.4','evo\n$\sigma$=0.6','evo\n$\sigma$=0.8','inv\n$\sigma$=0.4','inv\n$\sigma$=0.6','inv\n$\sigma$=0.8'])
    #plt.xticks([0,1,2,3,4,5],['evo\n1/25','evo\n1/100','evo\n1/1000','inv\n1/25','inv\n1/100','inv\n1/1000'])
    plt.xticks([0],['evo_abund.'])
    #plt.xticks([0,1,2],['evo\n$\Delta$=15','evo\n$\Delta$=25','evo\n$\Delta$=35'])
f.tight_layout()
