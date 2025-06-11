## ---------------------------
##
## Script name: Mut_Dataframe.py
##
## Purpose of script: This code takes the models' direct outputs and generates dataframes for
## each batch of simulations, naming them with a 'kind' number.
## These dataframes contain several variables calculated throughout the history of the simulation.
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
import pandas as pd
import pickle

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

def RebuildMatrix(snaps):
# Goal: rebuild interaction matrices from the compact form stored by the simulation code
# Inputs: an array with encoded interactions at a given time
# Outputs: the matrix of interactions
    new_snaps = []
    for snap in snaps:
        rm = np.zeros((int(snap[0][0]),int(snap[0][0])))

        for el in snap[1:]:
            rm[int(el[0])][int(el[1])]=el[2]

        new_snaps.append(rm)
    return new_snaps






def NetworkSC(G):
# Goal: calculate connectance and number of species
# Inputs: network (object from networkx)
# Outputs: connectance and number of nodes

    u_G = G.to_undirected()

    if(u_G.number_of_nodes()!=0):
        Cr = 2*u_G.number_of_edges()/u_G.number_of_nodes()**2
    else:
        Cr=0
    n_nodes = G.number_of_nodes()

    
    return Cr, n_nodes


def GenerateNetworks(aa):
# Goal: generate the newtorkx network
# Inputs: matrix with all interactions (both positive and negative)
# Outputs: the network object from networkx    

    aa = np.absolute(aa) # it doesn't care about interaction type

    RN = nx.from_numpy_array(aa)
    
    RN.remove_nodes_from(list(nx.isolates(RN))) # makes sure there are no isolated nodes left

    return RN


def InteractionTypes(a,species):
# Goal: counts the numbers of each interaction types and the total density
# Inputs: matrix of interactions, vector of abundances
# Outputs: # mutualisms, # consumer-resource, # competition, # total, total density
    


    # This code separates positive and negative interactions,
    # then separates each type and counts them
    pos = np.where(a<0, 0, a)
    neg = np.absolute(np.where(a>0, 0, a))
    
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0) # consumer
    MM = np.where(aux != 0, pos, 0) # mutualism
    
    bux = np.zeros_like(neg)
    bux[neg != 0] = 1
    bux = np.logical_and(bux, bux.T)
    pp = np.where(bux == 0, neg, 0) # resource
    CC = np.where(bux != 0, neg, 0) # competition
    
    # counts
    countM = np.count_nonzero(MM)/2 # each two entries are 1 mutualism interaction
    countP = np.count_nonzero(PP) # each consumer (or resource) entry is 1 consumer-resource interaction
    countC = np.count_nonzero(CC)/2 # each two entries are 1 competition interaction
    countT = countM+countP+countC
    

    
    total_density = np.sum(species)
    
    
           
    return countM,countP,countC,countT, total_density






def CalculateHistory(net,tvec):
# Goal: call InteractionTypes for all timesteps in the series and store the results
# Inputs: list with interaction matrix and vector of abundances, list of times
# Outputs: lists with values for all times for: original times, mutualism counts, consumer-resource counts, competition counts, total counts, total densities

    m_snaps,y_snaps,_ = net
    times = np.arange(len(y_snaps))
    
    
    M = []
    P = []
    C = []
    T = []

    total_d = []
    
    for t in range(len(tvec)):
        
        

        countM,countP,countC,countT, total_density = InteractionTypes(m_snaps[t],y_snaps[t])

        
        
        M.append(countM)
        P.append(countP)
        C.append(countC)
        T.append(countT)
        
        total_d.append(total_density)
            
        
    return times, M, P, C, T, total_d




def SaveDF():
    
    df = pd.DataFrame() # initialise the df
    evo_history=[] # this will store history of variables
    evo_historyX=[] # this stores the history to be saved together with the df
    rk_list=[] # this will store the list of additional variables defining the model
    tvec = [1,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,999] # chosen times for the series, at every 50 assembly events

    frange = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] # range of codes for each sample to be included in a dataframe
    
    if(isTest==1):
        frange=np.array([1,2,3,4,5,6,7,8,9,10])

    for j in frange: # These will run for all sample codes, the same as the num_sample in fevomodelM
    
        path = dir_path+'/data/'+folder+str(j)+'FEVO-O1.pkl'
        if(isTest==1):
            path = dir_path+'/'+folder+str(j)+'FEVO-O1.pkl'
        with open(path, 'rb') as file:
          

            sollist=[]
            sollist.append(pickle.load(file)) # receive the loaded file from the model
            
        for i in range(len(sollist[0])):

            # Storing variables
            numbers=sollist[0][i][1]
            rk_list.append(sollist[0][i][2])
            mat = RebuildMatrix(numbers[2]) # the compacted matrix is rebuilt
            numbers[2]=mat
            
            
            
              

       
            # Calculate the history of variables for all times in tvec
            evo_history.append([[numbers[2][ii] for ii in tvec],[numbers[3][ii] for ii in tvec],[]])
            times, M, P, C, T, total_d = CalculateHistory(evo_history[-1], tvec)
            
            
            dic = {} # Initialise the sample dic that'll entry the df

            
            
       
        
            # This block will compose the dic variables for each time in tvec
            # The final dic will have all variables for all chosen times, like a timeseries
            k=0
            for ss in range(len(tvec)):
                
                # calculate connectance and richness
                RN = GenerateNetworks(evo_history[-1][0][ss])
                Cr, n_nodes = NetworkSC(RN)
                
                
                dic['kind']=kind # store kind of simulations


                # store total density, richness, and connectance
                dic['tot_density'+str(tvec[ss])]=total_d[k]
                dic['n_species'+str(tvec[ss])]=n_nodes
                try:
                    dic['Con'+str(tvec[ss])] = 2*RN.number_of_edges()/n_nodes**2
                except:
                    dic['Con'+str(tvec[ss])] = 0
                
                
                 
                
                # store all proportions of interaction types
                jj=0
                for nn in [M,P,C,T]:

                    name = ['M','P','C','T']
                    
                    
                    dic[name[jj]+str(tvec[ss])] = nn[k]
                    
                    if(jj<3):
                        try:
                            dic['prop'+name[jj]+str(tvec[ss])] = nn[k]/T[k]
                        except:
                            dic['prop'+name[jj]+str(tvec[ss])] = 0
            
                    jj+=1
                
                
                
                k+=1
            
            

            
            # include the dic into the df
            df = pd.concat([df,pd.DataFrame([dic])], ignore_index = True)
            
            print(j) # log the processing of a sample
            
            
            evo_historyX.append(evo_history[-1])

       

                               
            
    # save the df with additional model info
    path = dir_path+'/df/'+saveFolder+str(kind)+'-df-0.pkl'
    if(isTest==1):
        path = dir_path+'/df/'+saveFolder+str(kind)+'-df-0.pkl'
    with open(path, 'wb') as file:
          
        pickle.dump([df,evo_historyX,rk_list], file)










################################################## CHANGE THESE VARIABLES TO RUN THE PRE-MADE SCENARIOS

isTest=1 # change this to 1 to run a custom test simulation, keep at 0 to run the pre-made simulations
community_type = 'main' # this control the type of simulation... the paper features 3 types: main, onoff, mixed

# TEST DATAFRAMES
if(isTest==1): # Choose these parameters according to the custom simulations you saved
    kind = 10
    folder = 'data/'
    saveFolder = '/'
    SaveDF()
    

# THESE ARE FOR DATAFRAMES OF PRE-MADE SIMULATIONS, EXACTLY AS THE PAPER
if(isTest==0):

    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE MAIN SCENARIOS        
    if(community_type == 'main'):
        
        modelType = ['evo','evo_noM','inv','inv_noM','inv_highM']
        
        for j in [0,1,2,3,4]:
            
            kind = j # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
            folder = 'mainData_1/main_'+modelType[j]+'_1/' # specific folder where the simulations are stored

            saveFolder='main_1/'

            print(kind,folder)
            
            SaveDF()
            
    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE ONOFF SCENARIOS
    if(community_type == 'onoff'):
        
        modelType = ['mut_intro','mut_stop']
        
        for j in [0,1]:
            
            kind = j # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
            folder = 'onoffData_1/onoff_'+modelType[j]+'_1/' # specific folder where the simulations are stored

            saveFolder='onoff_1/'

            print(kind,folder)
            
            SaveDF()
            
    # THIS BLOCK CHOOSES THE OPTIONS THAT RUN THE MIXED SCENARIOS
    if(community_type == 'mixed'):
        
        modelType = ['evo','mixed_02','mixed_05','mixed_08','inv']
        
        for j in [0,1,2,3,4]:
            
            kind = j # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
            folder = 'mixedData_1/mixed_'+modelType[j]+'_1/' # specific folder where the simulations are stored

            saveFolder='mixed_1/'

            print(kind,folder)
            
            SaveDF()



