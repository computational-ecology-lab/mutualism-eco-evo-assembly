## ---------------------------
##
## Script name: Mut_Eigenval.py
##
## Purpose of script: This code loads the dataframes of simulations, calculates the eigenvalues of the community matrix (jacobian) at the final equilibrium, and generates the figures.
## after dataframes are produced by Mut_Dataframe, they can be loaded in this code to extract the distribution of eigenvalues of each network and plot them.
## Each dataframe should contain a collection of simulations (samples).
## In the results, this is used only for the main scenarios
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
from scipy import linalg
from scipy.optimize import fsolve
from scipy.integrate import odeint
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import autograd.numpy as np
from autograd import grad, jacobian

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

# Set the parameters for figure
fsize=14
params = { 'figure.dpi': 600,
        'legend.fontsize': fsize-3,
          #'figure.figsize': (3.54, 3.54),
         'axes.labelsize': fsize,
         'axes.titlesize':fsize,
         'xtick.labelsize':fsize-2,
         'ytick.labelsize':fsize-2}
mpl.rcParams.update(params)



def PSystem(yy,t):
# This is the system of equations to run the simulation (the same one integrated to generate the simulations in fevomodelM).
# This function is used here to simulate the system from the point of the last assembly event until it reaches a new equilibrium
    
    p_load = 1/(1 + hp*np.dot(PP,yy))      
    mp_load = 1/(1+hm*np.dot(MM, yy))


    dy = yy*np.dot(PP,yy)*p_load - yy*np.dot(pp,np.multiply(yy,p_load)) - yy*np.dot(CC,yy) + yy*np.dot(MM,yy)*mp_load + (r-cl)*yy  -yy*np.dot(ins,yy)

    return dy

def PTSystem(yy):
# This function is the same as above, but is used to calculate the Jacobian
# The difference is only that it does not require a time variable as input
    
    p_load = 1/(1 + hp*np.dot(PP,yy))      
    mp_load = 1/(1+hm*np.dot(MM, yy))


    dy = yy*np.dot(PP,yy)*p_load - yy*np.dot(pp,np.multiply(yy,p_load)) - yy*np.dot(CC,yy) + yy*np.dot(MM,yy)*mp_load + (r-cl)*yy  -yy*np.dot(ins,yy)

    return dy

def PDynamics(y0,time):
# Goal: run the entire simulation to reach a new equilibrium
# Inputs: initial vector of abundances and number of time-steps 
# Outputs: final densities

    # the PDynamics function wraps the integration mainly so it's possible to access the system at each time-step of PDynamics. The system is integrated for several small durations
    for iii in range(time):
        
        # Call odeint to integrate the system for each time-step (each time-step here is actually the sum of probetime/timestep true odeint timesteps)
        sol = odeint(PSystem, y0, np.linspace(0,probetime,int(probetime/timestep)))
        
        y0 = sol[-1].copy() # the new initial densities take from the last step of the previous integration
        
        # species reaching below the threshold are made extinct
        for k in range(y0.shape[0]):
            if(y0[k]<exth):
                y0[k]=0
    

    return y0.copy()


def GenerateSample(aa,yy,rk):
# Prepare the sample to run and then calculate the eigenvalues
# Inputs: matrix of interactions, densities vactor, list with parameters and variables from the model
# Outputs: New community size, new densities vector, new matrices of positive and negative interactions, new matrix of intraspecific competitions, new vector of intrinsic growths
    
    
    # Initialise global variables
    global ins
    global r
    global S
    global y0
    global pos
    global neg
    global p_load
    global mp_load
    global PP
    global MM
    global cl
    
    # get the total size of the community
    S_ = yy.shape[0]
    
    
    # get the loaded functions to their respective variables
    ins=rk[0][1]
    r=rk[0][0]
    S=np.array(r).shape[0]
    y0=rk[0][2]
    pos=rk[0][3]
    neg=rk[0][4]
        
    # get interaction matrices for each type of positive interactions, used later to calculate the total cost of positive interactions
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0)
    MM = np.where(aux != 0, pos, 0)
    
    # set a total number of time-steps for PDynamics
    maxtime=5000
    
    # calculate costs of positive interactions
    cl = cost * ( np.count_nonzero(PP,axis=1)+np.count_nonzero(MM,axis=1) )

    # run the model until maxtime
    sol = PDynamics(y0.reshape(S,),maxtime)
    
    
    # calculate variables for the new equilibrium
    xx = []
    idx = []
    for i in range(S_):
        if(sol[i]>exth):
            xx.append(sol[i])
            idx.append(i)
    xx = np.array(xx)
    S = xx.shape[0]
    
    print('NEW SIZE: ',xx.shape[0])
    
    poss = np.zeros((S,S))
    negg = np.zeros((S,S))
    inss = np.zeros((S,S))
    rr = np.zeros(S)
    
    k=0
    kk=0
    for i in idx:
        kk=0
        rr[k]=r[i]
        for j in idx:
            poss[k][kk]=pos[i][j]
            negg[k][kk]=neg[i][j]
            inss[k][kk]=ins[i][j]
            
            
            
            kk+=1
        
        k+=1
    
    
    return S, xx, poss, negg, inss, rr



def CalculateJacobian(SS,poss,negg,inss,xx,rr,normalize=False):
# Goal: calculate the eigenvalues to be plotted
# Inputs: community size, positive interactions, negative interactions, intraspecific competition, densities, growth rates, (normalize divides the eigenvalues by the square root of community size)
# Outputs: max eigenvalue, list of eigenvalues, list of roots    

    # define and load global variables and matrices
    global S
    S=SS
    global pos
    pos=poss
    global neg
    neg=negg
    global ins
    ins=inss
    global y0
    y0=xx
    global r
    r=rr
    global cl
    global PP
    global MM
    
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0)
    MM = np.where(aux != 0, pos, 0)
    
    global p_load
    global mp_load
    

    # calculate costs
    cl = cost * ( np.count_nonzero(PP,axis=1)+np.count_nonzero(MM,axis=1) )
    
    # calculate roots of the system at equilibrium
    global root
    root = fsolve(PTSystem,y0)
    for i in range(root.shape[0]):
        if(root[i]<exth and i!=0):
            root[i]=0

    # check if everything is right
    print('ROOT REDUCTION: ',str(root.shape[0]-np.count_nonzero(root)))
    # calculate the jacobian and the eigenvalues
    ja=jacobian(PTSystem)
    if(normalize==False):
        eigenvs=linalg.eigvals(ja(root))
    else:
        eigenvs=linalg.eigvals(ja(root)/np.sqrt(S))
    maxeigen=max(eigenvs)
    return maxeigen,eigenvs,root




######################## Variables to be used globally in the calculations. Same meaning as in fevomodelM

probetime = 0.1
timestep = 0.01


r=0

p_load=0
mp_load=0
PP=0
MM=0
pp=0
CC=0

cl = 0
cost = 0.01

exth = 0.000001 # extinction threshold
pos=0
neg=0
ins=0
S=0
y0=0
hm=0.1
hp=0.1

root=0


# Codes based on the kinds... these are for the main scenarios
filecodes = [[0,0],[1,0],[2,0],[3,0],[4,0]]


# Load the dataframes of the main scenarios
evo_history = []
rk_history = []
for i in filecodes:
    with open(dir_path+'/df/main_1/'+str(i[0])+'-df-'+str(i[1])+'.pkl', 'rb') as file:
          
        sollist=pickle.load(file)
        
        evo_history.append([])
        rk_history.append([])
    
        evo_history[-1].extend(sollist[1])
        rk_history[-1].extend(sollist[2])


# For each scenario, call the functions to calculate the eigenvalues and store them
results=[]
for i in [0,1,2,3,4]:
    result=[]
    ii=0
    for net in evo_history[i]:
        tt=-1
        print('SIZE: ',str(net[1][tt].shape[0]))

        SS, xx, poss, negg, inss, rr = GenerateSample(net[0][tt],net[1][tt],rk=rk_history[i][ii])

        
        maxv,eigv,rootv = CalculateJacobian(SS,poss,negg,inss,xx,rr,normalize=False)
        
        re=[]
        im=[]
        for ev in eigv:
            re.append(np.real(ev))
            im.append(np.imag(ev))
        result.append((re,im))
        ii+=1
    results.append(result)
    

# Names of main scenarios
netype = ['Evo.', 'Evo. No-M', 'Inv.', 'Inv. No-M', 'Inv. High-M.' ]


# Plot the distributions of eigenvalues
i=0
for resulti in results:
    labeling=True
    for res in resulti:
        if(labeling==True):
            plt.scatter(res[0],res[1],s=5,alpha=0.3,label=netype[i])
            labeling=False
        else:
            plt.scatter(res[0],res[1],s=5,alpha=0.3)

    plt.xlim(right=0.01)

    plt.xlabel('Re($\lambda$)')
    plt.ylabel('Im($\lambda$)')
    plt.title(netype[i])
    plt.show()
    i+=1
plt.ylim(-0.4,0.4)
plt.xlim(-1,0.01)
