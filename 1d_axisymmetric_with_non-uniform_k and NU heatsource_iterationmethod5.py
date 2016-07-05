import numpy as np
from matplotlib import pyplot 
import time
#import numba
start_time = time.time()
x_file = open('C:/Users/lenovo/Desktop/M.Tech Project/captelli.txt')
data1= np.loadtxt(x_file)
T_k=data1[:,0]
k_T=data1[:,1]
sig=data1[:,2]



#This function is used to find the Thermal Conductivity values for specific Temperature
#def T2K(t,T_k,k_T):
#    if t>=20000:
#        k = 0.0003269*(t-20000) + 3.225
#    elif t<0:
#        print 'Error: Temperature below zero'       
#    else:
#        for i in range(0,34):
#            if T_k[i]<=t and t<T_k[i+1]:
#                k =(k_T[i+1]*(t-T_k[i]) + k_T[i]*(T_k[i+1]-t))/(T_k[i+1] - T_k[i])
#    return k


def T2K(t,T_k,k_T):
    #k = 0.0002*t + 0.0457
    k=2e-25*t**6 - 2e-20*t**5 + 7e-16*t**4 - 1e-11*t**3 + 1e-07*t**2 + 6e-6*t    
    return k


#This function is used to find Sigma values for specific Temperature
#def T2S(t):
#    if t>=4000:
#        s=0.6637*(t-4000)
#    elif t>=0 and t<4000:
#        s=0
#    else:
#        print 'Error: Temperature below zero'
#    return s   

#def T2S(t):
#    if t>=4000:
#        s=0.6638*(t-4000)
#    else:
#        s=0
#    return s

def T2S(t):
    #s=0.6638*(t)
    s=-9e-22*t**6 + 9e-17*t**5 - 3e-12*t**4 + 5e-08*t**3 - 0.0002*t**2 + 0.3405*t
    return s

# Geometry of cylinder and grids. Initializing the problem. 
R1=0.0
R2=0.04
dz=0.25
nr=500
R=R2-R1

E=127
#q0=100000
Tb=1500 # This is the boundary temperature

#spacing
dr = R/(nr)

#Creating radius
r=np.linspace(R1+0.5*dr,R2-0.5*dr,nr)
r=np.append(R1,r)
r=np.append(r,R2)

#Since the area is not constant like cartesian this step is done
An=[]
As=[]

for i in range(0,nr):
    temp1=2*np.pi*dz*(R1+dr*i)
    As=np.append(As,temp1)
    temp2=2*np.pi*dz*(R1+dr*(i+1))
    An=np.append(An,temp2)


#As[0]=2*np.pi*dz*(R1+dr/4)
#
#An[nr-1]=2*np.pi*dz*(R1+dr*(nr-0.25))


## Numerical Solution
# Initial assumption of T
#T_i=Tb*np.ones(nr+1)
#T_i=np.append(T_i,Tb)

T_i=np.linspace(2*Tb,Tb,nr)
T_i=np.append(2*Tb,T_i)
T_i=np.append(T_i,Tb)

# Iterate for new estimate
T=T_i.copy()
err=[]
conv=[]
k_temp=np.zeros(nr+1)


for i in range(0,500):
    #print i
    T_temp=(T[1:]+T[:-1])/2
    #pyplot.figure(4)    
    ##plotting the values
    #pyplot.plot(T_temp, 'bx')
    #pyplot.show()
    #pyplot.pause(0.03)
    for j in range(0,nr+1):
        k_temp[j]=T2K(T_temp[j],T_k,k_T)
    kn=k_temp[1:]
    ks=k_temp[:-1]

    kn[nr-1]=T2K(Tb,T_k,k_T)
    ks[0]=T2K(T[0],T_k,k_T)
       
    #Defining constant for interior nodes
    aS=np.multiply(ks,As)/dr
    aN=np.multiply(kn,An)/dr
    
    sP=np.zeros(nr)
    su=np.zeros(nr)
    
    for j in range(nr):
        su[j]=2*np.pi*(R1+j*dr/2)*dr*dz*T2S(T[j+1])*E**2
        #su[j]=2*np.pi*(R1+j*dr/2)*dr*dz*q0
   
    #Defining constant for first exterior node 
    aS[0]=0
    
    
    #Defining constant for last exterior node 
    aN[nr-1]=0
    sP[nr-1]=-kn[nr-1]*An[nr-1]/(dr/2)
    su[nr-1]=kn[nr-1]*An[nr-1]*Tb/(dr/2) + su[nr-1]
    
    aP=aS+aN-sP
    #print aP
    
    #Creating Matrix
    A = np.zeros((nr,nr))     # pre-allocate [A] array
    C=su
    
    #Assigning values to the matrix
    A[0, 0] = aP[0]
    A[0, 1] = -aN[0]
    for j in range(1, nr-1):
        A[j, j-1] = -aS[j]  
        A[j, j] = aP[j]            
        A[j, j+1] = -aN[j]
        
    A[nr-1, nr-1] = aP[nr-1]
    A[nr-1, nr-2] = -aS[nr-1]

    T_new = np.linalg.solve(A,C)
    #T_new = 0.5*T_new + 0.5*T[1:-1]    # SUR iteration
    #T_new = 1.5*T_new - 0.5*T[1:-1]   # SOR iteration   
    T_new=np.append(T_new[0],T_new)
    T_new=np.append(T_new,Tb)
    
    #pyplot.figure(1)    
    #plotting the values
    #pyplot.plot(r, T_new, 'kx')
    #pyplot.xlabel('Radius of the Disc')
    #pyplot.ylabel('Temperature');
    #pyplot.legend( loc='upper right',fontsize=10,  fancybox=True )
    #pyplot.title('Variation of Temperature for 1-D Steady state heat conduction')
    #pyplot.grid(True)
    #pyplot.pause(0.03)
    #pyplot.show()
    #
    temp1=np.linalg.norm(T_new-T)/np.linalg.norm(T_new)
    conv=np.append(conv,temp1)
    
    T=T_new.copy()
    
    
pyplot.figure(1)    
pyplot.plot(r, T, 'kx:') 
pyplot.show()
    
print("--- %s seconds ---" % (time.time() - start_time))
  
pyplot.figure(2)
pyplot.plot(conv,'r.-',label='Convergence')
pyplot.yscale('log')
#pyplot.xscale('log')
pyplot.show()
pyplot.legend( loc='upper right',fontsize=10,  fancybox=True )
