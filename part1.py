from scipy.sparse import diags
from scipy.sparse import linalg as splinalg
import scipy.sparse as sparse
import matplotlib.pyplot as plt
import numpy as np

#  00D#MMXX#

tol = float(1e-6)
one = int(1)
mone = int(-1)

L = float(1.0)
Re = float(1000)    # Reynolds number 
N = int(16)  		# mesh cells in x- and y-direction

u = np.zeros([2*N*(N+1),1], dtype = np.float)
p = np.zeros([N*N+4*N,1], dtype = np.float)
tx = np.zeros([N+1,1], dtype = np.float)     # grid points on primal grid
x = np.zeros([N+2,1], dtype = np.float)      # grid points on dual grid
th = np.zeros([N], dtype = np.float)       # mesh width primal grid
h = np.zeros([N+1], dtype = np.float)      # mesh width dual grid 

#Generation of a non-uniform grid
x[0] = 0
x[N+1] = 1
for i in range(N+1):
    xi = i*L/N
    tx[i] = 0.5*(1. - np.cos(np.pi*xi))     #  tx mesh point for primal mesh
    if i>0:
        th[i-1] = tx[i] - tx[i-1]           # th mesh width on primal mesh
        x[i] = 0.5*(tx[i-1]+tx[i])          # x mesh points for dual mesh
        
for i in range(N+1):
    h[i] = x[i+1]-x[i]                      # h mesh width on dual mesh
    
        
th_min = min(th)
h_min = min(h)
h_min = min(h_min,th_min)                   # determination smallest mesh size

dt = min(h_min,0.5*Re*h_min**2)             # dt for stable integration

dt = 5.*dt

x_array = np.zeros(len(x))
tx_array = np.zeros(len(tx))
plt.figure(num=1)
plt.scatter(x,x_array,s=1,c='r')
plt.scatter(tx,tx_array,s=1,c='b')
plt.legend(['Dual','Primal'])
plt.savefig('1.png')

