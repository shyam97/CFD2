# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 12:25:58 2020

@author: mgerritsma
"""


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
N = int(2)  		# mesh cells in x- and y-direction

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

#
#  Note that the time step is a bit conservative so it may be useful to see
#  if the time step can be slightly increased. This will speed up the
#  calculation.
#

#  Boundary conditions for the lid driven acvity test case
U_wall_top = -1
U_wall_bot = 0
U_wall_left = 0
U_wall_right = 0
V_wall_top = 0
V_wall_bot = 0
V_wall_left = 0
V_wall_right = 0


# Set up the sparse incidence matrix tE21. Use the orientations described
# in the assignment.
# Make sure to sue sparse matrices to avoid memory problems
tE21 = sparse.lil_matrix((N*N,2*N*(N+1)))
for y in range(N):
    for x in range(N):
        tE21[x+y*N,N*(N+1)+x+y*N]=1
        tE21[x+y*N,N*(N+1)+N+x+y*N]=-1
        tE21[x+y*N,x+y*(N+1)]=1
        tE21[x+y*N,x+y*(N+1)+1]=-1

#  Insert the normal boundary conditions and split of the vector u_norm
indices = []
indices.append([i*(N+1) for i in range(N)])
indices.append([(i+1)*(N+1)-1 for i in range(N)])
indices.append([N*(N+1)+i for i in range(N)])
indices.append([2*N*(N+1)-(i+1) for i in range(N)])
indices = np.sort(np.array(indices).flatten())
complement = np.array([i for i in range(2*N*(N+1)) if not i in indices])

tE21_norm = sparse.csc_matrix(tE21.toarray()[:,indices])
tE21 = sparse.csc_matrix(tE21.toarray()[:,complement])

u_norm = np.zeros(tE21_norm.shape[1])
for i in range(N):
    u_norm[i*2] = U_wall_left
    u_norm[i*2-1] = U_wall_right
    u_norm[2*N+i] = V_wall_bot
    u_norm[3*N+i] = V_wall_top

print(tE21_norm.shape)
print(u_norm.shape)
u_norm = tE21_norm * u_norm
print(u_norm.shape)
# Set up the outer-oriented incidence matrix tE10


#  Set up the sparse, inner-oriented  incidence matrix E10


#  Set up the (extended) sparse, inner-oriented incidence matrix E21


#  Split off the prescribed tangential velocity and store this in
#  the vector u_pres


#  Set up the Hodge matrices Ht11 and H1t1


#  Set up the Hodge matrix Ht02

'''
A = tE21@Ht11@E10
#print("A")
#print(A.toarray())
#input("Press Enter to continue...")

n = A.shape[0]
LU = splinalg.splu(A,diag_pivot_thresh=0) # sparse LU decomposition

u_pres_vort = Ht02@u_pres
temp = H1t1@tE10@Ht02@u_pres

u_pres = temp

VLaplace = H1t1@tE10@Ht02@E21
DIV = tE21@Ht11

ux_xi = np.zeros([(N+1)*(N+1),1], dtype = float)
uy_xi = np.zeros([(N+1)*(N+1),1], dtype = float)
convective = np.zeros([2*N*(N+1),1], dtype = float)

while (diff>tol):

    xi = Ht02@E21@u + u_pres_vort

    for i in range(N+1):
        for j in range(N+1):
            k = j + i*(N+1);
            if j==0:
                ux_xi[k] = U_wall_bot*xi[i+j*(N+1)]
                uy_xi[k] = V_wall_left*xi[j+i*(N+1)]
            elif j==N:
                ux_xi[k] = U_wall_top*xi[i+j*(N+1)]
                uy_xi[k] = V_wall_right*xi[j+i*(N+1)]
            else:
                ux_xi[k] = (u[i+j*(N+1)]+u[i+(j-1)*(N+1)])*xi[i+j*(N+1)]/(2.*h[i])                       # Klopt
                uy_xi[k] = (u[N*(N+1)+j+i*N] + u[N*(N+1)+j-1+i*N])*xi[j+i*(N+1)]/(2.*h[i])

    for  i in range(N):
        for j in range(N+1):
            convective[j+i*(N+1)] = -(uy_xi[j+i*(N+1)]+uy_xi[j+(i+1)*(N+1)])*h[j]/2.
            convective[N*(N+1)+i+j*N] = (ux_xi[j+i*(N+1)]+ux_xi[j+(i+1)*(N+1)])*h[j]/2.

    # Set up the right hand side for the equation for the pressure

    rhs_Pois = DIV@( u/dt - convective - VLaplace@u/Re - u_pres/Re) + u_norm/dt

    # Solve for the pressure

    p = LU.solve(rhs_Pois)

    # Store the velocity from the previous time level in the vector uold

    uold = u

    # Update the velocity field

    u = u - dt*(convective + E10@p + (VLaplace@u)/Re + u_pres/Re)

    # Every other 1000 iterations check whether you approach steady state and
    # check whether you satsify conservation of mass. The largest rate at whci
    # mass is created ot destroyed is denoted my 'maxdiv'. This number should
    # be close to machine precision.

    if ((iter%1000)==0):
        maxdiv = max(abs(DIV@u+u_norm))
        diff = max(abs(u-uold))/dt

        print("maxdiv : ",maxdiv)
        print("diff   : ", diff)


    iter += 1
'''
