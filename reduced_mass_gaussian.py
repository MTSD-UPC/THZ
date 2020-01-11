
#==============================================================================
# Module imports
import argparse
import numpy as np
from constants import *
from parse import *
# import scitools. as scinp
#==============================================================================
def normalize(v):
    N = 1.0 / np.sqrt(np.dot(v,v))
    return v * N, N

# parse arguments for input file names
parser = argparse.ArgumentParser()
parser.add_argument(
    "-gf", "--freqfile",
    default = "freq.log",
    help="G09 ground state vibration analysis logfile"
)
parser.add_argument(
    "-ef","--excitedfile", 
    default = "force.log",
    help="G09 excited state gradients logfile",
)
parser.add_argument(
    "-o", "--output",
    default = "disp.dat",
    help = "output file contains dimensionless displacements"
)

parser.add_argument(
    "-l","--linear",
    help="linear molecule (geometry)",
    action = "store_true"
)

args = parser.parse_args()

gs_fname = args.freqfile
ex_fname = args.excitedfile
ofname = args.output
linear = args.linear

print ("= using ground state vibration logfile %s"% gs_fname)
print ("= using excited state gradient logfile %s" % ex_fname)
print ("= linear molecule: ", linear)
# print
#=========================================================================

#=========================================================================
# get Number of atoms 
N = get_num_atoms(gs_fname)
coords = get_coordinates(gs_fname, N)
coords *= AUTOA
print ("Number of atoms: ", N)
if linear:
    nmodes = 3 * N - 5
else:
    nmodes = 3 * N - 6
print ("Number of vibration modes: ", nmodes)


# get atomic masses
atm_mass = get_atomic_masses(gs_fname, N)
# atm_mass in atomic mass unit
atm_mass *= AMUTOAU
M = np.zeros([3*N,1])
for i in range(3*N):
    M[i] = atm_mass[int(i/3)]


# calculate center of mass
sum_mass = np.sum(atm_mass)
sum_x = np.sum(atm_mass * coords[:,0])
sum_y = np.sum(atm_mass * coords[:,1])
sum_z = np.sum(atm_mass * coords[:,2])
R_com = np.array([sum_x, sum_y, sum_z]) / sum_mass
print ("center of mass:", R_com)
r_com = np.zeros([N,3])
for i in range(N):
    r_com[i] = coords[i] - R_com

# calculate moment of inertia tensor
I = np.zeros([3,3])
for i in range(3):
    for j in range(i,3):
        if i == j:
            p = (i+1) % 3
            q = (i+2) % 3
            I[i,j] = np.sum(atm_mass * (r_com[:,p]*r_com[:,p] + r_com[:,q]*r_com[:,q]))
        else:
            I[i,j] = -np.sum(atm_mass * (r_com[:,i]*r_com[:,j]))
    I[j,i] = I[i,j]
#print I

I_p, X = np.linalg.eig(I)

# generate transformation D
D = np.zeros([6,3*N])
mass_sqrt = np.sqrt(atm_mass)
for i in range(N):
    D[0,i*3] = mass_sqrt[i]
    D[1,i*3+1] = mass_sqrt[i]
    D[2,i*3+2] = mass_sqrt[i]

for i in range(N):
    for j in range(3):
        Pi = np.dot(X,r_com[i].reshape(3,-1))
        D[3,i*3+j] = (Pi[1]*X[j,2] - Pi[2]*X[j,1]) * mass_sqrt[i]
        D[4,i*3+j] = (Pi[2]*X[j,0] - Pi[0]*X[j,2]) * mass_sqrt[i]
        D[5,i*3+j] = (Pi[0]*X[j,1] - Pi[1]*X[j,0]) * mass_sqrt[i]

# normalize D vectors
new_D = []
small_tol = 1.0E-6
for i in range(6):
    summ = np.dot(D[i],D[i])
    if summ < small_tol:
        continue
    new_D.append(D[i] / np.sqrt(summ))

Ntr = len(new_D)
Nvib = 3 * N - Ntr

D = np.random.random([3*N,3*N])
for i in range(Ntr):
    D[i] = np.array(new_D[i])

# D = scinp.Gram_Schmidt(D,normalize=True,remove_null_vectors=False,
#                       remove_noise=True)

D = gram(D)
for i in range(3*N):
    for j in range(i+1):
        summ = np.dot(D[i],D[j])
        if i == j:
            if abs(summ - 1.0) > 1.0E-6:
                print ("%d-th row %f") % (i, summ)
        else:
            if abs(summ) > 1.0E-6:
                print ("%d-%d: %f") % (i, j, summ)

D = D[Ntr:].transpose()
# now D is the transformation matrix with dimension
# 3N x Nvib




# read Hessian Matrix from gs_fname (full matrix)
f_car = get_hessian(gs_fname, N)
f_mw = f_car / np.sqrt(np.outer(M,M.transpose()))
# f_mw = f_car / np.sqrt(np.outer(M,M.transpose()))
#w, L = np.linalg.eig(f_mw)
#idx = w.argsort()
#w = w[idx]
#L = L.T[idx].T
w2,L2 = np.linalg.eig(f_mw)
idx = w2.argsort()
w2 = w2[idx]
L2 = L2.T[idx].T

# calculate frequencies
f_int = np.dot(D.transpose(),np.dot(f_mw,D))
w,L = np.linalg.eig(f_int)
idx = w.argsort()
w = w[idx]
L = L.T[idx].T

# using "-" sign for imaginary frequencies
wave_num = np.zeros(Nvib)
for i in range(Nvib):
    if w[i] < 0:
        wave_num[i] = -np.sqrt(-w[i])
    else:
        wave_num[i] = np.sqrt(w[i])
wave_num = wave_num * FREQAU * 0.01 / C0 / (2*PI)
# print (wave_num)
# raise SystemExit
#print wave_num

# calculate Cartesian displacements
M = np.zeros([3*N, 3*N])
for i in range(3*N):
    M[i,i] = 1.0 / np.sqrt(atm_mass[int(i/3)])

l_cart = np.dot(M, np.dot(D,L))
# print (l_cart.shape)
# print (l_cart[:,0])
mu = np.zeros(Nvib)
L = np.zeros([3*N, Nvib])
for i in range(Nvib):
    L[:,i], mu[i] = normalize(l_cart[:,i])
mu = mu * mu / AMUTOAU

np.savetxt("reduced_mass_of_gaussian.txt",(wave_num, mu))