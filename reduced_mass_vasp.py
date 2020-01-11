'''
Please input the "filepath".
This program can give the reduced mass of the system.
By mhf.
Thanks very much for Mr Ren^_^
20191125

'''
#=========================================================================

#==============================================================================
# Module imports

import numpy as np

#=========================================================================
# some constants
AUTOA = 0.529177249 # Bohr to angstrom
PI = 3.141592653589793238
hbar = 1.054571628E-34
AMUTOAU = 1822.8889 # AMU to au
C0 = 299792458. # light speemd
WAVNUM2AU= 100.0*C0/6.57969E15
EVTOJ=1.60217733E-19
RYTOEV=13.605826
TIMEAU = hbar / (2.0 * RYTOEV * EVTOJ)
FREQAU = 1.0 / TIMEAU
VASP2GAUSSIAN=0.036749308136649*AUTOA*AUTOA
#==============================================================================
# filepath=("D:\mhfResut\OUTCAR-result\SiO2\SiO2vdWB88\Sio2-20190515vdwB88fre\OUTCAR")
filepath=("D:\mhfResut\C3N2H5-fre-result\OUTCAR")
# D:\mhfResut\C3N2H5-fre-result
def normalize(v):
    N = 1.0 / np.sqrt(np.dot(v,v))
    return v * N, N


def proj(u,v):
    return np.dot(v,u) / np.dot(u,u) * u

def mode(u):
    return np.sqrt(np.dot(u,u))

def gram(V):
    M,N = V.shape
    U = np.zeros([M,N])
    E = np.zeros([M,N])
    U[0] = V[0]
    E[0] = U[0] / mode(U[0])
    for k in range(1,M):
        tmp = np.zeros(N)
        for j in range(k):
            tmp += proj(U[j],V[k])
        U[k] = V[k] - tmp
        E[k] = U[k] / mode(U[k])
    return E





def normalize(v):
    N = 1.0 / np.sqrt(np.dot(v,v))
    return v * N, N


def getDOF(filepath):
    """
    Parse OUTCAR, get DOF: degrees of freedom.
    """
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith('   Degrees of freedom'):
                break

        dof = int(l.split()[-1])
        return dof

def getIBRION(param):

    pattern = "   {}".format(param)
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith(pattern):
                break
        value = l.split()[2]
    return value


def getNormalMode(imode):

    # if (int(IBRION) != 7) or (int(IBRION) != 8):
    #     raise ValueError("I need a VASP calculation with IBRION=7! exit..")

    pattern = '{:4d} f  ='.format(imode)  #
    pattern_i = '{:4d} f/i='.format(imode)  # imaginary frequency
    img_mode = False
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith(pattern):
                dump = l.split()
                Omega = dump[3:]
                break
            if l.startswith(pattern_i):
                img_mode = True
                dump = l.split()
                Omega = dump[2:]
                break
        # dump = l.split()
        # Omega = dump[3:]  # 4x2 datasets, \nu(THz), \omega(THz), wavenumber(cm-1), and energy(meV)
        for i in range(4):
            Omega[2 * i] = float(Omega[2 * i])

        l = fh.readline()
        mode = []
        for ii in range(NIons):
            l = fh.readline()
            mode.append([float(x) for x in l.split()[3:]])

    mode = np.array(mode)
    Mode = {}
    Mode['Omega'] = Omega
    Mode['mode'] = mode
    return Mode


def get_mass_weighted_NormalMode(NIons):
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if 'SQRT(mass)' in l:
                break
        nmodes = NIons * 3
        freqs = np.zeros([nmodes, 4])
        modes = np.zeros([nmodes, NIons, 3])
        for i in range(5):
            l = fh.readline()
        for imode in range(nmodes):
            l = fh.readline()
            sign = -1 if '/i' in l else 1
            freqs[imode] = [float(x) for x in l[9:].split()[:9:2]]
            freqs[imode] *= sign
            fh.readline()
            for i in range(NIons):
                l = fh.readline()
                modes[imode, i, :] = [float(x) for x in l.split()[3:]]
            fh.readline()
    return ( modes)




def getIonsPerType():
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith('   ions per type'):
                break
        IonsPerType = [int(x) for x in l.split()[4:]]
    return IonsPerType



def getIonMasses( ):
    M = []
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith("   POMASS"):
                this_m = float(l.split()[2][:-1])
                M.append(this_m)
                if len(M) == NTypes:
                    break
    masses = []
    for it in range(NTypes):
        masses.extend([M[it]]*IonsPerType[it])
    return np.array(masses)

def get_hessian_matrix():
    with open(filepath) as fh:
        while True:
            l = fh.readline()
            if l.startswith(' SECOND DERIVATIVES'):
                break

        hessian_matrix = []
        for ii in range(imode+2):
            l = fh.readline()
            hessian_matrix.append(l.split())

    hessian_matrix=np.array(hessian_matrix[2:])
    hessian_matrix=hessian_matrix[:,1:]
    return hessian_matrix


def getCoordinate( NIons):
    """
    Parse OUTCAR, get the equilibrium coordinates.
    Args:
        None
    Returns:
        coord: ndarray (Nions x 3), xyz coordinates
    Raises:
        None
    """
    coord = []
    with open(filepath, 'r') as fh:
        while True:
            l = fh.readline()
            if l.startswith(' position of ions in cartesian'):
                break
        for i in range(NIons):
            l = fh.readline()
            coord.append([float(x) for x in l.split()])
    return np.array(coord)

#=========================================================================

#=========================================================================
# get Number of atoms

imode=getDOF(filepath)
IonsPerType = getIonsPerType()
NTypes = len(IonsPerType)
NIons = sum(IonsPerType)
coord = getCoordinate(NIons )
N = NIons
# nmodes = 3 * N - 6
nmodes = 3 * N
masses=getIonMasses()

coords = coord
# coords *= AUTOA
print ("Number of atoms: "), N
# if linear:
#     nmodes = 3 * N - 5
# else:
#     nmodes = 3 * N - 6
# print ("Number of vibration modes: "), nmodes

# nmodes = 3 * N - 6
# get atomic masses
atm_mass = masses
# atm_mass in atomic mass unit
atm_mass = atm_mass*AMUTOAU
M = np.zeros([3*N,1])
for i in range(3*N):
    M[i] = atm_mass[i//3]


# calculate center of mass
sum_mass = np.sum(atm_mass)
sum_x = np.sum(atm_mass * coords[:,0])
sum_y = np.sum(atm_mass * coords[:,1])
sum_z = np.sum(atm_mass * coords[:,2])
R_com = np.array([sum_x, sum_y, sum_z]) / sum_mass
print ("center of mass:"), R_com
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

# Ntr = len(new_D)
# Nvib = 3 * N - Ntr
Nvib = 3 * N

D = np.random.random([3*N,3*N])
# for i in range(Ntr):
#     D[i] = np.array(new_D[i])

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

D = D[:].transpose()

# now D is the transformation matrix with dimension
# 3N x Nvib




# read Hessian Matrix from gs_fname (full matrix)
hessian_matrix=VASP2GAUSSIAN*get_hessian_matrix().astype(float)

hessian_matrix=(hessian_matrix+hessian_matrix.transpose())/2
f_car=hessian_matrix
# f_car = get_hessian(gs_fname, N)
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
        wave_num[i] = np.sqrt(-w2[i])
    else:
        wave_num[i] = -np.sqrt(w2[i])
wave_num = wave_num * FREQAU * 0.01 / C0 / (2*PI)
print (wave_num)
# raise SystemExit
#print wave_num

# calculate Cartesian displacements
M = np.zeros([3*N, 3*N])
for i in range(3*N):
    M[i,i] = 1.0 / np.sqrt(atm_mass[int(i/3)])

l_cart = np.dot(M, np.dot(D,L))
print (l_cart.shape)
print (l_cart[:,0])
Nvib=3*N
mu = np.zeros(Nvib)
L = np.zeros([3*N, Nvib])
for i in range(Nvib):
    L[:,i], mu[i] = normalize(l_cart[:,i])
mu = mu * mu / AMUTOAU
# print (L[:,0])

np.savetxt("reduced_mass_of_vasp.txt",(wave_num, mu))
