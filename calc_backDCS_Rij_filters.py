from scipy import optimize
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import readstruct as struct
import ADK
from os.path import exists
import os
import multiprocessing as mp
from functools import partial


nm = 18.897261328856432  # nm in a.u.
I0 = 3.509470*10**16  # W/cm^2 in a.u. ############double check this
c = 137.035999139  # speed of light in a.u.
eV = 0.03674932539796232  # eV in a.u.
atto = 0.0413413737896  # as in a.u.
fs = atto*1000.
angs = nm/10.    # Angstorm in a.u.

# ----------------------------------------------------------------------------
#             Input
# ----------------------------------------------------------------------------
ncore = 10
nrg = np.arange(20, 110, 2)*eV
nE = np.size(nrg)
k_r = np.sqrt(2.*nrg)
wavelen = 3030.*nm
omega = 2.*np.pi*c/wavelen
TS_case = 'eq'
# TS_case = 'ff'
# TS_case = 'Ey'
Target = 'camphor'
# Target = 'butane_anti'
# Target = 'butane_gauche'
dir = Target + '/'
if __name__ == '__main__':
    print('Target is', Target)
    if not exists(dir):
        os.makedirs(dir)


# ----------------------------------------------------------------------------
#             Interpolate the birth time
# ----------------------------------------------------------------------------
long_traj = np.loadtxt('birth_return_long.dat')
ion_func = interp1d(long_traj[:,0], long_traj[:,-1], kind='cubic')
def find_field(Enrg, U):
    func_time = lambda t: ion_func(t) - Enrg / U
    sol = optimize.root_scalar(func_time, bracket=[long_traj[0,0],long_traj[-1,0]],method='brentq')
    F = F0*np.cos(sol.root/180.*np.pi)
    return F

# Inten_max = 0.65
# # 3.17*Up = 40 eV --> calculate the minimum intensity.
# Inten_min = (nrg[0]/long_traj[-1,-1])*4.*omega**2*I0*1e-14
# # Add a small number to avoid error in the function find_field.
# Inten_min = np.round(Inten_min, decimals=2) + 0.01
# if __name__ == '__main__':
#     print('Inten_min=', Inten_min)
#     print('Inten_max=', Inten_max)
# dI = 0.05
# I_list = np.arange(Inten_min, Inten_max + dI, dI)

# Operate at constant intensity 1.7 E 13 , wavelength 4000 nm, Up 25.397
I = 0.41  # intensity in units of 1e14 W/cm^2 

##################double check this is correct intensity??

# ----------------------------------------------------------------------------
#             Read the Cartesian coordinate and Atom names from files
# ----------------------------------------------------------------------------

if TS_case == 'eq':
    inp = Target + '_cat_eq.xyz'
elif TS_case == 'ff':
    inp = Target + '_10fs_0E_cart.dat'
elif TS_case == 'Ey':
    inp = Target + '_10fs_Ey_cart.dat'


if __name__ == '__main__': print('Reading coordinates from', inp)

with open(inp) as file:
    target_cart = []
    atom = []
    for i, line in enumerate(file):
        if i < 2:
            continue
        cart = np.zeros(3)
        tmp = line.split()
        if len(tmp) < 4:
            continue
        atom.append(tmp[0])
        cart[0:3] = np.array(tmp[1:4], dtype=float)
        cart *= angs
        target_cart.append(cart)

# ----------------------------------------------------------------------------
#             Read the orbitals from files
# ----------------------------------------------------------------------------
para = struct.structure(Target)
Norb = len(para.Ip)
if __name__ == '__main__': print('Number of orbitals are', Norb)
# ----------------------------------------------------------------------------
#             Read the DCS from files at a range of energy
# ----------------------------------------------------------------------------
fold_H = '../../data/H/'
fold_C = '../../data/C/'
cm = 1e7*nm
# get theta grid
f_H = np.zeros( np.size(nrg), dtype=np.complex128 )
f_C = np.zeros( np.size(nrg), dtype=np.complex128 )
for i_E in range(50, np.size(nrg)):
    tmp_data = np.loadtxt(fold_H + 'E' + str(np.int32(np.rint(nrg[i_E]/eV)))  + '/scatamp.dat')
    tmp_data2 = np.loadtxt(fold_C + 'E' + str(np.int32(np.rint(nrg[i_E]/eV)))  + '/scatamp.dat')
    # Only take 180 degrees
    f_H[i_E] = (tmp_data[-1,2] + 1.j*tmp_data[-1,3])*cm
    f_C[i_E] = (tmp_data2[-1,2] + 1.j*tmp_data2[-1,3])*cm
I_A = 4.*np.abs(f_C)**2 + 10.*np.abs(f_H)**2

# ----------------------------------------------------------------------------
#             Rotation of the position vectors (z-y-z)
# ----------------------------------------------------------------------------
da = 0.02*np.pi
db = da
dc = da
alpha = np.arange(0., 2.*np.pi + da, da)
beta = np.arange(0., np.pi + db, db)
gamma = np.arange(0., 2.*np.pi + dc, dc)
if __name__ == '__main__':
    print('Size of alpha, beta, gamma', np.size(alpha), np.size(beta), np.size(gamma))

Na = np.size(alpha)
Nb = np.size(beta)
Nc = np.size(gamma)

# ----------------------------------------------------------------------------
#            Wrapping the ADK_rate function
# ----------------------------------------------------------------------------
def ADK_rate(F, m, Ip, Zeff, B):
    model = ADK.ADK(field=F, m=m, Ip=Ip, Zeff=Zeff, B=B)
    ADK_rate = model.ADK_rate()
    return ADK_rate

def R_z(alpha):
    R_z = np.zeros(([3, 3]))
    R_z[0, 0] = np.cos(alpha)
    R_z[0, 1] = -np.sin(alpha)
    R_z[1, 0] = np.sin(alpha)
    R_z[1, 1] = np.cos(alpha)
    R_z[2, 2] = 1.
    return R_z


def R_y(beta):
    R_y = np.zeros(([3, 3]))
    R_y[0, 0] = np.cos(beta)
    R_y[2, 2] = np.cos(beta)
    R_y[0, 2] = np.sin(beta)
    R_y[2, 0] = -np.sin(beta)
    R_y[1, 1] = 1.
    return R_y

def index_to_angs(index):
    ic = index % Nc
    tmp_ind = index // Nc
    ib = tmp_ind % Nb
    ia = tmp_ind // Nb
    return ia, ib, ic


def TI_IAM(index, i_E_max, target_w):
    ia, ib, ic = index_to_angs(index)
    c = gamma[ic]
    b = beta[ib]
    a = alpha[ia]
    target_cart_rot = []
    target_cart_rot = [row.copy() for row in target_cart]

    for i in range(0, len(atom)):
        target_cart_rot[i][:] = R_z(a) @ R_y(b) @ R_z(c) @ target_cart_rot[i][:]

    # Three cases based on Rij threshold
    thresholds = [4, 6, 8]
    results = {}
    
    for threshold in thresholds:
        MIT_target_tmp = np.zeros(i_E_max)
        for i in range(0, len(atom)):
            if atom[i] == 'C':
                f1 = f_C.copy()
            else:
                f1 = f_H.copy()
            for j in range(i+1, len(atom)):
                if atom[j] == 'C':
                    f2 = f_C.copy()
                else:
                    f2 = f_H.copy()
                Rij = target_cart_rot[i][:] - target_cart_rot[j][:]
                Rij_mag = np.linalg.norm(Rij)
                if Rij_mag < threshold:  # Filter based on threshold
                    q_dot = 2.*k_r[:i_E_max]*Rij[2]

                    # !!!!!!!!! The Jacobian for sin theta is done here!!!!!!!!!
                    MIT_target_tmp += 2. * (f1[:i_E_max]*np.conjugate(f2[:i_E_max])*np.exp(1.j*q_dot)).real * np.sin(b)* target_w[ib, ic]
        results[f'lt{threshold}'] = MIT_target_tmp
    
    return results

def iso_IAM(i_E_max):
    # ----------------------------------------------------------------------------
    #             Calculate the isotropic MIT from IAM
    # ----------------------------------------------------------------------------
    # Three cases based on R12 threshold
    thresholds = [4, 6, 8]
    
    for R12_threshold in thresholds:
        MIT_target_iso = np.zeros(nE)
        q_iso = 2.*k_r[:i_E_max]
        for i in range(0, len(atom)):
            if atom[i] == 'C':
                f1 = f_C.copy()
            else:
                f1 = f_H.copy()
            for j in range(i+1, len(atom)):
                R12 = np.linalg.norm(target_cart[i][:] - target_cart[j][:]) #computes distance between atoms
                if R12 < R12_threshold:  # Filter based on threshold
                    if atom[j] == 'C':
                        f2 = f_C.copy()
                    else:
                        f2 = f_H.copy()

                    MIT_target_iso[:i_E_max] += 2.*(f1[:i_E_max]*np.conjugate(f2[:i_E_max])).real * \
                        np.sin(q_iso*R12)/(q_iso*R12)
        data = np.zeros((nE,2))
        data[:, 0] = nrg/eV
        data[:i_E_max, 1] = MIT_target_iso[:i_E_max] #just MIT
        data[:i_E_max, 2] = MIT_target_iso[:i_E_max] + I_A[:i_E_max]
        np.savetxt(dir + Target + '_IAM_Rij_lt' + str(R12_threshold) + '_' + str(np.round(I,decimals=2)) +'.dat', data)
    return


def TI_rate(index, F0):
    rate = 0.
    i_c = index % Nc
    i_b = index // Nc
    for orb in range(0, Norb):
        B_m = para.calc_B(orb, beta[i_b], gamma[i_c])
        l_max = int((np.size(B_m)-1)/2)
        m_l = np.arange(-l_max, l_max+1, 1)
        Ip = para.Ip[orb]                
        wfat = para.calc_wfat(orb, beta[i_b], gamma[i_c])
    # sum over m for B_m
        for im in m_l:
            rate += ADK_rate(F0, im, Ip, 1., B_m[im+l_max])*wfat
    return rate

# ----------------------------------------------------------------------------
#            Run for single intensity
# ----------------------------------------------------------------------------
F0 = np.sqrt(I*1e14/I0)
Up = (F0/(2.*omega))**2

if __name__ == '__main__':
    print('Intensity is', I, 'e14 W/cm^2.')
    print('Up is (eV)', Up/eV)
    # ----------------------------------------------------------------------------
    #             Find the ionization time for a range of long trajectories
    # ----------------------------------------------------------------------------
    Field = np.zeros(nE)
    for i_E in range(0, nE):        
        i_E_max = i_E + 1
        if nrg[i_E]/Up > long_traj[-1,-1]:
            if __name__ == '__main__':
                print(i_E, nrg[i_E]/Up)
                print('E/Up > 3.17, breaking the loop for finding the classical ionization time.')
            break
        else:
            Field[i_E] = find_field(nrg[i_E], Up)
    # Field is 0 if the laser intensity cannot support higher return energy
    # -----------------------------------------------------------------------------------------------
    #             Obtain MO-ADK ionization rate
    #             The difference in the ionization rate for different ionization time is negligable!
    #             So only one set of ionization rate is stored for each intensity.
    # -----------------------------------------------------------------------------------------------
    print('Calculating MO-ADK ionization rate at that ionization time..')
    target_w = np.zeros((np.size(beta), np.size(gamma)))
    index_TI = np.arange(0,Nc*Nb,1, dtype=np.int_)
    iso_IAM(i_E_max)
    TIfunc = partial(TI_rate, F0 = Field[0])
    with mp.Pool(ncore) as p:
        target_w_raw = p.map(TIfunc, index_TI)
    for i in range(0, np.size(index_TI)):
        i_c = i % Nc
        i_b = i // Nc
        target_w[i_b, i_c] = target_w_raw[i]
    np.save( dir + 'TI_rate_' + str(np.round(I,decimals=2)) + 'e14_' + Target, target_w)
    norm1 = np.trapezoid(target_w, gamma, axis=1)
    norm_target = np.trapezoid(norm1*np.sin(beta), beta) * 2.*np.pi
    print('Angular intergal of the ADK rate is', norm_target)

    index = np.arange(0,Nc*Nb*Na,1, dtype=np.int_)

    # Create arrays for each threshold case
    MIT_target_lt4 = np.zeros((np.size(alpha), np.size(beta), np.size(gamma), i_E_max))
    MIT_target_lt6 = np.zeros((np.size(alpha), np.size(beta), np.size(gamma), i_E_max))
    MIT_target_lt8 = np.zeros((np.size(alpha), np.size(beta), np.size(gamma), i_E_max))
    
    func = partial(TI_IAM, i_E_max = i_E_max, target_w = target_w)
    with mp.Pool(ncore) as p:
        MIT_target_raw = p.map(func, index)
    for i in range(0, np.size(index)):
        ia, ib, ic = index_to_angs(i)
        MIT_target_lt4[ia, ib, ic, :] = MIT_target_raw[i]['lt4']
        MIT_target_lt6[ia, ib, ic, :] = MIT_target_raw[i]['lt6']
        MIT_target_lt8[ia, ib, ic, :] = MIT_target_raw[i]['lt8']

    # ----------------------------------------------------------------------------
    #             Orientation averaging
    # ----------------------------------------------------------------------------
    print('Orientation averaging...')
    
    # Process each threshold case
    for threshold_case, MIT_target in [('lt4', MIT_target_lt4), ('lt6', MIT_target_lt6), ('lt8', MIT_target_lt8)]:
        tmp = np.trapezoid(MIT_target, alpha, axis=0)
        tmp2 = np.trapezoid(tmp, beta, axis=0)
        MIT_target_avg = np.trapezoid(tmp2, gamma, axis=0) 
        
        data = np.zeros((nE,3))
        data[:, 0] = nrg/eV
        data[:i_E_max, 1] = MIT_target_avg/ norm_target + I_A[:i_E_max]
        data[:i_E_max, 2] = I_A[:i_E_max] * norm_target + MIT_target_avg #which one of these is correct?
        np.savetxt(dir + Target + '_TI-IAM_Rij_' + threshold_case + '_' + str(np.round(I,decimals=2)) + '.dat', data)
    
    print('Done.')

#PLOTTING RESULTS
fig, axs = plt.subplots(3, 1)  # 3 rows, 1 column

axs[0].plot(np.loadtxt('camphor_IAM_Rij_lt4_0.41.dat')[:, 0], np.loadtxt('camphor_IAM_Rij_lt4_0.41.dat')[:, 2], label= 'camphor')
axs[0].plot(np.loadtxt('fenchone_IAM_Rij_lt4_0.41.dat')[:, 0], np.loadtxt('fenchone_IAM_Rij_lt4_0.41.dat')[:, 2], label= 'fenchone')
axs[0].set_title("Rij < 4")

axs[1].plot(np.loadtxt('camphor_IAM_Rij_lt6_0.41.dat')[:, 0], np.loadtxt('camphor_IAM_Rij_lt6_0.41.dat')[:, 2], label= 'camphor')
axs[1].plot(np.loadtxt('fenchone_IAM_Rij_lt6_0.41.dat')[:, 0], np.loadtxt('fenchone_IAM_Rij_lt6_0.41.dat')[:, 2], label= 'fenchone')
axs[1].set_title("Rij < 6")

axs[2].plot(np.loadtxt('camphor_IAM_Rij_lt8_0.41.dat')[:, 0], np.loadtxt('camphor_IAM_Rij_lt8_0.41.dat')[:, 2], label= 'camphor')
axs[2].plot(np.loadtxt('fenchone_IAM_Rij_lt8_0.41.dat')[:, 0], np.loadtxt('fenchone_IAM_Rij_lt8_0.41.dat')[:, 2], label= 'fenchone')
axs[2].set_title("Rij < 8")

plt.tight_layout()
plt.savefig("Filtered_BL_analysis.pdf")
plt.show()
