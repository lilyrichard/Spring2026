from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 14           # default size for text
mpl.rcParams['axes.labelsize'] = 16      # axis labels
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

Target = 'camphor_cat'
path = '/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/Bond_Length_Analysis/data'

leq4 = np.loadtxt(f'{path}/DCS_R_{Target}_50_avg_leq4.dat')
leq6 = np.loadtxt(f'{path}/DCS_R_{Target}_50_avg_leq6.dat')
leq8 = np.loadtxt(f'{path}/DCS_R_{Target}_50_avg_leq8.dat')
all = np.loadtxt(f'{path}/DCS_R_{Target}_50_avg_all.dat')

leq4iso = np.loadtxt(f'{path}/DCS_{Target}_50_iso_leq4.dat')
leq6iso = np.loadtxt(f'{path}/DCS_{Target}_50_iso_leq6.dat')
leq8iso = np.loadtxt(f'{path}/DCS_{Target}_50_iso_leq8.dat')
alliso = np.loadtxt(f'{path}/DCS_{Target}_50_iso_all.dat')

Angle = leq4[:, 0]

#**************
plt.figure()
plt.plot(Angle, leq4iso[:,1], label='leq4')
plt.plot(Angle, leq6iso[:,1], label='leq6') 
plt.plot(Angle, leq8iso[:,1], label='leq8') 
plt.plot(Angle, alliso [:,1], label='all') 

plt.xlabel('Angle')
plt.ylabel('Molecular Interferance Term (iso)')

plt.grid(True)
plt.legend()

# Save figure

plt.savefig(f'{path}/DCS_iso_{Target}_AnalysisPlot.pdf')
plt.show()
#***********
plt.figure()
plt.plot(Angle, leq4[:,1], label='leq4')
plt.plot(Angle, leq6[:,1], label='leq6') 
plt.plot(Angle, leq8[:,1], label='leq8') 
plt.plot(Angle, all[:,1], label='all') 

plt.xlabel('Angle')
plt.ylabel('Molecular Interferance Term')

plt.grid(True)
plt.legend()

# Save figure

plt.savefig(f'{path}/DCS_R_{Target}_AnalysisPlot.pdf')
plt.show()


plt.close()