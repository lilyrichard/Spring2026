#plots BACKSCATTERED ICF as function of energy

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 16            # base font size
mpl.rcParams['axes.labelsize'] = 18        # x/y labels
mpl.rcParams['axes.titlesize'] = 20        # title
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 14

camphor_file = "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/DAMOP_2026/Results/Backscattered_DCS/camphor_cat_BackscatteredDCS.csv"
fenchone_file = "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/DAMOP_2026/Results/Backscattered_DCS/fenchone_cat_BackscatteredDCS.csv" 
vertonal_file = "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/DAMOP_2026/Results/Backscattered_DCS/vertonal_cat_BackscatteredDCS.csv"


camphor = np.loadtxt(camphor_file, delimiter=",", skiprows = 1)
fenchone = np.loadtxt(fenchone_file, delimiter=",", skiprows = 1)
vertonal = np.loadtxt(vertonal_file, delimiter=",", skiprows = 1)

Energy = camphor[:, 0]

camphor_IAM=camphor[:, 2]
fenchone_IAM=fenchone[:, 2]
vertonal_IAM=vertonal[:, 2]

camphor_TIIAM= camphor[:, 1]
fenchone_TIIAM= fenchone[:, 1]
vertonal_TIIAM= vertonal[:, 1]

ICF_TIIAM_CV = (camphor_TIIAM - vertonal_TIIAM) / (camphor_TIIAM + vertonal_TIIAM)
ICF_TIIAM_FV = (fenchone_TIIAM - vertonal_TIIAM) / (fenchone_TIIAM + vertonal_TIIAM)

ICF_IAM_CV = (camphor_IAM - vertonal_IAM) / (camphor_IAM + vertonal_IAM)
ICF_IAM_FV = (fenchone_IAM - vertonal_IAM) / (fenchone_IAM + vertonal_IAM)


##################DCS plots##########################
plt.plot(Energy, camphor_IAM, label = "Camphor")
plt.plot(Energy, fenchone_IAM, label = "Fenchone")
plt.plot(Energy, vertonal_IAM, label = "Vertonal")
plt.title("IAM Backscattered DCS")
plt.xlabel("Energy (eV)")
plt.ylabel("DCS")
plt.legend()
plt.show()

plt.plot(Energy, camphor_TIIAM, label = "Camphor")
plt.plot(Energy, fenchone_TIIAM, label = "Fenchone")
plt.plot(Energy, vertonal_TIIAM, label = "Vertonal")
plt.title("TI-IAM Backscattered DCS")
plt.xlabel("Energy (eV)")
plt.ylabel("DCS")
plt.legend()
plt.show()
##################ICF plots##########################

# Customize margins for the final ICF figures
left_margin = 0.16
right_margin = 0.96
top_margin = 0.94
bottom_margin = 0.14

fig1 = plt.figure()

plt.plot(Energy, ICF_TIIAM_CV, label='TI-IAM (Camphor–Vertonal)')
plt.plot(Energy, ICF_IAM_CV, label='IAM (Camphor–Vertonal)')

plt.xlabel("Energy (eV)")
plt.ylabel("ICF")
plt.legend()
plt.ylim(-0.12, 0.28)
fig1.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin)

fig1.savefig("/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/manuscript_figures/Backscattered CV ICF vs Energy.pdf", bbox_inches='tight', pad_inches=0.10)
plt.show()


fig2 = plt.figure()

plt.plot(Energy, ICF_TIIAM_FV, label='TI-IAM (Fenchone–Vertonal)')
plt.plot(Energy, ICF_IAM_FV, label='IAM (Fenchone–Vertonal)')

plt.xlabel("Energy (eV)")
plt.ylabel("ICF")
plt.legend()
plt.ylim(-0.12, 0.28)
fig2.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin)

fig2.savefig("/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/manuscript_figures/Backscattered FV ICF vs Energy.pdf", bbox_inches='tight', pad_inches=0.10)
plt.show()