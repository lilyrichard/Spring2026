import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

# -------- Set fonts BEFORE creating figure --------
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 16          # default size for text
mpl.rcParams['axes.labelsize'] = 16      # axis labels
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16

# Energy range
energies = np.arange(50, 82, 2)
base_path = "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research"
molecules = ["camphor", "fenchone", "vertonal"]

CvV_matrix_f1, FvV_matrix_f1 = [], []
CvV_matrix_f2, FvV_matrix_f2 = [], []

for energy in energies:
    # -------- filename1 --------
    datasets_f1 = {}
    for mol in molecules:
        folder = os.path.join(base_path, f"opt_{mol}_cat", "DCS_Ratio")
        filename1 = f"DCS_R_{mol}_cat_{energy}_avg.dat"
        datasets_f1[mol] = np.loadtxt(os.path.join(folder, filename1))

    dcs_camphor_f1  = datasets_f1["camphor"]
    dcs_fenchone_f1 = datasets_f1["fenchone"]
    dcs_vertonal_f1 = datasets_f1["vertonal"]

    CvV_matrix_f1.append((dcs_camphor_f1[:, 1] - dcs_vertonal_f1[:, 1]) / 
                         (dcs_camphor_f1[:, 1] + dcs_vertonal_f1[:, 1]))
    
    FvV_matrix_f1.append((dcs_fenchone_f1[:, 1] - dcs_vertonal_f1[:, 1]) / 
                         (dcs_vertonal_f1[:, 1] + dcs_fenchone_f1[:, 1]))

    # -------- filename2 --------
    datasets_f2 = {}
    for mol in molecules:
        folder = os.path.join(base_path, f"opt_{mol}_cat", "DCS_Ratio")
        filename2 = f"DCS_{mol}_cat_{energy}_iso.dat"
        datasets_f2[mol] = np.loadtxt(os.path.join(folder, filename2))

    dcs_camphor_f2  = datasets_f2["camphor"]
    dcs_fenchone_f2 = datasets_f2["fenchone"]
    dcs_vertonal_f2 = datasets_f2["vertonal"]

    CvV_matrix_f2.append((dcs_camphor_f2[:, 1] - dcs_vertonal_f2[:, 1]) / 
                         (dcs_camphor_f2[:, 1] + dcs_vertonal_f2[:, 1]))
    
    FvV_matrix_f2.append((-dcs_vertonal_f2[:, 1] + dcs_fenchone_f2[:, 1]) / 
                         (dcs_vertonal_f2[:, 1] + dcs_fenchone_f2[:, 1]))

CvV_matrix_f1 = np.array(CvV_matrix_f1)
FvV_matrix_f1 = np.array(FvV_matrix_f1)
CvV_matrix_f2 = np.array(CvV_matrix_f2)
FvV_matrix_f2 = np.array(FvV_matrix_f2)

angles = dcs_camphor_f1[:, 0]
A_grid, E_grid = np.meshgrid(angles, energies)

vmax = max(
    np.max(np.abs(CvV_matrix_f1)),
    np.max(np.abs(FvV_matrix_f1)),
    np.max(np.abs(CvV_matrix_f2)),
    np.max(np.abs(FvV_matrix_f2))
)

# -------- GridSpec layout --------
fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 3, width_ratios=[1, 1, 0.05], figure=fig, wspace=0.25)

ax00 = fig.add_subplot(gs[0, 0])
ax01 = fig.add_subplot(gs[0, 1])
cax0 = fig.add_subplot(gs[0, 2])  # colorbar top

ax10 = fig.add_subplot(gs[1, 0])
ax11 = fig.add_subplot(gs[1, 1])
cax1 = fig.add_subplot(gs[1, 2])  # colorbar bottom

# -------- Plots --------
cf1 = ax00.contourf(A_grid, E_grid, CvV_matrix_f1, levels=50, cmap="seismic", vmin=-vmax, vmax=vmax)
cf2 = ax01.contourf(A_grid, E_grid, FvV_matrix_f1, levels=50, cmap="seismic", vmin=-vmax, vmax=vmax)
cf3 = ax10.contourf(A_grid, E_grid, CvV_matrix_f2, levels=50, cmap="seismic", vmin=-vmax, vmax=vmax)
cf4 = ax11.contourf(A_grid, E_grid, FvV_matrix_f2, levels=50, cmap="seismic", vmin=-vmax, vmax=vmax)

# -------- Colorbars --------
cbar0 = fig.colorbar(cf2, cax=cax0)
cbar0.ax.tick_params(labelsize=14)  # colorbar font

cbar1 = fig.colorbar(cf4, cax=cax1)
cbar1.ax.tick_params(labelsize=14)

# -------- Labels --------
ax00.set_ylabel("Energy (eV)")
ax10.set_xlabel(r"Angle$^\circ$")
ax10.set_ylabel("Energy (eV)")
ax11.set_xlabel(r"Angle$^\circ$")

# -------- Tick cleanup --------
ax00.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
ax01.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
ax01.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
ax11.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

# -------- Subfigure labels --------
ax00.text(-0.15, 1.05, 'a)', transform=ax00.transAxes, fontsize=18, fontweight='bold', va='top')
ax01.text(-0.15, 1.05, 'b)', transform=ax01.transAxes, fontsize=18, fontweight='bold', va='top')
ax10.text(-0.15, 1.05, 'c)', transform=ax10.transAxes, fontsize=18, fontweight='bold', va='top')
ax11.text(-0.15, 1.05, 'd)', transform=ax11.transAxes, fontsize=18, fontweight='bold', va='top')

# -------- Set x-axis limits and ticks --------
for ax in [ax00, ax01, ax10, ax11]:
    ax.set_xlim(50, 175)
    ax.set_xticks(np.arange(50, 200, 25))

plt.subplots_adjust(left=0.08, right=0.93, top=0.95, bottom=0.08)
plt.savefig("/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/DAMOP_2026/Results/Combined_Asymmetry_Plots.pdf", dpi=150)
plt.show()