import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

# -------- Set fonts BEFORE creating figure --------
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16

# Energy range
energies = np.arange(50, 82, 2)
base_path = "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research"

molecules = ["camphor", "fenchone"]

CvF_matrix_f1 = []
CvF_matrix_f2 = []

for energy in energies:

    # -------- filename1 --------
    datasets_f1 = {}
    for mol in molecules:
        folder = os.path.join(base_path, f"opt_{mol}_cat", "DCS_Ratio")
        filename1 = f"DCS_R_{mol}_cat_{energy}_avg.dat"
        datasets_f1[mol] = np.loadtxt(os.path.join(folder, filename1))

    dcs_camphor_f1  = datasets_f1["camphor"]
    dcs_fenchone_f1 = datasets_f1["fenchone"]

    CvF_matrix_f1.append(
        (dcs_camphor_f1[:,1] - dcs_fenchone_f1[:,1]) /
        (dcs_camphor_f1[:,1] + dcs_fenchone_f1[:,1])
    )

    # -------- filename2 --------
    datasets_f2 = {}
    for mol in molecules:
        folder = os.path.join(base_path, f"opt_{mol}_cat", "DCS_Ratio")
        filename2 = f"DCS_{mol}_cat_{energy}_iso.dat"
        datasets_f2[mol] = np.loadtxt(os.path.join(folder, filename2))

    dcs_camphor_f2  = datasets_f2["camphor"]
    dcs_fenchone_f2 = datasets_f2["fenchone"]

    CvF_matrix_f2.append(
        (dcs_camphor_f2[:,1] - dcs_fenchone_f2[:,1]) /
        (dcs_camphor_f2[:,1] + dcs_fenchone_f2[:,1])
    )

CvF_matrix_f1 = np.array(CvF_matrix_f1)
CvF_matrix_f2 = np.array(CvF_matrix_f2)

angles = dcs_camphor_f1[:,0]
A_grid, E_grid = np.meshgrid(angles, energies)

vmax = max(
    np.max(np.abs(CvF_matrix_f1)),
    np.max(np.abs(CvF_matrix_f2))
)

# -------- Layout --------
fig = plt.figure(figsize=(13,5))
gs = fig.add_gridspec(1, 2, wspace=0.25)

ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[0,1])

# -------- Plots --------
cf1 = ax0.contourf(
    A_grid, E_grid, CvF_matrix_f1,
    levels=50, cmap="seismic",
    vmin=-vmax, vmax=vmax
)

cf2 = ax1.contourf(
    A_grid, E_grid, CvF_matrix_f2,
    levels=50, cmap="seismic",
    vmin=-vmax, vmax=vmax
)

# -------- Colorbars --------
cbar0 = fig.colorbar(cf1, ax=ax0, pad=0.02, fraction=0.046)
cbar0.ax.tick_params(labelsize=14)

cbar1 = fig.colorbar(cf2, ax=ax1, pad=0.02, fraction=0.046)
cbar1.ax.tick_params(labelsize=14)

# -------- Labels --------
ax0.set_xlabel(r"Angle$^\circ$")
ax0.set_ylabel("Energy (eV)")

ax1.set_xlabel(r"Angle$^\circ$")
ax1.tick_params(axis='y', which='both', left=False, labelleft=False)

# -------- Subfigure labels --------
ax0.text(-0.15, 1.05, 'a)', transform=ax0.transAxes,
         fontsize=18, fontweight='bold', va='top')
ax1.text(-0.15, 1.05, 'b)', transform=ax1.transAxes,
         fontsize=18, fontweight='bold', va='top')

# -------- Axis limits --------
for ax in [ax0, ax1]:
    ax.set_xlim(50, 175)
    ax.set_xticks(np.arange(50, 200, 25))

plt.subplots_adjust(left=0.06, right=0.94, top=0.95, bottom=0.12)

plt.savefig(
    "/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research/MDPI_Manuscript/updated_figures/CvF_ICF_heatmap.pdf",
    dpi=150
)

plt.show()