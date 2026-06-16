import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', family='Times New Roman')
plt.rc('axes', titlesize=18, labelsize=16)
plt.rc('legend', fontsize=14)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)


def load_plot_data(path):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    valid = (data[:, 0] > 0) & (data[:, 2] != 0)
    return data[valid]

fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6.5, 8.5))  # 3 rows, 1 column, publication size

camphor_lt4 = load_plot_data('./camphor/camphor_IAM_Rij_lt4_0.41.dat')
fenchone_lt4 = load_plot_data('./fenchone/fenchone_IAM_Rij_lt4_0.41.dat')
vertonal_lt4 = load_plot_data('./vertonal/vertonal_IAM_Rij_lt4_0.41.dat')
axs[0].plot(camphor_lt4[:, 0], camphor_lt4[:, 2], label='camphor')
axs[0].plot(fenchone_lt4[:, 0], fenchone_lt4[:, 2], label='fenchone')
axs[0].plot(vertonal_lt4[:, 0], vertonal_lt4[:, 2], label='vertonal')
axs[0].text(0.02, 0.95, 'Rij < 4', transform=axs[0].transAxes, va='top', ha='left', fontsize=16,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
#axs[0].set_xlabel('Energy (eV)')
axs[0].set_ylabel('DCS')
axs[0].set_xlim(50, 100)
axs[0].legend(loc='upper right')

camphor_lt6 = load_plot_data('./camphor/camphor_IAM_Rij_lt6_0.41.dat')
fenchone_lt6 = load_plot_data('./fenchone/fenchone_IAM_Rij_lt6_0.41.dat')
vertonal_lt6 = load_plot_data('./vertonal/vertonal_IAM_Rij_lt6_0.41.dat')
axs[1].plot(camphor_lt6[:, 0], camphor_lt6[:, 2], label='camphor')
axs[1].plot(fenchone_lt6[:, 0], fenchone_lt6[:, 2], label='fenchone')
axs[1].plot(vertonal_lt6[:, 0], vertonal_lt6[:, 2], label='vertonal')
axs[1].text(0.02, 0.95, 'Rij < 6', transform=axs[1].transAxes, va='top', ha='left', fontsize=16,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
#axs[1].set_xlabel('Energy (eV)')
axs[1].set_ylabel('DCS')

camphor_lt8 = load_plot_data('./camphor/camphor_IAM_Rij_lt8_0.41.dat')
fenchone_lt8 = load_plot_data('./fenchone/fenchone_IAM_Rij_lt8_0.41.dat')
vertonal_lt8 = load_plot_data('./vertonal/vertonal_IAM_Rij_lt8_0.41.dat')
axs[2].plot(camphor_lt8[:, 0], camphor_lt8[:, 2], label='camphor')
axs[2].plot(fenchone_lt8[:, 0], fenchone_lt8[:, 2], label='fenchone')
axs[2].plot(vertonal_lt8[:, 0], vertonal_lt8[:, 2], label='vertonal')
axs[2].text(0.02, 0.95, 'Rij < 8', transform=axs[2].transAxes, va='top', ha='left', fontsize=16,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
axs[2].set_xlabel('Energy (eV)')
axs[2].set_ylabel('DCS')

ymin, ymax = axs[0].get_ylim()
yticks = np.arange(np.floor(ymin), np.ceil(ymax) + 1, 1)
for ax in axs:
    ax.set_yticks(yticks)

plt.tight_layout()
plt.savefig("Filtered_BL_analysis.pdf")
plt.show()