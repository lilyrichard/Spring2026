import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

nm = 18.897261328856432  # nm in a.u.
I0 = 3.509470*10**16  # W/cm^2 in a.u.
c = 137.035999139  # speed of light in a.u.
eV = 0.03674932539796232  # eV in a.u.
angs = nm/10.    # Angstorm in a.u.

# ----------------------------------------------------------------------------
#             Input
# ----------------------------------------------------------------------------
TS_case = 'eq'
# TS_case = 'ff'
# TS_case = 'Ey'

Target_list = ['camphor_cat', 'fenchone_cat']

# ----------------------------------------------------------------------------
#             Set up the frame of the figure
# ----------------------------------------------------------------------------
plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['font.size'] = 28
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = "cm"
plt.rcParams['xtick.major.size'] = 8.
plt.rcParams['xtick.minor.size'] = 5.
plt.rcParams['ytick.major.size'] = 5.
plt.rcParams['ytick.minor.size'] = 3.
plt.rcParams['lines.linewidth'] = 2.5

# ----------------------------------------------------------------------------
#             Read the Cartesian coordinate and Atom names from files
# ----------------------------------------------------------------------------

# Bin settings are shared across targets so the x-axis limits can also be shared
bin_width = 0.2
bins = np.arange(2.0, 10, bin_width)

bond_histograms = {}
max_frequency = {'CC': 0, 'CH': 0, 'HH': 0}

for Target in Target_list:
    print('Target is', Target)

    if TS_case == 'eq':
        inp = Target + '_eq.xyz'
    elif TS_case == 'ff':
        inp = Target + '_6_7fs_0E_cart.dat'
    elif TS_case == 'Ey':
        inp = Target + '_10fs_Ey_cart.dat'

    print('Reading coordinates from', inp)
    with open(inp) as file:
        target_cart = []
        atom = []
        for line in file:
            tmp = line.rsplit()
            if len(tmp) < 4:
                continue
            cart = np.zeros(3)
            atom.append(tmp[0])
            cart[0:3] = tmp[1:4]
            cart *= angs
            target_cart.append(cart)

    CC_bin_values = []
    HH_bin_values = []
    CH_bin_values = []

    for i in range(len(atom)):
        for j in range(i+1, len(atom)):
            R12 = np.linalg.norm(target_cart[i] - target_cart[j])
            if atom[i] == 'H' and atom[j] == 'H':
                HH_bin_values.append(R12)
            elif atom[i] == 'C' and atom[j] == 'C':
                CC_bin_values.append(R12)
            else:
                CH_bin_values.append(R12)

    CC_frequency, _ = np.histogram(CC_bin_values, bins=bins)
    CH_frequency, _ = np.histogram(CH_bin_values, bins=bins)
    HH_frequency, _ = np.histogram(HH_bin_values, bins=bins)

    bond_histograms[Target] = {
        'CC': CC_frequency,
        'CH': CH_frequency,
        'HH': HH_frequency
    }

    max_frequency['CC'] = max(max_frequency['CC'], CC_frequency.max())
    max_frequency['CH'] = max(max_frequency['CH'], CH_frequency.max())
    max_frequency['HH'] = max(max_frequency['HH'], HH_frequency.max())

for Target in Target_list:
    CC_frequency = bond_histograms[Target]['CC']
    CH_frequency = bond_histograms[Target]['CH']
    HH_frequency = bond_histograms[Target]['HH']

    fig2, ax2 = plt.subplots(figsize=(6, 9), nrows=3, layout='constrained')
    ax2[0].barh(y=bins[:-1], width=CC_frequency, height=bin_width, align='edge', label='CC', color='b')
    ax2[1].barh(y=bins[:-1], width=CH_frequency, height=bin_width, align='edge', label='CH', color='g')
    ax2[2].barh(y=bins[:-1], width=HH_frequency, height=bin_width, align='edge', label='HH', color='r')

    ytick = np.arange(2, 11, dtype=int)
    ax2[0].set_yticks(ytick)
    ax2[1].set_yticks(ytick)
    ax2[2].set_yticks(ytick)
    ax2[0].set_yticklabels(ytick)
    ax2[1].set_yticklabels(ytick)
    ax2[2].set_yticklabels(ytick)

    xlim_cc = max_frequency['CC'] + 1
    xlim_ch = int(2 * np.ceil((max_frequency['CH'] + 1) / 2.0))
    xlim_hh = int(2 * np.ceil((max_frequency['HH'] + 1) / 2.0))

    ax2[0].set_xlim(0, xlim_cc)
    ax2[1].set_xlim(0, xlim_ch)
    ax2[2].set_xlim(0, xlim_hh)

    ax2[0].set_xticks(np.arange(0, xlim_cc + 1, 1, dtype=int))
    ax2[1].set_xticks(np.arange(0, xlim_ch + 1, 2, dtype=int))
    ax2[2].set_xticks(np.arange(0, xlim_hh + 1, 2, dtype=int))

    for axes in ax2:
        axes.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        axes.tick_params(axis='both', which='major', labelsize=16)
        axes.legend(fontsize=20, loc='upper right')
        axes.grid(True)

    ax2[1].set_ylabel(r'Bond length ($a_0$)')
    ax2[2].set_xlabel('Frequency')

    fig2.savefig(Target + '_bond_length_histo.pdf')
