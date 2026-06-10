#plots backscattered DCS from TI-IAM and IAM for three isomers 
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

# -----------------------------
# Settings
# -----------------------------
targets = ['camphor_cat', 'fenchone_cat', 'vertonal_cat']
Energy = np.arange(50, 81, 2)
base_path = Path("/Users/lilyrichard/Library/CloudStorage/OneDrive-KennesawStateUniversity/Research")

camphor_cat_title = "Camphor Cation"
fenchone_cat_title = "Fenchone Cation"
vertonal_cat_title = "Vertonal Cation"
# -----------------------------
# Loop over targets
# -----------------------------
for Target in targets:
    # Set the title based on target
    if Target == 'camphor_cat':
        title = camphor_cat_title
    elif Target == 'fenchone_cat':
        title = fenchone_cat_title
    elif Target == 'vertonal_cat':
        title = vertonal_cat_title

    # Lists to store values
    values_TI_IAM = []
    values_IAM = []

    for E in Energy:
        # -----------------------------
        # TI-IAM data
        # -----------------------------
        filename_TI_IAM = f'DCS_R_{Target}_{E}_avg.dat'
        full_path_TI_IAM = base_path / f"opt_{Target}" / "DCS_Ratio" / filename_TI_IAM

        if full_path_TI_IAM.exists():
            data = np.atleast_2d(np.loadtxt(full_path_TI_IAM))
            values_TI_IAM.append(data[-1, 1])
        else:
            print(f"Warning: TI-IAM file not found: {full_path_TI_IAM}")
            values_TI_IAM.append(np.nan)

        # -----------------------------
        # IAM data
        # -----------------------------
        filename_IAM = f'DCS_{Target}_{E}_iso.dat'
        full_path_IAM = base_path / f"opt_{Target}" / "DCS_Ratio" / filename_IAM

        if full_path_IAM.exists():
            data = np.atleast_2d(np.loadtxt(full_path_IAM))
            values_IAM.append(data[-1, 1])
        else:
            print(f"Warning: IAM file not found: {full_path_IAM}")
            values_IAM.append(np.nan)

        
        output_file = base_path / f"{Target}_BackscatteredDCS.csv"

        with open(output_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            
            # Header
            writer.writerow(["Energy (eV)", "TI-IAM DCS", "IAM DCS"])
            
            # Data rows
            for e, ti, iam in zip(Energy, values_TI_IAM, values_IAM):
                writer.writerow([e, ti, iam])

        print(f"Saved data file: {output_file}")



    # -----------------------------
    # Plot per target
    # -----------------------------
    plt.figure()
    plt.plot(Energy, values_TI_IAM, 'b-o', label='TI-IAM')
    plt.plot(Energy, values_IAM, 'r-s', label='IAM')  # different color/marker

    plt.xlabel('Rescattering Energy (eV)')
    plt.ylabel('Differential Cross Section (DCS)')
    plt.title(f'{title} Differential Cross Section')
    plt.grid(True)
    plt.legend()

    # Save figure
    
    plt.savefig(base_path/f'DCS_R_{Target}_comparison.pdf')
    plt.show()
    print(f"Saved figure: {base_path}")

    plt.close()