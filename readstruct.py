import numpy as np
import wigner
from scipy.special import factorial
import ast

# Conversion factor
I0 = 3.509470 * 10**16  # W/cm^2 in a.u.
nm = 18.897261328856432  # nm in a.u.
eV = 0.03674932539796232  # eV in a.u.
atto = 0.0413413737896  # as in a.u.
fs = atto * 1000.0
c = 137.035999139  # speed of light in a.u.


class structure:
    """Structure for MO-ADK."""

    def __init__(self, target):

        def read_Clm(self, inp):
            with open(inp, "r") as file:
                # Read Ip
                ip_line = next(
                    line for line in file if not line.startswith("#")
                ).strip()
                self.Ip.append(float(ip_line))

                # Read l
                ll_line = next(
                    line for line in file if not line.startswith("#")
                ).strip()
                self.ll.append(np.fromstring(ll_line[1:-1], dtype=int, sep=","))

                # Read m
                m_line = next(line for line in file if not line.startswith("#")).strip()
                self.m.append(ast.literal_eval(m_line[1:-1]))

                # Read C_lm
                clm_line = next(
                    line for line in file if not line.startswith("#")
                ).strip()
                self.C_lm.append(ast.literal_eval(clm_line[1:-1]))

                # Read dip
                dip_line = next(
                    line for line in file if not line.startswith("#")
                ).strip()
                self.dip.append(np.fromstring(dip_line[1:-1], sep=","))

                return

        self.m = []
        self.C_lm = []
        self.ll = []
        self.Ip = []
        self.dip = []

        # Start determining the structure factors
        if target == "butane_anti":
            Norb = 4
            HOMO_ind = 13
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self, "Clms/" + target + "/Clm_" + str(orb) + ".dat")

        elif target == "butane_gauche":
            Norb = 3
            HOMO_ind = 13
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self, "Clms/" + target + "/Clm_" + str(orb) + ".dat")

        elif target == "isobutane":
            Norb = 3
            HOMO_ind = 13
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self, "Clms/" + target + "/Clm_" + str(orb) + ".dat")

        elif target == "fenchone_cat":
            Norb = 2
            HOMO_ind = 31
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self,"/home/lricha95/LIED/Clms/" + target + "/Clm_" + str(orb) + ".dat")

        elif target == "camphor":
            Norb = 2 
            HOMO_ind = 31
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self,"/Users/lricha95/Library/CloudStorage/OneDrive-KennesawStateUniversity(2)/Research/opt_camphor_cat/Clm_optimized_camphor/Clm_" + str(orb) + ".dat")

        elif target == "fenchone":
            Norb = 2
            HOMO_ind = 31
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self,"/Users/lricha95/Library/CloudStorage/OneDrive-KennesawStateUniversity(2)/Research/opt_fenchone_cat/Clms/Clm_" + str(orb) + ".dat")
        elif target == "N2O4_1_8":
            Norb = 5
            HOMO_ind = 17
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self, "Clms/" + target + "/Clm_" + str(orb) + ".dat")
        elif target == "N2O4_2_6":
            Norb = 5
            HOMO_ind = 17
            for orb in range(HOMO_ind, HOMO_ind - Norb, -1):
                read_Clm(self, "Clms/" + target + "/Clm_" + str(orb) + ".dat")
        else:
            print("Structure parameters not found.")
            exit()

    def Q(self, l, m):
        Q = (-1) ** ((m + np.abs(m)) / 2.0) * np.sqrt(
            (2 * l + 1) * factorial(l + np.abs(m)) / (2.0 * factorial(l - np.abs(m)))
        )
        return Q

    def calc_B(self, orb, beta, gamma):
        l = self.ll[orb]
        m = self.m[orb]
        C_lm = self.C_lm[orb]
        # Determine the maximum m' in the summation
        l_max = l.max()
        mp = np.arange(-l_max, l_max + 1, 1)
        # B parameter when the field is negative (direct/backscattered electron to the right)
        B_right = np.zeros(np.size(mp), dtype=np.complex128)
        ang = wigner.D_matrix()
        for im in mp:
            # Sum over possible l
            for i in range(abs(im), l_max + 1):
                if i in l:
                    j = np.where(l == i)[0][0]
                    # Sum over all possible m
                    for jm in range(0, len(m[j])):
                        # im+l_max to map -mp to 0.
                        B_right[im + l_max] += (
                            C_lm[j][jm]
                            * self.Q(l[j], im)
                            * ang.Wigner_D(l[j], im, m[j][jm], 0.0, beta, gamma)
                        )
        return B_right

    def calc_wfat(self, orb, beta, gamma):
        muz = transform(beta, gamma, self.dip[orb]).dip[2].real
        kap = np.sqrt(2.0 * self.Ip[orb])
        # wfat parameter when the field is negative (direct/backscattered electron to the right)
        wfat_right = np.exp(2.0 * kap * muz)
        return wfat_right


class transform:
    """
    For a set of Euler angles (alpha, beta, 0), transform the vector in spherical basis
    Then perform the rotation, and finally transformed the vector back to Cartesian.
    """

    def __init__(self, beta, gamma, dip):
        dip_in = dip
        self.gamma = gamma
        self.beta = beta
        self.dip = self.cart2sph(dip_in)
        self.dip = self.rotate(self.dip)
        self.dip = self.sph2cart(self.dip)

    def cart2sph(self, dip_in):
        A = np.array(
            [
                [1.0 / np.sqrt(2), 1.0j / np.sqrt(2), 0],
                [0, 0, 1],
                [-1.0 / np.sqrt(2), 1.0j / np.sqrt(2), 0],
            ]
        )
        dip = np.zeros(3, dtype=np.complex128)
        dip = np.matmul(A, dip_in)
        return dip

    def sph2cart(self, dip_in):
        A = np.array(
            [
                [1.0 / np.sqrt(2), 0.0, -1.0 / np.sqrt(2)],
                [-1.0j / np.sqrt(2), 0.0, -1.0j / np.sqrt(2)],
                [0, 1, 0],
            ]
        )
        dip = np.zeros(3, dtype=np.complex128)
        dip = np.matmul(A, dip_in)
        return dip

    def rotate(self, dip_in):
        ang = wigner.D_matrix()
        D = np.zeros((3, 3), dtype=np.complex128)
        for mu in range(0, 3):
            for m in range(0, 3):
                D[mu, m] = ang.Wigner_D(1, mu - 1, m - 1, 0.0, self.beta, self.gamma)
        dip = np.zeros(3, dtype=np.complex128)
        dip = np.matmul(D, dip_in)
        return dip
