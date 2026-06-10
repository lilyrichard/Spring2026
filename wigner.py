# Copyright by C.H. Isaac Yuen, Kansas State University (2021).
# Please email me before redistributing this code.
# Email: iyuen@phys.ksu.edu

import numpy as np
from scipy.special import factorial

class D_matrix:
    '''Library for calculating Wigner D-matrix. This function calculate the
    D^j1_{m1p,m1}(alpha,beta,gamma), with
    ang1=alpha,ang2=beta,ang3=gamma.
    It adopt the z-y-z convention for the Euler rotation.
    Therefore, D^j1_{m1p,m1}(alpha,beta,gamma) =
    exp(-i*m1p*alpha) * d^j1_{m1p,m1}(beta) * exp(-i*m1*gamma)
    with REAL Wigner small d-matrix.
    The small d-matrix is defined using the Condon-Shortley phase convention.
    j1,m1p,m1 are ASSUMED to be integers!!!!'''

    def __init__(self):
        return

    def Wigner_small_d(self, j1, m1p, m1, ang2):
        tol = 1e-6
        Wigner_small_d = 0.
        # Condtion check
        if abs(m1p > j1):
            print('|m1p| > j1, program exit at Wigner_small_d.')
            exit()
        if abs(m1 > j1):
            print('|m1| > j1, program exit at Wigner_small_d.')
            exit()
        # If beta is 0, then d^j1_{m1p,m1}(0) = delta_{m1p,m1}
        if (abs(ang2) < tol):
            if (m1p == m1):
                Wigner_small_d = 1.
                return Wigner_small_d
            else:
                Wigner_small_d = 0.
                return Wigner_small_d
        # If beta is pi, then d^j1_{m1p,m1}(pi) = (-1)**(j1-m1) delta_{m1p,-m1}
        if (abs(np.pi - ang2) < tol):
            if (m1p == -m1):
                Wigner_small_d = (-1.)**(j1-m1)
                return Wigner_small_d
            else:
                Wigner_small_d = 0.
                return Wigner_small_d

        # If beta is -pi, then d^j1_{m1p,m1}(pi) = (-1)**(j1-m1) delta_{m1p,-m1}
        if (abs(-np.pi - ang2) < tol):
            if (m1p == -m1):
                Wigner_small_d = (-1.)**(j1-m1)
                return Wigner_small_d
            else:
                Wigner_small_d = 0.
                return Wigner_small_d

        # Calculate the summation over t

        tmin = max(0, m1-m1p)
        tmax = min(m1+j1, j1-m1p)

        for t in range(tmin, tmax+1):
            fac1 = factorial(j1-m1p-t)
            fac2 = factorial(j1+m1-t)
            fac3 = factorial(t+m1p-m1)
            fac4 = factorial(t)
            Wigner_small_d = Wigner_small_d + (-1.)**t * np.cos(ang2/2.)**(-2*t) \
                * np.sin(ang2/2.)**(2*t) / (fac1*fac2*fac3*fac4)

        # Gather the rest of the factors
        fac1 = float(factorial(j1+m1))
        fac2 = float(factorial(j1-m1))
        fac3 = float(factorial(j1+m1p))
        fac4 = float(factorial(j1-m1p))

        Wigner_small_d = (-1.)**(m1p-m1) * np.sqrt(fac1*fac2*fac3*fac4) \
            * np.cos(ang2/2.)**(2*j1 + m1 - m1p) * np.sin(ang2/2.)**(m1p - m1) \
            * Wigner_small_d

        return Wigner_small_d

    def Wigner_D(self, j1, m1p, m1, ang1, ang2, ang3):
        if(j1 == 0):
            Wigner_D = 0.
        fac1 = np.exp(-1j*m1p*ang1)
        fac2 = np.exp(-1j*m1*ang3)
        Wigner_D = fac1 * fac2 * self.Wigner_small_d(j1, m1p, m1, ang2)
        return Wigner_D
