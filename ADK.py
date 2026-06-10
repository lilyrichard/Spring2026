# Copyright by C.H. Isaac Yuen, Kansas State University (2021).
# Please email me before redistributing this code.
# Email: iyuen@phys.ksu.edu

import numpy as np
from scipy.special import factorial

class ADK:
    '''Class of functions for ADT rate.'''

    def __init__(self, field, m, Ip, Zeff, B):
        self.field = field + 10**-60  # Avoid division by zero
        self.m = m
        self.Ip = Ip
        self.Zeff = Zeff
        self.B = B
        self.kappa = np.sqrt(2.*self.Ip)

    def ADK_rate(self):
        fac2 = np.abs(self.B)**2/(2**(np.abs(self.m))
                                  * factorial(np.abs(self.m)))
        fac3 = 1./(self.kappa**(2.*self.Zeff/self.kappa - 1.))
        fac5 = (2.*self.kappa**3/np.abs(self.field)) ** \
            (2.*self.Zeff/self.kappa - np.abs(self.m) - 1.)
        fac7 = np.exp(-2./3.*self.kappa**3/np.abs(self.field))
        ADK_rate = fac2*fac3*fac5*fac7
        return ADK_rate