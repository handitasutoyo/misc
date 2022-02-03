import numpy as np
import math


# z factor calculator from New explicit correlation for the compressibility factor of natural gas: linearized z-factor isotherms
# DOI 10.1007/s13202-015-0209-3
class gas_properties:
    def __init__(self, pressure, temperature, sg):
        self.p = pressure * 14.5038  # psig
        self.t = (temperature * 9 / 5) + 491.67  # Rankine
        self.sg = sg

    def z(self):
        a1 = 0.317842
        a2 = 0.382216
        a3 = -7.768354
        a4 = 14.290531
        a5 = 0.000002
        a6 = -0.004693
        a7 = 0.096254
        a8 = 0.166720
        a9 = 0.966910
        a10 = 0.063069
        a11 = -1.966847
        a12 = 21.0581
        a13 = -27.0246
        a14 = 16.23
        a15 = 207.783
        a16 = -488.161
        a17 = 176.29
        a18 = 1.88453
        a19 = 3.05921
        Tpc = 169.2 + 349.5 * self.sg - 74 * self.sg ** 2
        Ppc = 756.8 - 131.07 * self.sg - 3.6 * self.sg ** 2
        Tpr = self.t / Tpc
        Ppr = self.p / Ppc

        # correlation uses equation nr. 14
        t = 1 / Tpr
        A = a1 * t * Ppr * (math.exp(a2 * (1 - t) ** 2))
        B = a3 * t + a4 * (t ** 2) + a5 * (Ppr ** 6) * (t ** 6)
        C = a9 + a8 * t * Ppr + a7 * (t ** 2) * (Ppr ** 2) + a6 * (t ** 3) * (Ppr ** 3)
        D = a10 * t * (math.exp(a11 * (1 - t) ** 2))
        E = a12 * t + a13 * (t ** 2) + a14 * (t ** 3)
        F = a15 * t + a16 * (t ** 2) + a17 * (t ** 3)
        G = a18 + a19 * t
        y = (D * Ppr) / (((1 + A ** 2) / C) - ((B * A ** 2) / C ** 3))
        z = (D * Ppr * (1 + y + (y ** 2) - (y ** 3))) / (D * Ppr + E * (y ** 2) - (F * (y ** G)) * (1 - y) ** 3)
        return z

    def gas_dens(self):
        z=self.z()
        dens = 28.97 * 0.001 * 100000 * self.sg * (self.p / 14.5038) / (z* 8.314 * (((self.t - 491.67) * (5 / 9)) + 273.15))
        return dens