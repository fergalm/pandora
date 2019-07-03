# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Thu Jun 27 15:32:38 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np

import obsstrategy as strat


def main():
    nTheta = 160
    nPhi = 160
    phi0 = 0
    phi1 = 180
#    nTheta = 40
#    nPhi = 40

    alpha_rad = np.radians(90)
    zeta_rad = np.radians(80)

    theta = np.linspace(-90, 90, nTheta)
    phi = np.linspace(phi0, phi1, nPhi)

    dutyCycle = np.zeros( (len(theta), len(phi)) )
    for i in range(nTheta):
        for j in range(nPhi):
            starVec = strat.computeStarUnitVector(phi[j], theta[i])

            try:
                dutyCycle[i, j] = strat.computeDutyCycle(alpha_rad, starVec, zeta_rad)
            except AssertionError:
                print("Assertion raised")
                dutyCycle[i, j] = -1
#        print(x1[i])
#        print( dutyCycle[i])

    plt.clf()
    cmap = plt.cm.RdBu
    cmap.set_under('w')
    extent=[phi0, phi1, -90, 90]

    plt.imshow(dutyCycle, origin='bottom', cmap=cmap, extent=extent)
    plt.clim(.3, .7)
    cb = plt.colorbar()
    cb.set_label("Duty Cycle")
    plt.title("Max Zenith Angle %.1f" %(np.degrees(zeta_rad)))
#    plt.clim(.2, .8)
    plt.xlabel("Ecliptic Longitude (degrees)")
    plt.ylabel("Ecliptic Latidutde (degrees)")