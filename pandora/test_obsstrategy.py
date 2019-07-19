# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Thu Jun 27 09:46:13 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import numpy as np

import obsstrategy as obs
import pytest

#def computeVisibilityV1(alpha_rad, starUnitVec, earthAvoidanceAngle_rad):
def isZeroAndPi(res):
    return np.allclose(res, [0, 180]) or np.allclose(res, [180, 360]), res


def test_CardinalPointsAlpha0Zeta90():
    """Simplest case with alpha=0 degrees  and eps=90 deg"""
    alpha_rad = 0
    maxZenithAngle_rad = np.radians(90)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = np.array([0,1,0])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(np.sort(res))
    assert np.allclose(res, [90, 270]), res

    starUnitVec = np.array([0,0,1])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(np.sort(res))
    assert isZeroAndPi(res)

    starUnitVec = np.array([0,0,-1])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(np.sort(res))
    assert isZeroAndPi(res)


def test_varyElngElatAlpha0Zeta90():
    """Increase epsilon, earth avoidance angle"""

    alpha_rad = 0
    maxZenithAngle_rad = np.radians(90)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = obs.computeStarUnitVector(0, 30)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert isZeroAndPi(res)

    starUnitVec = obs.computeStarUnitVector(0, -30)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert isZeroAndPi(res)

    starUnitVec = obs.computeStarUnitVector(30, 0)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert isZeroAndPi(res)

    starUnitVec = obs.computeStarUnitVector(-30, 0)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert isZeroAndPi(res)

    starUnitVec = obs.computeStarUnitVector(30, 30)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert np.allclose(res, [139.1066, 319.1066]), res

    starUnitVec = obs.computeStarUnitVector(90, 30)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.sort(np.degrees(res))
    assert np.allclose(res, [120, 300]), res

def test_varyZetaElng0():
    alpha_rad = 0
    func = obs.computeExtremaOfOrbitalPhase_rad
    starUnitVec = obs.computeStarUnitVector(0, 85)

    expected = [
                [81.329, 98.671],
                [60.381, 119.619],
                [30.12643986, 149.87356014],
                [10.03859332, 169.96140668],
               ]

    for i,z in enumerate([10, 30, 60, 80]):
        maxZenithAngle_rad = np.radians(z)
        res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        assert np.allclose(res, expected[i]), "%i %s" %(i, res)

    res = func(alpha_rad, starUnitVec, np.pi/2.)
    assert isZeroAndPi(res)


def test_varyElngZeta90():
    alpha_rad = 0
    maxZenithAngle_rad = np.radians(90)
    func = obs.computeExtremaOfOrbitalPhase_rad

    eLng = [30, 60, 90, 120, 150, 210, 240, 270, 300, 330]
    expected = [
                [153.435, 333.435],
                [139.106, 319.107],
                [135, 315],

                [139.106, 319.107],
                [153.435, 333.435],
                #180 is not tested in the loop

                [26.565, 206.565],
                [40.893, 220.893],
                [45, 225],

                [40.893, 220.893],
                [26.565, 206.565],
               ]

    starUnitVec = obs.computeStarUnitVector(0, 45)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert isZeroAndPi(res)

    starUnitVec = obs.computeStarUnitVector(180, 45)
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert isZeroAndPi(res)

    for i, lng in enumerate(eLng):
        starUnitVec = obs.computeStarUnitVector(lng, 45)
        res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        assert np.allclose(res, expected[i]), "%i %s" %(i, res)



#
# THIS TEST TURNED OFF
#

def test_varyElngZeta80():
    alpha_rad = 0
    elat_deg = 45
    maxZenithAngle_rad = np.radians(80)
    func = obs.computeExtremaOfOrbitalPhase_rad

    eLng = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    expected = [
                [14.2158, 165.7841],  #0
                [140.74650147, 346.12339618],  #30
                [128.40829078047892, 329.8049199212592],  #60
#
                [125, 325],  #90
                [128.40829078047892, 329.8049199212592],  #120
                [140.74650147, 346.12339618],  #150

                [14.2158, 165.7841],  #180
                [39.25349853102525, 193.87660382313072], #210
                [51.591709, 210.195], #240

                [55.000, 215],  #270
                [51.591709, 210.195], #300
                [39.25349853102525, 193.87660382313072], #330
               ]


    for i, lng in enumerate(eLng):
        starUnitVec = obs.computeStarUnitVector(lng, elat_deg)
        res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        assert np.allclose(res, expected[i]), "%i %s %s" %(i, res, expected[i])


def test_varyAlphaElng0Elat45Zeta90():
    """Test different earth orbital phase for zeta=90"""
    maxZenithAngle_rad = np.radians(90)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = obs.computeStarUnitVector(0, 45)

    alpha_deg = [0, 30, 60, 90, 120, 150, 210, 240, 270, 300, 330]
    alpha_rad = np.radians(alpha_deg)
    expected = [
                [0, 180],
                [26.565051, 206.565051],
                [40.893395, 220.893395],

                [45.000000, 225.000000], #90
                [40.893395, 220.893395],
                [26.565051, 206.565051],

                #180 is done outside the loop
                [153.434949, 333.434949],
                [139.106605, 319.106605],

                [135, 315],  #270
                [139.106605, 319.106605],
                [153.434949, 333.434949],
               ]

    assert len(alpha_rad) == len(expected)  #Sanity check
    #    debug()
    for i in range(len(alpha_rad)):
        res = func(alpha_rad[i], starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        msg = "alpha=%.1f Expected %s Computed %s" %(alpha_deg[i], expected[i], res)
        assert np.allclose(res, expected[i], atol=1e-3), msg

    res = func(np.pi, starUnitVec, maxZenithAngle_rad)
    assert isZeroAndPi(res)


def test_varyAlphaElng0Elat45Zeta60():
    """Test different earth orbital phase for zeta=60"""

    maxZenithAngle_rad = np.radians(60)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = obs.computeStarUnitVector(0, 45)

    alpha_deg = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    alpha_rad = np.radians(alpha_deg)
    expected = [
                [45.000000, 135.000000],
                [65.796572, 167.333531],
                [73.204928, 188.581861],

                [75.000000, 195.000000], #90
                [73.204928, 188.581861],
                [65.796572, 167.333531],

                [45.000000, 135.000000], #180
                [12.666469, 114.203428],
                [106.795072, 351.418139],

                [105.000000, 345.000000], #270
                [106.795072, 351.418139],
                [12.666469, 114.203428],
               ]

    assert len(alpha_rad) == len(expected)  #Sanity check
    #    debug()
    for i in range(len(alpha_rad)):
        res = func(alpha_rad[i], starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        msg = "alpha=%.1f Expected %s Computed %s" %(alpha_deg[i], expected[i], res)
        assert np.allclose(res, expected[i], atol=1e-3), msg

    res = func(np.pi, starUnitVec, maxZenithAngle_rad)
    assert isZeroAndPi(res)


def test_CardinalPointsAlpha0Zeta45():
    """Increase epsilon, earth avoidance angle"""

    alpha_rad = 0
    maxZenithAngle_rad = np.radians(45)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = np.array([0,1,0])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(np.sort(res), [45, 315]), res

    starUnitVec = np.array([0,-1,0])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(np.sort(res), [135, 225]), res

    starUnitVec = np.array([0,0,1])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(res, [45, 135]), res

    starUnitVec = np.array([0,0,-1])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(np.sort(res), [225, 315]), res

#
#def test_varyAlphaElng0Elat45Zeta80():
#    """Increase epsilon, earth avoidance angle"""
#
#    maxZenithAngle_rad = np.radians(80)
#    func = obs.computeExtremaOfOrbitalPhase_rad
#    starUnitVec = obs.computeStarUnitVector(0, 45)
#    print(starUnitVec)
#
#    nTest = 12
#    alpha_deg = 30* np.arange(nTest)
#    alpha_rad = np.radians(alpha_deg)
#
#    #These values are wrong, and shoudl be computed
#    expected = [ [14.215, 165.7841],
#                 [140.746, 346.123],
#                 [128.408, 329.804],
#
#                 [125, 325],
#                 [128.408, 329.804],
#                 [140.746, 346.123],
#
#                 [14.215, 165.7841],
#                 [39.253, 193.876],
#                 [51.591, 210.195],
#
#                 [55, 215],
#                 [51.591, 210.195],
#                 [39.253, 193.876],
#               ]
#
#    for i in range(nTest):
#        res = func(alpha_rad[i], starUnitVec, maxZenithAngle_rad)
#        res = np.degrees(res)
#        msg = "i=%i alpha=%.1f Expected %s Computed %s" %(i, alpha_deg[i], expected[i], res)
#        assert np.allclose(res, expected[i], atol=1e-3), msg
#
