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



def test_varyElngZeta80():
    alpha_rad = 0
    maxZenithAngle_rad = np.radians(80)
    func = obs.computeExtremaOfOrbitalPhase_rad

    eLng = [0, 30] #, 60, 90, 120, 150, 210, 240, 270, 300, 330]
    expected = [
                [14.2158, 165.7841],
                [140, 347],
#                [135, 315],
#
#                [139.106, 319.107],
#                [153.435, 333.435],
#                #180 is not tested in the loop
#
#                [26.565, 206.565],
#                [40.893, 220.893],
#                [45, 225],
#
#                [40.893, 220.893],
#                [26.565, 206.565],
               ]


    for i, lng in enumerate(eLng):
        starUnitVec = obs.computeStarUnitVector(lng, 45)
        res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
        res = np.degrees(res)
        assert np.allclose(res, expected[i]), "%i %s" %(i, res)


#def test_varyAlphaElng0Elat45Zeta90():
#    """Increase epsilon, earth avoidance angle"""
#
#    maxZenithAngle_rad = np.radians(90)
#    func = obs.computeExtremaOfOrbitalPhase_rad
#
#    starUnitVec = obs.computeStarUnitVector(0, 45)
#
#    nTest = 12
#    alpha_deg = 30* np.arange(nTest)
#    alpha_rad = np.radians(alpha_deg)
#
#
#    for i in range(nTest):
#        res = func(alpha_rad[i], starUnitVec, maxZenithAngle_rad)
#        res = np.degrees(res)
#        msg = "alpha=%.1f Expected %s Computed %s" %(alpha_deg[i], expected[i], res)
#        assert np.allclose(res, expected[i], atol=1e-3), msg



def test_CardinalPointsAlpha0Zeta45():
    """Increase epsilon, earth avoidance angle"""

    alpha_rad = 0
    maxZenithAngle_rad = np.radians(45)
    func = obs.computeExtremaOfOrbitalPhase_rad

    starUnitVec = np.array([0,1,0])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(np.sort(res), [135, 225]), res

    starUnitVec = np.array([0,-1,0])
    res = func(alpha_rad, starUnitVec, maxZenithAngle_rad)
    res = np.degrees(res)
    assert np.allclose(np.sort(res), [45, 315]), res

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
