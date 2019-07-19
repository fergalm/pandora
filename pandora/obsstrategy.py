# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Jun 25 20:55:47 2019

Everything is computed in eclipitic plane coordinates.

z-hat points toward ecliptic pole
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import numpy as np



def computeSolarAngleFromDoy(doy):
    pass

def computeEarthUnitVector(alpha_rad):
    """Compute vector from sun to Earth given earth orbit angle, alpha"""

    return np.array([np.cos(alpha_rad), np.sin(alpha_rad), 0])


def computeTelescopeUnitVector(alpha_rad, rho_rad):
    """

    alpha is earth orbital phase  angle,
    rho is telescope orbital phase angle
    """

    x = -np.sin(alpha_rad) * np.cos(rho_rad)
    y = np.cos(alpha_rad) * np.cos(rho_rad)
    z = np.sin(rho_rad)

    return np.array([x, y, z])



def computeStarUnitVector(ecliptic_lng_deg, ecliptic_lat_deg):

    lng_rad, lat_rad = np.radians([ecliptic_lng_deg, ecliptic_lat_deg])
    x = np.cos(lat_rad) * np.cos(lng_rad)
    y = np.cos(lat_rad) * np.sin(lng_rad)
    z = np.sin(lat_rad)

    return np.array([x, y, z])



def computeDutyCycle(alpha_rad, starUnitVec, maxZenithAngle_rad):

    #@TODO, make this a param
    maxSunEarthAngle_rad = np.radians(80)
    earthUnitVec = computeEarthUnitVector(alpha_rad)
#    debug()
#    if np.dot(earthUnitVec, starUnitVec) < np.cos(maxSunEarthAngle_rad):

    try:
        angles_rad = computeExtremaOfOrbitalPhase_rad(alpha_rad,
                                                      starUnitVec,
                                                      maxZenithAngle_rad)
    except ValueError:
        #Star is never observable
        return 0

    #Compute fraction of circle subtended by angles
    twopi = 2 * np.pi
    dutyCycle = np.diff(np.sort(angles_rad))[0] / twopi
    assert 0 <= dutyCycle and dutyCycle <= 1

    cos_zeta = np.cos(maxZenithAngle_rad)
    if dutyCycle > .5 and cos_zeta > 0:
        dutyCycle = 1 - dutyCycle

    return dutyCycle


def computeExtremaOfOrbitalPhase_rad(alpha_rad, starUnitVec, maxZenithAngle_rad):
    """Assume star is visible if t.s > 0


    TODO
    Write this up in latex

    http://mathforum.org/kb/servlet/JiveServlet/download/206-2697091-9755226-1202710/sin-cos-eqs.pdf
    Returns
    --------
    Orbital phases (0..2pi) where star becomes visible/unobservable
    """

    svec = starUnitVec
    term1 =  (svec[1] * np.cos(alpha_rad))
    term1 -= (svec[0] * np.sin(alpha_rad))
#    term1 = -(svec[1] * np.sin(alpha_rad))
#    term1 += (svec[0] * np.cos(alpha_rad))
    term2 = svec[2]
    term3 = np.cos(maxZenithAngle_rad)

    sin_rho = computeSineRho(term1, term2, term3)
    cos_rho = computeCosineRho(term1, term2, term3)

#    debug()
    method = 'default2pi'

    if method == 'default':
        rho_rad = getRhoFromSinCos(sin_rho, cos_rho, alpha_rad, svec, maxZenithAngle_rad)

        if np.all(rho_rad < 0):
            rho_rad += 2*np.pi
        elif np.any(rho_rad < 0):
            rho_rad += np.pi
    elif method == 'default2pi':
        rho_rad = getRhoFromSinCos(sin_rho, cos_rho, alpha_rad, svec, maxZenithAngle_rad)

        twopi = 2 * np.pi
        if rho_rad[0] < 0:
            rho_rad[0] += twopi
        if rho_rad[1] < 0:
            rho_rad[1] += twopi

    elif method == 'sineonly':
        rho_rad = np.zeros(2)
        for i in range(2):
            rho_rad[i] = np.arcsin(sin_rho[i])

            if not doesRhoSatisfyEqn(alpha_rad, rho_rad[i], svec, maxZenithAngle_rad):
                    rho_rad[i] = np.pi - rho_rad[i]

            if rho_rad[i] < 0:
                rho_rad[i] += 2 * np.pi
    else:
        assert False, "Unknown trig resolution method"


    return np.sort(rho_rad)


def getRhoFromSinCos(sin_rho, cos_rho, alpha_rad, svec, maxZenithAngle_rad):
    """The ugly, messy part of the code"""
    rho_rad = np.arctan2(sin_rho, cos_rho)

    if not doesRhoSatisfyEqn(alpha_rad, rho_rad[0], svec, maxZenithAngle_rad):
        rho_rad = np.arctan2(sin_rho, cos_rho[::-1])

#    if not doesRhoSatisfyEqn(alpha_rad, rho_rad[0], svec, maxZenithAngle_rad):
#            rho_rad = np.arctan2(sin_rho, -cos_rho[::-1])
##
#    if not doesRhoSatisfyEqn(alpha_rad, rho_rad[0], svec, maxZenithAngle_rad):
#            rho_rad = np.arctan2(-cos_rho, sin_rho[::-1])

    assert doesRhoSatisfyEqn(alpha_rad, rho_rad[0], svec, maxZenithAngle_rad), rho_rad

    return rho_rad


def doesRhoSatisfyEqn(alpha_rad, rho_rad, starUnitVec, maxZenithAngle_rad):
    tvec = computeTelescopeUnitVector(alpha_rad, rho_rad)

#    print( np.dot(tvec, starUnitVec))
#    print(np.cos(maxZenithAngle_rad))
    opt1 = np.allclose(np.dot(tvec, starUnitVec),
                       np.cos(maxZenithAngle_rad), 1e-6)

    opt2 = np.allclose(np.dot(tvec, starUnitVec),
                       np.cos(np.pi + maxZenithAngle_rad), 1e-6)

    return opt1 or opt2

def computeSineRho(term1, term2, term3):
    """Compute sin of orbital phase angle corresponding to
    star entering, or leaving, active part of duty cycle
    (where it can be observed

    terms are computed elsewhere
    """
    A = term1**2 + term2**2
    B = -2 * term2 * term3
    C = term3**2 - term1**2

    #Compute roots of quadratic equation. Result sometimes has some
    #non-zero imaginary component due to round off error. We just trucate
    #this away
    sin_rho = np.roots([A,B,C])
    sin_rho = np.real(sin_rho) #Debugging code

    assert len(sin_rho) == 2
    assert np.all(np.isreal(sin_rho)), sin_rho
    assert np.all(np.fabs(sin_rho) <= 1), sin_rho
    return sin_rho


def computeCosineRho(term1, term2, term3):
    """Compute sin of orbital phase angle corresponding to
    star entering, or leaving, active part of duty cycle
    (where it can be observed

    terms are computed elsewhere
    """
    A = term1**2 + term2**2
    B = -2 * term1 * term3
    C = term3**2 - term2**2

    #Compute roots of quadratic equation. Result sometimes has some
    #non-zero imaginary component due to round off error. We just trucate
    #this away
    cos_rho = np.roots([A,B,C])
    cos_rho = np.real(cos_rho)

    assert len(cos_rho) == 2
    assert np.all(np.isreal(cos_rho)), cos_rho
    assert np.all(np.fabs(cos_rho) <= 1), cos_rho
    return cos_rho


def plotEarthAngleCurve(alpha_rad, eLng_deg, eLat_deg, maxOffZenith_deg):
    """

    maxOffZenith_deg
        How far away from zenith can we still see the star
        """

    twopi = 2 * np.pi
#    earthAvoidanceAngle_rad = np.radians(earthAvoidanceAngle_deg)
    maxOffZenith_rad = np.radians(maxOffZenith_deg)
    svec = computeStarUnitVector(eLng_deg, eLat_deg)

    dotprod = []
    phase = np.linspace(0, twopi, 37)
    for rho_rad in phase:
        tvec = computeTelescopeUnitVector(alpha_rad, rho_rad)
        dotprod.append(np.dot(tvec, svec))

    angle_from_earth_deg = np.degrees( np.arccos(dotprod) )
    phase_deg = phase * 180 / np.pi
    label = "PLACEHOLDER"
    plt.plot(phase_deg, angle_from_earth_deg, label=label)

    plt.fill_between(phase_deg, 0, maxOffZenith_deg, alpha=.3, color='g')
    plt.fill_between(phase_deg, maxOffZenith_deg, 180, alpha=.3, color='r')
    plt.pause(.01)

    x1, x2 = computeExtremaOfOrbitalPhase_rad(alpha_rad, svec, maxOffZenith_rad)
    if x1 < 0:
        x1 += twopi
    if x2 < 0:
        x2 += twopi

    x1, x2 = np.degrees([x1, x2])
    plt.axvline(x1)
    plt.axvline(x2, ls=":")

    alpha_rad = 0
    dotprod=[]
    for rho_rad in phase:
        tvec = computeTelescopeUnitVector(alpha_rad, rho_rad)
        dotprod.append(np.dot(tvec, svec))

    angle_from_earth_deg = np.degrees( np.arccos(dotprod) )
    phase_deg = phase * 180 / np.pi
    plt.plot(phase_deg, angle_from_earth_deg, alpha=.4, label="$\alpha=0$", color='C1')
    x1, x2 = computeExtremaOfOrbitalPhase_rad(alpha_rad, svec, maxOffZenith_rad)
    if x1 < 0:
        x1 += np.pi
    x1, x2 = np.degrees([x1, x2])
    plt.axvline(x1, alpha=.4, color='C1')
    plt.axvline(x2, alpha=.4, color='C1')
    print(x1, x2)


#    plt.axvline(x1 + 90, ls="--")
#    plt.axvline( (x2 + 90) % 360, ls="--")

#    print(x1, x2)
#    if np.allclose([x2], [360]):
#        plt.axvline(0)


    plt.ylim(0, 180)
    plt.ylabel("Zenith Angle (deg) \nlower is better")
    plt.xlabel("Telescope orbital phase")

def main():
    alpha_rad = 0

    plt.clf()
    twopi = 2 * np.pi
    phase = np.linspace(0, twopi, 64)
    phase_deg = phase * 180 / np.pi
    plt.fill_between(phase_deg, 0, 105, alpha=.3)

    elng_deg = 0
    for elat in [0, 16, 30, 45, 60, 75, 90]:
        svec = computeStarUnitVector(elng_deg, elat)

        dotprod = []
        for rho_rad in phase:
            tvec = computeTelescopeUnitVector(alpha_rad, rho_rad)
#            print(np.degrees(rho_rad), np.dot(tvec, svec))
            dotprod.append(np.dot(tvec, svec))

        angle_from_earth_deg = np.degrees( np.arccos(dotprod) )
        plt.plot(phase_deg, angle_from_earth_deg, label="Eclip lat=%.0f" %(elat))

#        print(elng_deg, elat, computeVisibilityV1(alpha_rad, svec, 94))
    plt.legend()

    plt.xlabel("Telescope orbital phase")
    plt.ylabel("Angle between target and Earth")
    plt.title("For ecliptic longitude of %.1f deg" %(elng_deg))



def isind(sine):
    """Compute inverse sine of argument in degrees"""
    return np.degrees(np.arcsin(sine))

def icosd(cosine):
    """Compute inverse sine of argument in degrees"""
    return np.degrees(np.arcsin(cosine))