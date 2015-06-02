#!/usr/bin/env python

# ------------------------------------------------------------------------------
#
# Authors: Jonathan Zrake, Bez Laderman: NYU CCPP
# Date: May 7th, 2012
#
# This piece of code implements the left and right eigenvectors of the ideal
# special relativistic hydrodynamics equations. The formulas are an exact
# translation of those given in the literature:
#
# R. Donat, J.A. Font, J.M. Ibanez, & A. Marquina
# JCP, 1998, 146, 58
#
# http://www.sciencedirect.com/science/article/pii/S0021999198959551
#
#
# Having these eigenvectors in a hydrodynamics code is good. They can be used
# for any scheme which requires flux splitting with characteristic
# decomposition, such as high order ENO or WENO schemes.
#
# ------------------------------------------------------------------------------

from math import sqrt
import numpy as np

Gamma = 1.4 # adiabatic index

def run_evec():

    D = 1.0 # rest mass density
    p = 1.0 # pressure
    u = 0.2 # vx (three velocity)
    v = 0.3 # vy
    w = 0.4 # vz

    sie = (p/D) / (Gamma - 1) # specific internal energy
    h = 1 + sie + p/D         # specific enthalpy
    cs2 = Gamma * p / (D*h)   # sound speed squared
    V2 = u*u + v*v + w*w
    W = 1.0 / sqrt(1 - V2)    # Lorentz factor
    W2 = W*W
    K = h                     # for gamma-law only, K = h
    hW = h*W

    # equations (14) and (15)
    lp = (u*(1-cs2) + sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2)
    lm = (u*(1-cs2) - sqrt(cs2*(1-V2)*(1-V2*cs2-u*u*(1-cs2))))/(1-V2*cs2)

    Ap = (1 - u*u) / (1 - u*lp)
    Am = (1 - u*u) / (1 - u*lm)

    # Equations (17) through (20)
    # --------------------------------------------------------------------------
    RR = [[K/hW, u, v, w, 1-K/hW],
          [W*v, 2*h*W2*u*v, h*(1+2*W2*v*v), 2*h*W2*v*w, 2*h*W2*v - W*v],
          [W*w, 2*h*W2*u*w, 2*h*W2*v*w, h*(1+2*W2*w*w), 2*h*W2*w - W*w],
          [1, hW*Ap*lp, hW*v, hW*w, hW*Ap - 1],
          [1, hW*Am*lm, hW*v, hW*w, hW*Am - 1]]

    # NOTES
    # --------------------------------------------------------------------------
    # (1) Font writes the columns of the left eigenvector matrix
    # horizontally, which is how they are written below. So we take the
    # transpose at the end of the day.
    #
    # (2) Font's notation uses L_{-/+} for the last left eigenvectors, but that
    # naming is weird, since L_{+} contains lm and Am and vice-versa.
    # --------------------------------------------------------------------------

    Delta = h*h*h*W*(K-1)*(1-u*u)*(Ap*lp - Am*lm) # equation (21)
    a = W / (K-1)
    b = 1 / (h*(1 - u*u))
    c = 1 / (h*(1 - u*u))
    d = -h*h / Delta
    e = +h*h / Delta

    LL = [[a*(h-W), a*W*u, a*W*v, a*W*w, -a*W],
          [-b*v, b*u*v, b*(1-u*u), 0, -b*v],
          [-c*w, c*u*w, 0, c*(1-u*u), -c*w],
          [d*(hW*Am*(u-lm) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm),
           d*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Am) - K*Am),
           d*(W2*v*(2*K - 1)*Am*(u - lm)),
           d*(W2*w*(2*K - 1)*Am*(u - lm)),
           d*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Am*lm) + K*Am*lm)],
          [e*(hW*Ap*(u-lp) - u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp),
           e*(1 + W2*(V2 - u*u)*(2*K - 1)*(1 - Ap) - K*Ap),
           e*(W2*v*(2*K - 1)*Ap*(u - lp)),
           e*(W2*w*(2*K - 1)*Ap*(u - lp)),
           e*(-u - W2*(V2 - u*u)*(2*K - 1)*(u - Ap*lp) + K*Ap*lp)]]

    R = np.matrix(RR).T
    L = np.matrix(LL)

    Ry = np.matrix(np.zeros((5,5)))
    Ly = np.matrix(np.zeros((5,5)))

    Rz = np.matrix(np.zeros((5,5)))
    Lz = np.matrix(np.zeros((5,5)))

    Ry[:,0] = R[:,1]
    Ry[:,1] = R[:,2]
    Ry[:,2] = R[:,0]
    Ry[:,3] = R[:,3]
    Ry[:,4] = R[:,4]

    Ly[0,:] = L[1,:]
    Ly[1,:] = L[2,:]
    Ly[2,:] = L[0,:]
    Ly[3,:] = L[3,:]
    Ly[4,:] = L[4,:]

    print np.around(Ry*Ly, 14)

    print np.around(R*L, 14) # ignore values near machine precision


if __name__ == "__main__":
    run_evec()
