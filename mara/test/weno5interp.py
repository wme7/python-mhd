#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq
from matplotlib import pyplot as plt



CeesA2C = [ [23./24.,  1./12.,  -1./24.],
            [-1./24., 13./12.,  -1./24.],
            [-1./24.,  1./12.,  23./24.] ]
CeesC2A = [ [25./24., -1./12.,   1./24.],
            [ 1./24., 11./12.,   1./24.],
            [ 1./24., -1./12.,  25./24.] ]
CeesC2L = [ [15./8., -5./4.,  3./8.],
            [ 3./8.,  3./4., -1./8.],
            [-1./8.,  3./4.,  3./8.] ]
CeesC2R = [ [ 3./8., 3./4.,  -1./8.],
            [-1./8., 3./4.,   3./8.],
            [ 3./8.,-5./4.,  15./8.] ]

CeesC2L_A = [ [11./6., -7./6.,  1./3. ],
              [ 1./3.,  5./6., -1./6. ],
              [-1./6.,  5./6.,  1./3. ] ]

CeesC2R_A = [ [ 1./3.,  5./6., -1./6. ],
              [-1./6.,  5./6.,  1./3. ],
              [ 1./3., -7./6., 11./6. ] ]

DeesC2L = [   1./ 16.,   5./  8.,   5./ 16. ]
DeesC2R = [   5./ 16.,   5./  8.,   1./ 16. ]
DeesA2C = [  -9./ 80.,  49./ 40.,  -9./ 80. ]
DeesC2A = [ -17./240., 137./120., -17./240. ]

DeesC2L_A = [ 0.1, 0.6, 0.3 ]
DeesC2R_A = [ 0.3, 0.6, 0.1 ]


def get_log_slope(x, y):
    def errfunc(v):
        return (v[0] - v[1]*np.log10(x)) - np.log10(y)
    v0 = [0.0, -2.0]
    v, success = leastsq(errfunc, v0)
    return v[1]

def SQU(x):
    return x**2

def weno5(v, c, d):
    eps = 1e-16
    B = [(13.0/12.0)*np.power(  v[ 2] - 2*v[ 3] +   v[ 4], 2.0) +
         ( 1.0/ 4.0)*np.power(3*v[ 2] - 4*v[ 3] +   v[ 4], 2.0),

         (13.0/12.0)*np.power(  v[ 1] - 2*v[ 2] +   v[ 3], 2.0) +
         ( 1.0/ 4.0)*np.power(3*v[ 1] - 0*v[ 2] -   v[ 3], 2.0),

         (13.0/12.0)*np.power(  v[ 0] - 2*v[ 1] +   v[ 2], 2.0) +
         ( 1.0/ 4.0)*np.power(  v[ 0] - 4*v[ 1] + 3*v[ 2], 2.0)]

    vs = [c[0][0]*v[ 2] + c[0][1]*v[ 3] + c[0][2]*v[4],
          c[1][0]*v[ 1] + c[1][1]*v[ 2] + c[1][2]*v[3],
          c[2][0]*v[ 0] + c[2][1]*v[ 1] + c[2][2]*v[2]]

    w = [d[0] / SQU(eps + B[0]),
         d[1] / SQU(eps + B[1]),
         d[2] / SQU(eps + B[2])]

    wtot = w[0] + w[1] + w[2]
    return (w[0]*vs[0] + w[1]*vs[1] + w[2]*vs[2])/wtot



def run_weno_C2L(N):
    def thefunc(x):
        return np.sin(8*np.pi*np.array(x))
    x = np.linspace(0.0, 1.0, N)
    y = thefunc(x)

    # interpolate to the point x[i+1/2]
    # choose i = N/2 for fun
    x_bar = [ ]
    y_bar = [ ]
    for i in range(N/8, 7*N/8):
        x_bar.append(0.5*(x[i-1] + x[i]))
        y_bar.append(weno5(y[i-2:i+3], CeesC2L, DeesC2L))
    return sum(abs(y_bar - thefunc(x_bar))) / N


def run_weno_deriv(N):
    def thefunc(x):
        return np.sin(8*np.pi*np.array(x))
    def thefunc_deriv(x):
        return np.cos(8*np.pi*np.array(x))*8*np.pi

    x, dx = np.linspace(0.0, 1.0, N, retstep=True)
    f = thefunc(x)

    # interpolate to the point x[i+1/2]
    # choose i = N/2 for fun
    x_bar = [ ]
    y_bar = [ ] # is f'(x)
    for i in range(N/8, 7*N/8):

        fiph = weno5(f[i-2:i+3], CeesC2R_A, DeesC2R_A)
        fimh = weno5(f[i-2:i+3], CeesC2L_A, DeesC2L_A)

        dfdx = (fiph - fimh) / dx

        x_bar.append(x[i])
        y_bar.append(dfdx)

    return sum(abs(y_bar - thefunc_deriv(x_bar))) * dx


def run_weno_C2R(N):
    def thefunc(x):
        return np.sin(8*np.pi*np.array(x))
    x = np.linspace(0.0, 1.0, N)
    y = thefunc(x)

    # interpolate to the point x[i+1/2]
    # choose i = N/2 for fun
    x_bar = [ ]
    y_bar = [ ]
    for i in range(N/8, 7*N/8):
        x_bar.append(0.5*(x[i+1] + x[i]))
        y_bar.append(weno5(y[i-2:i+3], CeesC2R_A, DeesC2R_A))
    return sum(abs(y_bar - thefunc(x_bar))) / N



def run_weno_C2A(N):
    def thefunc(x):
        return np.sin(8*np.pi*np.array(x))
    def thefunc_avg(x0, x1):
        return (np.cos(8*np.pi*x0) - np.cos(8*np.pi*x1)) / ((x1 - x0) * (8*np.pi))
    x = np.linspace(0.0, 1.0, N)
    y = thefunc(x)

    # interpolate to the point x[i+1/2]
    # choose i = N/2 for fun
    x_bar = [ ]
    y_bar = [ ]
    y_avg = [ ]
    for i in range(N/8, 7*N/8):
        x_bar.append(x[i])
        y_bar.append(weno5(y[i-2:i+3], CeesC2A, DeesC2A))
        y_avg.append(thefunc_avg(0.5*(x[i] + x[i-1]), 0.5*(x[i] + x[i+1])))
    return sum(abs(np.array(y_bar) - np.array(y_avg))) / N



def run_weno_A2C(N):
    """
    NOT ACTUALLY MADE!
    """
    def thefunc(x):
        return np.sin(8*np.pi*np.array(x))
    def thefunc_avg(x0, x1):
        return (np.cos(8*np.pi*x0) - np.cos(8*np.pi*x1)) / ((x1 - x0) * (8*np.pi))
    x = np.linspace(0.0, 1.0, N)
    y_avg = thefunc_avg(0.5*(x[i]+x[i-1]),0.5*(x[i]+x[i+1]))

    # interpolate to the point x[i+1/2]
    # choose i = N/2 for fun
    x_bar = [ ]
    y_bar = [ ]
    y = [ ]
    for i in range(N/8, 7*N/8):
        x_bar.append(x[i])
        y_bar.append(weno5(y_avg[i-2:i+3], CeesA2C, DeesA2C))
        y.append(thefunc(0.5*x))
    return sum(abs(np.array(y_bar) - np.array(y))) / N



if __name__ == "__main__":
    Ns = [2**n for n in range(4, 10)]
    res = [ ]
    for N in Ns:
        res.append(run_weno_deriv(N))

    order = get_log_slope(Ns, res)

    plt.loglog(Ns, res, '-o', label=r"order=$%3.2f$" % order)
    plt.legend()
    plt.show()

