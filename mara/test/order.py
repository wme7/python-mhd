#!/usr/bin/env python


if __name__ == "__main__":
    import numpy as np
    from scipy.optimize import leastsq
    from matplotlib import pyplot as plt
    from optparse import OptionParser

    parser = OptionParser()
    opts, args = parser.parse_args()
    data = np.loadtxt(args[0])

    def errfunc(v, x, y):
        return (v[0] - v[1]*np.log10(x)) - np.log10(y)

    x = data[:,0]
    y = data[:,1]
    v0 = [0.0, -2.0]
    v, success = leastsq(errfunc, v0, args=(x,y), maxfev=10000)

    plt.loglog(x, y, '-o', label=r"%s: order=$%3.2f$" % (args[0], v[1]))
    plt.legend()
    plt.show()

