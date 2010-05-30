#!/usr/bin/env python



if __name__ == "__main__":

    from numpy.distutils.core import setup, Extension
    from os import listdir

    librmhd1 = Extension('librmhd-1', ['src/rmhd-1.c'])
    librmhd2 = Extension('librmhd-2', ['src/rmhd-2.c', 'src/quartic.c'])
    sr_riemann = Extension('sr_riemann', ['src/sr_riemann.f', 'src/sr_riemann-vt.f'])

    setup(name='rmhd', ext_modules=[librmhd1, librmhd2, sr_riemann])

