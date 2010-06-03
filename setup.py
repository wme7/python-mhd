#!/usr/bin/env python



if __name__ == "__main__":

    from numpy.distutils.core import setup, Extension
    from os import listdir

    srcrmhd = ['src/'+c for c in listdir('src') if c.endswith('.c')]
    librmhd = Extension('librmhd', srcrmhd)
    sr_riemann = Extension('sr_riemann', ['src/sr_riemann.f', 'src/sr_riemann-vt.f'])

    setup(name='rmhd', ext_modules=[librmhd])
