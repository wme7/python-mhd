#!/usr/bin/env python



if __name__ == "__main__":

    from numpy.distutils.core import setup, Extension

    srceuler = ['integrate.c', 'hll.c', 'euler.c']
    srcrmhd  = ['integrate.c', 'hll.c',  'rmhd.c', 'quartic.c']

    ext_euler = Extension(name               = 'euler',
                          sources            = ['src/'+s for s in srceuler],
                          extra_compile_args = ["-Wall", "-O3"])

    ext_rmhd  = Extension(name               = 'rmhd',
                          sources            = ['src/'+s for s in srcrmhd],
                          extra_compile_args = ["-Wall", "-O3"])

    setup(name='python-mhd', ext_modules=[ext_euler, ext_rmhd])

