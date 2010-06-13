#!/usr/bin/env python



if __name__ == "__main__":

    from numpy.distutils.core import setup, Extension

    srcscalar = ['integrate.c', 'hll.c', 'scalar.c']
    srceuler  = ['integrate.c', 'hll.c', 'euler.c']
    srcsrhd   = ['integrate.c', 'hll.c',  'srhd.c']
    srcrmhd   = ['integrate.c', 'hll.c',  'rmhd.c', 'quartic.c']

    ext_scalar = Extension(name               = 'scalar',
                           sources            = ['src/'+s for s in srcscalar],
                           extra_compile_args = ["-Wall", "-O3"])

    ext_euler  = Extension(name               = 'euler',
                           sources            = ['src/'+s for s in srceuler],
                           extra_compile_args = ["-Wall", "-O3"])

    ext_srhd   = Extension(name               = 'srhd',
                           sources            = ['src/'+s for s in srcsrhd],
                           extra_compile_args = ["-Wall", "-O3"])

    ext_rmhd   = Extension(name               = 'rmhd',
                           sources            = ['src/'+s for s in srcrmhd],
                           extra_compile_args = ["-Wall", "-O3"])

    setup(name='python-mhd', ext_modules=[ext_scalar, ext_euler, ext_rmhd, ext_srhd])

