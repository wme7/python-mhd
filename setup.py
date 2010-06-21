#!/usr/bin/env python

if __name__ == "__main__":
    from numpy.distutils.core import setup, Extension
    from hydro.config import host_config

    cflags = host_config['cflags']

    src_scalar = ['integrate.c', 'hll.c', 'scalar.c']
    src_euler  = ['integrate.c', 'hll.c', 'euler.c']
    src_srhd   = ['integrate.c', 'hll.c', 'srhd.c']
    src_rmhd   = ['integrate.c', 'hll.c', 'hllc_rmhd.c', 'rmhd.c', 'quartic.c']

    ext_scalar = Extension(name               = 'scalar',
                           sources            = ['src/'+s for s in src_scalar],
                           extra_compile_args = cflags)
    ext_euler  = Extension(name               = 'euler',
                           sources            = ['src/'+s for s in src_euler],
                           extra_compile_args = cflags)
    ext_srhd   = Extension(name               = 'srhd',
                           sources            = ['src/'+s for s in src_srhd],
                           extra_compile_args = cflags)
    ext_rmhd   = Extension(name               = 'rmhd',
                           sources            = ['src/'+s for s in src_rmhd],
                           extra_compile_args = cflags)
    setup(name='python-mhd', ext_modules=[ext_scalar, ext_euler, ext_rmhd, ext_srhd])
