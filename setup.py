#!/usr/bin/env python



if __name__ == "__main__":

    from numpy.distutils.core import setup, Extension

    srceulr = ['src/integrate.c', 'src/euler.c', 'src/hll.c']

    ext_euler = Extension(name               = 'euler',
                          sources            = srceulr,
                          define_macros      = [('NQ',5)],
                          extra_compile_args = ["-Wall", "-O3"])

    setup(name='python-mhd', ext_modules=[ext_euler])
