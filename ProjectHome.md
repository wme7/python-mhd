# Relativistic Magnetohydrodynamics in Python #

## Summary ##
Python MHD is a collection of fast C-libraries for solving the equations of Relativistic MHD on fixed grids using
different classes of finite volume high-order Godunov techniques. The libraries are written as standalone C-code, and may be driven by the provided Python modules, or extracted for used in other applications.

![http://python-mhd.googlecode.com/files/KH2d-noB.png](http://python-mhd.googlecode.com/files/KH2d-noB.png)

## Requirements ##
The back-end C-library compiles without any external dependencies at all, just a standard C compiler. An exact riemann solver for Special Relativistic hydrodynamics is provided for comparison testing, and may be compiled optionally if a Fortran compiler is available. The Python bindings require the use of [numpy](http://numpy.scipy.org/) and [pylab](http://matplotlib.sourceforge.net/), and optionally [mpi4py](http://code.google.com/p/mpi4py/). These libraries are all provided in the [Enthought Python Distribution](http://www.enthought.com/products/epd.php), but may be installed individually.


---

## Usage ##
The following snippet of code illustrates how to run a simple 2 dimensional test problem from Python. The code will generate the image shown below, using the [matplotlib](http://matplotlib.sourceforge.net/) libraries.
```
def blast_wave_example():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(132,132), L=(2,2))
    problem = RMHDCylindricalA(pre=100.0)
    state = LibraryState(plm_theta=2.0, mode_reconstruct=2, mode_riemann_solver=1)

    run_args = {'name': "cylindrical blast wave",
                'RK_order': 3, 'CFL': 0.6, 'tfinal': 0.2}

    P = driver.run(_lib, state, problem, **run_args)

    visual.four_pane_2d(P, extent=[-1,1,-1,1])
    visual.show()
```


---

## Test Cases ##
#### Cylindrical blast wave, strong magnetic field ####
![http://python-mhd.googlecode.com/files/cyl-blast.png](http://python-mhd.googlecode.com/files/cyl-blast.png)

#### Quadrant problems ####
These test setups are meant to indicate the presence of asymmetries in the integration algorithms. The exact symmetry across the diagonal of the contour lines indicates that the code is behaving properly.
![http://python-mhd.googlecode.com/files/QuadrantA.png](http://python-mhd.googlecode.com/files/QuadrantA.png)
![http://python-mhd.googlecode.com/files/QuadrantB.png](http://python-mhd.googlecode.com/files/QuadrantB.png)