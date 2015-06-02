
host = require 'host'

TestFile = "/Users/jzrake/Work/drvtrb/data/test/chkpt.0001.h5"

h5_open_file(TestFile, "r+")
Nx = h5_get_ndims("prim/rho")
Nq = h5_get_nsets("prim")
set_domain({0,0,0}, {1,1,1}, Nx, Nq, 3)
h5_close_file()

set_fluid("rmhd")
print_mara()
read_prim(TestFile, host.CheckpointOptions)


print("mach:", measure_mean_max_sonic_mach())
print("divB:", measure_mean_max_divB())
