
local N = 16

set_domain({0,0,0}, {1,1,1}, {N,N,N}, 5, 2)
init_prim(function(x,y,z)
             return { x*x + y*y + z*z, 1.0, x, y, z }
          end)

local prim = get_prim()
local rho_0 = prim.rho

local start = os.time()
local rho_k = fft_forward(rho_0)
local rho_x = fft_reverse(rho_k)

print("length is", #rho_x)


local vx = prim.vx
local vy = prim.vy
local vz = prim.vz

local start = os.time()
prim.vx, prim.vy, prim.vz = fft_helmholtz(vx, vy, vz)
init_prim(prim)




local hid = 0
if mpi_get_rank() == 0 then
   h5_open_file("pspec.h5", "w")
   hid = h5_open_group("pspec", "w")
end

fft_power_vector_field(vx, vy, vz, hid, "velocity")
fft_power_scalar_field(prim.rho, hid, "density")

if mpi_get_rank() == 0 then
   h5_close_file()
end

local host = require 'host'
write_prim("pspec.h5", host.CheckpointOptions)
