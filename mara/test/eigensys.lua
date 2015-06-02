

local matrix = require 'matrix'
local util = require 'util'


local Direction = 3
local PrimState = { 1.0, 4.0, 0.1, 0.1, 0.1 }
local FluidType = "srhd"



-- .............................................................................
-- Flux Jacobian for SRHD computed numerically by differencing the fluxes.
-- .............................................................................
local function JacobianSrhd(U, n)
   local eps = 1e-12
   local J = { }

   for i=0,5-1 do
      J[i+1] = { }
      for j=0,5-1 do

         local U0 = U:copy()
         local U1 = U:copy()

         U0[j] = U0[j] - eps
         U1[j] = U1[j] + eps

         local P0 = fluid.ConsToPrim(U0)
         local P1 = fluid.ConsToPrim(U1)

         local F0 = fluid.FluxFunction(P0, 1)
         local F1 = fluid.FluxFunction(P1, 1)

         local G0 = fluid.FluxFunction(P0, 2)
         local G1 = fluid.FluxFunction(P1, 2)

         local H0 = fluid.FluxFunction(P0, 3)
         local H1 = fluid.FluxFunction(P1, 3)

	 local Fn0 = n[1]*F0 + n[2]*G0 + n[3]*H0
	 local Fn1 = n[1]*F1 + n[2]*G1 + n[3]*H1

         J[i+1][j+1] = (Fn1[i] - Fn0[i]) / (U1[j] - U0[j])
      end
   end
   return matrix(J)
end


-- .............................................................................
-- Flux Jacobian as given by http://www.microcfd.com/download/pdf/aiaa2001-2609.pdf
-- output is A := dF{D,E,px,py,pz}/d{D,E,px,py,pz}, Mara's convention
-- .............................................................................
local function JacobianEuler(U, n) -- n must be a unit vector

   local gm = 1.4
   local g1 = gm - 1.0
   local g2 = gm - 2.0

   local D = U[0] -- inputs are Mara convention: {D,E,px,py,pz}
   local E = U[1]
   local u = U[2]/D
   local v = U[3]/D
   local w = U[4]/D

   local nx, ny, nz = n[1], n[2], n[3]

   local vn = u*nx + v*ny + w*nz
   local ek = 0.5*(u^2 + v^2 + w^2)
   local p0 = (E-D*ek)*g1 -- pressure
   local a2 = gm*p0/D     -- sound speed
   local h0 = a2 / g1 + ek

   local A = { { 0, 0, nx, ny, nz },
               { (g1*ek-h0)*vn,
                 gm*vn,
                 h0*nx - g1*u*vn,
                 h0*ny - g1*v*vn,
                 h0*nz - g1*w*vn },
               { g1*ek*nx - u*vn,
                 g1*nx,
                 1*vn - g2*u*nx,
                 u*ny - g1*v*nx,
                 u*nz - g1*w*nx },
               { g1*ek*ny - v*vn,
                 g1*ny,
                 v*nx - g1*u*ny,
                 1*vn - g2*v*ny,
                 v*nz - g1*w*ny },
               { g1*ek*nz - w*vn,
                 g1*nz,
                 w*nx - g1*u*nz,
                 w*ny - g1*v*nz,
                 1*vn - g2*w*nz } }
   return matrix(A)
end


-- .............................................................................
-- Flux Jacobian as given by Toro (3.79)
-- output is A := dF{D,E,px,py,pz}/d{D,E,px,py,pz}, Mara's convention
-- .............................................................................
local function JacobianToro(U, n) -- n must be {1,0,0}
   local gm = 1.4
   local g1 = gm - 1.0
   local g2 = gm - 2.0

   local D = U[0] -- inputs are Mara convention: {D,E,px,py,pz}
   local E = U[1]
   local u = U[2]/D
   local v = U[3]/D
   local w = U[4]/D

   local ek = 0.5*(u^2 + v^2 + w^2)
   local p0 = (E-D*ek)*g1 -- pressure
   local a2 = gm*p0/D     -- sound speed
   local H = a2 / g1 + ek

   local A = { { 0, 0, 1, 0, 0 },
               { 0.5*u*(H*(gm-3) + g1*ek - a2), gm*u, H - g1*u^2, -g1*u*v, -g1*u*w },
               { g1*H - u^2 - a2, g1, (3-gm)*u, -g1*v, -g1*w },
               { -u*v, 0, v, u, 0 },
               { -u*w, 0, w, 0, u } }
   return matrix(A)
end


local function TestState(P, dim)
   local sep = "---------------------------------"
   local thresh = {euler=1e-12, srhd=1e-3}
   local function clip(x)
      if math.abs(x) < thresh[FluidType] then
         return 0.0
      else
         return x
      end
   end

   print("\n\n")
   print(sep)
   print("Testing state on dim:", dim)
   print(sep)

   Jacobians = { euler=JacobianEuler, srhd=JacobianSrhd }
   set_fluid(FluidType)

   local n = { 0,0,0 }
   n[dim] = 1.0

   local U = fluid.PrimToCons(P)
   local A = Jacobians[FluidType](U, n)
   local Lv, Rv, lm = fluid.Eigensystem(P, dim)

   local L = matrix(Lv)
   local R = matrix(Rv)

   print("\nR.L (expect I):")
   print(sep)
   print((R*L):replace(clip))

   print("\neigenvalues:")
   print(sep)
   print(lm)

   print("\nL.A.R (expect eigenvalues):")
   print(sep)
   print((L*A*R):replace(clip))
end


TestState(PrimState, Direction)
