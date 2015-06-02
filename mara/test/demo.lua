

local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'

local RunArgs = {
   N         = 64,
   dim       = 1,
   id        = "test",
   cpi       = -1.0,      -- interval between checkpoints (cpi < 0 => none)
   CFL       = 0.6,
   fixdt     = -1.0,      -- value for uniform time stepping
   tmax      = 0.2,       -- run simulation until
   eosfile   = "none",    -- a tabulated equation of state tabeos.h5
   fluid     = "euler",
   bound     = "outflow",
   advance   = "rk3",
   riemann   = "hllc",
   godunov   = "weno-split",

   -- --------------------------------------------------------------------------
   -- flags sent to config_solver
   -- --------------------------------------------------------------------------
   fsplit    = "llf",  -- one of [llf, marq]       ... flux splitting mode
   extrap    = "weno5",-- one of [pcm, plm, weno5] ... reconstruction type
   theta     = 2.0,    -- must be [0,2]            ... theta value for PLM/minmod
   IS        = "js96", -- one of [js96, b08, sz10] ... smoothness indicator
   sz10A     = 100.0,  -- should be in [0,100]     ... used by sz10 (see weno.c)

   -- --------------------------------------------------------------------------
   -- flags for particular initial conditions setups
   -- --------------------------------------------------------------------------
   eos       = "gamma-law",
   adgam     = 1.4,          -- adiabatic index if using gamma-law EOS
   quiet     = false,
   ic        = "Shocktube1", -- name of test a problem
   problem   = "shocktube",
   angle     = "{1,0,0}",    -- normal used for shock tubes and implosion tests

   -- --------------------------------------------------------------------------
   -- inline plotting options
   -- --------------------------------------------------------------------------
   pdf       = false,        -- write 1d figure to a pdf file instead of window
   noplot    = false,        -- suppress opening a plot window
   pltvar    = "rho,pre,vx,vy,vz",
}
math.randomseed(mpi_get_rank())
util.parse_args(RunArgs)
util.RunArgs = RunArgs
tests.RunArgs = RunArgs

-- This function may be used by problem setups if the domain is more complicated
-- than square:
local custom_domain = nil

local function HandleErrors(Status, attempt)
   return 0
end

local function plot_prim()
   if not RunArgs.noplot then
      local P = get_prim()
      --[[
      local pre = P.pre
      local rho = P.rho
      local W2 = 1.0 / (1.0 - (P.vx^2 + P.vy^2))
      local e = (pre/rho) * 1.0/(RunArgs.adgam - 1.0)
      local h = 1 + e + pre/rho
      local sx = rho * h * W2 * P.vx
      util.plot({sx=sx})
       ]]--
      local pltdict =  { }
      for k,v in pairs(util.string_split(RunArgs.pltvar, ",")) do
	 pltdict[v] = P[v]
      end
      if RunArgs.dim == 1 then
	 util.plot(pltdict)
      end
   end
end

local function cfg_mara()
   local N = RunArgs.N
   local NumberOfConserved = { rmhd=8, srhd=5, euler=5 }
   local Nq = NumberOfConserved[RunArgs.fluid]

   if custom_domain then
      custom_domain(Nq, 3)
   elseif RunArgs.dim == 1 then
      set_domain({0}, {1}, {N}, Nq, 3)
   elseif RunArgs.dim == 2 then
      set_domain({0,0}, {1,1}, {N,N}, Nq, 3)
   elseif RunArgs.dim == 3 then
      set_domain({0,0,0}, {1,1,1}, {N,N,N}, Nq, 3)
   else
      error("Invalid Dimension")
   end
   
   if type(RunArgs.bound) == "string" then
      set_boundary(RunArgs.bound)
   elseif type(RunArgs.bound) == "table" then
      set_boundary(unpack(RunArgs.bound))
   end
   config_solver(RunArgs, RunArgs.quiet)
   set_fluid(RunArgs.fluid)
   set_advance(RunArgs.advance)
   set_riemann(RunArgs.riemann)
   set_godunov(RunArgs.godunov)
   set_eos(RunArgs.eos,RunArgs.adgam)
end


local ProblemList = { }

function ProblemList.DensityWaveConvergenceRate()
   RunArgs.bound = "periodic"
   local res_values
   local L1_values = { }

   if RunArgs.dim == 1 then
      res_values = { 64, 128, 256, 512, 1024 }
   else
      res_values = { 16, 32, 64, 128}
   end

   print("------------------------------")
   print("N\t log10(L1)\t order")
   print("------------------------------")
   for n,N in pairs(res_values) do
      RunArgs.N = N
      local problem = tests.DensityWave
      problem.eps = 3.2e-1

      problem.velocity = { 0.1, 0.0, 0.0 }

      local status = util.run_simulation(problem:get_pinit(), cfg_mara, RunArgs)
      local P_comp = get_prim()

      init_prim(problem:get_pinit(status.CurrentTime))

      local P_true = get_prim()
      local diff = P_comp.rho - P_true.rho

      local L1 = 0.0
      for i=0,diff:size()-1 do
         L1 = L1 + math.abs(diff[i]) / N
      end

      L1_values[n] = math.log10(L1)
      local order = n == 1 and 0.0 or (
	 (L1_values[n] - L1_values[n-1]) /
	 (math.log10(res_values[n]) - math.log10(res_values[n-1])))

      print(string.format("%d\t %+f\t %+f", N, math.log10(L1), order))
   end
end


function ProblemList.IsentopicConvergenceRate()
   RunArgs.bound = "periodic"
   local res_values
   local L1_values = { }

   if RunArgs.dim == 1 then
      res_values = { 64, 128, 256, 512, 1024 }
   else
      res_values = { 16, 32, 64, 128}
   end

   print("------------------------------")
   print("N\t log10(L1)\t order")
   print("------------------------------")
   for n,N in pairs(res_values) do
      RunArgs.N = N
      local problem = tests.IsentropicPulse
      problem.mode = 2
      util.run_simulation(problem:get_pinit(), cfg_mara, RunArgs)

      local P = get_prim()
      local dS = problem:entropy(P.rho, P.pre) - problem.entropy_ref

      local L1 = 0.0
      for i=0,dS:size()-1 do
         L1 = L1 + math.abs(dS[i]) / N
      end

      L1_values[n] = math.log10(L1)
      local order = n == 1 and 0.0 or (
	 (L1_values[n] - L1_values[n-1]) /
	 (math.log10(res_values[n]) - math.log10(res_values[n-1])))

      print(string.format("%d\t %+f\t %+f", N, math.log10(L1), order))
   end
end



function ProblemList.ImplosionProblem()
   RunArgs.fluid = "euler"
   RunArgs.bound = {"reflect2d", 2, 3}
   util.run_simulation(tests[RunArgs.ic]:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaShocktube()
   local problem = tests[RunArgs.ic]
   if not problem then
      error("ic="..RunArgs.ic.." does not exist")
   end
   util.run_simulation(tests[RunArgs.ic]:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaExplosion()
   if RunArgs.fluid == "rmhd" then
      RunArgs.riemann = "hlld"
      RunArgs.advance = "single"
      RunArgs.godunov = "plm-muscl"
   end
   util.run_simulation(tests.Explosion:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaKelvinHelmholtz()
   RunArgs.dim = 2
   RunArgs.bound = "periodic"
   custom_domain = function(Nq, Ng)
      set_domain({0,0}, {2,1}, {2*RunArgs.N, RunArgs.N}, Nq, Ng)
   end
   util.run_simulation(tests.KelvinHelmholtz:get_pinit(), cfg_mara , RunArgs)
end
function ProblemList.SingleLayerKelvinHelmholtz()
   RunArgs.dim = 2
   RunArgs.bound = "perxouty"
--   custom_domain = function(Nq, Ng)
--      set_domain({0,0}, {2,1}, {2*RunArgs.N, RunArgs.N}, Nq, Ng)
--   end
   tests.KelvinHelmholtz.rho_profile = function(self, x, y, z)
      return math.exp(math.abs(y/0.2))
   end
   tests.KelvinHelmholtz.vx_profile = function(self, x, y, z)
      return y > 0.0 and -.5 or 0.5
   end
   util.run_simulation(tests.KelvinHelmholtz:get_pinit(), cfg_mara , RunArgs)
end
function ProblemList.SmoothKelvinHelmholtz()
   RunArgs.dim = 2
   RunArgs.bound = "periodic"
   util.run_simulation(tests.SmoothKelvinHelmholtz:get_pinit(), cfg_mara , RunArgs)
end
function ProblemList.VanillaDensityWave()
   util.run_simulation(tests.DensityWave:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaIsentropicPulse()
   util.run_simulation(tests.IsentropicPulse:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.TangentialVelocity()
   util.run_simulation(tests.TangentialVelocity:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end


-- -----------------------------------------------------------------------------
-- Some shortcuts
-- -----------------------------------------------------------------------------
ProblemList["explode"] = ProblemList["VanillaExplosion"]
ProblemList["implode"] = ProblemList["ImplosionProblem"]
ProblemList["denswave"] = ProblemList["VanillaDensityWave"]
ProblemList["isenwave"] = ProblemList["VanillaIsentropicPulse"]
ProblemList["kh"] = ProblemList["VanillaKelvinHelmholtz"]
ProblemList["kh1L"] = ProblemList["SingleLayerKelvinHelmholtz"]
ProblemList["shocktube"] = ProblemList["VanillaShocktube"]

local prob_to_run = ProblemList[RunArgs.problem]
if not prob_to_run then
   error("problem="..RunArgs.problem.." not found")
else
   prob_to_run()
end
