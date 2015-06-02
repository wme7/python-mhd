



function RunSimulation(CFL, EndTime)

   local dt = 0.0
   local Iteration = 0
   local CurrentTime = 0.0
   local piter = true
   local start = os.clock()

   print_mara()
   init_prim(ShockTube1)

   while CurrentTime < EndTime do

      local zps = advance(dt)

      if piter then
         print(string.format("%05d: t=%5.4f, dt=%5.4e %5.4ez/s, %3.2fus/(z*Nq)",
                             Iteration, CurrentTime, dt, zps, 1e6/(zps*5)))
      end

      dt = get_timestep(CFL)
      CurrentTime = CurrentTime + dt
      Iteration = Iteration + 1

   end
   print(string.format("run took %3.2f seconds", os.clock() - start))
   return get_prim()

end


function ShockTube1(x,y,z)
   if x < 0.5 then
      return { 1.000, 1.0, 0, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0 }
   end
end


set_domain({0.0}, {1.0}, {256}, 5, 2)
set_fluid("srhd")
set_eos("gamma-law", 1.4)
set_units(1.0, 1.0, 1.0)
set_boundary("outflow")
set_riemann("hllc")


runs = { }

set_advance("rk2")
set_godunov("plm-split")
runs["plm-split"] = RunSimulation(0.8, 0.2)

set_advance("single")
set_godunov("plm-muscl")
runs["plm-muscl"] = RunSimulation(0.8, 0.2)

function plot(runs, var, gp)
   for k,v in pairs(runs) do
      gp {
         f = v[var],
         post = string.format("with points title '%s'", k),
      }
   end
end

require 'gnuplot_t'

gp1 = gnuplot_t.new()
gp2 = gnuplot_t.new()
plot(runs, "pre", gp1)
plot(runs, "rho", gp1)
plot(runs, "vx", gp1)
