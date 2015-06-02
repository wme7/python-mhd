local Nzone = 16

set_fluid("euler")
set_domain({ -0.5,-0.5,-0.5 }, { 0.5,0.5,0.5 }, {Nzone,Nzone,Nzone}, 5, 2)

visual.open_window({clear_color={0.7,0.8,0.9}})

local Ntime = 100

for i=0,Ntime do
   local t = i / Ntime
   local w = 0.1*math.sin(4*math.pi*t)

   local function Pinit(x,y,z)
      return {1,1, -10*y, 10*x, w}
   end

   init_prim(Pinit)

   set_boundary("outflow")
   boundary.ApplyBoundaries()
   local S = streamline({0.1, 0.1, 0.0}, 6.6, 1e-2, "velocity")
   visual.draw_lines3d(S)
   collectgarbage()
end
