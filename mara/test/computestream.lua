
local host = require 'host'
local ReadFile = cmdline.args[1]

h5_open_file(ReadFile, "r")
set_domain({ -0.5,-0.5,-0.5 }, { 0.5,0.5,0.5 }, h5_get_ndims("prim/rho"), 8, 2)
h5_close_file()

set_fluid("rmhd")
read_prim(ReadFile, host.CheckpointOptions)

set_boundary("periodic")
boundary.ApplyBoundaries()

visual.open_window({clear_color={0.05,0.025,0.025}})

local nframes = 1000
local nlines  = 40

for i=1,nframes do
   local lines = { }
   for n=1,nlines do
      local x = -0.25 + i/nframes + 0.05*n/nlines
      lines[n] = streamline({x, x, x}, 1.0, 1e-2, "velocity")
   end
   visual.draw_lines3d(lines)
end


