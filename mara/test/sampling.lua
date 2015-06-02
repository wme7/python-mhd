

-- *****************************************************************************
--
-- Test of parallel sampling routines.
--
-- *****************************************************************************

local Nsamp = 10000
local Nzone = 32


local function TestSamplingInternal()
   set_domain({0,0,0}, {1,1,1}, {Nzone,Nzone,Nzone}, 8, 2)
   init_prim(function(x,y,z)
		return {1,1, 0,0,0, x,y,z}
	     end)
   set_boundary("outflow")
   boundary.ApplyBoundaries()
   local trun = test_sampling(Nsamp)
   print(string.format("took %d samples in %f seconds", Nsamp, trun))
end


local function TestSamplingInternalMany()
   set_domain({0,0,0}, {1,1,1}, {Nzone,Nzone,Nzone}, 8, 2)
   init_prim(function(x,y,z)
		return {1,1, 0,0,0, x,y,z}
	     end)
   set_boundary("outflow")
   boundary.ApplyBoundaries()
   local trun = test_sampling_many(Nsamp)
   print(string.format("took %d samples in %f seconds", Nsamp, trun))
end


local function TestSamplingNd(dims, mode, verbose)

   if     dims == 1 then
      set_domain({0}, {1}, {Nzone}, 8, 2)

   elseif dims == 2 then
      set_domain({0,0}, {1,1}, {Nzone,Nzone}, 8, 2)

   elseif dims == 3 then
      set_domain({0,0,0}, {1,1,1}, {Nzone,Nzone,Nzone}, 8, 2)
   end


   local coords = { } -- buffer the grid points sampled by init_prim

   init_prim(function(x,y,z)
--		table.insert(coords, lunum.array({x,y,z}))
		return {1,1, 0,0,0, x,y,z}
	     end)

   set_boundary("outflow")
   boundary.ApplyBoundaries()

   if mode == 'random' then
      -- ***********************************************************************
      -- Choose random points to be sampled throughout the domain

      local start = os.clock()
      math.randomseed(12345)

      for i=0,Nsamp do
	 local x, y, z =  math.random(), math.random(), math.random()
	 local P = prim_at_point{x,y,z}
	 if verbose then print(lunum.array{x,y,z}, P) end	
	 collectgarbage()
      end
      
      print(string.format("took %d samples in %f seconds",
			  Nsamp, os.clock() - start))


   elseif mode == 'grid' then
      -- ***********************************************************************
      -- Run the samples over the same points sampled by the init_prim function

      for k,v in pairs(coords) do
	 local P = prim_at_point(v)
	 if verbose then print(v,P) end
	 collectgarbage()
      end

      visual.open_window()
      visual.draw_texture(get_prim().By)


   elseif mode == 'prolong' then
      -- ***********************************************************************
      -- Sample the grid one level finer to test the quality of interpolation

      local Bx = lunum.zeros{2*Nzone, 2*Nzone}

      for i,j in Bx:indices() do
	 local x = (i + 0.5) / (2*Nzone)
	 local y = (j + 0.5) / (2*Nzone)
	 Bx[{i,j}] = prim_at_point{x,y}[6]
	 collectgarbage()
      end

      visual.open_window()
      visual.draw_texture(Bx)
   end
end


set_fluid("rmhd")
TestSamplingInternal()
TestSamplingInternalMany()

--TestSamplingNd(3, 'random', true)
--TestSamplingNd(2, 'grid', true)
--TestSamplingNd(2, 'prolong', false)
