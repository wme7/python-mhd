

local function TestTexture()

   local Nx = 100
   local Ny = 240

   visual.open_window()

   local x = lunum.zeros{Nx,Ny}
   local y = lunum.zeros{Nx,Ny}

   for i,j in x:indices() do
      x[{i,j}] = 0.08 * i * math.pi
   end
   for i,j in y:indices() do
      y[{i,j}] = 0.04 * j * math.pi
   end

   
   for iter=0,1000 do

      local A = lunum.sin(x+math.cos(0.1*iter)) + lunum.cos(y-math.sin(0.1*iter))
      print(iter)
      collectgarbage()
      visual.draw_texture(A)

   end
end


local function TestLines3d()

   local N = 100
   local points1 = lunum.zeros({N,4})
   local points2 = lunum.zeros({N,4})

   for i=0,N-1 do
      local t = i/N

      local x =  0.5*(math.sin(20*t) - 0.5)
      local y =  0.5*(math.tanh(12*t) - 0.5)
      local z =  0.5*(math.cos(5*t) - 0.5)
      local w =       math.cos(5*t)

      points1[{i,0}] = x
      points1[{i,1}] = y
      points1[{i,2}] = z
      points1[{i,3}] = w

      points2[{i,0}] = x/2
      points2[{i,1}] = y/2
      points2[{i,2}] = z/4
      points2[{i,3}] = w
   end

   visual.open_window({clear_color={0.2,0.8,0.2}, window_size={1024,768}})
   visual.draw_lines3d({points1, points2})

end


--TestTexture()
TestLines3d()
