



set_eos("gamma-law", 1.4)
set_fluid("rmhd")
math.randomseed(12345)

local cos = math.cos
local sin = math.sin
local random = math.random
local acos   = math.acos
local pi     = math.pi
local sqrt   = math.sqrt
local log10  = math.log10
local pow    = math.pow


local function random_state(Vel, Mag)

   local rho,pre,vx,vy,vz,Bx,By,Bz = 1,2,3,4,5,6,7,8

   local thtv = acos(random() * 2.0 - 1.0)
   local thtb = acos(random() * 2.0 - 1.0)
   local phiv = random() * 2*pi
   local phib = random() * 2*pi

   local P = { }

   P[rho] = 1.0
   P[pre] = 1.0

   P[vx ] = Vel*sin(thtv)*cos(phiv);
   P[vy ] = Vel*sin(thtv)*sin(phiv);
   P[vz ] = Vel*cos(thtv);

   P[Bx ] = Mag*sin(thtb)*cos(phib);
   P[By ] = Mag*sin(thtb)*sin(phib);
   P[Bz ] = Mag*cos(thtb);

   return P
end


function run_test(ntrials, Vel, Mag, args)

   if args.print_loud then
      print("\nTesting conserved to primitive solvers for gamma-law EOS.\n")
   end

   if args.print_iter then
      print("(code iterations error time(s))\n")
   end

   if args.print_loud then
      print(string.format("%22s %24s %24s %24s", "noble2dzt", "duffell3d",
			  "noble1dw", "anton2dzw"))
      print(string.format("--------------------------------------------------------"..
			  "-----------------------------------------"))
   end

   local codes_av = {0,0,0,0}
   local iters_av = {0,0,0,0}
   local error_av = {0,0,0,0}
   local times_av = {0,0,0,0}
   local npass    = {0,0,0,0}

   for m=1,ntrials do

      P = random_state(Vel, Mag)
      codes, iters, error, times = test_rmhd_c2p(P)

      if args.print_iter then
         print(string.format("(%d %2d %2.1e %2.1e)   (%d %2d %2.1e %2.1e)   " ..
                             "(%d %2d %2.1e %2.1e)   (%d %2d %2.1e %2.1e) ",
                          codes[1], iters[1], error[1], times[1],
                          codes[2], iters[2], error[2], times[2],
                          codes[3], iters[3], error[3], times[3],
                          codes[4], iters[4], error[4], times[4]))
      end

      for i=0,3 do
         if codes[i] == 0 then npass[i] = npass[i] + 1 end
         codes_av[i] = codes_av[i] + codes[i] / ntrials
         iters_av[i] = iters_av[i] + iters[i] / ntrials
         error_av[i] = error_av[i] + error[i] / ntrials
         times_av[i] = times_av[i] + times[i] / ntrials
      end
   end

   local pass = { }
   for i=0,3 do
      pass[i] = npass[i]/ntrials
   end

   if args.print_loud then
      print(string.format("--------------------------------------------------------"..
			  "-----------------------------------------\n"))
      print("percent success, average time(s):\n")
      print(string.format("    %3.2f%%   %3.2e       %3.2f%%   %3.2e   "..
			  "    %3.2f%%   %3.2e       %3.2f%%   %3.2e",
		       100*pass[1], times_av[1],
		       100*pass[2], times_av[2],
		       100*pass[3], times_av[3],
		       100*pass[4], times_av[4]))
   end
   return pass, times_av

end


local pass = { }
local time = { }

local Nsamp_G = 100
local Nsamp_B = 100


for i=1, Nsamp_G do

   pass[i] = { }
   time[i] = { }

   for j=1, Nsamp_B do

      local G = 1.0 + 1000.0 * (i-1) / Nsamp_G
      local B = pow(10.0, 3.0 * (j-1) / Nsamp_B)
      local V = sqrt(1.0 - 1.0/(G*G))

      p, t = run_test(100, V, B, {print_iter=false, print_loud=false})

      pass[i][j] = {G, log10(B), p[1], p[2], p[3], p[4]}
      time[i][j] = {G, log10(B), t[1], t[2], t[3], t[4]}

      print(p)
   end
end


local json = require 'json'

local f = io.open("c2p_pass.json", "w")
f:write(json.encode(pass))

local f = io.open("c2p_time.json", "w")
f:write(json.encode(time))
