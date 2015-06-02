

local json = require 'json'


local CheckpointOptions = {
   input_function="H5MPI",
   output_function="H5MPI",
   disk_align_threshold=1,
   stripe_size_mb=4,
   enable_chunking=1,
   enable_alignment=1
}

function RunSimulation(Status, Howlong)

   local t0 = Status.CurrentTime
   while Status.CurrentTime - t0 < Howlong do

      local dt = Status.Timestep
      local zps = advance(dt)

      driving.Advance(dt)
      driving.Resample()


      print(string.format("%05d: t=%5.4f, dt=%5.4e %4.3ez/s",
			  Status.Iteration, Status.CurrentTime, dt, zps))


      Status.Timestep = get_timestep(Status.CFL)
      Status.CurrentTime = Status.CurrentTime + Status.Timestep
      Status.Iteration = Status.Iteration + 1

   end
end



function InitSimulation()


   function pinit(x,y,z)
      r2 = x*x + y*y + z*z
      if r2 < 0.2 then
	 return { 1.000, 1.0, 0, 0, 0 }
      else
	 return { 0.125, 0.1, 0, 0, 0 }
      end
   end


   set_domain({-0.5,-0.5,-0.5}, {0.5,0.5,0.5}, {16,16,16}, 5, 2)
   set_fluid("euler")
   set_eos("gamma-law", 1.4)
   set_units(1.0, 1.0, 1.0)
   set_boundary("periodic")
   set_riemann("hllc")
   set_advance("single")
   set_godunov("plm-muscl")
   set_driving(new_ou_field(3, 1.0, 1.0, 2, 12345))

   init_prim(pinit)

   local Status = {
      CurrentTime = 0.0,
      Iteration = 0,
      Checkpoint = 0,
      Timestep = 0.0,
      CFL = 0.4,
   }
   return Status
end


function CheckpointWrite(Status)

   Status.Checkpoint = Status.Checkpoint + 1
   local chkpt = string.format("chkpt.%04d.h5", Status.Checkpoint)
   h5_open_file(chkpt, "w")
   h5_write_numeric_table("status", Status)
   h5_write_string("driving", json.encode(driving.Serialize()))
   h5_close_file()

   write_prim(chkpt, CheckpointOptions)
end


function CheckpointRead(chkpt)

   h5_open_file(chkpt, "r")
   local Status = h5_read_numeric_table("status")
   local field = json.decode(h5_read_string("driving"))
   h5_close_file()

   set_driving(field)
   read_prim(chkpt, CheckpointOptions)

   return Status
end



function Through()
   local Status = InitSimulation()

   for n=1,10 do
      RunSimulation(Status, 0.1)
      CheckpointWrite(Status)
   end
   os.execute("mv chkpt.0010.h5 through.h5")
end

function Restart()
   local Status = InitSimulation()

   for n=1,5 do
      RunSimulation(Status, 0.1)
      CheckpointWrite(Status)
   end

   local Status = CheckpointRead("chkpt.0005.h5")

   for n=1,5 do
      RunSimulation(Status, 0.1)
      CheckpointWrite(Status)
   end
   os.execute("mv chkpt.0010.h5 restart.h5")
end


Through()
Restart()
os.execute("rm chkpt.*.h5")
print("checking for differences in files: [expect no output]")
os.execute("h5diff through.h5 restart.h5")
