

function TestSerialize()
   local json = require 'json'
   local field = new_ou_field(2, 1.0, 1.0, 2, 12345)
   print("json encoding works?", json.decode(json.encode(field)) == field)


   set_domain({-0.5, -0.5}, {0.5, 0.5}, {128, 128}, 5, 2)
   set_driving(field)

   driving.Resample()
   print("driving reconstruction works?", driving.Serialize() == field)

   driving.Advance(0.1)
   print("changes after advancing?", not (driving.Serialize() == field))


   h5_open_file("driving.h5", "w")
   h5_write_string("driving", json.encode(driving.Serialize()))
end


function TestRestart()

   set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {32, 32, 32}, 5, 2)
   local initial_field = new_ou_field(3, 1.0, 1.0, 2, 12345)

   -- Going through to 20s with a stop in between
   -- ..........................................................................
   set_driving(initial_field)

   local t = 0.0
   local dt = 0.1

   while t < 10.0 do
      driving.Advance(dt)
      t = t + dt
   end

   local intermed_field = driving.Serialize()
   set_driving(intermed_field)

   while t < 20.0 do
      driving.Advance(dt)
      t = t + dt
   end

   local final_field = driving.Serialize()

   -- Going through to 20s continuously
   -- ..........................................................................
   set_driving(initial_field)

   t = 0.0
   dt = 0.1

   while t < 20.0 do
      driving.Advance(dt)
      t = t + dt
   end

   CheckpointOptions = {
      input_function="H5SER",
      output_function="H5SER",
      disk_align_threshold=1,
      stripe_size_mb=4,
      enable_chunking=0,
      enable_alignment=0
   }

   print("restarted field agrees?", final_field == driving.Serialize())
end



TestSerialize()
TestRestart()

