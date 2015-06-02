
local system = { }

-- special options for the mara io library, geared toward Lustre file systems
system.CheckpointOptions = {
   input_function="H5SER",
   output_function="H5SER",
   disk_align_threshold=4096,
   stripe_size_mb=4,
   enable_chunking=0,
   enable_alignment=0
}

-- command to run after 'mkdir -p' on special filesystems
system.Filesystem = function(d) return "" end

-- location of the file 'eos3.tab' on this system
system.ShenFile = ""

return system
