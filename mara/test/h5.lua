


local json = require('json')

local TestFile = "unit-test.h5"

h5_open_file(TestFile, "w")
h5_open_group("group1", "w")
h5_open_group("group2", "w")


h5_open_file(TestFile, "r+")
h5_open_group("group3", "w")
h5_open_group("group4", "w")


test_data = { }

for i=1, 2000 do
   test_data[i] = {num=3.14, thing="a word"}
end

h5_write_string("group1/message", "message for group1")
h5_write_string("group2/message", "message for group2")
h5_write_string("group3/thing", json.encode(test_data))
h5_write_string("group4/thing", json.encode(test_data))
h5_close_file()


h5_open_file(TestFile, "r")

message1 = h5_read_string("group1/message")
message2 = h5_read_string("group2/message")

print("decoding...")
data3 = json.decode(h5_read_string("group3/thing"))
data4 = json.decode(h5_read_string("group4/thing"))

print(message1)
print(message2)
print(data3[1]["num"])
print(data4[1]["thing"])


h5_close_file()

CheckpointOptions = {
   input_function="H5SER",
   output_function="H5SER",
   disk_align_threshold=1,
   stripe_size_mb=4,
   enable_chunking=0,
   enable_alignment=0
}

set_domain({0.0, 0.0}, {1.0, 1.0}, {64, 128}, 5, 2)
init_prim(function(x,y,z) return {1,2,3,4,5} end)
set_fluid("euler")

sec = write_prim(TestFile, CheckpointOptions)
print("write_prim took " .. sec .. " seconds")

sec = read_prim(TestFile, CheckpointOptions)
print("read_prim took " .. sec .. " seconds")



local Status = {
   CurrentTime = 1.3,
   Iteration = 12,
   Checkpoint = 3,
   Timestep = 0.0123,
   CFL = 0.4,
}


h5_open_file(TestFile, "r+")
h5_write_numeric_table("status", Status)

for k,v in pairs(h5_read_numeric_table("status")) do
   print(string.format("%s = %f (%f)", k, v, Status[k]))
end

h5_close_file()

h5_open_file(TestFile, "r+")
print("number of data sets:", h5_get_nsets("prim"))
print("dimension of data:", h5_get_ndims("prim/rho"))
h5_close_file()



h5_open_file(TestFile, "r+")

local A = lunum.range(20):reshape{4,5} * 2.0
print(A:shape'array', A:dtype())
h5_write_array("nd-array", A)

h5_close_file()


h5_open_file(TestFile, "r")

local B = h5_read_array("nd-array")
print(B:shape'array', B:dtype())
print(B:eq(A))

h5_close_file()



h5_open_file(TestFile, "r")
for k,v in pairs(h5_get_setnames(".")) do
   print(v)
end
h5_close_file()
