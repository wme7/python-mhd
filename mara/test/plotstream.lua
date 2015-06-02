
h5_open_file(cmdline.args[1], "r")

visual.open_window()

for k,v in pairs(h5_get_setnames('.')) do
   local S = h5_read_array(v)
   visual.draw_lines3d(S)
   os.execute('sleep 0.05')
end


h5_close_file()
