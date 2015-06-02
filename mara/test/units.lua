
host = require 'host'

unitsys = { }

local function neutron_star()
   local LIGHT_SPEED = 2.99792458000e+10 -- cm/s

   local Density = 1e13;        -- gm/cm^3
   local V       = LIGHT_SPEED; -- cm/s
   local Length  = 1e2;         -- cm
   local Mass    = Density * Length^3.0
   local Time    = Length / V;
   set_units(Length, Mass, Time)
   units.Print()
end

unitsys = { neutron_star = neutron_star }


unitsys.neutron_star()

set_eos("gamma-law", 1.4)
print(eos.TemperatureMeV(1.0, 0.05))

local tab = load_shen(host.ShenFile, 0.08, {1e12, 1e14}, {0.1, 200.0})
set_eos("shen", tab)
print(eos.TemperatureMeV(1.0, 0.05))
