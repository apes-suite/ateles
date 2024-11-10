logging = {level = 10}
--debug = {logging = {level=10}}
-------------------------------------------

order = 16
degree = order - 1
nElems = 128
------------------------------------------
tmax= 5 * 1e-06
-- global simulation options
simulation_name='ateles'
--wallclock = 2*60*60 - 4*60
sim_control = {
             time_control = {
                  min = 0,
                  max = tmax,
                  --max = {sim=tmax, clock = wallclock},
                  interval = {iter = 1},
                }
}

check = {interval = 1}
---- Restart settings
NOrestart = {
 -- read = './restart/ateles_lastHeader.lua',
  write = './restart/',
  time_control = {
    min = 0,
    max = tmax,
    interval = tmax/20,
    align_trigger = {sim = true},
  },
}

-- Mesh definition --
mesh = {
  predefined = 'line_bounded',
  origin = {0.0},
  length = 1.0,
  element_count = nElems
}

-- Segments for tracking--
segments = math.ceil(mesh.length*nElems) * (degree + 1) * 4
eps = mesh.length / segments * 1.e-6

-- tracking --
tracking = {
  { label = 'line',
    variable = {'density', 'pressure', 'velocity'},
    shape = {kind='canoND',
    object = { origin = {0.0+eps, 0.0, 0.0},
               vec = {1.0-2*eps, 0.0, 0.0},
               segments = {segments}},
     },
   time_control = {
   min = 0,
   max = tmax,
   interval = tmax
  },
    folder = './',
    output = {format = 'asciiSpatial', use_get_point = true}
  },
}

--physical data
gamma = 1.4
velocityX = 150.0
dens = 1.0
press = 1.0e5
-- Define a piston inside the domain --
xmin = 0.4
xmax = 0.44

temp = press/(dens*(1.0/1.4))
vel_init = 0.0
c = math.sqrt((press*gamma)/dens)
u = velocityX
a = - (gamma/(2*c^2))*((4*c^2/gamma)+(gamma + 1) * u^2)
b = (gamma/(2*c^2))* ((2*c^2/gamma) - u^2*(gamma - 1))
p_ratio1 = -a/2 + math.sqrt((a/2)^2 - b)
p_ratio2 = -a/2 - math.sqrt((a/2)^2 - b)
density_ratio = (1 + ((gamma + 1)/(gamma - 1 ))*
  math.max(p_ratio1,p_ratio2))/(((gamma + 1)/(gamma -1))+
  math.max(p_ratio1,p_ratio2))

densL = density_ratio * dens
pressL = math.max(p_ratio1,p_ratio2) * press

temperature_ratio = (1 +(pressL - press)/press ) *
  ((2*gamma + (gamma-1)*((pressL-press)/press))/
  (2*gamma + (gamma+1)*((pressL - press)/press)))
tempL = temperature_ratio * temp

function inside_piston(x,y,z,t)
  xa = xmin + velocityX * t
  xb = xmax + velocityX * t
  if (xa <= x and x <= xb ) then
    return true
  else
    return false
  end
end

function velocityRelax(x,y,z,t)
    return velocityX
end

function temperature(x,y,z,t)
  xa = xmin + velocityX * t
  xb = xmax + velocityX * t
  diff_half = (xb -xa)/2
  if (x >= xa and  x <=(xa + diff_half)) then
    return temp
  elseif (x >(xa + diff_half) and x<=xb) then
    return tempL
  else
    return temp
  end
end

function characteristic(x,y,z,t)
  if inside_piston(x,y,z,t) then
    return 1.0
  else
    return 0.0
  end
end
dx = 0.00000001
variable = {
  {
    name = 'Xi',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      { const = 0.0 },
      {
         fun   = characteristic,
         shape = {
           kind = 'canoND',
           object = {
             origin = {0.3 ,0.0, 0.0},
             vec = {
               { 0.6, 0.0, 0.0 },
               { 0.0, dx, 0.0 },
               { 0.0, 0.0, dx },
             },
           }
         }
      }
    }
  },
  {
    name = 'relax_velocity',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      { const = 0.0 },
      {
         fun   = velocityRelax,
         shape = {
           kind = 'canoND',
           object = {
             origin = {0.3 ,0.0,0.0},
             vec = {
               { 0.6, 0.0, 0.0 },
               { 0.0, dx, 0.0 },
               { 0.0, 0.0, dx },
             },
           }
         }
      }
    }
  },
  {
    name = 'relax_temperatur',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      { const = 0.0 },
      {
         fun   = temperature,
         shape = {
           kind = 'canoND',
           object = {
             origin = {0.3 ,0.0,0.0},
             vec = {
               { 0.6, 0.0, 0.0 },
               { 0.0, dx, 0.0 },
               { 0.0, 0.0, dx },
             },
           }
         }
      }
    }
  },
}
-- timing settings (i.e. output for performance measurements, this table is otional)
timing_file = 'timing.res'         -- the filename of the timing results

phi = 1.0
beta = 1e-6
eta_v = phi^2 * beta^2
eta_t = 0.4 * phi * beta
-- Equation definitions --
equation = {
  name = 'euler_1d',
  numflux = 'hll',
  isen_coef = 1.4,
  r = 1.0/1.4,
  porosity             = phi,
  viscous_permeability = eta_v,
  thermal_permeability = eta_t,
  material = {
    characteristic = 'Xi',
    relax_velocity = 'relax_velocity',
    relax_temperature = 'relax_temperatur',
    mode_reduction = true,
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_1d',
    modg_space = 'Q',
    m = degree,
  },
  stabilization = {

  {
    name = 'spectral_viscosity',
    alpha = 36,
    order = 24,
    isAdaptive = true,
    --recovery_order = 1.0e-4,
    recovery_density = 1.0e-1,
    recovery_pressure = 1.0e-1
  },
  --{
  --   name = 'covolume',
  --   alpha = 36,
  --   order = 24,
  --   beta = 1.0 - 2.0/degree,
  -- },
},
  -- the temporal discretization scheme
  temporal = {
    name = 'imexRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl = 0.2,
    },
  },
}
-- ...the general projection table
projection = {
  kind = 'l2p',
  material = {
    factor = 3.0
    },
}
vel_init = 0.0
c = math.sqrt((press*gamma)/dens)
u = velocityX
a = - (gamma/(2*c^2))*((4*c^2/gamma)+(gamma + 1) * u^2)
b = (gamma/(2*c^2))* ((2*c^2/gamma) - u^2*(gamma - 1))
p_ratio1 = -a/2 + math.sqrt((a/2)^2 - b)
p_ratio2 = -a/2 - math.sqrt((a/2)^2 - b)
density_ratio = (1 + ((gamma + 1)/(gamma - 1 ))* math.max(p_ratio1,p_ratio2))/
                 (((gamma + 1)/(gamma -1))+ math.max(p_ratio1,p_ratio2))

densL = density_ratio * dens
pressL = math.max(p_ratio1,p_ratio2) * press
ushock= (u/c)/(1 - density_ratio^(-1))* c
velocity_2 = ushock - u
velocity_1 = ushock
x_position = 0.44 + ushock*tmax

-- Initial condition
function iniVel(x,y,z,t)
  if inside_piston(x,y,z,0.0) then
    return velocityX
  else
    return 0.0
  end
end

function inipress(x,y,z,t)
  if (x > 0.42 and x < 0.44) then
    return pressL
  else
    return press
  end
end

function inidens(x,y,z,t)
  if (x > 0.42 and x < 0.44) then
    return densL
  else
    return dens
  end
end

initial_condition = {
  density = inidens,
  pressure = inipress,
  velocity = iniVel,
}

 -- Boundary definitions
boundary_condition = {
  {
    label = 'west',
    kind = 'outflow',
    pressure = press,
  }
  ,
  {
    label = 'east',
    kind = 'outflow',
    pressure = press,
  }
  ,
}
