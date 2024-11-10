require('wedge')
logging = {level=10}
--debug = {logging = {level=10}}
-------------------------------------------
--dt = 1.5e-5
--degree = 5
degree = 27
--degree = 3
-------------------------------------------
-- global simulation options
simulation_name='ateles' 
tmax = 0.2
--wallclock = 24*60*60 - 5*60
sim_control = {
             time_control = {
                  min = 0,
                  max = {iter=200},
                  --max = tmax,
                  --max = {sim=tmax, clock = wallclock},                 
                  interval = {iter = 1},
                }
}

check = {interval = 1}
---- Restart settings
NOrestart = { 
  --read = './restart/ateles_lastHeader.lua',
  write = './restart/',
  time_control = {
    min = 0, 
    max ={iter = 1000} ,
    --interval = {iter=1},
    align_trigger = {sim = true},
  },
}
segments = 200 * (degree +1)

-- tracking --
tracking = {
  { label = 'point',
    variable = {'density', 'pressure', 'velocity','polygon'},
    shape = {kind='canoND',
    object = { origin = {1.6, 0.0, 0.0},
               --segments = {segments, segments, segments, segments}
             },
     },
   time_control = {
   min = 0, 
   max = tmax,
   interval = {iter=1}
  },
    folder = './',
    output = {format = 'ascii', use_get_point = true}
  },
}

--physical data
gamma = 1.4
velocityX = 0.0
velocityY = 0.0
dens = 1.4
press = 400 
wedge_vel = 40.0
profile = {}
for i, v in ipairs(wedge_vertex) do
  point_X = wedge_vertex[i][1]
  point_Y = wedge_vertex[i][2]
  table.insert(profile, {point_X, point_Y})
end

function velocityRelax(x,y,z,t)
    return {-wedge_vel, 0.0}
end

-- Mesh definitions --
mesh = { predefined = 'slice',
         origin = {0, -1.0, 0},
         refinementLevel = 3,
         length = 2.0
       }

eps = 0.00001
variable = {
 { 
  name = 'polygon',
  ncomponents = 1,
  vartype = 'st_fun',
  st_fun = {
    {const = 0.0},
  {
    predefined = 'polygon_body_2d',
    movement = {movement_kind = 'lin_movement_2d'},
    lin_parameter = {-wedge_vel, 0.0},
    inval = {1.0},
    outval = {0.0},
    vertices = {profile},
    shape = {
      kind = 'canoND',
      object = {
           origin = {
             0.0,-0.6,0.0
           },
           vec = {
             { 2.0, 0.0, 0.0 },
             { 0.0, 1.2 , 0.0 },
             { 0.0, 0.0, eps },
           },
         },
       },
     },
    },
   },
  {
    name = 'relax_velocity',
    ncomponents = 2,
    vartype = "st_fun",
    st_fun = {
      { const = {0.0, 0.0} },
      {
        fun   = velocityRelax,
        shape = {
          kind = 'canoND',
          object = {
           origin = {
             0.0,-0.6,0.0
           },
           vec = {
             { 2.0, 0.0, 0.0 },
             { 0.0, 1.2 , 0.0 },
             { 0.0, 0.0, eps },
           },
          }
        }
      }
    }
  },
}

-- timing settings (i.e. output for performance measurements, this table is otional)
timing_file = 'timing.res'         -- the filename of the timing results

eps = 1e-5

phi = 1.0
beta = 1e-12
eta_v = phi^2 * beta^2 
eta_t = 0.4 * phi * beta
-- Equation definitions --
equation = {
  name = 'euler_2d',
  numflux = 'hll',
  isen_coef = 1.4,
  r = 1.0/1.4,
  --ensure_positivity = true,
  porosity             = phi,
  viscous_permeability = eta_v,
  thermal_permeability = eta_t,
  material = {
    characteristic = 'polygon',
    relax_velocity = 'relax_velocity',
    relax_temperature = press/(dens*(1.0/1.4)),
  --  mode_reduction = true,
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
  spatial =  {
    name = 'modg_2d',
    modg_space = 'Q',
    m = degree,
  },
  -- the spatial discretization scheme
  stabilization = {
   {
    name = 'spectral_viscosity',
    alpha = 36,
    order = 20,
    isAdaptive = true,
    --recovery_order = 1.0e-2,
    recovery_density = 1.0e-1,
    recovery_pressure = 1.0e-1
   },
  {
    name = 'covolume',
    alpha = 36,
    order = 20,
    beta = 1.0 - 2.0/(degree-1),
    isAdaptive = true,
   ----recovery_order = 1.0e-2,
    recovery_density = 1.0e-1,
    recovery_pressure = 1.0e-1
   },
  -- {
  --   name = 'cons_positivity_preserv',
  --   eps = 1.0e-07,
  -- }
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'imexRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl = 0.20, 
    },
  },
}
-- ...the general projection table
projection = {
  kind = 'l2p',
  --material = {factor = 3.0}
 -- lobattoPoints = true,  -- if lobatto points should be used,
 -- factor = 3.0,
}

-- Initial condition

initial_condition = {
  density = dens, 
  pressure = press, 
  velocityX = 0.0,
  velocityY = 0.0,
}

 -- Boundary definitions
--boundary_condition = {
--  
--  {
--    label = 'west',
--    kind = 'inflow',
--    density = dens,
--    velocityX = 0.0,
--    velocityY = 0.0,
--    pressure = press,
--  }
--  ,
--  {
--    label = 'wall',
--    kind = 'wall',
--  },
--  {
--    label = 'east',
--    kind = 'supersonic_outflow',
--  },
-- }
