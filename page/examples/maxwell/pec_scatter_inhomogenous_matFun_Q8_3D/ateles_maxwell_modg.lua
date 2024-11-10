-- Configuration file for Ateles --
require('seeder')

timestep_info = 100

logging = {
  level = 10
}

-- global simulation options
simulation_name = 'pec_scatter_maxwell_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0.0,
    max = 1.0e-02 -- final simulation time
  }
}

--commpattern = 'gathered_type'

-- Mesh definitions --
mesh = './mesh/'

-- Tracking
tracking = {
  label = 'probe_displacementField_Q8',
  folder = './',
  variable = {'displacement_field'},
  shape = {kind = 'canoND', object= { origin ={-0.25,-0.25,0.0} } },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/2.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- Restart settings
--restart = { 
--            -- file to restart from
--            read = './restart/pec_scatter_maxwell_modg_lastHeader.lua',
--            -- folder to write restart data to
--            write = './restart/',
--            -- temporal definition of restart write
--            time_control = {
--                             min = 0,
--                             max = sim_control.time_control.max, 
--                             interval = sim_control.time_control.max/4.0,
--                          },
--          }

-- timing settings (i.e. output for performance measurements, this table is otional)
timing_file = 'timing.res'         -- the filename of the timing results

cube_permea = 2.0
cube_permit = 2.0

variable = {
  {
    name = 'permeability',
    ncomponents = 1,
    vartype = 'st_fun',
    evaltype = 'first',
    st_fun = {
      {
        const = { 1.0 }
      },
      {
        const = { cube_permea },
        shape={
          kind = 'canoND',
          object = {
            origin = { -0.5, -0.5, -0.25 },
            vec = {
              { 0.5-scatter_eps, 0.0, 0.0 },
              { 0.0, 0.5-scatter_eps, 0.0 },
              { 0.0, 0.0, 0.5-scatter_eps },
            },
            segments = { 100, 100, 100 }
          }
        }
      }
    }
  },
  {
    name = 'permittivity',
    ncomponents = 1,
    vartype = 'st_fun',
    evaltype = 'first',
    st_fun = {
      {
        const = { 1.0 }
      },
      {
        const = { cube_permit },
        shape = {
          kind = 'canoND',
          object = {
            origin = { -0.5, -0.5, -0.25 },
            vec = {
              { 0.5-scatter_eps, 0.0, 0.0 },
              { 0.0, 0.5-scatter_eps, 0.0 },
              { 0.0, 0.0, 0.5-scatter_eps },
            },
            segments = { 100, 100, 100 }
          }
        }
      }
    }
  },
  {
    name = 'conductivity',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 0.0
  }
}

-- Equation definitions --
equation = {
  name   = 'maxwell', -- we solve maxwell's equations
  material = {
    permeability = 'permeability',
    permittivity = 'permittivity',
    conductivity = 'conductivity'
  }
}

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',            -- we use the weno scheme for reconstruction
    m =  7,                   -- the reconstructed polynomial is of degree 1
    modg_space = 'Q'
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',     -- the name of the timestep control mechanism
      cfl  = 0.5,     -- Courant-Friedrichs-Lewy number
    }
  }
}

-- the postion and width of the pulse for the scatter test
pulseWidth = 0.025
pulseCenter = -0.25

function iniMagneticZ(x,y,z)
  r = math.sqrt( (x-pulseCenter)^2 + (y-pulseCenter)^2 )
  return math.exp((-1.0)*r/pulseWidth)
end

-- .. the general projection table
projection = {
  kind = 'fpt', -- 'fpt' or 'l2p', default 'l2p'
                -- for fpt the  nodes are automatically 'chebyshev'
                -- for lep the  nodes are automatically 'gauss-legendre'
  -- lobattoPoints = false  -- if lobatto points should be used, default = false
  factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
  -- blocksize = 32,        -- for fpt, default -1
  -- fftMultiThread = false -- for fpt, logical, default false
}
-- ...the initial condition table
initial_condition = {
  displacement_fieldX = 0.0,  -- electric field , x component
  displacement_fieldY = 0.0,  -- electric field , y component
  displacement_fieldZ = 0.0,  -- electric field , z component
  magnetic_fieldX = 0.0,  -- magnetic induction , x component
  magnetic_fieldY = 0.0,  -- magnetic induction , y component
  magnetic_fieldZ = iniMagneticZ,  -- magnetic induction , z component
}

-- Boundary definitions
boundary_condition = {
  {
    label = 'pecEast',
    kind = 'pec',
  },
  {
    label = 'pecWest',
    kind = 'pec',
  },
  {
    label = 'pecSouth',
    kind = 'pec',
  },
  {
    label = 'pecNorth',
    kind = 'pec',
  },
  {
    label = 'pecScatter',
    kind = 'pec',
  }
}
