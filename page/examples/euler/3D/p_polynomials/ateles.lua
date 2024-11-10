-- Euler 3D setup that employs the P-polynomial strategy --
-- This setup simulates a plane of higher velocity fluid
-- in an otherwise still fluid within a rectangular tube.
-- It's main purpose is to show the use of P-polynomials
-- for the construction of the multidimensional polynomials.
require 'hyperfun'

-- ...the length of the cube
cubeLength = 2.0

-- global simulation options
simulation_name = 'shear_layer_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0.0,
    max = 1.0e-05
  }
}

-- Mesh definitions --
mesh = 'mesh/'

-- Tracking
tracking = {
  label = 'probe_momentum_P8',
  folder = './',
  variable = {'momentum'},
  shape = {kind = 'canoND', object= { origin ={0., 0., 0.} } },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/2.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- the filename of the timing results i.e. output for performance measurements
timing_file = 'timing.res'


-- Equation definitions --
equation = {
    name       = 'euler',
    isen_coef  = 1.4,
    r          = 296.0,
    material = {
      characteristic = 0,
      relax_velocity = {0,0,0},
      relax_temperature = 0
    }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    modg_space = 'P',
    m =  7
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
       name = 'cfl',
       cfl  = 0.6
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0
}


width = 0.07
xDist = 0.1
vel =  50.0
periods = 2

initial_condition = {
  density = 1.225,
  pressure = 100000,
  --velocityX = iniVelX,
  velocityX = function(x,y,z)
    center = math.sin(periods*math.pi*x)*xDist
    return vel*tanh((y-center)/width)
  end,
  velocityY = 0.0,
  velocityZ = 0.0
}

-- Boundary definitions
boundary_condition = {
  {
    label = 'slipSouth',
    kind = 'slipwall'
  },
  {
    label = 'slipNorth',
    kind = 'slipwall',
  },
  {
    label = 'slipTop',
    kind = 'slipwall',
  },
  {
    label = 'slipBottom',
    kind = 'slipwall',
  }
}
