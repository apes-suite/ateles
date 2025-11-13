-- Setup of a velocity layer in a channel with 3D Euler --
-- This setup implements a basic Euler 3D simulation with a velocity layer
-- in rectangular channel, discretized with fourth order Q polynomials.

require 'hyperfun'

logging = { level = 10 }

-- global simulation options
simulation_name = 'shear_layer_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0.0,
    max = 1.0e-05,
  }
}

-- Mesh definitions --
mesh = 'mesh/'

-- Tracking
tracking = {
  label = 'probe_momentum_Q4',
  folder = '',
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
timing = { filename = 'timing.res' }

-- Equation definitions --
equation = {
  name = 'euler',
  isen_coef = 1.4,
  r = 296.0,
  material = {
    characteristic = 0,
    relax_velocity = {0, 0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
degree = 3
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m =  degree
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.6*(3*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

width = 0.07
xDist = 0.1
vel =  50.0
periods=2

initial_condition = {
  density = 1.225,
  pressure = 100000,

  velocityX = function(x,y,z)
    center = math.sin(periods*math.pi*x)*xDist
    return vel*tanh((y-center)/width)
  end,

  velocityY = 0.0,
  velocityZ = 0.0
}

boundary_condition = {
  {
    label = 'slipSouth',
    kind = 'slipwall'
  },
  {
    label = 'slipNorth',
    kind = 'slipwall'
  },
  {
    label = 'slipTop',
    kind = 'slipwall'
  },
  {
    label = 'slipBottom',
    kind = 'slipwall'
  }
}
