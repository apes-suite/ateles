-- Euler 3D setup with FPT projection --
logging = { level = 10 }

-- ...the length of the cube
cubeLength = 2.0

-- the refinement level of the octree
level = 1

-- Transport velocity of the pulse in x direction.
velocityX = 100

-- global simulation options
simulation_name = 'gPulseDens_euler_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0,
    max = cubeLength/velocityX/4 -- final simulation time
  }
}

-- Mesh definitions --
mesh = {
  predefined = 'cube',
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0
  },
  length = cubeLength,
  refinementLevel = level
}

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
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m = 7
  },
  -- the temporal discretization scheme
  temporal = {
    -- Explicit Runge Kutta
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.8
    }
  }
}

-- ...the general projection table
projection = {
  kind = 'fpt',
  factor = 1.0
}

-- This is a very simple example to define constant boundary condtions.
initial_condition = {
  density = {
    predefined = 'gausspulse',
    center = { 0.0, 0.0, 0.0 },
    halfwidth = 0.20,
    amplitude = 2.0,
    background = 1.225
  },
  pressure = 100000,
  velocityX = velocityX,
  velocityY = 0.0,
  velocityZ = 0.0
}

-- Tracking
tracking = {
  label = 'track_momentum',
  folder = '',
  variable = {'momentum'},
  shape = {kind = 'canoND', object= { origin ={0., 0., 0.} } },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/8.0
  },
  output = { format = 'ascii' , ndofs = 1}
}
