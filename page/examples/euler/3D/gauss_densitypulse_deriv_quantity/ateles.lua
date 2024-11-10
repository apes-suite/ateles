-- Setup illustrating tracking of derived quantities for Euler 3D --
-- This is a simple setup with a convected pulse in density that is moved
-- through a domain with periodicity in all directions.
--
-- The solver operates on the conservative quantities density, momentum and
-- energy, but it is possible to compute other quantities out of those. Here
-- we track kinetic energy and velocity as derived quantities in a single
-- element.

logging = {level = 10}

-- ...the length of the cube
cubeLength = 2.0

-- the refinement level of the octree
level = 2

-- Transport velocity of the pulse in x direction.
velocityX = 100

-- global simulation options
simulation_name = 'gPulseDens_euler_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0,
    max = 0.0007 -- final simulation time
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
    m = 7   -- the maximal polynomial degree for each spatial direction
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.9,
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

initial_condition = {
  density = {
    predefined = 'gausspulse',
    center     = {0.0, 0.0, 0.0},
    halfwidth  = 0.20,
    amplitude  = 2.0,
    background = 1.225
  },
  pressure = 100000,
  velocityX = velocityX,
  velocityY = 0,
  velocityZ = 0
}

-- Tracking the integral mean of the given quantities in the element
-- containing the given point (the origin).
tracking = {
  label = 'track_ke',
  folder = '',
  variable = { 'momentum', 'velocity', 'density', 'kinetic_energy' },
  shape = {
    kind = 'canoND',
    object= {
      origin ={ 0., 0., 0. }
    }
  },
  time_control = {
    min = 0,
    max = 0.005,
    interval = sim_control.time_control.max/10
  },
  output = { format = 'ascii', ndofs = 1 }
}
