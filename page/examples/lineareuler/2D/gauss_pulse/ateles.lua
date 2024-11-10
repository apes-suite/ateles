-- Density pulse in 2D linearized Euler equations
-- This setup simulates the transport of a density pulse through a periodic
-- domain in 2D linearized Euler.
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Parameters to vary --
degree = 12
poly_space = 'Q'
-- ...the uniform refinement level for the periodic cube
level = 2

logging = { level = 4 }
-- Check for Nans and unphysical values
check =  { interval = 1 }

-- ...the general projection table
projection = {
  kind = 'fpt',
  factor = 1.0
}

--...Configuration of simulation time
sim_control = {
  time_control = {
    max = 0.02,
    min = 0.0,
    interval = {iter = 1}
  }
}

-- End  Parameters to vary --
--------------------------------------------------------------------------------
-- Definition of the test-case.

-- Mesh definitions --
cubeLength = 2.0
level = 2.0
mesh = {
  predefined = 'slice',
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    0.0
  },
  length = cubeLength,
  refinementLevel = level
}

-- Tracking
eps=cubeLength/(2^(level+1))
tracking = {
  label = 'track_2d_density',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= { origin = {0., 0., 0.} }
  },
  time_control = {
    max = sim_control.time_control.max,
    min = 0,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

background_dens = 1.225
background_velX = 100.0
background_press = 100000.0

-- Equation definitions --
equation = {
  name   = 'LinearEuler_2d',
  isen_coef = 1.4,
  background = {
    density =  background_dens,
    velocityX = background_velX,
    velocityY =  0.0,
    pressure = background_press
  }
}

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_2d',
    m =  degree,
    modg_space = poly_space
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.95
    }
  }
}

-- variables for gaussian pluse
c = math.sqrt(equation.isen_coef* background_press / background_dens)
ampl_density= background_dens/c
ampl_pressure= background_press/c

function gaus_dens(x,y,z)
  d= (x*x)+(y*y)
  return( ampl_density* math.exp(-d/0.01*math.log(2)) )
end

-- Initial Condition definitions --
initial_condition = {
  density = gaus_dens,
  velocityX = 0.0,
  velocityY = 0.0,
  pressure = 0.0
}
