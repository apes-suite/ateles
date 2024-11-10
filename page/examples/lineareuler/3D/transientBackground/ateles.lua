-- Linearized Euler with timedependent background state --
--
-- This setup illustrates the definition of a background state that
-- varies in time.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Parameters to vary --
degree = 13
nelems = 4
poly_space = 'Q'
-- End  Parameters to vary --

-- The variation of the background state in the linearized
-- equations over time.
background_dens = 1.0
background_velX = 100.0
background_press = 2.0

function sinus_dens(t)
  return (background_dens + 0.5*math.sin(t*2*math.pi))
end

function sinus_velX(t)
  return (background_velX + 0.1*math.sin(t*2*math.pi) )
end

function sinus_press(t)
  return (background_press + 0.5*math.sin(t*2*math.pi))
end

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
    max = 0.01,
    min = 0.0,
    interval = {iter = 100}
  }
}

--------------------------------------------------------------------------------
-- Details.
simulation_name = 'lineareuler'

-- Mesh definitions --
cubeLength = 2.0
mesh = {
  predefined = 'line',
  origin = { -2.0 },
  length = cubeLength,
  element_count = nelems
}

-- Tracking
eps=cubeLength*2^(-21)
tracking = {
  label = 'transientBackground',
  folder = './',
  variable = {'density', 'completeState'},
  shape = {
    kind = 'canoND',
    object= { origin = {0.0, 0.0, 0.0} }
  },
  time_control = {
    max = sim_control.time_control.max,
    min = 0,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}


-- Equation definitions --
equation = {
  name   = 'LinearEuler',
  isen_coef = 1.4,
  background = {
    density   = sinus_dens,
    velocityX = sinus_velX,
    velocityY = 0.0,
    velocityZ = 0.0,
    pressure  = sinus_press,
  }
}

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m =  degree,
    modg_space = poly_space
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.95
    }
  }
}

-- Initial Condition definitions --
initial_condition = {
  density   = 0.0,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
  pressure  = 0.0
}
