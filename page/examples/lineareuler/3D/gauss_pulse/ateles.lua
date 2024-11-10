-- Pulse in density transported in linearized 3D Euler equations
-- Parameters to vary --
degree = 12
poly_space = 'Q'

logging = { level = 4 }
-- Check for Nans and unphysical values
check = { interval = 1 }

--...Configuration of simulation time
sim_control = {
  time_control = {
    max = 0.01,
    min = 0.0,
    interval = {iter = 1}
  }
}

-- End  Parameters to vary --
--------------------------------------------------------------------------------
-- Definition of the test-case.

-- Mesh definitions --
mesh = {
  predefined = 'line_bounded', -- use the predefined line with boundaries
  origin = { -2.0 },           -- origin of the line
  length = 4,                  -- length of the line
  element_count = 4            -- number of elements
}

-- Equation definitions --
equation = {
  name   = 'linearEuler',
  isen_coef = 1.4,
  background = {
    density = 1.225,
    velocityX = 100.0,
    velocityY = 0.0,
    velocityZ = 0.0,
    pressure = 100000.0
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
    control = {
      name = 'cfl',
      cfl  = 0.95
    }
  }
}

-- ...the general projection table
projection = {
  kind = 'l2p',
  nodes_kind = 'chebyshev',
  factor = 1.0
}

-- variables for gaussian pluse
c = math.sqrt( equation.isen_coef * equation.background.pressure
                                  / equation.background.density )
ampl_density= equation.background.density/c
ampl_pressure= equation.background.pressure/c

function ic_gauss_density(x,y,z)
  d= x*x+y*y+z*z
  return( ampl_density * math.exp(-d/0.01*math.log(2)) )
end

initial_condition = {
  density = function(x,y,z)
    d= x*x+y*y+z*z
    return( ampl_density * math.exp(-d/0.01*math.log(2)) )
  end,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
  pressure = 0.0
}

boundary_condition = {
  {
    label = 'west',
    kind = 'primitives',
    density = 0.0,
    velocityX = 0.0,
    velocityY = 0.0,
    velocityZ = 0.0,
    pressure = 0.0
  },
  {
    label = 'east',
    kind = 'primitives',
    density = 0.0,
    velocityX = 0.0,
    velocityY = 0.0,
    velocityZ = 0.0,
    pressure = 0.0
  }
}

tracking = {
  label = 'track_density',
  folder = '',
  variable = {'density'},
  shape = { kind = 'canoND', object= { origin = {0., 0., 0.} } },
  time_control = {
    max = sim_control.time_control.max,
    min = 0,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
