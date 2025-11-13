-- Simulation of an acoustic pulse with linearized 3D Euler equations
-- This setup illustrates the definition of derived quantities to track
-- in the domain.
-- See the 'variable' table.
simulation_name = 'linearEuler_gradients'

-- Parameters to vary --
degree = 11
poly_space = 'Q'
-- ...the uniform refinement level for the periodic cube
level = 2

logging = { level = 4 }
-- Check for Nans and unphysical values
check =  { interval = 1 }

-- ...the general projection table
projection = {
  kind = 'l2p',
  factor = 1.0
}

--...Configuration of simulation time
sim_control = {
  time_control = {
    max = {iter = 20},
    min = 0.0,
    interval = {iter = 1}
  }
}

-- Equation definitions --
bg_dens = 1.225
bg_velX = 100.0
bg_velY = 0.0
bg_velZ = 0.0
bg_press = 100000.0

equation = {
  name   = 'linearEuler',
  numflux = 'godunov',
  isen_coef = 1.4,
  background = {
    density = bg_dens,
    velocityX = bg_velX,
    velocityY = bg_velY,
    velocityZ = bg_velZ,
    pressure = bg_press
  }
}

-- Mesh definitions --
cubeLength = 4.0
mesh = {
  predefined = 'cube',
  origin = { -2.0, 0.0, 0.0 },
  length = 4.0,
  refinementLevel = 2
}

variable = {
  {
    -- Arbitrary space-time function for a variable,
    -- here the constant background density.
    name = 'bg_density',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun =  bg_dens
  },
  -- Operations allow us to combine any of the available
  -- variables to new derived ones:
  {
    name = 'full_density',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'addition',
      input_varname = {'density', 'bg_density'}
    }
  },
  {
    name = 'grad_fulldensity',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'gradient',
      input_varname = 'full_density',
    }
  },
  {
    name = 'gradX_fulldensity',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'grad_fulldensity',
      input_varindex = {1}
    }
  },
  {
    name = 'grad_density',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'gradient',
      input_varname = 'density',
    }
  },
  {
    name = 'gradX_density',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'grad_density',
      input_varindex = {1}
    }
  }
}

-- Tracking
tracking = {
  label = 'track_grads',
  folder = './',
  variable = {
    'density',
    'full_density',
    'grad_density',
    'gradX_density',
    'grad_fulldensity',
    'gradX_fulldensity'
  },
  shape = {
    kind = 'canoND',
    object= {
      origin = {0.0, 0.0, 0.0}
    }
  },
  time_control = {
    max = sim_control.time_control.max,  -- final Simulated time
    min = 0,
    interval = {iter = 10}
  },
  output = {
    format = 'ascii', use_get_point = true
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
    name = 'explicitRungeKuttaTaylor',
    steps = 8,
    control = {
      name = 'cfl',
      cfl = 2*(3*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

-- variables for gaussian pluse of the initial condition
c = math.sqrt( equation.isen_coef * equation.background.pressure
                                  / equation.background.density )
ampl_density= equation.background.density/c
ampl_pressure= equation.background.pressure/c

function ic_gauss_density(x,y,z)
  d= x*x+y*y+z*z
  return( ampl_density * math.exp(-d/0.01*math.log(2)) )
end

function ic_gauss_pressure(x,y,z)
  d= x*x+y*y+z*z
  return( ampl_pressure * math.exp(-d/0.01*math.log(2)) )
end

-- Initial Condition definitions --
initial_condition = {
  density = ic_gauss_density,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
  pressure = ic_gauss_pressure
}
