-- Setup of a constant state in Navier-Stokes 2D --
--
-- This very basic setup can be used to check, whether something
-- is fundamentally wrong. Nothing should happen.
--
-- It makes use of the explicit Euler time integration scheme, which
-- actually is not stable for higher order discretizations.
-- Though, we use it when looking into single steps.

sim_name = 'ns_const_state'

-- Polynomial degree used in spatial approximation
degree = 2

-- Isothermal coefficient
isen = 1.4
-- Gas constant
R    = 287.0
-- Density
rho  = 1.4
-- Pressure
p    = 1.0
--Mach number
mach = 0.1

-- Fluid velocity
vel  = mach * math.sqrt(isen * p / rho )

-- logging, with higher level, the solver writes out
-- more information regarding the settings.
logging = { level = 8 }

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraition
-- information should be written out to the output file.
sim_control = {
  time_control = {
    max = {iter=20000},
    min = 0.0,
    interval = { iter = 1000 },
  }
}
-- Check for Nans and unphysical values
check =  { interval = 1 }


-- Mesh configuration --
-- We use just a single element here and achieve this by using the internal
-- mesh definition for a line with an element_count of 1.
cube_length = 2.0
mesh = {
  predefined = 'line',
  origin = {
    (-1.0)*cube_length/2.0,
    (-1.0)*cube_length/2.0,
    (-1.0)*cube_length/2.0
  },
  length = cube_length,
  element_count = 1
}

-- Equation definitions --
-- To solve the 2D compressible Navier Stokes equations we need to
-- define the physical fluid parameters, and additionally an internal
-- penalization factor, which is used in the implementation of the viscous
-- fluxes.
equation = {
  name      = 'navier_stokes_2d',
  isen_coef = isen,
  r         = R,
  -- Viscous parameters
  therm_cond = 0.5,
  mu         = 2.0,
  ip_param   = 4.0,
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0},
    relax_temperature = 0.0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
-- For 2D simulations we need to use the modg_2d spatial scheme.
-- In the time discretization we use an explicit Euler scheme and
-- control the timestep length by the CFL condition for which we
-- need to define a factor for both, convective and viscous parts.
scheme = {
  spatial =  {
    name = 'modg_2d',
    m = degree,
  },
  temporal = {
    name = 'explicitEuler',
    control = {
      name     = 'cfl',
      cfl      = 0.1*(2*degree+1)^2/(2*(degree+1)^2),
      cfl_visc = 0.4*((2*degree+1)^2/(2*(degree+1)^2))^2
    }
  }
}

-- Projection type --
-- We employ the fast polynomial transformation (FPT)
-- with an oversampling factor of 2, so we use twice as many points as
-- spectral modes to evaluate nonlinear terms.
-- This helps to minimize aliasing errors.
projection = {
  kind = 'fpt',
  factor = 2.0
}

-- Define the inital conditions --
-- We need to set density, pressure and
-- the velocity in x and y direction
initial_condition = {
  density  = rho,
  pressure = p,
  velocityX = vel,
  velocityY = 0.0
}

-- Tracking --
-- We the the field variables momentum, density and energy in one
-- point.
-- The data is written every 1000 iteration, as defined in the
-- time_control part.
tracking = {
  label = 'track_const_state_l2p',
  variable = {'momentum','density','energy'},
  shape = {
    kind = 'canoND',
    object= {
      origin ={0, 0, 0}
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = {iter = 1000}
  },
  output = {
    format = 'ascii',
    ndofs = 1
  }
}
