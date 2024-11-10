-- Configuration for a viscous vortex in 2D --
-- This setup implements a vortex in 2D that is simulated by the
-- compressible Navier-Stokes equations with an additional right-hand side
-- term that allows us to compare numerical results to a known analytical
-- solution.

sim_name = 'tgv_nvrstk_modg_2d'

-- Variables to be set for the simulation --

-- Define density, velocity and pressure as
-- lua functions used in the initial condition
function ini_dens(x,y,z)
  return 4.0
end
function ini_vel_x(x,y,z)
  return math.sin(x)*math.cos(y)
end
function ini_vel_y(x,y,z)
  return math.cos(x) * math.sin(y)
end
function ini_press(x,y,z)
  return (1.0/4.0)*(math.cos(2.0*x)+math.cos(2.0*y)) + 8.0
end

-- Define the source term to create the analytical solution
-- for validation. First we define the density, then the momentum
-- and lastly the energy.
function source_fun_dens(x,y)
  return 8*math.cos(x)*math.cos(y)
end
function source_fun_mom_x(x,y)
  return ( 16*math.cos(x)*math.sin(x)*(math.cos(y)^2)
           + (16.0/3.0)*math.sin(x)*math.cos(y) - 2.5*math.sin(2*x) )
end
function source_fun_mom_y(x,y)
  return ( 16*math.cos(y)*math.sin(y)*(math.cos(x)^2)
           + (16.0/3.0)*math.sin(y)*math.cos(x) - 2.5*math.sin(2.0*y) )
end
function source_fun_en(x,y)
  return ( 0.0004464285714*math.cos(2.0*x)
           + 0.0004464285714*math.cos(2.0*y)
           + 48*math.cos(x)*math.cos(y)
           + (40.0/3.0)*(math.cos(x)^2)
           + (40.0/3.0)*(math.cos(y)^2)
           + 16*math.cos(x)*(math.cos(y)^3)
           + 16*(math.cos(x)^3)*math.cos(y)
           - (64.0/3.0)*(math.cos(x)^2)*(math.cos(y)^2)
           - 24*(math.cos(x)^3)*(math.cos(y)^3)
           + 1.75*math.cos(2.0*x)*math.cos(x)*math.cos(y)
           + 1.75*math.cos(2.0*y)*math.cos(x)*math.cos(y)
           - 1.75*math.sin(2.0*x)*math.cos(y)*math.sin(x)
           - 1.75*math.sin(2.0*y)*math.cos(x)*math.sin(y)
           - 8.0 )
end
-- Bind all source function into one only, to
-- use just one call for all of them
function source_fun(x,y,z,t)
  src_dens = source_fun_dens(x,y)
  src_mom_x = source_fun_mom_x(x,y)
  src_mom_y = source_fun_mom_y(x,y)
  src_ener = source_fun_en(x,y)
  return {src_dens, src_mom_x, src_mom_y, src_ener}
end
-- logging, with higher level, the solver writes out
-- more information regarding the settings.
logging = {level=10}

-- Variables as space time functions --
-- In this table we can define all variables, as
-- space time function. We define here the lua
-- function defined previously as space time function
variable = {
  {
    name = "arbitrary_var",
    ncomponents = 4,
    vartype = "st_fun",
    st_fun = source_fun
  }
}

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraition
-- information should be written out to the output file.
sim_control = {
  time_control = {
    max =  {iter = 10},
  }
}
-- Check for Nans and unphysical values
check = { interval = 1 }

-- Mesh configuration --
-- length of bounding cube
cubeLength = 2.0*math.pi
-- the refinement level of the octree
level = 2
mesh = {
  predefined = 'slice',
  origin = { 0.0, 0.0, 0.0 },
  length = cubeLength,
  refinementLevel = level
}

-- Equation definitions --
-- We solve the Navier-Stokes 2D equations
-- Besides the fluid parameters, we need to set an internal penalization
-- parameter that is used for the implementation of the viscous fluxes.
equation = {
  name       = 'navier_stokes_2d',
  isen_coef  = 1.4,
  r          = 280.0,
  -- viscous parameters
  therm_cond = 0.5,
  mu         = 2.0,
  ip_param   = 4.0,
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0},
    relax_temperature = 0.0
  }
}
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
-- modg_2d results in a 2D simulation
-- the temporal table defines which time stepping
-- scheme is to be used.
-- Here we use the classical explicit RK4 time integration scheme.
-- We adapt the timestep length according to the CFL condition, for
-- which we set the Courant factor.
-- For the viscous terms we use a second factor (cfl_visc), the timestep
-- is taken such, that both conditions, convective and viscous are satisfied.
scheme = {
  spatial =  {
    name = 'modg_2d',
    modg_2d_space = 'P',
    m = 7
  },
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl = 0.1,
      cfl_visc = 0.4,
    }
  }
}

-- Projection type --
-- We consider here fast polynomial
-- transformation with an oversampling
-- factor of 2, which means that we
-- use two times more points for the
-- approximation
projection = {
  kind = 'fpt',
  factor = 1.0
}

-- Source terms --
-- Activate the source term, this function will be added to the right hand
-- side of the equation.
-- Here we use it to achieve a known analytical solution.
source = {
  arbitrary = "arbitrary_var"
}

-- Define the inital conditions --
-- We need to set density, pressure and the velocity in x and y.
initial_condition = {
  density = ini_dens,
  pressure = ini_press,
  velocityX = ini_vel_x,
  velocityY = ini_vel_y,
  useFpt = true
}

-- Tracking --
-- We track here a the momentum in one element (selected by the defined point).
-- After each iteration, the current value of the tracked variable will be
-- written on a new line in the corresponding file.
tracking = {
  label = 'track_momentum',
  variable = { 'momentum' },
  shape = {
    kind = 'canoND',
    object= {
      origin ={ 1e-5, 1e-5, 1e-5 }
    }
  },
  time_control = {
    -- We write the tracked variable every iteration over the complete
    -- simulation time.
    min = { iter = 0 },
    max = { sim_control.time_control.max.iter },
    interval = { iter = 1 }
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- Write out a timing file which includes
-- all timings of the solver
timing = { filename = 'timing.res' }
