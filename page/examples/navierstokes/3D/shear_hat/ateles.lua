-- Setup of a shear hat in 3D Navier-Stokes --
--
-- A simple setup with an initial y-velocity profile that forms a hat
-- around the y-z plane in a fully periodic domain.
sim_name = 'shear-hat'

-- Variables to be set for the simulation --

-- Polynomial degree used for approximation
degree = 3
-- Isothermal coefficient
isen = 1.4
--Boltzman constant
boltz = 287.0
-- Density
rho = 1.4
-- Pressure
p = 1.0
--Mach number
mach = 0.1
-- Viscousity parameters
conducT = 0.5
mu = 2.0
c2 = 8./3.
-- Spatial function to describe the y-velocity profile
vel  = function (x,y,z)
  if x > 0 then
    return mach * math.sqrt(isen * p / rho ) * (1. - (x-0.5*dx)/dx)
  else
    return mach * math.sqrt(isen * p / rho ) * (1. + (x+0.5*dx)/dx)
  end
end

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraition
-- information should be written out to the output file.
sim_control = {
  time_control = {
    max      = { iter = 400 },
    min      = 0.0,
    interval = { iter = 1 },
  }
}

-- Mesh configuration --
-- We just use two elements along the x-axis.
-- Element size for the mesh
dx = 1.0e-4
mesh = {
    predefined = 'line',
    origin = { -dx, -0.5*dx, -0.5*dx },
    length = 2*dx,
    element_count = 2
}

-- Equation definitions --
-- For the 3D Navier-Stokes equations we need to define the
-- properties of the fluid and an internal penalization material
-- that is used in the implementation of the viscous fluxes.
equation = {
  name      = 'navier_stokes',
  isen_coef = isen,
  r         = boltz,
  -- Viscous parameters
  therm_cond = conducT,
  mu         = mu,
  ip_param   = c2*3*(degree+2)/(2*(degree+3)),
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0, 0.0},
    relax_temperature = 0.0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
-- In the spatial discretization we have to use the modg scheme for 3D.
-- In time we use the classical Runge-Kutta 4 stage scheme with adaptive
-- timesteps that are chosen according to the CFL condition, given by
-- the courant factor as 'cfl'.
-- The viscous timestep limit is considered by the factor in 'cfl_visc'.
scheme = {
  spatial =  {
    name = 'modg',
    m = degree
  },
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name     = 'cfl',
      cfl      = 0.8,
      cfl_visc = 0.4
    }
  }
}

-- For the projection to nodal values in the evaluation of nonlinear terms
-- we use the FPT method.
projection = {
  kind = 'fpt',
  factor = 1.0
}

-- Define the inital conditions for all primitive state variables.
-- The velocity profile in y direction is defined by a Lua function that
-- describes the linear profile left and right of the yz plane.
initial_condition = {
  density   = rho,
  pressure  = p,
  velocityX = 0.0,
  velocityY = vel,
  velocityZ = 0.0
}

-- Tracking --
-- We track the conservative state variables in a single point over time.
-- The values are written every 10 iterations into the given file.
tracking = {
  label = 'track_shearhat3D',
  folder = './',
  variable = {'momentum','density','energy'},
  shape = {
    kind = 'canoND',
    object= {
      origin ={0.01*dx, 0., 0.}
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = {iter = 10}
  },
  output = {
    format = 'ascii',
    use_get_point = true
  }
}
