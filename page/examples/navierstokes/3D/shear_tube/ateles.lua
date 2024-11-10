-- Shear-Tube example for Navier-Stokes 3D --
-- This example illustrates the use of a positiviy preserving filter.
-- It utilizes a cubical fully periodic domain with a cylindrical jet
-- tube along the x-axis surrounded by a fluid at rest.
simulation_name = 'shear_tube'

-- Variables to be set for the simulation --

-- Polynomial degree used for approximation
degree = 3
-- Pressure
press = 1
-- Density
dens = 1.4
--Mach number
mach = 0.8

-- Control the Simulation times --
-- Set the simulation time and when information should be written out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraitions
-- information should be written to the screen.
sim_control = {
  time_control = {
    min = 0,
    max = {iter = 15},
    interval = {iter = 1},
    check_iter = 1
  }
}

-- Mesh configuration -- 
-- This setup uses a fully periodic cubical mesh.
mesh = {
    predefined = 'cube',
    refinementLevel = 2,
    length = 0.125,
    origin = {0.0, -0.0625, -0.0625}
}

-- Equation definitions --
-- For the 3D Navier-Stokes equations we need to define the
-- properties of the fluid and an internal penalization material
-- that is used in the implementation of the viscous fluxes.
equation = {
  name       = 'navier_stokes',
  isen_coef  = 1.4,
  r          = 296,
  therm_cond = 5.92e-3,
  mu         = 4e-6, 
  ip_param   = 4, 
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0, 0.0},
    relax_temperature = 1.0
  }
}
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
-- In the spatial discretization we have to use the modg scheme for 3D.
-- In time we use a 2 stage strong stability preserving explicit
-- Runge-Kutta scheme with fixed timesteps given by 'dt'.
-- The simulation is stabilized by a positivity filter that ensures that
-- pressure and density remain positive everywhere.
-- To guarantee this, the strong stability preserving scheme and Lobatto
-- integration points are required.
scheme = {
  spatial = {
    name = 'modg',
    m    = degree
  },
  temporal = {
    name = 'explicitSSPRungeKutta',
    steps = 2,
    control = {
      name = 'fixed',
      dt = 1.25e-4
    }
  },
  stabilization = {
    name = 'cons_positivity_preserv',
    eps  = 1.0e-8
  }
}

-- For the projection to nodal values in the evaluation of nonlinear terms
-- we use the FPT method.
-- Instead of the internal only Chebyshev integration nodes we use
-- Lobatto points for the nodal representation.
-- This is required for the positivity preserving filter.
projection = {
  kind = 'fpt',
  factor = 2.0,
  lobattoPoints = true
}

-- Description of the Jet-Tube
jet_radius = mesh.length * 2^(-mesh.refinementLevel)
velAmpl = mach*math.sqrt(equation.isen_coef*press/dens)

velX = function(x,y,z)
  if (math.abs(y) < jet_radius and math.abs(z) < jet_radius) then
    return velAmpl
  else
    return 0.0
  end
end

-- Define the inital conditions --
-- We need to set density, pressure and 
-- the velocity in x, y and z direction
initial_condition = { 
  density = dens,
  pressure = press,
  velocityX = velX,
  velocityY = 0.0,
  velocityZ = 0.0
}

-- Tracking --
-- We track here a point (just origin is given)
-- and the quantities momentum, density and 
-- energy. The interval defines after how 
-- many iterations the quantity information 
-- should be writen out.
tracking = {
  label = 'track_momentum',
  folder = './',
  variable = { 'momentum' },
  shape = {
    kind = 'canoND',
    object= { 
      origin ={ jet_radius/2, jet_radius*3/2, jet_radius*3/2 } 
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = {iter = 1}
  },
  output = { 
    format = 'ascii', 
    use_get_point=true 
  },
}
