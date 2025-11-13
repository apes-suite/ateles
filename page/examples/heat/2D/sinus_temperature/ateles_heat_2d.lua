-- Configuration file for Ateles --
sim_name = 'heat_2d_modg'

-- Variables to be set for the simulation --

-- Polynomial degree for the approximation
degree = 8
-- Constant for the sinus function
const = 2.0*math.pi
-- function to define the sinus changing temperature
function tempY(x,y,z)
  return( math.sin(const*y))
end

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraition
-- information should be written out to the output file.
sim_control = {
  time_control = {
    min = 0,
    max = 0.001,
    interval = {iter = 100}
  }
}

-- Mesh configuration --
-- The mesh is here a line, with starting
-- points origin and a length of linelength
-- and the number of elements of element_count

-- Refinement level for the simulation
-- domain.
level = 2
-- Length of the bounding cube
lineLength = 1.0
mesh = {
  predefined = 'slice',
  origin = {
    (-1.0)*lineLength/2.0,
    (-1.0)*lineLength/2.0,
    (-1.0)*lineLength/2.0
  },
  length = lineLength,
  refinementLevel = level
}

-- Equation definitions --
-- We use here Heat equation 1D
-- therefore we need to set the
-- thermal diffusity.
equation = {
  name   = 'heat_2d',
  k = 1,
}

-- Scheme definitions --
-- modg_2d results in a 2D simulation
-- the temporal table defines which time stepping
-- scheme is used here. In this test case we consider
-- the four step explite Runge-Kutta scheme.
-- The Cfl defines the timestep width.
scheme = {
  spatial =  {
    name = 'modg_2d',
    m = degree,
  },
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.07*0.4*((2*degree+1)/(degree+1))^4,
    },
  },
}

-- Projection type --
-- We consider here fast polynomial
-- transformation
projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32,
  fftMultiThread = false
}

-- Define the inital conditions --
-- We need to set the temperature
initial_condition = {
  temperature = tempY,
}

-- Tracking --
-- We track here a point (just origin is given)
-- and the quantity temperature. The interval defines
-- after how many iterations the quantity information
-- should be writen out.
tracking = {
  label = 'track_temp',
  variable = {'temperature'},
  shape={
    kind = 'canoND',
    object = {
      origin ={0.0, 0.0, -0.5}
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/15.0
  },
  output = {
    format = 'ascii',
    ndofs = 1 },
    folder = './'
  }

