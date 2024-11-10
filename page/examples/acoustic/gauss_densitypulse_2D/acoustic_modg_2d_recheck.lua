-- Configuration file for Ateles --

-- Variables to be set for the simulation --

-- Polynomial degree for the approximation
degree = 12

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to 
-- start the simulation and the interval, after how many iteraition 
-- information should be written out to the output file.
sim_control = { 
  time_control = { 
    max = 0.004,
    min = 0,
    interval = {iter = 100} 
  }
}
-- Check for Nans and unphysical values
check =  { interval = 100 }

-- Mesh configuration -- 
-- The mesh is here a line, with starting 
-- points origin and a length of cube_length
-- and the number of elements of element_count

-- Length of the bounding cube
cubeLength = 2.0
-- Refinement level for the simulation 
-- domain.
level = 3 
mesh = {
  predefined = 'slice',
  origin = { 
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
  },
  length = cubeLength,
  refinementLevel = level
}

-- Equation definitions --
-- We use here Acoustics 2D 
-- therefore we need to set the
-- background information, which are
-- density, velocity and pressure
-- Equation definitions --
equation = {
  name   = 'acoustic_2d',
  background = {
    density = 1.225, 
    velocityX = 0.0,
    velocityY = 0.0,
    pressure = 100000.0
  }
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
    m =  degree,    
    modg_space = poly_space
  }, 
  temporal = {
    name = 'explicitRungeKutta', 
    steps = 4,
    control = {
      name = 'cfl',   
      cfl  = 0.95,
    },
  },
}

-- Projection type --
-- We consider here fast polynomial 
-- transformation
projection = {
  kind = 'fpt',  
  factor = 1.0,   
}

-- variables for gaussian pluse             
c = math.sqrt(equation.background.pressure / equation.background.density)
ampl_density= equation.background.density/c 

-- Gaussian pulse, prescribed at the initail condition
function ic_gauss_density(x,y)
d= x*x+y*y 
return( ampl_density * math.exp(-d/0.002*math.log(2)) )
end

-- Define the inital conditions --
-- We need to set density, pressure and 
-- the velocity in x and y direction
initial_condition = { 
  density = ic_gauss_density,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
}

-- Tracking --
-- We track here a point (just origin is given)
-- and the quantity density. The interval defines
-- after how many iterations the quantity information 
-- should be writen out.

-- A small fraction to shift the point inside 
-- the domain, if it might lay outside of the domain
eps=cubeLength/(2^(level+1))
tracking = {
  label = 'track_density',
  folder = './',
  variable = {'density'},
  shape = {
    kind = 'canoND', 
    object= { 
      origin = {0.5,0.5,0.0+eps}, 
    }
  },
  time_control = { 
    max = sim_control.time_control.max,
    min = 0,
    interval = sim_control.time_control.max/20.0
  },
  output = { 
    format = 'ascii',
    ndofs = 1 
  }
}
