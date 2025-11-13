-- Configuration of a convected Gaussian pulse in density --
-- This simple setup illustrates the use of a sponge.
simulation_name = 'sponge_2d_modg'

-- ...the length of the cube
cubeLength = 2.0

logging = {level = 10}

-- the refinement level of the octree
level = 3

-- Transport velocity of the pulse in x direction.
velocityX = 100

dx = cubeLength/(2^(level))

-- smallness parameter
eps = cubeLength/(2.0^(level+8))


-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = {sim = 0.008}
  }
}

-- Provide the definition of the Sponge in the form of a variable, that
-- will be added to the variable system in the solver.
-- A predefined spongelayer is available for this, it makes use of a plane
-- to describe where the sponge is to start and how long it should be (given
-- by the length of the normal vector).
variable = {
  {
    name = "spongelayer_var",
    ncomponents = 5,
    vartype = "st_fun",

    st_fun = {
      predefined = 'combined',
      spatial = {
        predefined = 'spongelayer_plane_2d',
        origin = { 0.25, 0.0, 0,0 },
        normal = { 0.6, 0, 0 },
        damp_profile = 'exponential',
        damp_factor = 1500,
        target_state = {
          density = 1.225,
          velocityX = velocityX,
          velocityY = 0.0,
          pressure = 100000
        }
      },
      temporal = 1.0,
      shape = {
        kind = 'canoND',
        object= {
          origin = { 0.2, -cubeLength/2.0, -cubeLength/2.0 },
          vec = {
            { 0.8, 0, 0 },
            { 0, cubeLength, 0 },
            { 0, 0, dx }
          },
          segments = { 50, 50, 3 }
        }
      }
    } -- st_fun

  }
}

-- Mesh definitions --
mesh = {
  predefined = 'slice',
  -- The slice describes a plane of elements in the complete octree.
  -- There are 2^refinementLevel elements in x and y direction, and 1 element
  -- in z direction.
  -- The mesh is periodic in each direction.
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0
  },
  length = cubeLength,
  refinementLevel = level
}

timing = { filename = 'timing.res' }

-- Equation definitions --
equation = {
  name   = 'euler_2d',
  therm_cond = 2.555e-02,
  isen_coef = 1.4,
  r      = 296.0,
  material = {
    characteristic = 0.0,
    relax_velocity = {0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
degree = 3
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_2d',
    m = degree,
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.9*(2*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32,
  fftMultiThread = false
}

function ic_gauss(x,y,z)
  -- A circular Gaussian pulse
  d = x*x + y*y
  return( 1.225 + 2* math.exp(-d/0.1*math.log(2.0)) )
end

-- This is a very simple example to define constant boundary condtions.
initial_condition = {
  density = ic_gauss,
  pressure = 100000,
  velocityX = velocityX,
  velocityY = 0,
}

-- Tracking
-- This tracks a single element and reports the average density of that element.
-- The element is selected by the defined point in the shape subtable.
tracking = {
  label = 'sponge_2d',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= { origin = {0.55 , 0, -0.875} }
  },
  time_control = {
    min = 0,
    interval =  sim_control.time_control.max.sim/10.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- Sponge definition
-- The sponge is implemented as a source term and its definition is given in the
-- form of a variable, defined above in the `variable` table.
source = {
  spongelayer = "spongelayer_var"
}
