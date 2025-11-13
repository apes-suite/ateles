-- Euler 3D setup of a pulse in density with a sponge layer --
-- This example shows the use of a sponge to dampen out the state
-- and reduce reflections at boundaries.
-- The setup does not have boundaries, but a sponge layer is used
-- to dampen out the Gaussian pulse in density that is convected
-- into it.

-- ...the length of the cube
cubeLength = 2.0

logging = {level = 10}

-- the refinement level of the octree
level = 2

-- smallness parameter
eps = cubeLength/(2.0^(level+8))

-- Transport velocity of the pulse in x direction.
velocityX = 100

simulation_name = 'sponge_layer_modg'

-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = 0.01
  }
}

iniVel = velocityX
iniDens = 1.225

variable = {
  {
     -- Variable to describe the source term acting as a sponge.
     name = "spongelayer_var",
     ncomponents = 6,
     vartype = "st_fun",
     st_fun = {
       predefined = 'combined',
       spatial = {
         -- The predefined function spongelayer employs a plane to
         -- describe the area where the sponge is to be applied.
         -- The dampening increases along the normal of the plane
         -- until the final value is reached at the length of the
         -- normal vector.
         predefined = 'spongelayer_plane',
         origin = {0.4, 0, 0},
         normal = {0.6, 0, 0},
         damp_factor = 800,
         damp_profile = 'exponential',
         target_state = {
           density = iniDens,
           velocityX = iniVel,
           velocityY = 0.0,
           velocityZ = 0.0,
           pressure =100000
         }
       },
       temporal = 1.0, -- the sponge is constant over time.
       shape = {
         kind = 'canoND',
         object= {
           origin = {
             0.2+eps,
             -cubeLength/2.0 + eps,
             -cubeLength/2.0 + eps,
           },
           vec = {
             {0.8,0,0},
             {0,cubeLength,0},
             {0,0,cubeLength}
           },
           segments = {50,50,50}
         }
       }
     }
  }
}

-- Mesh definitions --
mesh = {
  predefined = 'cube',
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0
  },
  length = cubeLength,
  refinementLevel = level
}

-- Equation definitions --
equation = {
  name = 'euler',
  isen_coef = 1.4,
  r = 296.0,
  material = {
    characteristic = 0,
    relax_velocity = {0, 0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Sponge definition
source = {
  -- Use the variable defined above as sponge.
  spongelayer = "spongelayer_var"
}

-- Scheme definitions --
degree = 3
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m = 3
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.9*(3*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

-- This is a very simple example to define constant boundary condtions.
initial_condition = {
  density = {
    predefined = 'gausspulse',
    center     = {-0.5, 0.0, 0.0},
    halfwidth  = 0.20,
    amplitude  = 1.0,
    background = 1.225
  },
  pressure = 100000,
  velocityX = velocityX,
  velocityY = 0,
  velocityZ = 0
}

-- Tracking
tracking = {
  label = 'sponge',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object = {
      origin = {0.25, 0., 0.}
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/8.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
