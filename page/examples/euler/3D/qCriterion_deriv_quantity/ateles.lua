-- Tracking of q-criterion and lambda 2 in Euler 3D --
-- This setup uses some initial polynomial velocity profiles in all
-- directions and tracks the q-criterion and lambda 2 at a single
-- point in the domain.

simulation_name = 'q_crit_3d'

-- ...the length of the cube
cubeLength = 2.0

logging = {level = 10}

-- the refinement level of the octree
level = 4

-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = {iter = 2},
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
  name   = 'euler',
  isen_coef = 1.4,
  r      = 296.0,
  material = {
    characteristic = 0,
    relax_velocity = {0,0,0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m = 4
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl = 0.9,
      use_modal_estimate = true
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

function velX(x,y,z,t)
  vel = y*y*y + z*z*z + 2.0*x*x*x
  --print(vel)
  return vel
end

function velY(x,y,z,t)
  vel = x*x*x + z*z*z + 2.0*y*y*y
  --print(vel)
  return vel
end

function velZ(x,y,z,t)
  vel = x*x*x + y*y*y + 2.0*z*z*z
  --print(vel)
  return vel
end

initial_condition = {
  density = 1.4,
  pressure = 1,
  velocityX = velX,
  velocityY = velY,
  velocityZ = velZ,
}

-- Tracking
tracking = {
  label = 'vort',
  folder = '',
  variable = { 'q_criterion', 'lambda2' },
  shape = {
    kind = 'canoND',
    object= { origin = { 0.6875, 0.3125, -0.4375 } }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max.sim,
    interval = { iter = 1 }
  },
  output = { format = 'ascii', use_get_point = true }
}
