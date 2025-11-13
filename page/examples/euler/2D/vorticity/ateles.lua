-- Configuration file for an example to track vorticity --
-- This setup shows the tracking of derived quantities in a single point.
simulation_name = 'vorticity_2d_modg'

-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = {iter =1, sim=0.1},
  }
}

-- Mesh definitions --
-- ...the length of the cube
cubeLength = 2.0
-- the refinement level of the octree
level = 5
-- smallness parameter
eps = cubeLength/(2^(level+8))
mesh = {
  predefined = 'slice',
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    0
  },
  length = cubeLength,
  refinementLevel = level
}

-- Equation definitions --
equation = {
  name   = 'euler_2d',
  therm_cond = 2.555e-02,
  isen_coef = 1.4,
  r      = 296.0,
  material = {
    characteristic = 0,
    relax_velocity = {0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
degree = 5
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_2d',
    modg_2d_space = 'Q',
    m = degree,
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.7*(2*degree+1)^2/(2*(degree+1)^2)
    },
  },
}

function velX(x,y,z,t)
  vel = 3*y*y*y
  return vel
end

function velY(x,y,z,t)
  vel = 4*x*x*x
  return vel
end

projection = {
  kind = 'l2p',
  factor = 1.0
}

initial_condition = {
  density = 1.4,
  pressure = 1,
  velocityX = velX,
  velocityY = velY,
}

-- Tracking of derived quantities like velocity and vorticity
tracking = {
  label = 'vort',
  folder = '',
  variable = {'density', 'momentum', 'velocity', 'vorticity'},
  shape = {
    kind = 'canoND',
    object= { origin = {0.65625, -0.65625, 0.03125} }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max.sim,
    interval = {iter = 1}
  },
  output = { format = 'ascii', use_get_point = true }
}
