-- Tracking Vorticity in 3D Euler equations --
simulation_name = 'vorticity_modg'

-- ...the length of the cube
cubeLength = 2.0

logging = {level = 3}

-- the refinement level of the octree
level = 4

-- smallness parameter
eps = cubeLength/(2^(level+8))


-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = {iter = 1, sim = 0.1}
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
  name      = 'euler',
  isen_coef = 1.4,
  r         = 296.0,
  material  = {
    characteristic = 0,
    relax_velocity = {0, 0, 0},
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
    name = 'explicitRungeKutta',  --'explicitEuler',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      cfl  = 0.9
    }
  }
}

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

function velX(x,y,z,t)
  vel = y*y*y + z*z*z
  return vel
end

function velY(x,y,z,t)
  vel = x*x*x + z*z*z
  return vel
end

function velZ(x,y,z,t)
  vel = x*x*x + y*y*y
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
  folder = './',
  variable = {'momentum', 'velocity', 'vorticity'},
  shape = { kind = 'canoND', object = { origin = {0.6875, 0.3125, -0.4375} } },
  time_control = {
    min = 0,
    max = sim_control.time_control.max.sim,
    interval = sim_control.time_control.max.sim/10.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
