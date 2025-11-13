-- This setup simulates a 1D Riemann problem in X direction in the 2D Euler
-- equations.
simulation_name = 'toro1_x_euler_modg_2d'

-- global simulation options
sim_control = {
  time_control = {
    min = 0,
    max = 0.025
  }
}

-- Mesh definitions --
channel_length = 1.0
mesh = {
  predefined = 'line_bounded',
  origin = {0, 0, 0},
  length = channel_length,
  element_count = 100
}

-- Equation definitions --
equation = {
  name   = 'euler_2d',
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


-- Tracking
tracking = {
  label = 'probe_density_Q4_toro_x',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = {
        (channel_length/2.0) + 0.0001,
        0.0001,
        0.0001
      }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- The initial conditions for the Riemann problem
-- ... left state
rho_l = 1.0
u_l = 0.0
p_l = 1.0
-- ... right state
rho_r = 0.125
u_r = 0.0
p_r = 0.1

function rho(x,y,z)
  if ( x < channel_length/2.0 ) then
    return rho_l
  else
    return rho_r
  end
end

function p(x,y,z)
  if ( x < channel_length/2.0 ) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if ( x < channel_length/2.0 ) then
    return u_l
  else
    return u_r
  end
end

function velY(x,y,z)
    return 0.0
end

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = u,
  velocityY = 0
}

-- Scheme definitions --
degree = 3
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_2d',
    m =  degree
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.6*(2*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

-- Boundary conditions
boundary_condition = {
  {
    label = 'west',
    kind = 'inflow_normal',
    density = rho_l,
    v_norm = u_l,
    v_tan = 0.0
  },
  {
    label = 'east',
    kind = 'outflow',
    pressure = p_r
  }
}
