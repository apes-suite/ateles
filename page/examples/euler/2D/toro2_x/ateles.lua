-- Toro2 setup solving a Riemann problem in the Euler 2D equations.
-- This setup uses a Finite Volume discretization (polynomial degree=0) and the
-- Godunov numerical flux.
-- For the timeintegration the two-stage, strong stability preserving
-- Runge-Kutta scheme is used.
--
-- The Riemann problem prescribed here makes use of velocities with opposite
-- signs, while all other quantities remain the same.

-- global simulation options
simulation_name = 'toro2_euler_modg_2d'
sim_control = {
  time_control = {
    min = 0,
    max = 0.0015
  }
}

-- Mesh definitions --
channel_length = 1.0
epsx = channel_length*1.0e-6
mesh = {
  predefined = 'line_bounded',
  origin = {0, 0, 0},
  length = channel_length,
  element_count = 1000
}

-- Equation definitions --
equation = {
  name   = 'euler_2d',
  isen_coef = 1.4,
  r      = 296.0,
  numflux = 'godunov',
  material = {
    characteristic = 0,
    relax_velocity = {0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)


-- The initial conditions for the Riemann problem
-- ... left state
rho_l = 1.0
u_l = -2.0
p_l = 0.4
-- ... right state
rho_r = 1.0
u_r = 2.0
p_r = 0.4

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

projection = {
  kind = 'l2p',
  factor = 1.0
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = u,
  velocityY = 0
}

-- Scheme definitions --
degree = 0
scheme = {
  -- the spatial discretization scheme
  spatial = {
    name = 'modg_2d',
    m = degree
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitSSPRungeKutta',
    steps = 2,
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

-- Tracking
tracking = {
  label = 'probe_density_Q1_toro',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = {
        (channel_length/2.0) + epsx,
        epsx,
        epsx
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
