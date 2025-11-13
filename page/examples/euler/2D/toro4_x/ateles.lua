-- Toro4 setup solving a Riemann problem in the Euler 2D equations.
-- This setup uses a higher order and the positivity preserving filter and
-- oversampling to compute a strong shock.
-- For the timeintegration the two-stage, strong stability preserving
-- Runge-Kutta scheme is used.
--

-- The initial conditions for the Riemann problem
-- ... left state
rho_l = 0.125
u_l = 0.0
p_l = 0.1
-- ... right state
rho_r = 2.0
u_r = 0.0
p_r = 2.0

-- global simulation options
simulation_name = 'toro4_euler_modg_2d'
sim_control = {
  time_control = {
    min = 0,
    max = 0.012,
    interval = {iter = 1}
  }
}

check = {interval = 1}

-- Mesh definitions --
channel_length = 1.0
epsx = channel_length*1.0e-6
mesh = {
  predefined = 'line_bounded',
  origin = {0, 0, 0},
  length = channel_length,
  element_count = 80
}

-- Equation definitions --
equation = {
  name = 'euler_2d',
  isen_coef = 1.4,
  r = 296.0,
  numflux = 'godunov',
  material = {
    characteristic = 0,
    relax_velocity = {0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)


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
    name = 'explicitSSPRungeKutta',
    steps = 2,
    control = {
      name = 'cfl',
      cfl = 0.9
    }
  },
  stabilization = {
    { name = 'spectral_viscosity',
      alpha = 36,
      order = 30
    },
    { name = 'cons_positivity_preserv',
      eps = 1.0e-6
    }
  },
}

projection = {
  kind = 'l2p',
  factor = 3.0,
  nodes_kind = 'chebyshev',
  lobattoPoints = true
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = u,
  velocityY = 0
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
  label = 'probe_density_toro',
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
