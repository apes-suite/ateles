-- Riemann problem Toro1 in z direction for 3D Euler equations
require('seeder')

timing_file = 'timing.res'

-- global simulation options
simulation_name = 'toro1_z_euler_modg'
sim_control = {
  time_control = {
    min = 0,
    max = 0.004
  }
}

-- Mesh definitions --
mesh = 'mesh/'

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

-- Tracking
tracking = {
  label = 'probe_density_Q4_toro_z',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object = {
      origin = {
        epsx,
        epsx,
        (channel_length/2.0) + epsx,
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
  if ( z < channel_length/2.0 ) then
    return rho_l
  else
    return rho_r
  end
end

function p(x,y,z)
  if ( z < channel_length/2.0 ) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if ( z < channel_length/2.0 ) then
    return u_l
  else
    return u_r
  end
end

projection = {
  kind = 'fpt',
  factor = 1.0,
  blocksize = 32
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = 0,
  velocityY = 0,
  velocityZ = u
}

-- Scheme definitions --
degree = 3
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m =  degree
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.6*(3*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

-- Boundary conditions
boundary_condition = {
  {
    label = 'inlet',
    kind = 'inflow_normal',
    density = rho_l,
    v_norm = u_l,
    v_tan = 0.0
  },
  {
    label = 'outlet',
    kind = 'outflow',
    pressure = p_r
  },
  {
    label = 'bottom',
    kind = 'slipwall'
  },
  {
    label = 'top',
    kind = 'slipwall'
  },
  {
    label = 'south',
    kind = 'slipwall'
  },
  {
    label = 'north',
    kind = 'slipwall'
  }
}
