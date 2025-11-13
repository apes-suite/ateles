-- Stabilized shock setup for 2D Euler equations.
-- This setup simulates a planar shock that travels along the y axis.
-- The simulation is stabilized by a covolume filter with spectral viscosity.

require('seeder')
-- seeder.lua sets a different timing file name, set it to the usual one again.
timing_file = 'timing.res'

-- global simulation options
simulation_name = 'euler_2d'

sim_control = {
  time_control = {
    min = 0,
    max = 0.015
  }
}

check = { interval = 1 }


-- Mesh definitions --
mesh = 'mesh/'

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

-- The state right of the shock
rho_r = 1.0
u_r = 0.0
p_r = 1.0
mach_r = u_r/math.sqrt( equation.isen_coef * p_r / rho_r )

-- Shock properties
shockMach = 2.0
shockXCoord = -1.2
shockSpeed = shockMach * math.sqrt(equation.isen_coef * p_r / rho_r )

-- The state left of the shock (evaluated by Rankine-Huginoit condition)
gm1 = equation.isen_coef-1
gp1 = equation.isen_coef+1

chi = ( u_r - shockSpeed ) / math.sqrt(equation.isen_coef * p_r / rho_r )
rho_l = rho_r * ( (gp1*chi*chi) / (gm1*chi*chi+2) )
u_l = shockSpeed + ( u_r - shockSpeed ) * (rho_r/rho_l)
p_l = p_r * ( (2*equation.isen_coef*chi*chi-gm1) / gp1 )
mach_l = u_l/math.sqrt( equation.isen_coef * p_l / rho_l )

function rho(x,y,z)
  if ( y < channel_length/3.0 ) then
    return rho_l
  else
    return rho_r
  end
end

function p(x,y,z)
  if ( y < channel_length/3.0 ) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if ( y < channel_length/3.0 ) then
    return u_l
  else
    return u_r
  end
end

projection = {
  kind = 'fpt',
  factor = 2.0,
  blocksize = 32
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = 0.0,
  velocityY = u
}

-- Scheme definitions --
filter_order = 14
degree = 4
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_2d',
    m = degree
  },
  ---- the stabilzation of the scheme
  stabilization = {
    {
      name = 'spectral_viscosity',
      alpha = 36,
      order = filter_order
    },
    {
      name = 'covolume',
      alpha = 36,
      order = filter_order,
      beta = 1.0
    }
  },
  -- temporal discretization
  temporal = {
    name = 'explicitRungeKuttaTaylor',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.3*(2*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

-- Boundary conditions
boundary_condition = {
  {
    label = 'inlet',
    density = rho_l,
    v_norm = u_l,
    v_tan = 0.0,
    pressure = p_l
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

if (mach_l > 1)  then
  boundary_condition[1].kind = 'supersonic_inflow_normal'
else
  boundary_condition[1].kind = 'inflow_normal'
end

-- Tracking
tracking = {
  label = 'probe_density_Q4_covolume_rktaylor_y',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = { epsx, (channel_length/2.0) + epsx, epsx }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
