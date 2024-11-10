-- Setup of a travelling shock in a tube along the z axis.
-- This example illustrates the use of stabilization filters in
-- Euler 3D.
-- It simulates a shock that is travelling along the z-axis in
-- a domain with a local refinement in some part.
-- Slipwalls enclose the tube and inflow and outflow boundary
-- conditions are used at its ends.
--
require('seeder')

-- global simulation options
simulation_name = 'euler_3d'
timing_file = 'timing.res'
sim_control = {
  time_control = {
    min = 0,
    max = {iter = 10},
    interval = {iter = 1}
  }
}

check = { interval = 1 }

-- Mesh definitions --
mesh = 'mesh/'

-- Equation definitions --
equation = {
  name      = 'euler',
  isen_coef = 1.4,
  r         = 296.0,
  numflux   = 'godunov',
  material = {
    characteristic = 0,
    relax_velocity = {0, 0, 0},
    relax_temperature = 0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Tracking
tracking = {
  label = 'probe_density_covolume_multilevel',
  folder = './',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = {
        channel_width/2.,
        channel_width/2.,
        channel_length/3. + channel_width/20.
      }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.interval,
  },
  output = { format = 'ascii', ndofs = 1 }
}


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
gp1 = equation.isen_coef+1
gm1 = equation.isen_coef-1

chi = ( u_r - shockSpeed ) / math.sqrt(equation.isen_coef * p_r / rho_r )
rho_l = rho_r * ( (gp1*chi*chi) / (gm1*chi*chi+2) )
u_l = shockSpeed + ( u_r - shockSpeed ) * (rho_r/rho_l)
p_l = p_r * ( (2*equation.isen_coef*chi*chi-gm1) / gp1 )
mach_l = u_l/math.sqrt( equation.isen_coef * p_l / rho_l )

if (mach_l > 1)  then
  bkindinlet = 'supersonic_inflow_normal'
else
  bkindinlet = 'inflow_normal'
end

function rho(x,y,z)
  if ( z < channel_length/3.0 ) then
    return rho_l + y
  else
    return rho_r + x
  end
end

function p(x,y,z)
  if ( z < channel_length/3.0 ) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if ( z < channel_length/3.0 ) then
    return u_l
  else
    return u_r
  end
end

projection = {
  kind = 'l2p',
  nodes_kind = 'chebyshev',
  factor = 2.0
}

initial_condition = {
  density = rho,
  pressure = p,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = u
}

-- Scheme definitions --
filter_order = 14
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m = 3
  },
  temporal ={
    name = 'explicitSSPRungeKutta',
    steps = 2,
    control = {
      name = 'cfl',
      cfl  = 0.3
    }
  },

  ---- the stabilzation of the scheme
  stabilization = {
    {
      name = 'spectral_viscosity',
      alpha = 36,
      order = filter_order,
    },
    {
      name = 'covolume',
      beta = 1.0,
    }
  }
}


-- Boundary conditions
boundary_condition = {
  {
    label = 'inlet',
    kind = bkindinlet,
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
