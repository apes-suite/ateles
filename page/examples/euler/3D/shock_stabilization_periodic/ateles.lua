-- Simulation of a planar shock in 3D Euler equations
-- through a domain that is periodic in all directions.
-- global simulation options
simulation_name = 'euler_3d'
sim_control = {
  time_control = {
    min = 0,
    max = { iter = 5 }
  }
}

check = { interval = 1 }

-- Mesh definitions --
mesh = 'mesh/'

-- Equation definitions --
equation = {
  name   = 'euler',
  isen_coef = 1.4,
  r      = 296.0,
  material = {
    characteristic = 0,
    relax_velocity = {0, 0, 0},
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
gp1 = equation.isen_coef + 1
gm1 = equation.isen_coef - 1

chi = ( u_r - shockSpeed ) / math.sqrt(equation.isen_coef * p_r / rho_r )
rho_l = rho_r * ( (gp1*chi*chi) / (gm1*chi*chi+2) )
u_l = shockSpeed + ( u_r - shockSpeed ) * (rho_r/rho_l)
p_l = p_r * ( (2*equation.isen_coef*chi*chi-gm1) / gp1 )
mach_l = u_l/math.sqrt( equation.isen_coef * p_l / rho_l )

function rho(x,y,z)
  if ( z < 1/3.0 ) then
    return rho_l
  else
    return rho_r
  end
end

function p(x,y,z)
  if ( z < 1/3.0 ) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if ( z < 1/3.0 ) then
    return u_l
  else
    return u_r
  end
end

projection = {
  kind = 'l2p',
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
degree = 4
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
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
    --name = 'explicitRungeKutta',
    --steps = 4,
    name = 'explicitSSPRungeKutta',
    steps = 2,
    control = {
      name = 'cfl',
      --dt = 1e-5
      cfl  = 0.3*(3*degree+1)^2/(2*(degree+1)^2)
    }
  }
}

epsx = 1.e-3
tracking = {
  label = 'probe_density_Q4_periodic_covolume_z',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = { (1/2) + epsx, (1/2) + epsx, (1/2) + epsx }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max
  },
  output = { format = 'ascii', ndofs = 1 }
}
