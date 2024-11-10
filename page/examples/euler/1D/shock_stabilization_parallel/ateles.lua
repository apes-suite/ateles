-- This setup illustrates the use of covolume with spectral filtering for
-- a single shock in 1D Euler equations moving from left to right through the
-- domain.
-- The simulation uses the 1D modal DG scheme in space and a Taylor
-- Runge-Kutta scheme in time.
simulation_name = 'euler_1d'

sim_control = {
  time_control = {
    min = 0,
    max = 0.015
  }
}

check = { interval = 1 }

-- Mesh definitions --
channel_length = 1
mesh = {
  predefined = 'line_bounded',
  origin = {0, 0, 0},
  length = channel_length,
  element_count = 8
}

-- Equation definitions --
equation = {
  name   = 'euler_1d',
  isen_coef = 1.4,
  r      = 296.0,
  material = {
      characteristic = 0.0,
      relax_velocity = 0.0,
      relax_temperature = 0.0
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- The state right of the shock (fluid at rest)
rho_r = 1.0
u_r = 0.0
p_r = 1.0

-- Shock properties
shockMach = 2.0 -- Shock moves with Mach 2 into the fluid at rest.
shockSpeed = shockMach * math.sqrt(equation.isen_coef * p_r / rho_r )

-- The state left of the shock (given by the Rankine-Huginoit condition for
-- the chosen shock Mach number)
gm1 = equation.isen_coef-1
gp1 = equation.isen_coef+1
chi = ( u_r - shockSpeed ) / math.sqrt(equation.isen_coef * p_r / rho_r )

rho_l = rho_r * ( (gp1*chi*chi) / (gm1*chi*chi+2) )
u_l = shockSpeed + ( u_r - shockSpeed ) * (rho_r/rho_l)
p_l = p_r * ( (2*equation.isen_coef*chi*chi-gm1) / gp1 )
mach_l = u_l/math.sqrt( equation.isen_coef * p_l / rho_l )

function rho(x,y,z)
  if (x < channel_length/3.0) then
    return rho_l
  else
    return rho_r
  end
end

function p(x,y,z)
  if (x < channel_length/3.0) then
    return p_l
  else
    return p_r
  end
end

function u(x,y,z)
  if (x < channel_length/3.0) then
    return u_l
  else
    return u_r
  end
end

projection = {
  kind = 'fpt',
  factor = 2.0,
  blocksize = 32,
}

initial_condition = {
  density = rho,
  pressure = p,
  velocity = u,
}

-- Scheme definitions --
filter_order = 14
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_1d',
    m = 4,
  },
  ---- the stabilization of the scheme
  stabilization = {
    {
      name = 'spectral_viscosity',
      alpha = 36,
      order = filter_order,
    },
    {
      name = 'covolume',
      alpha = 36,
      order = filter_order,
      beta = 1.0,
    }
  },
  -- temporal discretization
  temporal = {
    name = 'explicitRungeKuttaTaylor',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.3
    }
  }
}

-- Boundary conditions
boundary_condition = {
  {
    label = 'west',
    density = rho_l,
    v_norm = u_l,
    pressure = p_l
  },
  {
    label = 'east',
    kind = 'outflow',
    pressure = p_r
  }
}

if (mach_l > 1)  then
  boundary_condition[1].kind = 'supersonic_inflow_normal'
else
  boundary_condition[1].kind = 'inflow_normal'
end

-- Tracking
tracking = {
  label = 'probe_density_Q4_covolume_rktaylor_x',
  folder = './',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin = { (channel_length/2.0) + 0.001, 0.001, 0.001 }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
