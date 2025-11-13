-- Reflection of an acoustic pulse in a multilevel mesh and wall by penalization
--
-- This is a basic setup of an acoustic pulse with a 3D Euler simulation
-- in a mesh with different resolutions.
-- As we use penalization terms here, we employ an IMEX time integration scheme
-- that allows us to solve the stiff penalization terms implicitly.
-- The domain is of size = 1 and the material occupies half
-- of the domain. i.e in the interval (0.5,1), which is resolved coarsely
-- with just one element.
-- The left half is built by four refined elements and contains the fluid.
-- We define an acoustic pulse that runs in this part of the domain
-- against the wall and is to be reflected there.

dx = 0.5
eps = 1.e-10

wall_x = 0.5
wall_thickness = 0.5

-- ... the penalization parameters
porosity = 1.0
alpha_v = 1.e-18
alpha_t = 1.e-6

--Simulation
spatial_degree = 7
cfl = 0.66*(3*spatial_degree+1)^2/(2*(spatial_degree+1)^2)

--  parameters for IC
ini_vel_x = 0.0
ini_dens  = 1.4
ini_press = 1.0

-- equation parameters
isen  = 1.4
boltz = 1.0
-- background state
back_press = ini_press/isen                   --equation.background.pressure
back_dens = ini_dens                          --equation.background.density
ini_temp = back_press/back_dens/boltz
c = math.sqrt((isen*back_press)/back_dens)    --speed of sound

--GaussPulse parameter
--Amplitude
amp_dens  = 0.001
amp_press = 0.001
--Center
axis_x = 0.25
--Width Control
width_con = 0.01

-- global simulation options
simulation_name = 'matml_reflected_pulse'

sim_control = {
  time_control = {
    max =  { sim = 0.25 },
    interval = {iter=10},
  }
}

-- Mesh definitions --
mesh = 'mesh_'

function p_prime(x, amp)
  d = ((x-axis_x)^2)
  arg = (-d/width_con)* math.log(2)
  res =  amp*math.exp(arg)
  return res
end

function pulse_press (x,y,z)
  res = back_press + p_prime(x,amp_press)
  return res
end

function pulse_dens (x,y,z)
  res = back_dens + p_prime(x,amp_dens)
  return res
end

function pulse_velx (x,y,z)
  res =  p_prime(x,amp_press/(back_dens*c))
  return res
end

-- Now, define the initial conditions
initial_condition = {
  density = pulse_dens,
  pressure = pulse_press,
  velocityX = pulse_velx,
  velocityY = 0.0,
  velocityZ = 0.0
}

obstacle_shape = {
  kind = 'canoND',
  object = {
    origin = { wall_x + eps, -eps, -eps },
    vec = {
      { wall_thickness - 2*eps,       0.0,       0.0 },
      {                    0.0,  dx+2*eps,       0.0 },
      {                    0.0,       0.0,  dx+2*eps },
    },
    segments = { 20, 20, 20 }
  }
}

variable = {
  {
     name = 'characteristic',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = {
       -- The masking function is 0 except in the part of the domain defined
       -- by the shape above.
       { const = { 0.0 } },
       {
         const = {1.0},
         shape = obstacle_shape
       }
     }
  },

  {
     name = 'bk_press',
     ncomponents = 1,
     vartype = 'st_fun',
     st_fun = back_press
  },

  {
     name = 'pert_press',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'difference',
       input_varname = {'pressure', 'bk_press'}
     }
  }
}

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m    = spatial_degree
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'imexRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = cfl
    }
  }
}

-- Equation definitions --
equation = {
  name       = 'euler',
  -- general fluid's parameter
  isen_coef  = isen,
  r          = boltz,
  -- viscous parameters
  -- Parameters of the penalization
  porosity             = porosity,
  viscous_permeability = alpha_v*porosity,
  thermal_permeability = alpha_t*porosity,
  material = {
    characteristic    = 'characteristic',
    relax_velocity    = {0.0, 0.0, 0.0},
    relax_temperature = ini_temp
  }
}
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

projection = { kind = 'fpt' }

-- Boundary definitions
boundary_condition = {
  {
    label = 'east',
    kind = 'wall'
  },
  {
    label = 'west',
    kind = 'outflow',
    pressure = back_press,
  }
}

tracking = {
  label = 'microphone',
  folder = '',
  variable = {'density'},
  shape = {
    kind = 'canoND',
    object= {
      origin ={ 0.475, 0.25, 0.25 }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = {iter=1}
  },
  output = { format = 'ascii', use_get_point = true }
}
