-- This setup illustrates the use of a modal estimation for the adaptive
-- timestep computation for the 1D Euler equations.
-- ----------------------------------------------- --
-- --------------- General options --------------- --
simulation_name = 'modalest_1d'
sim_control = {
  time_control = {
    min = 0,
    max = {iter=20},
    interval = {iter=1}
  }
}

p_ref   = 101325 -- Reference pressure in Pascal
T_ref   = 288.15 -- Reference temperature in Kelvin
rho_ref = 1.225  -- Reference density in kg/m^3
c_ref   = 340    -- Reference speed of sound
-- --------------- General options --------------- --
-- ----------------------------------------------- --


-- ------------------- --
-- ------ Mesh ------- --
mesh = {
  predefined = 'line',
  origin = { 0, -0.125, -0.125 },
  length = 1,
  refinementLevel = 2
}
-- ------ Mesh ------- --
-- ------------------- --


-- --------------------------------------------------------- --
-- ----------------------- Equation ------------------------ --
equation = {
  name = 'euler_1d',
  isen_coef = c_ref^2*rho_ref/p_ref,
  r = p_ref/(rho_ref*T_ref),
  material = {
    characteristic = 0,
    relax_velocity = 0,
    relax_temperature = 1
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)
-- ----------------------- Equation ------------------------ --
-- --------------------------------------------------------- --


-- ----------------------- Scheme -------------------------- --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg_1d',
    m = 63
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',
      use_modal_estimate = true,
      cfl  = 0.9
    }
  }
}

projection = {
  kind = 'l2p',
  factor = 1.0
}
-- ----------------------- Scheme -------------------------- --
-- --------------------------------------------------------- --


-- ------ Initial conditions ------- --
initial_condition = {
  density = {
    predefined = 'gausspulse',
    center = { 0.5, 0, 0 },
    halfwidth = 0.1,
    amplitude = 0.1,
    background = 1
  },
  pressure = {
    predefined = 'gausspulse',
    center = { 0.5, 0, 0 },
    halfwidth = 0.1,
    amplitude = 0.1,
    background = 1
  },
  velocity = 0.0
}
-- ------ Initial conditions ------- --
-- --------------------------------- --


-- -------------------- Tracking ---------------------- --
tracking = {
  label = 'point_series',
  folder = './',
  variable = { 'density', 'momentum', 'energy' },
  shape = {
    kind = 'canoND',
    object= { origin = { 0.7, 0, 0} }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = {iter=1}
  },
  output = { format = 'ascii', ndofs = 1 }
}
-- ---------------------------------------------------- --
