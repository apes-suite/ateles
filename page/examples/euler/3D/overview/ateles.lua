-- Example Ateles configuration for simulation of Euler3D equations --
--
-- This configuration file illustrates the various settings and     --
-- options relevant for Euler 3D simulations.                       --
-- You may want to use it as a template for your configurations.    --

-- NOTE: This run is configured to run only 2 iterations, for       --
--       use you'll want to change that. See sim_control!           --

-- -------------------------------------------------------------------------- --
gamma = 1.4 -- isentropic expansion coefficient
R     = 1.0 -- gas constant

equation = { -- >>>>>>>>>>>>>>>
  name              = 'euler',
  isen_coef         = gamma,
  r                 = R,
  cv                = R/(gamma-1),

  numflux           = 'hll', -- Numerical flux to use between elements.
                             -- Available are 'hll', 'godunov' and
                             -- 'lax_friedrich'. 'lax_friedrich' is not
                             -- recommended, as it may cause instabilities,
                             -- especially at boundaries.

  linear_limit      = 0.001, -- Linearization limit, if this is larger than 0
                             -- elements, where the indicator is smaller than
                             -- this limit will be computed with a linearized
                             -- flux.
  linearization_indicator = 'error',
                             -- The estimator to use to decide, whether
                             -- linearization should be used. Only relevant if
                             -- linear_limit>0.
                             -- Available options are:
                             -- 'error': An estimation for the terms neglected
                             --          by linearization, involves all
                             --          variables. This is the default.
                             -- 'energy': Use deviation of energy from its
                             --           mean value to decide the
                             --           linearization.
                             -- 'density': Use deviation of density from its
                             --            mean to decide the linearization.

  ensure_positivity = false, -- Ensure that density and energy remain positive
                             -- by only considering higher modes up to the point
                             -- where positivity is guaranteed.

  porosity          = 1.0,   -- Porosity to use in material modelling for
                             -- wall representation in elements.
  viscous_permeability = 1.0e-6, -- Viscous permeability for the porous medium
                                 -- to represent wall geometries in elements.
  thermal_permeability = 1.0e-3, -- Thermal permeability for the material to
                                 -- represent walls.
  material = {
    -- Description of the material distribution to define obstacles inside the
    -- domain.
    characteristic = 0.0, -- Masking function (may be a variable) that describes
                          -- where Material is to be found (Chi(x,y,z)), should
                          -- be 1 inside material and 0 everywhere else.
    relax_velocity = {0, 0, 0}, -- Velocity of the obstacle.
    relax_temperature = 0.0     -- Temperature of the obstacle.
  }
} -- equation <<<<<<<<<<<<
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
hours = 3600
jobtime = 2 * hours -- Time available for the computation.
                    -- This variable is not a setting for Ateles, but
                    -- used below. You might want to get this from an
                    -- environment variable in your job script.
                    -- use os.getenv('varname') to get environment variables.

sim_control = { -- >>>>>>>>>>>>>>>
  -- Controlling the execution of the simulation.
  time_control = {
    -- Time definitions are provided in terms of:
    -- * sim:   simulation time, this is the time measurement for the
    --          transient phenomenon that is simulated.
    -- * iter:  the number of iterations (time steps) that were done
    -- * clock: the real time that has passed since the beginning of
    --          the simulation in seconds.
    -- All three measures might be set, whatever is encountered first,
    -- will trigger the setting:
    -- {sim = 1.23, iter = 123, clock=12.3}
    -- If instead of a table a single number is provided, this is
    -- interpreted as the setting for sim.
    -- All not-given times are set to never.
    max = { -- Point in time when to stop the simulation.
     sim   = 1.0,
     clock = jobtime - 5*60,
     iter  = 10  -- This is just for a quick check that the setup works
                -- for real runs, you'd take this away.
    },
    -- Providing a clock setting for the max time allows you to ensure that
    -- a restart file is written before the computation is ended by the
    -- scheduler. We stop here five minutes early to allow for the restart to
    -- be written. This also implies that single iterations should take less
    -- than 5 minutes. If they take longer, you might run into the end without
    -- writing a restart.
    min = 0.0,
    interval = {iter = 10}, -- This controls the output frequency to
                            -- report the current time step.
  },

  abort_criteria = {
    -- Criteria upon which the simulation should be stopped.
    stop_file = 'STOP' -- Stop file, can be used to signal the simulation to
                       -- gracefully end.
                       -- If this is empty (the default), this is disabled.
                       -- Non-empty settings here, will cause the simulation to
                       -- come to an end when a file of this name exists in the
                       -- working directory.
                       -- Thus, with stop_file = 'STOP', you can cause the
                       -- simulation to stop and write a restart by doing
                       -- "touch STOP" in the working directory of the run.
                       -- Empty files like the one created by touch, will be
                       -- deleted by Ateles. If you want to keep the file, it
                       -- needs to have some content. You can achieve this for
                       -- example by "echo keep > STOP".
  }

} -- sim_control <<<<<<<<<<<<
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
mesh = { -- >>>>>>>>>>>>
  -- Simple predefined cubical mesh, periodic in all directions.
  predefined = 'cube',
  origin = {0, 0, 0},
  length = 1,
  refinementLevel = 3
} -- mesh <<<<<<<<<<<<
-- Alternatively, for meshes from seeder you set the prefix (may be a path), to
-- where the mesh is found:
-- mesh = 'mesh/'
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
scheme = { -- >>>>>>>>>
  spatial =  {
    -- the spatial discretization scheme (needs to match dimensionality of
    -- equation)
    name = 'modg',
    modg_space = 'Q', -- How to build multidimensional polynomials.
    -- 'Q' means we include all modes up to the maximal polynomial degree
    -- in each direction, no matter the others, resulting in modes with
    -- three times the maximal polynomial degree.
    -- 'P' means that only modes are considered for which the sum of
    -- the polynomial degree in each direction does not exceed the maximal
    -- polynomial degree.
    m = 11  -- the maximal polynomial degree for each spatial direction
  },
  temporal = {
    -- the temporal discretization scheme
    name = 'explicitRungeKutta',
    steps = 4,
    -- Strong stability preserving explicit Runge Kutta:
    -- For stability preserving filters
    -- name = 'explicitSSPRungeKutta',
    -- steps = 2,
    -- Implicit-Explicit Runge Kutta (DIRK)
    -- If Material is used, this time scheme should be employed!
    -- name = 'imexRungeKutta',
    -- steps = 4,
    -- Explicit Taylor expanded Runge Kutta
    -- Higher order for linear autonomous equations.
    -- name = 'explicitRungeKuttaTaylor'
    -- steps = 10, -- arbitrary number of steps
    -- Explicit Euler (unstable for higher spatial orders):
    -- only for testing!
    -- name = 'explicitEuler',
    control = {
      -- Choice of time step width.
      name = 'cfl',  -- the name of the timestep control mechanism
                     -- (alternatively a fixed timestep may be used).
      cfl  = 0.9, -- Courant number
      -- Alternatively a fixed time step may be configured:
      -- name = 'fixed',
      -- dt = 0.123,
      -- Modal to nodal conversion for the timestep width can be avoided
      -- by using a modal estimate. This is especially relevant when using
      -- linearization.
      -- Less accurate (may result in smaller time steps).
      use_modal_estimate = false
    }
  }
} -- scheme <<<<<<<<<<<<<<
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
initial_condition = { -- >>>>>>>>>>>>>>
  -- As initial conditions all primitive variables (density, velocity and
  -- pressure) must be provided as spatial functions (f(x,y,z)).
  -- There are some predefined spatial functions, like the gaussian pulse that
  -- may be used. Constants and Lua functions are also fine.
  density = {
    predefined = 'gausspulse',
    center = { 0.5, 0.5, 0.5 },
    halfwidth = 0.20,
    amplitude = 2.0,
    background = 1.225
  },
  pressure = 1,
  velocityX = 1,
  velocityY = 0.0,
  velocityZ = 0.0
} -- initial_condition <<<<<<<<<<<<
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
boundary_condition = { -- >>>>>>>>>>
  -- This minimal example is completely periodic and has no boundary conditions.
  -- Each boundary condition is described by a kind, and space-time functions
  -- for all required variables in the that boundary.
  -- Space time functions are either predefined functions, Lua functions,
  -- constants or combined functions, where temporal and spatial parts are
  -- superimposed.
  -- If there are values that are to be extrapolated (Neumann boundary
  -- condition), you can set enforce_zero_grad to true, to use an extrapolation
  -- of a polynomial with zero gradient at the boundary.
  -- This is achieved by computing the last mode to fulfill this condition.
  -- If you set neumann_mode_fraction to a smaller value than 1, then only
  -- this fraction of lower modes will be used in the enforce_zero_grad
  -- procedure and higher modes will be set to 0.

  --  { -- SLIPWALL
  --    --   Velocity in normal direction 0, other values extrapolated.
  --    label = 'cylinder',
  --    kind  = 'slipwall', -- or 'wall'
  --    enforce_zero_grad = true,
  --    neumann_mode_fraction = 1.0
  --  },

  --  { -- PRIMITIVES
  --    --   Prescribe all primitive variables.
  --    label = 'outside',
  --    kind = 'primitives',
  --    density = 1.23,
  --    velocityX = 0.2,
  --    velocityY = 0.3,
  --    velocityZ = 0.4,
  --    pressure = 1
  --  },

  --  { -- CONSERVATIVES
  --    --   Prescribe all conservative variables.
  --    label = 'left',
  --    kind = 'conservatives',
  --    density = 1.23,
  --    momentumX = 0.2,
  --    momentumY = 0.3,
  --    momentumZ = 0.4,
  --    energy = 4
  --  },

  --  { -- INFLOW
  --    --   Prescribe density and velocity, extrapolate pressure.
  --    label = 'outside',
  --    kind = 'conservatives',
  --    density = 1.23,
  --    velocityX = 0.2,
  --    velocityY = 0.3,
  --    velocityZ = 0.4
  --    enforce_zero_grad = true,
  --    neumann_mode_fraction = 1.0
  --  },

  --  { -- INFLOW_NORMAL
  --    --   Prescribe density and velocity normal to boundary,
  --    --   extrapolate pressure.
  --    label = 'inlet',
  --    kind = 'inflow_normal',
  --    density = 1.23,
  --    v_norm = 0.5,
  --    enforce_zero_grad = true,
  --    neumann_mode_fraction = 1.0
  --  },

  --  { -- SUPERSONIC_INFLOW_NORMAL
  --    --   Prescribe all primitive variables with velocity normal to boundary.
  --    label = 'superin',
  --    kind = 'inflow_normal',
  --    density = 1.23,
  --    v_norm = 2.5,
  --    pressure = 1
  --  },

  --  { -- OUTFLOW
  --    --   Prescribe pressure, extrapolate all other variables.
  --    label = 'east',
  --    kind = 'outflow',
  --    pressure = 1
  --    enforce_zero_grad = true,
  --    neumann_mode_fraction = 1.0
  --  },

  --  { -- SUPERSONIC_OUTFLOW
  --    --   Extrapolate all variables.
  --    label = 'superout',
  --    kind = 'supersonic_outflow',
  --    enforce_zero_grad = true,
  --    neumann_mode_fraction = 1.0
  --  },

} -- boundary_condition <<<<<<<<<<
-- -------------------------------------------------------------------------- --


-- -------------------------------------------------------------------------- --
restart = { -- >>>>>>>>>>>>>>>
  -- Header of restart, to start simulation from.
  -- example:
  read = './restart_simulation_lastHeader.lua',

  -- If an initialization should be done, when the file given in the
  -- read setting above is missing, set the init_on_missing to true.
  -- With an init_on_missing = false, the computation will stop, if
  -- the file in read is not found. This option has no meaning when no
  -- read is provided.
  -- default: false
  init_on_missing = true,

  -- Prefix for the files, that are to be written during the
  -- simulation. If this ends with a path separator, the restart files
  -- will be written into the specified directory and that directory
  -- has to exist
  write = 'restart_', -- 'restart/' for a directory

  time_control = {
    -- Starting point after which restart files should
    -- be written.
    -- Setting iter to 0 here, results in restart files
    -- being written from the initial condition onwards.
    min = {iter = 0},

    -- The maximal point in time, up to which, restarts
    -- should be written.
    -- Note, that if this is not defined at all, it will
    -- be set to never, resulting in doing restarts for
    -- all times, including a final restart after
    -- reaching the termination of the time loop.
    max      = sim_control.time_control.max.sim,

    -- Frequency at which restart files are to be
    -- written between min and max.
    interval = 0.1,

    -- If a restart is read, it may be that the restart
    -- writing did not happen at the desired interval,
    -- nevertheless you might want to write the next
    -- restart file according to the original interval
    -- rythm. To achieve this you can tell the restart
    -- to align the trigger to last interval before
    -- the restart time with the align_trigger option:
    align_trigger = { sim = true },
    -- For each time component you can define whether
    -- an alignment should be done or not.
    -- If you only set one flag like:
    -- align_trigger = false
    -- It will be applied to the simtime, the others
    -- will be set to false.
    -- Default is false.

    -- Controlling if a restart file should be written
    -- or not involves communication, if a clock setting
    -- is used.
    -- This might have a negative impact on the
    -- performance, if it is done every iteration and
    -- a single iteration is too short.
    -- By setting check_iter the interval at which
    -- these checks are done, can be increased, and
    -- thus performance impacts reduced.
    -- CAUTION: be aware that increasing this to values
    --          larger than 1, decreases the accuracy
    --          of the points in time at which restarts
    --          are written!
    -- default:
    -- check_iter = 1
  }
} -- restart <<<<<<<<<<<<
-- -------------------------------------------------------------------------- --

-- -------------------------------------------------------------------------- --
tracking = { -- >>>>>>>>>>>>
  {
    label = 'QuantityOfInterest',
    variable = { 'density' },
    -- Available variables in Euler:
    -- * density
    -- * momentum
    -- * energy
    -- * pressure
    -- * velocity
    -- * speedOfsound:   local speed of sound
    -- * temperature:    temperature of the fluid
    -- * mach_number:    local Mach number
    -- * mach_vector:    local velocity vector scaled by the speed of sound
    -- * kinetic_energy: the kinetic energy of the fluid
    -- * gradv:          gradient of the velocity field
    -- * vorticity:      vorticity of the flow field
    -- * q_criterion:    Q-Criterion (positive second invariant of velocity
    --                                gradient tensor)
    -- * lambda2:        Lambda 2 criterion: largest eigenvalue of shear and
    --                   rotational contributions to the velocity gradient
    --                   tensor.
    -- * linindicator:   Indicator that is used to decide whether to just
    --                   compute the linearized Euler flux in the element.
    --                   This depends on the chosen linearization_indicator.
    --                   It will be 1 in elements that are computed nonlinearly
    --                   and 0 in elements where the linearized flux is used.
    shape = { kind='global' },
    time_control = {
      min = { iter = 0 },
      max = sim_control.time_control.max.sim,
      interval = sim_control.time_control.max.sim/10
    },
    folder = 'track_',
    output = { format = 'vtk' }
    --  output = { format = 'asciispatial',
    --             -- All degrees of freedom will be written, specify,
    --             -- ndofs = 1 for dumping the average
    --             ndofs=1 }
    --  -- output = {format = 'harvester'} -- write restart files
    -- Trackings can be reduced in space to single values:
    -- reduction = {'sum','average'},  -- (default = 'none')
  },
  {
    label = 'PointProbe',
    variable = { 'pressure' },
    shape = { kind='canoND', object = {origin={0.5,0.5,0.5}} },
    folder = 'pp_',
    output = { format = 'ascii',
               -- track the value exactly at the given point (instead of
               -- getting the complete element):
               use_get_point = true },
    time_control = {
      min = { iter = 0 },
      max = sim_control.time_control.max.sim,
      interval = {iter = 1}
    },
  }
} -- tracking <<<<<<<<<<<<<<<
-- -------------------------------------------------------------------------- --

-- -------------------------------------------------------------------------- --
ply_sampling = { -- >>>>>>>>>>>>>>
  -- Subsampling for tracking, define a ply_sampling table to activate
  -- subsampling for all tracking objects (except those with use_get_point).
  nlevels = 1,    -- maximal level to use in subsampling
                  -- defaults to 0, which deactivates subsampling

  --method  = 'fixed', -- method to use for subsampling
                       -- 'adaptive': (recommended default) adaptive refinement
                       --             of the mesh based on solution
                       -- 'fixed': will refine all elements by nlevels and
                       --          evaluate the polynomials at the barycenter of
                       --          each fine element

  -- Parameters for adaptive sampling:
  --tolerance = 0,                -- threshold for ignoring higher modes,
                                  -- if the sum of absolute values of higher
                                  -- modes in relation to the first mode is
                                  -- below the tolerance, they will be cut off.
                                  -- Default: 0 (never ignore modes).
                                  -- Recommended: 0.01 - 0.05
  --dof_reduction = 0.5,          -- Factor to multiply the 1D number of degrees
                                  -- of freedom with in each refinement.
                                  -- (0.5 implies to only keep half the degrees
                                  -- of freedom after the refinement).
                                  -- Values smaller 0.5 are not useful here and
                                  -- 0.5 already yields quite poor polynomial
                                  -- representations in the sampling.
                                  -- Default: 0.5
                                  -- Recommended: 2/3 - 3/4, 0.7 generally looks
                                  -- quite good.
  --adaptiveDofReduction = false, -- Increase the number of modes to keep after
                                  -- refinement above the configured
                                  -- dof_reduction factor if no more than the
                                  -- original memory will be used.
                                  -- Default: false
                                  -- Recommended: true (generally advisable, no
                                  --                    big drawbacks)
  --absUpperBoundLevel = 0        -- Maximal absolute level to refine to.
                                  -- 0 means no absolute level, all elements
                                  -- will be limited in the refinement by the
                                  -- nlevels set above
                                  -- This is useful for multilevel meshes, where
                                  -- each element might start out on a different
                                  -- level.
                                  -- For a complete mesh overview you might want
                                  -- to use a very high nlevels setting but
                                  -- limit the refinement by absUpperBoundLevel
                                  -- to allow for a more uniform resolution in
                                  -- the sampled data.
                                  -- Default: 0 (no absolute upper bound)
} -- ply_sampling <<<<<<<<<<<<
-- -------------------------------------------------------------------------- --

-- -------------------------------------------------------------------------- --
projection = { -- >>>>>>>>>>>
  -- Projection between modal to nodal space
  -- for BC, IC and non-linear problems is given
  -- kind = 'l2p',  -- 'fpt', 'fxt' or 'l2p',
  --                --  default: 'l2p'

  -- for l2p the nodes are by default 'gauss-legendre',
  -- can be changed by setting the nodes_kind:
  -- nodes_kind = 'chebyshev',

  -- for fpt the nodes are automatically 'chebyshev'
  -- for fxt the nodes are automatically 'gauss-legendre'

  -- lobattoPoints = false, -- if lobatto points should be used,
                            -- default = false,
                            -- only working for Chebyshev points!

  factor = 1.0,             -- dealising factor for fpt
                            -- oversampling factor to remove aliasing
                            -- effects by padding, default: 1
                            -- Note that for FXT an evenly oversampled
                            -- order is required, if this is not the
                            -- case, the next higher even order will be
                            -- used and the actual factor might
                            -- accordingly be higher.

  -- FXT settings:
  -- prec = 1.5e-8,         -- precision to use for the fast multipole
                            -- computation during initialization.
                            -- Defaults to sqrt of epsilon for double
                            -- precision numbers (1.5e-8).

  -- FPT settings:
  -- approx_terms = 18,     -- Number of terms used to approximate the
                            -- matrix multiplication for blocks, that
                            -- are detached from the diagonal.
                            -- The default of 18 is recommended for
                            -- double precision.
  -- blocksize = 64,        -- for FPT, default 64. The blocksize
                            -- defines how big the minimal block
                            -- should be that is approximated in
                            -- fast algorithm.
                            -- The smaller it is, the more operations
                            -- are merely approximated.
                            -- Recommended for double precision is a
                            -- setting of 64.
                            -- The fast algorithm will only be used
                            -- for m >= blocksize.
                            -- Note, that this has to be larger than
                            -- 2*approx_terms to provide any
                            -- reduction in operation counts.
  -- striplen = 256,        -- This provides the length for arrays to
                            -- apply the matrix operation to
                            -- simultaneously.
                            -- Default is the vlen from the
                            -- tem_compileconf_module.
  -- subblockingWidth = 8,  -- The subblockingWidth is used during the
                            -- unrolling of the diagonal multiplication
                            -- during the projection. By setting this
                            -- value to an appropriate value a better
                            -- cache usage can be achieved.
                            -- Default is 8
  -- adapt_factor_pow2 = true, -- for FPT, default false. Should the
                               -- oversampling factor be adjusted to
                               -- obtain a power of 2 in the
                               -- oversampled order?
  -- fftMultiThread = false -- for FPT, default false. Should nested
                            -- multithreading be activated for FFTW?

  -- There can be individual projection settings for various parts
  -- of the computation.
  -- These are useful to allow for a higher oversampling in specific
  -- parts.

  -- individual projection methods for source terms
  -- if none is provided, the general projection method is used
  source_terms = {
    -- the configuration parameter are similar to the
    -- general projection method
    -- kind = 'fpt',
    factor = 1.0
  },
  -- individual projection methods for initial condition
  -- if none is provided the general projection method is used
  initial_condition = {
    -- the configuration parameter are similar to the
    -- general projection method
    factor = 2.0
  },
  -- individual projection methods for boundary condition
  -- if none is provided the general projection method is used
  boundary_condition = {
    -- the configuration parameter are similar to the
    -- general projection method
    -- kind = 'fpt',
    factor = 2.0
  },
  -- individual projection methods for boundary condition
  -- if none is provided the general projection method is used
  --material= {
  --  -- the configuration parameter are similar to the
  --  -- general projection method
  --  -- kind = 'fpt',
  --  factor = 2.0
  --}
} -- projection <<<<<<<<<<<<
