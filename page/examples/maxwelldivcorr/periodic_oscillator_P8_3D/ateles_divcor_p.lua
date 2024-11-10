-- Configuration file for Ateles (Periodic Oscillator) --

-- This is a configuration file for the Finite Volume / Discontinuous Galerkin
-- Solver ATELES.
-- It provides a testcase for the simulation of Maxwell equations in a
-- homogenous media. The simulation domain is a periodic cube with edge length
-- 2.0. Therefore this is a very good way to verify your algorithmic
-- implementations, since this testcase does not involve any boundary
-- conditions. The testcase simulates the temporal development of standing waves
-- for displacement and magnetic fields. Since we are considering a very simple
-- domain, an analytic solution is well known and given as Lua functions in this
-- script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- To verify diffusion and dispersion relations this testcases allows you to
-- modify the spatial harmonics by varying the integer mode number in x and y
-- direction by modifying the lua variables m and n. Please notice,
-- this testcase is correct only for homogenous media with epsi = mu = 1 
-- (see equations table).
-- This testcase can be run in serial (only one execution unit) or in parallel
-- (with multiple mpi ranks).
-- To calculate a grid convergence behavior please modify the level variable.
-- An increment of one will half the radius of your elements.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

--...Configuration of simulation time
simulation_name = 'maxwell'
sim_control = {
  time_control = {
    max =  math.sqrt(2), -- Number of iterations / Simulated time
    min = 0.0
  }
}

variable = {
  -- This is the global material for Maxwell. It consists of two different 
  -- components, permittivity, and conductivity.
  -- As this is the global fallback material, we define each material to be a 
  -- neutral term, which in this case is 0.
  {
    name = "global_permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 1.0 }
  },
  {
    name = "global_permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 1.0 }
  },
  {
    name = "global_gam",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 1.0 }
  },
  {
    name = "global_chi",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 1.0 }
  },
}

--... Mesh definitions --
cubeLength=2.0
mesh = {
  predefined = 'cube',
  origin = {
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0,
    (-1.0)*cubeLength/2.0
  },
  length = 2,
  refinementLevel = 2
}

--...Equation definitions --

equation = {
  name = 'maxwellDivCorrection',    -- we solve maxwell’s equations with divergence correction
  material = {
    permeability = "global_permeability",
    permittivity = "global_permittivity",
    gam = "global_gam",
    chi = "global_chi"
  }
}

--equation = {
--   		name   = 'maxwell_DivCorr',   -- we solve maxwell's equations
--    		mu  = 1.0,            -- the magnetic permeability (vacuum has 4.0*math.pi*(10.0^(-7.0)))
--    		epsi = 1.0,           -- the displacement permitivity (vacuum has 8.85418781762*(10.0^(-12.0)))
--	   }
--
--------------------------------------------------------------------------------
-- Definition of the Periodic oscillator test-case.
-- Some global parameters for the T_{nm} mode testcase

-- ...the integer number of the mode in x direction
amplX = 1.0
-- ...the integer number of the mode in y direction
amplY = 1.0

-- The analytic solution for this testcase is given by the following functions
-- (only correct for epsi = mu = 1):
-- ... definition of temporal angular frequency
w = math.sqrt(amplX^2 + amplY^2)
-- ... E_x = 0.0
function displacementX(x,y,z,t)
  return 0.0
end 
-- ... E_y = 0.0
function displacementY(x,y,z,t)
  return 0.0 --math.sin(amplX*math.pi*x)*math.sin(amplY*math.pi*z)*math.cos(w*t)
end 
-- ... E_z(x,y,z,t) = sin(amplX \pi x) sin(amplY \pi y) cos(w t)
function displacementZ(x,y,z,t)
  return math.sin(amplX*math.pi*x)*math.sin(amplY*math.pi*y)*math.cos(w*t)
end 
-- ... B_x(x,y,z,t) = -\frac{\pi n}{w} sin(m \pi x) cos(n \pi y) sin(w t)
function magneticX(x,y,z,t)
  return (-1.0)*(math.pi*amplY/w)*math.sin(amplX*math.pi*x)*math.cos(amplY*math.pi*y)*math.sin(w*t)
end 
-- ... B_y(x,y,z,t) = \frac{\pi m}{w} cos(m \pi x) sin(n \pi y) sin(w t)
function magneticY(x,y,z,t)
   return (math.pi*amplX/w)*math.cos(amplX*math.pi*x)*math.sin(amplY*math.pi*y)*math.sin(w*t)
end 
-- ... B_z = 0.0
function magneticZ(x,y,z,t)
  return 0.0 --(math.pi*amplX/w)*math.cos(amplX*math.pi*x)*math.sin(amplY*math.pi*z)*math.sin(w*t )
end

-- ...the initial condition table. 
-- ...initial condition function for displacement field (z component)
function ic_displacementZ(x,y,z)
  return displacementZ(x,y,z,0.0)
end

projection = {
  kind = 'l2p',  -- 'fpt' or 'l2p', default 'l2p'
                 -- for fpt the  nodes are automatically 'chebyshev'
                 -- for lep the  nodes are automatically 'gauss-legendre'
  -- lobattoPoints = false  -- if lobatto points should be used, default = false
  factor = 2.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
  -- blocksize = 32,        -- for fpt, default -1
  -- fftMultiThread = false -- for fpt, logical, default false
}

-- ...Initial condition
initial_condition = {
  displacement_fieldX = 0.0,           -- displacement field , x component
  displacement_fieldY = 0.0,           -- displacement field , z component
  displacement_fieldZ = ic_displacementZ,  -- displacement field , z component
  magnetic_fieldX = 0.0,           -- magnetic induction , x component
  magnetic_fieldY = 0.0,           -- magnetic induction , y component
  magnetic_fieldZ = 0.0,           -- magnetic induction , z component
  electric_correction = 0.0,           -- magnetic induction , z component
  magnetic_correction = 0.0,           -- magnetic induction , z component
}

--... Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',            -- we use the modal discontinuous Galerkin scheme 
    m =  7,                   -- the maximal polynomial degree for each spatial direction
    modg_space = 'P',
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta', --'explicitEuler',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',   -- the name of the timestep control mechanism
      cfl  = 0.95,     -- CourantÐFriedrichsÐLewy number
    }
  }
}

--...Configuration for the restart file
estart = {
--            -- file to restart from
--            read = './restart/maxwell/per_osci_maxwell_modg_lastHeader.lua',
  -- folder to write restart data to
  write = 'restart/DivCor/p/',
  -- temporal definition of restart write
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max
  }
}
          -- Tracking used for validation.

tracking = {
  label = 'divcor_probe_electricField_P8',
  folder = './',
  variable = { 'displacement_field' },
  shape = {
    kind = 'canoND', 
    object= { 
      origin ={ 0., 0., 0. }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/5.0,
  },
  output = { format = 'ascii', ndofs = 1 }
}
