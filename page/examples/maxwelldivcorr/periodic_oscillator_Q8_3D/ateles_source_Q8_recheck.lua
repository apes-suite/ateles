-- Configuration file for Ateles (Periodic Oscillator) --

--------------------------------------------------------------------------------
--...Configuration of simulation time
simulation_name = 'maxwell_source'  
sim_control = {
  time_control = { 
    max =  math.sqrt(2), -- Number of iterations / Simulated time
    min = 0.0
  }
}

--... Mesh definitions --
mesh = 'mesh/'

-- Source term definition, i.e. in Maxwell equations we are talking about space charges and 
-- current densities. In general they can depend on spatial coordinates and time.
function currentDensitySpaceTime(x, y, z, t)
  d = math.sqrt(x^2.0+y^2.0+z^2.0)
  if d <= r then
    jx=Q*r*freq*math.sin(freq*t)
    jy=Q*r*freq*math.sin(freq*t)
    jz=Q*r*freq*math.sin(freq*t)
    return {jx,jy,jz}
  else
    return {0.0, 0.0, 0.0}
  end
end

--...Equation definitions --
permea = 1.0
permit = 1.0
gam = 1.0
chi = 1.0

variable = {
  {
    name = "var_current_density",
    ncomponents = 3,
    vartype = "st_fun",
    st_fun = {
      fun = currentDensitySpaceTime,
      shape = {
        kind = 'canoND',
        object= {
          origin = { -0.1, -0.1, -0.1 },
          vec = { { 0.3, 0.0, 0.0 }, { 0.0, 0.3, 0.0 }, { 0.0, 0.0, 0.3 } },
          -- length = 0.6,
          segments = { 100, 100, 100 }
        }
      }
    }
  },
  -- This is the global material for Maxwell. It consists of two different 
  -- components, permittivity, and conductivity.
  -- As this is the global fallback material, we define each material to be a 
  -- neutral term, which in this case is 0.
  {
    name = "permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { permea }
    }
  },
  {
    name = "permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { permit }
    }
  },
  {
    name = "gam",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { gam }
    }
  },
  {
    name = "chi",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { chi }
    }
  }
}

source = {
  currentDensity = 'var_current_density'
}

equation = {
  -- we solve maxwell’s equations with divergence correction
  name = 'maxwellDivCorrection',
  material = { 
    permeability = 'permeability',
    permittivity = 'permittivity',
    gam = 'gam',
    chi = 'chi'
  }
}
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

---- ...the initial condition table. 
---- ...initial condition function for displacement field (z component)
--function ic_displacementZ(x,y,z)
--  return displacementZ(x,y,z,0.0)
--end
 
-- ...Initial condition
initial_condition = {
  displacement_fieldX = 0.0,           -- displacement field , x component
  displacement_fieldY = 0.0,           -- displacement field , z component
  displacement_fieldZ = 0.0,  -- displacement field , z component
  magnetic_fieldX = 0.0,           -- magnetic induction , x component
  magnetic_fieldY = 0.0,           -- magnetic induction , y component
  magnetic_fieldZ = 0.0,           -- magnetic induction , z component
  magnetic_correction = 0.0, -- magnetic div correction
  electric_correction = 0.0, -- displacement div correction
}

--... Definition of the projection method
projection = {
  kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
                 -- for fpt the  nodes are automatically 'chebyshev'
                 -- for lep the  nodes are automatically 'gauss-legendre'
  -- lobattoPoints = false  -- if lobatto points should be used, default = false
  factor = 1.0          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
  -- blocksize = 32,        -- for fpt, default -1
  -- fftMultiThread = false -- for fpt, logical, default false
}

--... Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',           -- we use the modal discontinuous Galerkin scheme 
    m = 7,                   -- the maximal polynomial degree for each spatial direction
    modg_space = 'Q'
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta', --'explicitEuler',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',   -- the name of the timestep control mechanism
      cfl  = 0.95     -- CourantÐFriedrichsÐLewy number
    }
  }
}

--...Configuration for the restart file
estart = { 
  write = 'restart/',
  -- temporal definition of restart write
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/5.0
  }
}

 -- Tracking used for validation.    
tracking = {
  label = 'divcor_source_probe_electricField_Q8',
  folder = './',
  variable = {'displacement_field'},
  shape = {kind = 'canoND', object= { origin ={1., 1., 1.} } },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/5.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

-- Source term
-- ... charge of the source 
 Q = 1.0 
-- ... radius of sphere source
 r = 0.4
-- ... parameters for the analytic solution
freq = ( 2.0*math.pi/math.sqrt(permit*permea) ) *10
-- ... the temporal period of the waveguide
T = 2.0*math.pi/freq

 
-- Boundary definitions
boundary_condition = {
  {
    label = 'wall_1',   
    kind = 'pec'
  },
  {
    label = 'wall_2',   
    kind = 'pec'
  },
  {
    label = 'wall_3',   
    kind = 'pec'
  },
  {
    label = 'wall_4',   
    kind = 'pec'
  },
  {
    label = 'wall_5',   
    kind = 'pec'
  },
  {
    label = 'wall_6',   
    kind = 'pec'
  }
}

-- DEBUG OPTIONS
--debug = { 
--         debugMode = true,        -- default= false
--         verbose = 99,             -- default= 0
--         debugFiles = true,       -- default= false
--         -- What to dump into debugFiles? --
--           dumpTreeIDs = true,      -- default= false
--           dumpPropBits = true,     -- default= false
--           dumpAuxLists = true,     -- default= false
--           dumpDependencies = true, -- default= false
--           dumpState = true,        -- default= false
--           dumpHaloState = true,    -- default= false
--         --  end debugFiles  --
--         debugDependencies = true, -- default= false
--         checkDependencies = true, -- default= false
--         checkNans = true,         -- default= false
--         checkSteps = true,        -- default= false
--         debugMesh = 'dbg/mesh_',  -- default= ''
--         debugSource = true,       -- default= false
--         debugRestart = true,      -- default= false
--         traceMemory = true,       -- default= false
--        }
