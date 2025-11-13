-- Configuration file for Ateles --

-- The analytic solution of this testcase.
-- ... length of the waveguide
lengthX = 1.0
lengthY = 1.0
lengthZ = 2.0
-- ... shift of the coordinate system in x and y direction
xOffset = 0.5
yOffset = 0.5
-- ... the magnetic permeability
permea = 2
-- ... the displacement permitivity
permit = 1
-- .. speed of light
-- ... the physical mode indices
m = 1.0 -- x mode (this can be any integer)
n = (lengthY/lengthX)*m -- y mode (fixed by integer for m)
l = 2.0 -- z mode (only even integers are allowed for the analytic solution)
-- ... parameters for the analytic solution
gammaSq = 2.0 * (m*math.pi/lengthX)^2.0
w = math.sqrt( (gammaSq + (l*math.pi/lengthZ)^2.0) / (permea*permit) )
-- ... the temporal period of the waveguide
T = 2.0*math.pi/w

-- Definition of the analytical solution for all components of the Maxwell equation
-- ... displacement field - x component
B = 1.0
function displacementX(x,y,z,t)
  return B * (n*math.pi/lengthY) * math.cos(m*math.pi*(x+xOffset)/lengthX) * math.sin(n*math.pi*(y+yOffset)/lengthY) * ((-1.0)*w/gammaSq) * math.sin(l*math.pi*z/lengthZ) * math.sin(w*t)
end
-- ... displacement field - y component
function displacementY(x,y,z,t)
  return (-1.0) * B * (n*math.pi/lengthY) * math.sin(m*math.pi*(x+xOffset)/lengthX) * math.cos(n*math.pi*(y+yOffset)/lengthY) * ((-1.0)*w/gammaSq) * math.sin(l*math.pi*z/lengthZ) * math.sin(w*t)
end
-- ... displacement field - z component
function displacementZ(x,y,z,t)
  return 0.0
end
-- ... magnetic field - x component
function magneticX(x,y,z,t)
  return (1.0/gammaSq) * (l*math.pi/lengthZ) * math.cos(l*math.pi*z/lengthZ) * (-1.0) * B * (n*math.pi/lengthY) * math.sin(m*math.pi*(x+xOffset)/lengthX) * math.cos(n*math.pi*(y+yOffset)/lengthY) * math.cos(w*t)
end
-- ... magnetic field - y component
function magneticY(x,y,z,t)
  return (1.0/gammaSq) * (l*math.pi/lengthZ) * math.cos(l*math.pi*z/lengthZ) * (-1.0) * B * (n*math.pi/lengthY) * math.cos(m*math.pi*(x+xOffset)/lengthX) * math.sin(n*math.pi*(y+yOffset)/lengthY) * math.cos(w*t)
end
-- ... magnetic field - z component
function magneticZ(x,y,z,t)
  return B * math.cos(m*math.pi*(x+xOffset)/lengthX) * math.cos(n*math.pi*(y+yOffset)/lengthY) * math.sin(l*math.pi*z/lengthZ) * math.cos(w*t)
end

-- global simulation options
simulation_name = 'rectang_waveguide_maxwell_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    min = 0.0,
    max = 1.0*T/1000 -- final simulation time
  }
}

--commpattern = 'gathered_type'

-- Mesh definitions --
mesh = './mesh/'


-- timing settings (i.e. output for performance measurements, this table is otional)
timing = {
  folder = './',                  -- the folder for the timing results
  filename = 'timing.res'         -- the filename of the timing results
}

-- Tracking
tracking = {
  label = 'probe_electricField_Q8',
  folder = './',
  variable = {'displacement_field'},
  shape = {
    kind = 'canoND',
    object= {
      origin ={ 0., 0., 0. }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/2.0
  },
  output = { format = 'ascii', ndofs = 1 }
}

variable = {
  -- This is the global material for Maxwell. It consists of three different
  -- components, permeability, permittivity, and conductivity, each a scalar, so
  -- that we need three scalars
  -- for this equation system.
  -- As this is the global fallback material, we define each material to be a
  -- neutral term, which in this case is 0.
  {
    name = "global_permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { permea }
    }
  },
  {
    name = "global_permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { permit }
    }
  },
  {
    name = "global_conductivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = { 0.0 }
    }
  }
}

-- Equation definitions --
equation = {
  name = 'maxwell',  -- we solve maxwell's equations
  material = {
    permeability = 'global_permeability',
    permittivity = 'global_permittivity',
    conductivity = 'global_conductivity'
  }
}

-- Scheme definitions --
degree = 7
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    m =  degree,
    modg_space = 'Q',
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',     -- the name of the timestep control mechanism
      cfl  = 0.064*(3*degree+1)^2/((degree+1)^2),     -- Courant-Friedrichs-Lewy number
    }
  }
}

-- ...initial condition function for displacement field
function ic_displacementX(x,y,z)
  return displacementX(x,y,z,0)
end
function ic_displacementY(x,y,z)
  return displacementY(x,y,z,0)
end
function ic_displacementZ(x,y,z)
  return displacementZ(x,y,z,0)
end
-- ...initial condition function for magnetic field
function ic_magneticX(x,y,z)
  return magneticX(x,y,z,0)
end
function ic_magneticY(x,y,z)
  return magneticY(x,y,z,0)
end
function ic_magneticZ(x,y,z)
  return magneticZ(x,y,z,0)
end

projection = {
  kind = 'fpt', -- 'fpt' or 'l2p', default 'l2p'
                -- for fpt the  nodes are automatically 'chebyshev'
                -- for lep the  nodes are automatically 'gauss-legendre'
  -- lobattoPoints = false  -- if lobatto points should be used, default = false
  factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
  blocksize = 32,        -- for fpt, default -1
  -- fftMultiThread = false -- for fpt, logical, default false
}

-- ...the initial condition table
initial_condition = {
  displacement_fieldX = ic_displacementX,  -- displacement field , x component
  displacement_fieldY = ic_displacementY,  -- displacement field , y component
  displacement_fieldZ = ic_displacementZ,  -- displacement field , z component
  magnetic_fieldX = ic_magneticX,  -- magnetic induction , x component
  magnetic_fieldY = ic_magneticY,  -- magnetic induction , y component
  magnetic_fieldZ = ic_magneticZ,  -- magnetic induction , z component
}


-- Boundary definitions
boundary_condition = {
  {
    label = 'pecEast',   -- boundary for the inner cylinder
    kind = 'pec',              -- perfectly electric conductor
  },
  {
    label = 'pecWest',   -- boundary for the outer cylinder
    kind = 'pec',              -- perfectly electric conductor
  },
  {
    label = 'pecSouth',       -- boundary for the inner cylinder
    kind = 'pec',         -- perfectly electric conductor
  },
  {
    label = 'pecNorth',       -- boundary for the inner cylinder
    kind = 'pec',         -- perfectly electric conductor
  }
}

