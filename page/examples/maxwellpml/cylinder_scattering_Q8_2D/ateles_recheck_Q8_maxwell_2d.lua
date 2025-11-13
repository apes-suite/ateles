-- Configuration file for Ateles --

require('seeder')
require('hyperfun')

logging = {
  level = 10
}
timing_file = 'timing.res'         -- the filename of the timing results

-- global simulation options
sim_control = {
  name = 'pml_maxwell',
  time_control = {
    min = 0,
    max = 0.1,
    --interval = {iter = 1000},
  }
}
             --
-- The definition of the current density
x_center = -elemsize/2.0
y_center = -elemsize/2.0
-- ... pulse width in the x and y direction
width = 0.00125
-- ... peak-value of the current density
curConst = 1.0
-- ... the temporal mode
w = 10
-- ... width of the temporal ramping
widthTime = 0.05
-- ... the starting point of the temporal ramping
t_center = -0.5 -- TODO -2.0
function current(x,y,z,t)
  r =  (x-x_center)^2.0 + (y-y_center)^2.0
  tRamp = ( tanh( (t + t_center ) / widthTime ) + 1 ) / 2.0
  cur = curConst * math.exp( -r/width ) * math.sin( 2.0*math.pi*w*t ) * tRamp
  return { 0.0, cur }
end

-- ... material parameters inside te circle
cube_permea_circ = 2.0
cube_permit_circ = 2.0
cube_conduc_circ = 10.0 --TODO order 128: 160000.0 --TODO order 32: 10000.0
-- ... material parameters of the background material
cube_permea_back = 1.0
cube_permit_back = 1.0
-- ... the center of the circle
x_circle_center = elemsize
y_circle_center = elemsize
function permeability(x,y,z,t)
  -- Calculate the distance from the circle center
  dist = math.sqrt( (x-x_circle_center)^2 + (y-y_circle_center)^2 )
  -- Check if the point is inside the circle
  if dist <= elemsize  then
    return cube_permea_circ
  else
    return cube_permea_back
  end
end
function permittivity(x,y,z,t)
  -- Calculate the distance from the circle center
  dist = math.sqrt( (x-x_circle_center)^2 + (y-y_circle_center)^2 )
  -- Check if the point is inside the circle
  if dist <= elemsize  then
    return cube_permit_circ
  else
    return cube_permit_back
  end
end
function conductivity(x,y,z,t)
  -- Calculate the distance from the circle center
  dist = math.sqrt( (x-x_circle_center)^2 + (y-y_circle_center)^2 )
  -- Check if the point is inside the circle
  if dist <= elemsize  then
    return cube_conduc_circ
  else
    return 0.0
  end
end
pml_damp_exp = 2
pml_damp_factor = 100
variable = {
  -- PML for the left domain boundary
  {
    name = "pml_var",
    ncomponents = 4,
    vartype = "st_fun",
    st_fun = {
      -- PML for the left domain boundary
      {
        shape = {
          kind = 'canoND',
          object= {
            origin = { -cubeLength / 2 + elemsize - eps, -cubeLength / 2, 0.0 },
            vec = {
              { -elemsize + 2 * eps, 0, 0 },
              { 0, cubeLength, 0 },
              { 0, 0, eps }
            },
            segments = { 50, 50, 50 }
          }
        },
        predefined = 'combined',
        temporal = 1.0,
        spatial = {
          predefined = 'pml',
          plane_origin = { -cubeLength / 2 + elemsize, 0, 0 },
          plane_normal = { -1.0, 0, 0 },
          damp_factor = pml_damp_factor,
          damp_exponent = pml_damp_exp,
        }
      },
      -- PML for the right domain boundary
      {
        shape = {
          kind = 'canoND',
          object= {
            origin = { cubeLength / 2 - elemsize + eps, -cubeLength / 2, 0.0 },
            vec = {
              { elemsize - 2 * eps, 0, 0 },
              { 0, cubeLength, 0 },
              { 0, 0, eps }
            },
            segments = { 50, 50, 50 }
          }
        },
        predefined = 'combined',
        temporal = 1.0,
        spatial = {
          predefined = 'pml',
          plane_origin = { cubeLength / 2 - elemsize, 0, 0 },
          plane_normal = { 1.0, 0, 0 },
          damp_factor = pml_damp_factor,
          damp_exponent = pml_damp_exp,
        }
      },
      -- PML for the north domain boundary
      {
        shape = {
          kind = 'canoND',
          object= {
            origin = { -cubeLength / 2, cubeLength / 2 - elemsize + eps, 0.0, },
            vec = {
              { cubeLength, 0, 0 },
              { 0, elemsize - 2 * eps, 0 },
              { 0, 0, eps }
            },
            segments = { 50, 50, 50 },
          }
        },
        predefined = 'combined',
        temporal = 1.0,
        spatial = {
          predefined = 'pml',
          plane_origin = { 0, cubeLength / 2 - elemsize, 0 },
          plane_normal = { 0, 1.0, 0 },
          damp_factor = pml_damp_factor,
          damp_exponent = pml_damp_exp,
        }
      },
      -- PML for the south domain boundary
      {
        shape = {
          kind = 'canoND',
          object= {
            origin = { -cubeLength / 2, -cubeLength / 2 + elemsize - eps, 0.0 },
            vec = {
              { cubeLength, 0, 0 },
              { 0.0, -elemsize + 2 * eps, 0 },
              { 0, 0, eps }
            },
            segments = {50,50,50}
          }
        },
        predefined = 'combined',
        temporal = 1.0,
        spatial = {
          predefined = 'pml',
          plane_origin = { 0, -cubeLength / 2 + elemsize, 0 },
          plane_normal = { 0, -1.0, 0 },
          damp_factor = pml_damp_factor,
          damp_exponent = pml_damp_exp,
        }
      }
    }
  },
  -- The source of the waves
  {
    name = "current_density_var",
    ncomponents = 2,
    vartype = "st_fun",
    st_fun = current
  },
  {
    name = "permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      {
        const = 1.0
      },
      {
        fun = permeability,
        shape = {
          kind = "canoND",
          object = {
            origin = { eps, eps, eps },
            vec = {
              { 2*elemsize-2*eps, 0.0, 0.0 },
              { 0.0, 2*elemsize-2*eps, 0.0 },
              { 0.0, 0.0, 2*elemsize-2*eps }
            },
            segments = { 50, 50, 50 }
          }
        }
      }
    }
  },
  {
    name = "permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      {
        const = 1.0
      },
      {
        fun = permittivity,
        shape = {
          kind = "canoND",
          object = {
            origin = { eps, eps, eps },
            vec = {
              { 2*elemsize-2*eps, 0.0, 0.0 },
              { 0.0, 2*elemsize-2*eps, 0.0 },
              { 0.0, 0.0, 2*elemsize-2*eps }
            },
            segments = { 50, 50, 50 }
          }
        }
      }
    }
  },
  {
    name = "conductivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      {
        const = 0.0
      },
      {
        fun = conductivity,
        shape = {
          kind = "canoND",
          object = {
            origin = { eps, eps, eps },
            vec = {
              { 2*elemsize-2*eps, 0.0, 0.0 },
              { 0.0, 2*elemsize-2*eps, 0.0 },
              { 0.0, 0.0, 2*elemsize-2*eps }
            },
            segments = { 50, 50, 50 }
          }
        }
      }
    }
  }
}

-- Mesh definitions --
mesh = 'mesh/'

---- Restart settings
--restart = {
--            --read = './restart/inhomogeneous_material_lastHeader.lua',
--            write = './restart/',
--            time_control = {
--                             min = 0,
--                             max = sim_control.time_control.max,
--                             interval = sim_control.time_control.max/800,
--                   },
--          }

-- the filename of the timing results i.e. output for performance measurements
timing = { filename = 'timing.res' }

-- Equation definitions --
equation = {
  name   = 'maxwell_2d',
  material = {
    permeability = 'permeability',
    permittivity = 'permittivity',
    conductivity = 'conductivity'
  }
}

-- Check for Nans and unphysical values
check =  {
  interval = 1
}

-- Scheme definitions --
degree = 8
scheme = {
  spatial =  {
    name = 'modg_2d',
    m = degree,
    modg_space = 'Q'
  },
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name = 'cfl',
      cfl  = 0.6*(2*degree+1)^2/((degree+1)^2)
    }
  }
}

function ini_displ_x(x,y,z)
  return math.sin(x)
end

-- ...the initial condition table.
initial_condition = {
  displacement_fieldX = ini_displ_x, -- electric field , x component
  displacement_fieldY = 0.0, -- electric field , z component
  magnetic_field     = 0.0, -- magnetic induction , z component
  pml_PX = 0.0,
  pml_PY = 0.0,
  pml_QX = 0.0,
  pml_QY = 0.0,
}

-- ...the general projection table
projection = {
  kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
  -- for fpt the  nodes are automatically 'chebyshev'
  -- for lep the  nodes are automatically 'gauss-legendre'
  -- lobattopoints = false  -- if lobatto points should be used, default = false,
                            -- only working for chebyshev points --> fpt
  factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
  -- blocksize = 32,        -- for fpt, default -1
  -- fftmultithread = false -- for fpt, logical, default false
}


-- Source terms
source = {
  -- PML for the left domain boundary
  pml = "pml_var",
  -- The source of the waves
  current_density = "current_density_var"
}

tracking = {
  label = 'probe_electricField_Q8_pml_circMaterial',
  folder = './',
  variable = {'displacement_field'},
  shape = {
    kind = 'canoND',
    object= {
      origin = { 0.0, 0.0, 0.0 }
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/20.0
  },
  output = { format = 'ascii', ndofs = 1 }
}
