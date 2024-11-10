-- Configuration file for Ateles (Periodic Oscillator) --

-- This is a configuration file for the Finite Volume / Discontinuous Galerkin
-- Solver ATELES. 
-- It provides a testcase for the simulation of Acoustic equations in a
-- homogenous media. The simulation domain is a periodic cube with edge length
-- 2.0. Therefore this is a very good way to verify your algorithmic
-- implementations, since this testcase does not involve any boundary
-- conditions.
-- The testcase simulates the temporal development of standing waves for
-- acoustic equation. Since we are considering a very simple domain, an
-- analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- To verify diffusion and dispersion relations this testcases allows you to
-- modify the spatial harmonics by varying the integer mode number in x and y
-- direction by modifying the lua variables m and n.
-- This testcase can be run in serial (only one execution unit) or in parallel
-- (with multiple mpi ranks).
-- To calculate a grid convergence behavior please modify the level variable.
-- An increment of one will half the radius of your elements.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Parameters to vary --
degree = 7
poly_space = 'Q'


-- Check for Nans and unphysical values
check =  {
           interval = 100,
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

--...Configuration of simulation time
simulation_name = 'plane_wave_XQ8' -- the name of the simualtion
timing_file = 'timing.res'         -- the filename of the timing results
sim_control = { 
                time_control = { max = 0.01*1,  -- final Simulated time
                                 min = 0,
                                 interval = {iter = 100} 
                                }
              }

-- End  Parameters to vary --
--------------------------------------------------------------------------------
-- Definition of the test-case.

-- Mesh definitions --
-- ...the uniform refinement level for the periodic cube
level = 2
-- ...the length of the cube
cubeLength = 2.0
mesh = { predefined = 'cube',
         origin = { 
                    (-1.0)*cubeLength/2.0,
                    (-1.0)*cubeLength/2.0,
                    (-1.0)*cubeLength/2.0
                  },
         length = cubeLength,
         refinementLevel = level
       }

-- Tracking              
eps=cubeLength/(2^(level+1))
tracking = {
             label = 'track_line_density_m7',
             folder = './',
             variable = {'density'},
             shape = {kind = 'canoND', object= { origin = {0.0+eps,0.0,0.0}, 
                                               }
                     },
             time_control= { min = 0,
                           --  max = sim_control.time_control.max,
                             interval = sim_control.time_control.max/8.0,
                           },
             output = { format = 'ascii', ndofs = 1 },
           }


-- Equation definitions --
equation = {
             name   = 'acoustic',
             background = {
                 density = 1.225, 
                 velocityX = 0.0,
                 velocityY = 0.0,
                 velocityZ = 0.0,
                 pressure = 100000.0
                 }
           }

-- Scheme definitions --
scheme = {
    -- the spatial discretization scheme
    spatial =  {
               name = 'modg',            -- we use the modal discontinuous Galerkin scheme 
               m =  degree,                   -- the maximal polynomial degree for each spatial direction
               modg_space = poly_space
               }, 
    -- the temporal discretization scheme
    temporal = {
                 name = 'explicitRungeKutta', 
                 steps = 4,
              -- how to control the timestep
                 control = {
                          name = 'cfl',   -- the name of the timestep control mechanism
                          cfl  = 0.95,     -- CourantÐFriedrichsÐLewy number
                         },
               },
             }


c = math.sqrt(equation.background.pressure / equation.background.density)
ampl_density= (equation.background.density/c) 
ampl_X= ampl_density /equation.background.density * ( equation.background.velocityX - c)
ampl_Y = 1.0--0.001225
ampl_Z = 1.0 --0.001225

function velocityX(x,y,z,t)
  return ( (ampl_X * math.cos(math.pi*x-c*t)) )
end 
function velocityY(x,y,z,t)
  return ( (ampl_Y * math.cos(math.pi*y-c*t)) )
end 
function velocityZ(x,y,z,t)
  return ( (ampl_Z * math.cos(math.pi*z-c*t)) )
end 
function density(x,y,z,t)
  return ( ampl_density * math.cos(math.pi*x-c*t) )
end 

function init_velocityX(x,y,z)
  return (velocityX(x,y,z,0.0)  )
end 
function init_velocityY(x,y,z)
  return (velocityY(x,y,z,0.0)  )
end 
function init_velocityZ(x,y,z)
  return (velocityZ(x,y,z,0.0)  )
end 
function init_density(x,y,z)
  return ( density(x,y,z,0.0) )
end 

-- Initial Condition definitions --
initial_condition = { 
                      density = init_density, 
                      velocityX = init_velocityX,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                    }
