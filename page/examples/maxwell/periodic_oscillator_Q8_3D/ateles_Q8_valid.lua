-- Configuration file for Ateles (Periodic Oscillator) --

-- This is a configuration file for the Finite Volume / Discontinuous Galerkin Solver ATELES. 
-- It provides a testcase for the simulation of Maxwell equations in a homogenous media. The simulation domain
-- is a periodic cube with edge length 2.0. Therefore this is a very good way to verify your algorithmic implementations, 
-- since this testcase does not involve any boundary conditions. 
-- The testcase simulates the temporal development of standing waves for electric and magnetic fields. Since we 
-- are considering a very simple domain, an analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- To verify diffusion and dispersion relations this testcases allows you to modify the spatial harmonics by
-- varying the integer mode number in x and y direction by modifying the lua variables m and n. Please notice,
-- this testcase is correct only for homogenous media with epsi = mu = 1 (see equations table).
-- This testcase can be run in serial (only one execution unit) or in parallel (with multiple mpi ranks).
-- To calculate a grid convergence behavior please modify the level variable. An increment of one will half the radius of your elements.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Parameters to vary --
degree = 7
poly_space = 'Q'

-- ...the uniform refinement level for the periodic cube
level = 2

--...Configuration of simulation time
--simtime = {
--            max = 1.0,            -- Number of iterations / Simulated time
--            min = 0
--          }

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
           -- lobattoPoints = false  -- if lobatto points should be used, default = false
           -- factor = 2.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
           -- blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }
--------------------------------------------------------------------------------
-- General settings
-- The io_buffer_size specifies the amount of memory reserved for blocks of
-- temporary memory allocated for copies of data in MB.
-- (You do not have to specify it at all).
-- io_buffer_size = 1 -- default is 8 MB

-- Communication pattern to use in parallel executions
commpattern = 'isend_irecv' -- default is simple isend-irecv exchange
--          = 'isend_irecv_overlap' allow overlapping of isends and irecvs
--          = 'typed_isend_irecv' use a MPI datatype to describe the data
--          = 'gathered_type' minimized number of blocks

-- Print memory upon finalize?
-- This is obtained from the Linux pseudo file /proc/self/status.
-- It is only printed by the first process in a prallel run, and only provides
-- data on this first process.
printRuntimeInfo = true -- default is true
--------------------------------------------------------------------------------

-- END Parameters to vary --
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
require("valid_tracking")
require("posci_common")
