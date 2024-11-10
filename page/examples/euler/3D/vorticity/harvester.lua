name = 'gPulseDens_euler_modg'

logging = {level=10} 
-- define the input
input = {
         --read = './tracking/track_pres_gPulseDens_euler_modg_header_0.000E+00.lua',
         read = './restart/simulation_header_0.000E+00.lua',
         --read = './tracking/simulation_euler_conservative_vort_simulation_header_0.000E+00.lua',
                      
         --read = './restart/custom.lua',
         --mesh = './mesh/',

         -- define the subsampling parameters
         subsampling = {
                         levels = 1,  -- the number of subsampling steps, i.e. 1 means we increase the number of elements by 8
                         --projection = 'QLegendre', -- the type of projection we use for the subsampling of the data.
                         --dofReductionFactor = 1.5,  -- the reduction factor for the dofs of each subsampling step (per spatial dir)
                         projection = 'QLegendrePoint', -- the type of projection we use for the subsampling of the data.
                       }
        }


-- define the output
output = {  -- some general information
            folder = './harvest/',     -- Output location 

           {

            --requestedData = { variable = {{'velocity'},{'density'},{'kinetic_energy'} }},    
            --requestedData = {variable = {'vorticity'}},    
            --requestedData = {variable = {{'density'},{'momentum'},{'velocity'},{'vorticity'}}},    
            --requestedData = {variable = {{'density'},{'momentum'},{'velocity'}}},    
            --requestedData = { variable = {{'energy'},{'density'},{'momentum'}}},    
            format = 'VTU',   -- Output format 

            binary = true,
            vrtx = {           -- when this table is defined set use_vrtx = .true.
--                     updated_vertices  = 'prefix_of_new_vertices',            -- when this is given set  new_vrtx = .true.
--                     old_vertices  = 'prefix_of_old_vertices',                -- when this is given set  old_vrtx = .true.
--                     dump_vertices = 'prefix_of_where_to_dump_vertices'       -- when this is given set dump_vrtx = .true.
                   } 
           }    

          }

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
           -- lobattoPoints = false  -- if lobatto points should be used, default = false
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
              blocksize = 32,        -- for fpt, default -1
              fftMultiThread = false -- for fpt, logical, default false
             }         

