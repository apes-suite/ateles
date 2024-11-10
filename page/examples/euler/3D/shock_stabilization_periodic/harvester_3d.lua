name = 'euler_modg_3d'
 
-- define the input
input = {
         read = './restart/euler_3d_lastHeader.lua',
         --mesh = './mesh/',

         ---- define the subsampling parameters
         subsampling = {
                         levels = 2,  
                         projection = 'QLegendrePoint', 
                       }
  }


-- define the output
output = {  -- some general information
            folder = './harvest/',     -- Output location 

           {

    
            format = 'VTU',   -- Output format 

            binary = true,
            --requestedData = {variable = {{'density'},{'momentum'},{'energy'},{'pressure'},{'mach_number'},{'mach_vector'}}},
            vrtx = {           -- when this table is defined set use_vrtx = .true.
--                     updated_vertices  = 'prefix_of_new_vertices',            -- when this is given set  new_vrtx = .true.
--                     old_vertices  = 'prefix_of_old_vertices',                -- when this is given set  old_vrtx = .true.
--                     dump_vertices = 'prefix_of_where_to_dump_vertices'       -- when this is given set dump_vrtx = .true.
                   } 
           }    

         }

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
              blocksize = 32,        -- for fpt, default -1
              fftMultiThread = false -- for fpt, logical, default false
             }
