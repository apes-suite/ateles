name = 'mesh_modg_3d'
 
-- define the input
input = {
         mesh = './mesh/',

  }


-- define the output
output = {  -- some general information
            folder = './harvest/',     -- Output location 
           {
            format = 'VTU',   -- Output format 
            binary = true,
            requestedData = {variable = {{'treeID'},{'process'}}},
            vrtx = {           -- when this table is defined set use_vrtx = .true.
                   } 
           }    

         }

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
              blocksize = 32,        -- for fpt, default -1
              fftMultiThread = false -- for fpt, logical, default false
             }
