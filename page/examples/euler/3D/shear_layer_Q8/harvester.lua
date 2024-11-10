name = 'shear_layer'
 
-- define the input
input = {
         read = './restart/shear_layer_modg_lastHeader.lua',
         --read = os.getenv('RESTARTHEADER'),
         --mesh = './mesh/'

         ---- define the subsampling parameters
         subsampling = {
                         levels = 3,  -- the number of subsampling steps, i.e. 1 means we increase the number of elements by 8
                         --projection = 'QLegendre', -- the type of projection we use for the subsampling of the data.
                         --dofReductionFactor = 1.5,  -- the reduction factor for the dofs of each subsampling step (per spatial dir)
                         projection = 'QLegendrePoint', -- the type of projection we use for the subsampling of the data.
                       }
         
  }

-- define the output
output = {  -- some general information
            folder = './harvest/',     -- Output location 

           {

    
            output_format = 'VTK',   -- Output format 

            binary = true,
            vrtx = {           -- when this table is defined set use_vrtx = .true.
--                     updated_vertices  = 'prefix_of_new_vertices',            -- when this is given set  new_vrtx = .true.
--                     old_vertices  = 'prefix_of_old_vertices',                -- when this is given set  old_vrtx = .true.
--                     dump_vertices = 'prefix_of_where_to_dump_vertices'       -- when this is given set dump_vrtx = .true.
                   } 
           }    

         }

---- output settings
--output = {   folder = './harvest/maxwell/', 
--
--       { 
--
--
--          dumpAll = true,  -- only read if requested data table is not present else default false
--
--          output_format = 'VTK',   -- Output format 
--          binary = true,
--          vrtx = {      
--                },
--           
--
--          label = 'First_subset',
--          shape ={ { kind = 'all' } } ,
--        
--          subsampling = {
--              levels = 4,  -- the number of subsampling steps, i.e. 1 means we increase the number of elements by 8
----              --projection = 'QLegendre', -- the type of projection we use for the subsampling of the data.
----              --dofReductionFactor = 1.5,  -- the reduction factor for the dofs of each subsampling step (per spatial dir)
--              projection = 'QLegendrePoint', -- the type of projection we use for the subsampling of the data.
--                      }
--       },
--
--}
