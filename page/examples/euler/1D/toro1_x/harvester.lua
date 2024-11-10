name = 'toro'
 
input = {
         --read = './restart/simulation_lastHeader.lua',
         read = '$!filename!$',
         
         subsampling = {
                         levels = 12, 
                         projection = 'QLegendrePoint',
                       },

  }


output = {  
            --folder = './harvest/', 
            folder = '$!folder!$', 

           {

    
            format = 'VTU',

            binary = true,
            vrtx = {
                   },
           }   ,
           {

            label = '_track',
            format = 'ASCII',

            vrtx = {
                   },
           }    

         }

