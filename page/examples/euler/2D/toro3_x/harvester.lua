name = 'toro'
 
input = {
         --read = './restart/simulation_lastHeader.lua',
         read = '$!filename!$',
         
         subsampling = {
                         levels = 8, 
                         projection = 'QLegendrePoint',
                       },

  }


output = {  
            --folder = './harvest/', 
            folder = '$!folder!$', 

           {

    
            format = 'VTU',

            binary = true,
         --TODO   requestedData = {
         --TODO                    variable = {
         --TODO                                 { name = 'density', },
         --TODO                                 { name = 'momentum',},
         --TODO                                 { name = 'energy',  },
         --TODO                                 { name = 'velocity' },
         --TODO                                 { name = 'vorticity' },
         --TODO                                 --{ name = 'pressure', ncomponents = 1 },
         --TODO                                 --{ name = 'temperature', ncomponents = 1 },
         --TODO                               },
         --TODO                   },
            vrtx = {
                   },
           }    

         }

