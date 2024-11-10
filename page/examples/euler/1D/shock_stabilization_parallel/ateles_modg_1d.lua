require('common_1d')

scheme.temporal ={
                      name = 'explicitRungeKutta',  
                      steps = 4,
                      control = {
                                 name = 'cfl',   
                                 cfl  = 0.3,     
                                },
                 }
