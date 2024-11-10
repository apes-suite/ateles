require('common_1d')

scheme.temporal ={
                      name = 'explicitRungeKuttaTaylor',  
                      steps = 4,
                      control = {
                                 name = 'cfl',   
                                 cfl  = 0.3*27, --TODO 0.3*9, 0.3*3,
                                },
                 }
