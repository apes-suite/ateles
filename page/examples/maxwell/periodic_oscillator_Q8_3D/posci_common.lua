-- Definition of the P:wqeriodic oscillator test-case.
-- TO BE included in actual config scripts, setting parameters of interest.

-- Some global parameters for the T_{nm} mode testcase

-- ...the length of the cube
cubeLength = 2.0

-- ...the integer number of the mode in x direction
amplX = 1.0
-- ...the integer number of the mode in y direction
amplY = 1.0

logging = {level = 10}
-- global simulation options
simulation_name = 'posci_modg'..poly_space..(degree+1) -- the name of the simualtion
sim_control = {
                time_control = {
                  min = simtime.min,
                  max = simtime.max
                }
              }

variable = {
  {  
     name = 'permeability',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = 1.0

  },
  
  {  
     name = 'permittivity',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = 1.0

  },
  
  {  
     name = 'conductivity',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = 0.0

  },
}

-- Equation definitions --
equation = {
  name = 'maxwell',                   -- we solve maxwell's equations
    material = {
                 permeability = 'permeability',
                 permittivity = 'permittivity',
                 conductivity = 'conductivity'

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
               name = --'explicitEuler',
               'explicitRungeKutta', 
               steps = 4,
               -- how to control the timestep
               control = {
                          name = 'cfl',   -- the name of the timestep control mechanism
                          cfl  = 0.95*(3*degree+1)^2/((degree+1)^2),     
                         },
               },
}

-- Mesh definitions --
mesh = { predefined = 'cube',
         origin = { 
                    (-1.0)*cubeLength/2.0,
                    (-1.0)*cubeLength/2.0,
                    (-1.0)*cubeLength/2.0
                  },
         length = cubeLength,
         refinementLevel = level
       }

-- The analytic solution for this testcase is given by the following functions
-- (only correct for epsi = mu = 1):
-- ... definition of temporal angular frequency
w = math.sqrt(amplX^2 + amplY^2)
-- ... E_x = 0.0
function displacementX(x,y,z,t)
  return 0.0
end 
-- ... E_y = 0.0
function displacementY(x,y,z,t)
  return 0.0 --math.sin(amplX*math.pi*x)*math.sin(amplY*math.pi*z)*math.cos(w*t)
end 
-- ... E_z(x,y,z,t) = sin(amplX \pi x) sin(amplY \pi y) cos(w t)
function displacementZ(x,y,z,t)
  return math.sin(amplX*math.pi*x)*math.sin(amplY*math.pi*y)*math.cos(w*t)
end 
-- ... B_x(x,y,z,t) = -\frac{\pi n}{w} sin(m \pi x) cos(n \pi y) sin(w t)
function magneticX(x,y,z,t)
  return (-1.0)*(math.pi*amplY/w)*math.sin(amplX*math.pi*x)*math.cos(amplY*math.pi*y)*math.sin(w*t)
end 
-- ... B_y(x,y,z,t) = \frac{\pi m}{w} cos(m \pi x) sin(n \pi y) sin(w t)
function magneticY(x,y,z,t)
   return (math.pi*amplX/w)*math.cos(amplX*math.pi*x)*math.sin(amplY*math.pi*y)*math.sin(w*t)
end 
-- ... B_z = 0.0
function magneticZ(x,y,z,t)
  return 0.0 --(math.pi*amplX/w)*math.cos(amplX*math.pi*x)*math.sin(amplY*math.pi*z)*math.sin(w*t )
end

-- ...the initial condition table. 
-- ...initial condition function for displacement field (z component)
function ic_displacementZ(x,y,z)
  return displacementZ(x,y,z,0.0)
end

initial_condition = { 
                      displacement_fieldX = 0.0,           -- displacement field , x component
                      displacement_fieldY = 0.0,  -- displacement field , z component
                      displacement_fieldZ = ic_displacementZ,  -- displacement field , z component
                      magnetic_fieldX = 0.0,  -- magnetic induction , x component
                      magnetic_fieldY = 0.0,  -- magnetic induction , y component
                      magnetic_fieldZ = 0.0,           -- magnetic induction , z component

                    --  projOverSampling = 2.0,
                    --  useFpt = true
                    }

