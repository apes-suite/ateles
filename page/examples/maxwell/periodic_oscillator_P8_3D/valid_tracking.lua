--...Configuration of simulation time
-- Run one complete period.
simtime = { 
            max = math.sqrt(2), -- Number of iterations / Simulated time
            min = 0
          }

-- Tracking used for validation.
tracking = { 
             label = 'probe_electricField_'..poly_space..(degree+1),
             folder = './',
             variable = {'displacement_field'}, 
             shape = {kind = 'canoND', object= { origin ={0., 0., 0.} } },
             time_control = {
               min = simtime.min, 
               max = simtime.max, 
               interval = simtime.max/2.0
             },  
             output = { format = 'ascii', ndofs = 1 }, -- 'asciiSpatial', 'harvester', 'convergence'
}
