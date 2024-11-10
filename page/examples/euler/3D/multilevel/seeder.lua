-- Minimalistic mesh describing a channel with left half being refined once,
-- boundaries in the west and east and all other directions periodic.

-- Size of coarse element
dx = 0.5

-- Length of the channel 
Lch = 1.0

folder = "mesh_"
timing_file = 'sdr_timing.res'

ebug = {debugMode = true, debugMesh = 'debug/'}

minlevel = math.ceil(math.log((Lch/dx)+2)/math.log(2))
cubeLength = dx * 2^minlevel

bounding_cube = {
                 origin = {-dx,-dx,-dx},
                 length = cubeLength
                }

smoothbounds = false

eps = cubeLength/2^(20)
spatial_object = { --> spatial objects
                  
  {--> western inflow bnd
     attribute = {
       kind = 'boundary',
       label = 'west',
       level = minlevel+1
     },
     geometry = {
        kind = 'canoND',
        object = {
          vec= {
                 {0.0, dx+eps,    0.0},
                 {0.0,    0.0, dx+eps}
               },
          origin = { -eps, -eps, -eps }
        }
     }
  }, --< western inflow bnd


  {--> eastern outflow bnd
     attribute = {
       kind = 'boundary',
       label = 'east',
       level = minlevel
     },
     geometry = {
        kind = 'canoND',
        object = {
          vec= {
                 {0.0, dx+eps,    0.0},
                 {0.0,    0.0, dx+eps}
               },
          origin = { Lch+eps, -eps, -eps }
        }
     }
  }, --< eastern outflow bnd

  {--> Refinement of left half of the channel
     attribute = {
       kind = 'refinement',
       level = minlevel+1
     },
     geometry = {
        kind = 'canoND',
        object = {
          vec= {
                 {dx,     0.0,    0.0},
                 {0.0, dx+eps,    0.0},
                 {0.0,    0.0, dx+eps}
               },
          origin = { -eps, -eps, -eps }
        }
     }
  }, --< Refinement of left half of the channel


  { --> periodic in Z direction
    attribute = {
      kind = 'periodic',
      level = minlevel
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          vec = {
            {      0.0, dx+2*eps, 0.0},
            {Lch+2*eps,      0.0, 0.0}
          },
          origin = {-eps, -eps, -eps}
        },
        plane2 = {
          vec = {
            {Lch+2*eps,      0.0, 0.0},
            {      0.0, dx+2*eps, 0.0}
          },
          origin={ -eps, -eps, dx+eps}
        }
      }
    }
  }, --< periodic in X direction


  { --> periodic in Y direction
    attribute = {
      kind = 'periodic',
      level = minlevel
    },
    geometry = {
       kind = 'periodic',
       object = {
         plane1 = {
           vec= {
                  {Lch+2*eps, 0.0,      0.0},
                  {      0.0, 0.0, dx+2*eps}
                },
           origin = {-eps, -eps, -eps}
         },
         plane2 = {
           vec= {
                  {0.0,       0.0, dx+2*eps},
                  {Lch+2*eps, 0.0,      0.0}
                },
           origin={-eps, dx+eps, -eps}
         }
       }
    }
  }, --< periodic in Y direction
 
  { --> seed 
    attribute = { kind = 'seed' },
    geometry = {
      kind = 'canoND',
      object = { { origin = {Lch/4, dx/2, dx/2} } }
    }
  } --< seed

} --< spatial objects
