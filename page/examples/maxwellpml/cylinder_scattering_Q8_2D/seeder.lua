-- Use this file as template. Do not modify this file for running some testcases

outputname = 'inhomMax'

folder = 'mesh/' --default: 'mesh_'
timing_file = 'sdr_timing.res'

-- ATTENTION: minlevel has to be larger than 1
minlevel = 3

cubeLength = 4.0

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
                 origin = {cubeLength/(-2.0),cubeLength/(-2.0),cubeLength/(-2.0)},
                 length = cubeLength
                }

-- smallness parameter
eps = bounding_cube.length/(2^(minlevel+5))

-- the size of an element
elemsize = bounding_cube.length/(2^(minlevel))

spatial_object = 
                {
                  -- periodic boundary in z direction
                  {
                    attribute = {
                                 kind = 'periodic',
                                 level = minlevel,
                                },
                    geometry = {
                                 kind = 'periodic',
                                 object = { 
                                            plane1={
                                                      vec= {
                                                             {0.0,bounding_cube.length,0.0},
                                                             {bounding_cube.length,0.0,0.0},
                                                           },
                                                      origin={
                                                               -bounding_cube.length/2+eps,
                                                               -bounding_cube.length/2+eps,
                                                               -eps,
                                                             },
                                                   },
                                            plane2={
                                                      vec= {
                                                             {bounding_cube.length,0.0,0.0},
                                                             {0.0,bounding_cube.length,0.0},
                                                           },
                                                      origin={
                                                               -bounding_cube.length/2+eps,
                                                               -bounding_cube.length/2+eps,
                                                               elemsize+eps,
                                                             },
                                                   },
                                          },
                               },
                  },
                  -- seed 
                  {
                  attribute = {
                               kind = 'seed'
                              },
                  geometry = {
                               kind = 'canoND',
                               object =  
                                        {
                                           origin = {0,0,0},

                                        },
                             },
                  },
                }
