-- Use this file as template. Do not modify this file for running some testcases

outputname = 'rectangular_wave_guide'
outputpreview = true 
folder = 'mesh/'
timing_file = 'sdr_timing.res'

-- The global refinement level of the mesh
level = 3
minlevel = level

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
                  origin = {-1.0, -1.0, 0.0},
                  length = 2.0
                }

eps = bounding_cube.length/2^(level+1)
geom_level = level+1
spatial_object = {
              -- eastern pec bnd
              {
               attribute = {
                             kind = 'boundary',
                             label = 'pecEast',
                             level = geom_level,
                           },
               geometry = {
                             kind = 'canoND',
                             object = 
                                      {
                                       vec= {
                                              {0.0,2.0,0.0},
                                              {0.0,0.0,2.0},
                                            },
                                       origin={-0.5-eps,-1.0, 0.0},
                                      },
                          }
              },
              -- western pec bnd
              {
               attribute = {
                             kind = 'boundary',
                             label = 'pecWest',
                             level = geom_level,
                           },
               geometry = {
                             kind = 'canoND',
                             object = 
                                      {
                                       vec= {
                                              {0.0,2.0,0.0},
                                              {0.0,0.0,2.0},
                                            },
                                       origin={0.5+eps, -1.0, 0.0},
                                      },
                          }
              },
              -- southern pec bnd
              {
               attribute = {
                             kind = 'boundary',
                             label = 'pecSouth',
                             level = geom_level,
                           },
               geometry = {
                             kind = 'canoND',
                             object = 
                                      {
                                       vec= {
                                              {2.0,0.0,0.0},
                                              {0.0,0.0,2.0},
                                            },
                                       origin={-1.0,-0.5-eps, 0.0},
                                      },
                          }
              },
              -- northern pec bnd
              {
               attribute = {
                             kind = 'boundary',
                             label = 'pecNorth',
                             level = geom_level,
                           },
               geometry = {
                             kind = 'canoND',
                             object = 
                                      {
                                       vec= {
                                              {2.0,0.0,0.0},
                                              {0.0,0.0,2.0},
                                            },
                                       origin={-1.0,0.5+eps, 0.0},
                                      },
                          }
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
                                       origin = {0.0, 0.0, 0.0 },
                                    },
                         },
              },
        }

