-- Use this file as template. Do not modify this file for running some testcases

outputname = 'pec_scatter'
outputpreview = true 
folder = 'mesh/'
timing_file = 'sdr_timing.res'

-- The global refinement level of the mesh
level = 3
minlevel = level

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
                  origin = {-1.0, -1.0, -1.0},
                  length = 2.0
                }

-- the refinement level of the PEC boundaries at the outside
geom_level = level

-- the refinemente level of the PEC scatterer in the domain
scatter_level = level+2

-- the physical length of the scattering cube 
scatterLength = 0.01

-- now we define the geometrical objects in our domain.
eps = bounding_cube.length/2^(level+1)
scatter_eps = bounding_cube.length/2^(scatter_level+1)
scatter_dx = bounding_cube.length/2^(scatter_level)
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
                                       origin={-0.5-eps,-1.0+eps,-1.0+eps},
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
                                       origin={0.5+eps,-1.0+eps,-1.0+eps},
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
                                       origin={-1.0+eps,-0.5-eps,-1.0+eps},
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
                                       origin={-1.0+eps,0.5+eps,-1.0+eps},
                                      },
                          }
              },
              -- periodic bnd at top and bottom
              {
               attribute = {
                             kind = 'periodic',
                             label = 'periodic',
                             level = level,
                           },
               geometry = {
                             kind = 'periodic',
                             object = 
                                      {
                                       plane1 = {
                                                  vec= {
                                                         {2.0,0.0,0.0},
                                                         {0.0,2.0,0.0},
                                                       },
                                                  origin={-1.0+eps/2.0,-1.0+eps/2.0,2.0*eps+eps/4.0},
                                                },
                                       plane2 = {
                                                  vec= {
                                                         {0.0,2.0,0.0},
                                                         {2.0,0.0,0.0},
                                                       },
                                                  origin={-1.0+eps/2.0,-1.0+eps/2.0,-2.0*eps-eps/4.0},
                                                },
                                      }
                          }
              },
              -- pec scatterer
              {
               attribute = {
                             kind = 'boundary',
                             label = 'pecScatter',
                             level = scatter_level,
                             distance_refine = {
                               radius = 3*scatter_dx,
                               level_offset = 0
                             }
                           },
               geometry = {
                             kind = 'canoND',
                             object = 
                                      {
                                       vec= {
                                              {scatterLength,0.0,0.0},
                                              {0.0,scatterLength,0.0},
                                              {0.0,0.0,scatterLength},
                                            },
                                       origin={
                                               --scatter_eps,
                                               --scatter_eps,
                                               --scatter_eps,
                                               (-1.0)*scatterLength/2.0,
                                               (-1.0)*scatterLength/2.0,
                                               (-1.0)*scatterLength/2.0,
                                              },
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
                                       origin = {scatterLength+eps, 0.0,0.0 },
                                    },
                         },
              },
        }

