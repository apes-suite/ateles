-- Use this file as template. Do not modify this file for running some testcases

outputname = 'shearLayer'
outputpreview = true 
folder = 'mesh/'
timing_file = 'sdr_timing.res'

-- ATTENTION: the minimum is 3, otherwise you won't get any fluid element.
minlevel = 4 -- 3

cubeLength = 2.0

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
                 origin = {cubeLength/(-2.0),cubeLength/(-2.0),cubeLength/(-2.0)},
                 length = cubeLength
                }


geom_level = minlevel -- +3
eps = bounding_cube.length/2^(geom_level+3)
elemSize = bounding_cube.length/2^(minlevel)
spatial_object = 
                {
                  -- southern slip bnd
                  {
                   attribute = {
                                 kind = 'boundary',
                                 label = 'slipSouth',
                                 level = minlevel,
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
                  -- northern slip bnd
                  {
                   attribute = {
                                 kind = 'boundary',
                                 label = 'slipNorth',
                                 level = minlevel,
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
                  -- bottom slip bnd
                  {
                   attribute = {
                                 kind = 'boundary',
                                 label = 'slipBottom',
                                 level = minlevel,
                               },
                   geometry = {
                                 kind = 'canoND',
                                 object = 
                                          {
                                           vec= {
                                                  {2.0,0.0,0.0},
                                                  {0.0,2.0,0.0},
                                                },
                                           origin={-1.0+eps,-1.0+eps,-eps},
                                          },
                              }
                  },
                  -- top slip bnd
                  {
                   attribute = {
                                 kind = 'boundary',
                                 label = 'slipTop',
                                 level = minlevel,
                               },
                   geometry = {
                                 kind = 'canoND',
                                 object = 
                                          {
                                           vec= {
                                                  {2.0,0.0,0.0},
                                                  {0.0,2.0,0.0},
                                                },
                                           origin={-1.0+eps,-1.0+eps,elemSize+eps},
                                          },
                              }
                  },
                  -- refinement box around the transition zone
                  {
                   attribute = {
                                 kind = 'refinement',
                                 level = geom_level,
                               },
                   geometry = {
                                 kind = 'canoND',
                                 object = 
                                          {
                                           vec= {
                                                  {2.0-2*eps,0.0,0.0},
                                                  {0.0,eps,0.0},
                                                  {0.0,0.0,elemSize-eps},
                                                },
                                           origin={-1.0+eps,-0.5*eps,0.5*eps},
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
                                           origin = {0.0,0.0,0.0 },
                                        },
                             },
                  },
                }
	    
