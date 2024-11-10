folder = 'mesh/'
timing_file = 'sdr_timing.res'
minlevel = 3

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = { origin = { 0.0, 0.0, 0.0 },
                  length = 1 }

spatial_object = {
  { attribute = { kind = 'seed',  },
    geometry = { kind = 'canoND', object = { origin = { 
                                                        1/2,
                                                        1/2,
                                                        1/2-1e-5,
                                                      },
                                           }
               }
  },
  { attribute = { kind = 'refinement', level = minlevel+1, },
    geometry = { kind = 'canoND', object = { origin = { 
                                                        1/2-1e-5,
                                                        1/2-1e-5,
                                                        1/2-1e-5,
                                                      },
                                           }
               }
  },
}
