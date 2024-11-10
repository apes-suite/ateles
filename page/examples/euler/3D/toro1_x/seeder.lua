-- CHANGE VALUES ONLY HERE --
nelems = 100 --os.getenv('nelem')         -- number of elements through the channel
nelems_side = 1 --os.getenv('sidelem')  -- number of elements on side directions
channel_length = 1.0                -- physical length of the channel

-- DO NOT CHANGE VALUES BELOW (UNLESS NECESSARY) --
epsx = channel_length / (2^21) -- an epsilon to offset geometry
dx = channel_length / nelems
channel_width = nelems_side*dx    -- width of the channel
reflevel = math.floor(math.log(2*nelems+4) / math.log(2)) + 1
bounding_length = (2^reflevel)*dx


folder = 'mesh/'
timing_file = 'sdr_timing.res'
minlevel = reflevel

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = { origin = { -dx, -dx, -dx },
                  length = bounding_length }

spatial_object = {
  { attribute = { kind = 'seed' },
    geometry = { kind = 'canoND', object = { origin = { channel_length/2.0,
                                                        channel_width/2.0,
                                                        channel_width/2.0
                                                      }
                                           }
               }
  },
  { attribute = { kind = 'boundary', label = 'inlet', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { -epsx, -epsx, -epsx },
                            vec = { { 0.0, channel_width+2*epsx, 0.0 },
                                    { 0.0, 0.0, channel_width+2*epsx }
                                  }
                          }
               }
  },
  { attribute = { kind = 'boundary', label = 'outlet', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { channel_length+epsx, -epsx, -epsx },
                            vec = { { 0.0, channel_width+2*epsx, 0.0 },
                                    { 0.0, 0.0, channel_width+2*epsx }
                                  }
                          }
               }
  },
  { attribute = { kind = 'boundary', label = 'bottom', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { -epsx, -epsx, -epsx },
                            vec = { { channel_length+2*epsx, 0.0, 0.0 },
                                    { 0.0, channel_width+2*epsx, 0.0 }
                                  }
                          }
               }
  },
  { attribute = { kind = 'boundary', label = 'top', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { -epsx, -epsx, channel_width+epsx },
                            vec = { { channel_length+2*epsx, 0.0, 0.0 },
                                    { 0.0, channel_width+2*epsx, 0.0 }
                                  }
                          }
               }
  },
  { attribute = { kind = 'boundary', label = 'south', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { -epsx, -epsx, -epsx },
                            vec = { { channel_length+2*epsx, 0.0, 0.0 },
                                    { 0.0, 0.0, channel_width+2*epsx }
                                  }
                          }
               }
  },
  { attribute = { kind = 'boundary', label = 'north', level = reflevel },
    geometry = { kind = 'canoND',
                 object = { origin = { -epsx, channel_width+epsx, -epsx },
                            vec = { { channel_length+2*epsx, 0.0, 0.0 },
                                    { 0.0, 0.0, channel_width+2*epsx }
                                  }
                          }
               }
  }
}
