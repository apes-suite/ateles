--           2 
--          /
--       -------    
--      /      /|
--     /  5   / |
--     -------  | 4
--    |  |   |  |
-- 3  |  |   |  | 
--    |   ---|---
--    | /  6 | /
--    |/     |/
--     ------
--     /
--    1      
--       
printRuntimeInfo = false
outputname= 'box_4'
comment = 'box'
minlevel = 4
folder = 'mesh/'
level = 4
timing_file = 'sdr_timing.res'

--debug = {debugMode=true, debugFiles=true, debugMesh='debug/'}

bounding_cube = { origin = {-4.0, -4.0, -4.0}, length = 8.0 }

eps=bounding_cube.length*0.5^(level+1)

spatial_object = {
  { attribute = { kind = 'seed', label = 'seed', },
    geometry = {
      kind = 'canoND',
      object = { origin = { 0.0, 0.0, 0.0 },
      }
    }
  }, -- seed
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_1',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { (-2.0-eps), (-2.0-eps), (2.0+eps) },
        vec = { { (4.0+2*eps), 0.0, 0.0 },
                { 0.0, (4.0+2*eps), 0.0 },
        },
        only_surface = true,
      } -- object
    },
  },
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_2',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { (-2.0-eps), (-2.0-eps), (-2.0-eps) },
        vec = { { (4.0+2*eps), 0.0, 0.0 },
                { 0.0, (4.0+2*eps), 0.0 },
        },
        only_surface = true,
      } -- object
    },
  },
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_3',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { (-2.0-eps), (-2.0-eps), (2.0+eps) },
        vec = { { 0.0, 0.0, (-4.0-2*eps) },
                { 0.0, (4.0+2*eps), 0.0 },
        },
        only_surface = true,
      } -- object
    },
  }, 
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_4',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { (2.0+eps), (-2.0-eps), (2.0+eps) },
        vec = { { 0.0, 0.0, (-4.0-2*eps) },
                { 0.0, (4.0+2*eps), 0.0 },
        },
        only_surface = true,
      } -- object
    },
  }, 
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_5',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {(-2.0-eps), (2.0+eps),(2.0+eps) },
        vec = { { (4.0+2*eps), 0.0, 0.0 },
                { 0.0, 0.0, (-4.0-2*eps) },
        },
        only_surface = true,
      } -- object
    },
  }, 
  --------------------------------------------
  { attribute = {
      kind = 'boundary', label = 'wall_6',
      level = level, calc_dist = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { (-2.0-eps), (-2.0-eps), (2.0+eps) },
        vec = { { (4.0+2*eps), 0.0, 0.0 },
                { 0.0, 0.0, (-4.0-2*eps) },
        },
        only_surface = true,
      } -- object
    },
  }, 
} -- spatial object
