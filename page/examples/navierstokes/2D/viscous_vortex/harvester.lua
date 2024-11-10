require('ateles')

-- define the input
name = 'gPulseDens_euler_modg'
input = {
  read = './restart/simulation_lastHeader.lua',

  subsampling = {
    levels = 9-level,
    projection = 'QLegendrePoint',
  },
  add_variable = {
    {
       name = 'density_analytic',
       ncomponents = 1,
       density_analytic = {
         kind = 'lua_fun',
         fun = ini_dens,
       },
    },
    {
       name = 'velocity_x_analytic',
       ncomponents = 1,
       velocity_x_analytic = {
         kind = 'lua_fun',
         fun = ini_vel_x,
       },
    },
    {
       name = 'velocity_y_analytic',
       ncomponents = 1,
       velocity_y_analytic = {
         kind = 'lua_fun',
         fun = ini_vel_y,
       },
    },
    {
       name = 'pressure_analytic',
       ncomponents = 1,
       pressure_analytic = {
         kind = 'lua_fun',
         fun = ini_press,
       }
    }
  }
}

-- define the output
eps = 1e-4
nElems_perDir = 2^(input.subsampling.levels+level)
print('Number of elements after post-processing (2D): ',nElems_perDir^2)
post_elemsize = cubeLength/nElems_perDir
print('Element size after post-processing: ', post_elemsize)
output = {
  folder = './harvest/',
  {
     format = 'VTU',
     binary = true,
     vrtx = { },
     requestedData = {
       variable = {
         { name = 'density', ncomponents = 1 },
         { name = 'momentum', ncomponents = 2 },
         { name = 'energy', ncomponents = 1 },
         { name = 'pressure', ncomponents = 1 },
         { name = 'temperature', ncomponents = 1 },
         { name = 'velocity', ncomponents = 2 },
         { name = 'density_analytic', ncomponents = 1 },
         { name = 'velocity_x_analytic', ncomponents = 1 },
         { name = 'velocity_y_analytic', ncomponents = 1 },
         { name = 'pressure_analytic', ncomponents = 1 },
         { name = 'difference', ncomponents = 1,
           dep = {'pressure','pressure_analytic'} },
       },
     },
     shape = {
       kind = 'canoND',
       object = {
         origin = { 0,0,eps },
         vec ={
           {cubeLength,0,0},
           {0,cubeLength,0},
         },
         segments = {2*cubeLength/post_elemsize,2*cubeLength/post_elemsize},
       }
     }
  }
}
