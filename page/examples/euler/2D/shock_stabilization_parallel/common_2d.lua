require('seeder')

-- global simulation options
simulation_name = 'euler_2d'
sim_control = {
                time_control = {   
                                min = 0, 
                                max = 0.15,
                                interval = {iter = 100},
                               }
              }
check = {
          interval = 1,
        }


-- Mesh definitions --
mesh = './mesh/'

-- Equation definitions --
equation = {
    name   = 'euler_2d',
    therm_cond = 2.555e-02,
    isen_coef = 1.4,
    r      = 296.0,
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

--restart = {
--           write = './restart/',
--           time_control = {   
--                            min = 0.0, 
--                            max = sim_control.time_control.max,
--                            interval = sim_control.time_control.max,
--                          }
--          }


-- The state right of the shock
rho_r = 1.0
u_r = 0.0
p_r = 1.0
mach_r = u_r/math.sqrt( equation.isen_coef * p_r / rho_r )

-- Shock properties
shockMach = 2.0
shockXCoord = -1.2
shockSpeed = shockMach * math.sqrt(equation.isen_coef * p_r / rho_r )

-- The state left of the shock (evaluated by Rankine-Huginoit condition)
chi = ( u_r - shockSpeed ) / math.sqrt(equation.isen_coef * p_r / rho_r )
rho_l = rho_r * ( ((equation.isen_coef+1)*chi*chi) / ((equation.isen_coef-1)*chi*chi+2) )  
u_l = shockSpeed + ( u_r - shockSpeed ) * (rho_r/rho_l) 
p_l = p_r * ( (2*equation.isen_coef*chi*chi-(equation.isen_coef-1)) / (equation.isen_coef+1) ) 
mach_l = u_l/math.sqrt( equation.isen_coef * p_l / rho_l )

print('Mach number relevant for inflow: ', mach_l)
print('Mach number relevant for outflow: ', mach_r)

function rho(x,y,z)
    if( y < channel_length/3.0 ) then
        return rho_l
    else
        return rho_r
    end
end

function p(x,y,z)
    if( y < channel_length/3.0 ) then
        return p_l
    else
        return p_r
    end
end

function u(x,y,z)
    if( y < channel_length/3.0 ) then
        return u_l
    else
        return u_r
    end
end

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
           -- lobattoPoints = false  -- if lobatto points should be used, default = false
              factor = 2.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
              blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }

initial_condition = { 
                      density = rho,
                      pressure = p,
                      velocityX = 0.0,
                      velocityY = u,
             }

-- Scheme definitions --
filter_order = 14
scheme = {
           -- the spatial discretization scheme
           spatial =  {
                        name = 'modg_2d',             
                        m = 31, --TODO 127, --TODO 63, --TODO 31, --TODO 15,                
                      }, 
           ---- the stabilzation of the scheme
           stabilization = {
                             {
                              name = 'spectral_viscosity',
                              alpha = 36,
                              order = filter_order,
                             },
                             {
                              name = 'covolume',
                              alpha = 36,
                              order = filter_order,
                              beta = 1.0,
                             },
                           },
}

-- Boundary conditions
boundary_condition = {
        { 
         label = 'inlet', 
         density = rho_l,
         v_norm = u_l,
         v_tan = 0.0,
         pressure = p_l,
	}
         ,
        { 
         label = 'outlet', 
         kind = 'outflow',
         pressure = p_r,
         }
         ,
        { 
         label = 'bottom', 
         kind = 'slipwall', 
         }
         ,
        { 
         label = 'top', 
         kind = 'slipwall', 
         }
         ,
        { 
         label = 'south', 
         kind = 'slipwall', 
         }
         ,
        { 
         label = 'north', 
         kind = 'slipwall', 
         }
}
if mach_l > 1  then
  boundary_condition[1].kind = 'supersonic_inflow_normal'
else
  boundary_condition[1].kind = 'inflow_normal'
end
