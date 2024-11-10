 format = 'ascii'
 solver = 'ATELES_v0.9'
 simname = 'shear_layer_modg'
 basename = 'shear_layer_modg_probe_momentum_Q4'
 glob_rank = 0
 glob_nprocs = 1
 sub_rank = 0
 sub_nprocs = 1
 resultfile = 'shear_layer_modg_probe_momentum_Q4_p*'
 nDofs = 1
 nElems = 1
 time_control = {
    min = {
        sim =    0.000000000000000E+00 
    },
    max = {
        sim =   10.000000000000001E-06 
    },
    interval = {
        sim =    5.000000000000000E-06 
    },
    check_iter = 1 
}
 shape = {
    {
        kind = 'canoND',
        object = {
            {
                origin = {    0.000000000000000E+00,    0.000000000000000E+00,    0.000000000000000E+00 },
                distribution = 'equal' 
            } 
        } 
    } 
}
 varsys = {
    systemname = 'euler_conservative',
    variable = {
        {
            name = 'momentum',
            ncomponents = 3,
            state_varpos = { 2, 3, 4 } 
        } 
    },
    nScalars = 3,
    nStateVars = 1,
    nAuxScalars = 0,
    nAuxVars = 0 
}
