 format = 'ascii'
 solver = 'ATELES_v0.9'
 simname = 'lineareuler'
 basename = './lineareuler_transientBackground'
 glob_rank = 0
 glob_nprocs = 1
 sub_rank = 0
 sub_nprocs = 1
 resultfile = './lineareuler_transientBackground_p*'
 nDofs = 1
 nElems = 1
 time_control = {
    min = {
        sim =    0.000000000000000E+00 
    },
    max = {
        sim =   10.000000000000000E-03 
    },
    interval = {
        sim =  500.000000000000010E-06 
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
    systemname = 'LinearEuler',
    variable = {
        {
            name = 'density',
            ncomponents = 1,
            state_varpos = { 1 } 
        },
        {
            name = 'completeState',
            ncomponents = 5 
        } 
    },
    nScalars = 6,
    nStateVars = 2,
    nAuxScalars = 0,
    nAuxVars = 0 
}
