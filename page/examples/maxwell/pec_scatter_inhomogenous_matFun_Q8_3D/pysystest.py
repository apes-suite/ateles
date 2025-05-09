__pysys_title__   = r""" Maxwell: Q8 PEC scatter inhomogenous material function """ 
#                        ==========================
__pysys_purpose__ = r""" Check for the Maxwell equations with PEC scattering material function and Q8 polynomials """ 
    
__pysys_created__ = "2025-05-07"
#__pysys_skipped_reason__   = "Skipped until Bug-1234 is fixed"

#__pysys_traceability_ids__ = "Bug-1234, UserStory-456" 
#__pysys_groups__           = "myGroup, disableCoverage, performance"
#__pysys_modes__            = lambda helper: helper.inheritedModes + [ {'mode':'MyMode', 'myModeParam':123}, ]
#__pysys_parameterized_test_modes__ = {'MyParameterizedSubtestModeA':{'myModeParam':123}, 'MyParameterizedSubtestModeB':{'myModeParam':456}, }

import pysys.basetest, pysys.mappers
from pysys.constants import *

trackfile = 'pec_scatter_maxwell_modg_probe_displacementField_Q8_p00000.res'

from apes.apeshelper import ApesHelper
class PySysTest(ApesHelper, pysys.basetest.BaseTest):
    def setup(self):
        self.apes.setupAteles(atlfile='ateles_maxwell_modg.lua')
        self.deleteFile(trackfile)

    def execute(self):
        atlrun = self.apes.runAteles(np = 4)

    def validate(self):
        self.apes.checkAtlLog()
        self.assertPathExists(trackfile, abortOnError = True)
        self.apes.assertIsClose(trackfile, ref_file = 'ref_' + trackfile)

