__pysys_title__   = r""" Euler: Toro setup in 3D in X """
#                        ==========================
__pysys_purpose__ = r""" Check for the Euler equations with a Riemann problem in 3D in X """
    
__pysys_created__ = "2025-05-07"
#__pysys_skipped_reason__   = "Skipped until Bug-1234 is fixed"

#__pysys_traceability_ids__ = "Bug-1234, UserStory-456" 
#__pysys_groups__           = "myGroup, disableCoverage, performance"
#__pysys_modes__            = lambda helper: helper.inheritedModes + [ {'mode':'MyMode', 'myModeParam':123}, ]
#__pysys_parameterized_test_modes__ = {'MyParameterizedSubtestModeA':{'myModeParam':123}, 'MyParameterizedSubtestModeB':{'myModeParam':456}, }

import pysys.basetest, pysys.mappers
from pysys.constants import *

trackfile = 'toro1_x_euler_modg_probe_density_Q4_toro_x_p00000.res'

from apes.apeshelper import ApesHelper
class PySysTest(ApesHelper, pysys.basetest.BaseTest):
    def setup(self):
        self.apes.setupAteles()
        self.deleteFile(trackfile)

    def execute(self):
        atlrun = self.apes.runAteles(np = 1)

    def validate(self):
        self.apes.checkAtlLog()
        self.assertPathExists(trackfile, abortOnError = True)
        self.apes.assertIsClose(trackfile, ref_file = 'ref_' + trackfile)


