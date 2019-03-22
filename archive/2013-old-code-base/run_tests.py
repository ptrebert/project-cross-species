#!/usr/bin/python3.1

import unittest as unittest


from tests.test_wigfix_to_bed import TestWigFixToBed
from tests.test_datastructs import TestRegionCluster
from tests.test_mapping import TestMappingProcessBlank
from tests.test_mapping import TestMapFiles
from tests.test_prediction import TestPrepareDatasets

if __name__ == "__main__":
	suite_wigfix = unittest.TestLoader().loadTestsFromTestCase(TestWigFixToBed)
	suite_regclust = unittest.TestLoader().loadTestsFromTestCase(TestRegionCluster)
	suite_mapprocblank = unittest.TestLoader().loadTestsFromTestCase(TestMappingProcessBlank)
	suite_mapfiles = unittest.TestLoader().loadTestsFromTestCase(TestMapFiles)
	suite_prepdata = unittest.TestLoader().loadTestsFromTestCase(TestPrepareDatasets)
	unittest.TextTestRunner().run(suite_wigfix)
	unittest.TextTestRunner().run(suite_regclust)
	unittest.TextTestRunner().run(suite_mapprocblank)
	unittest.TextTestRunner().run(suite_mapfiles)
	unittest.TextTestRunner().run(suite_prepdata)
