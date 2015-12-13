'''
Created on 2015-12-12

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca>
'''
from nose.tools import assert_is_instance
import unittest

from pybold.specimen import SpecimensClient
import pybold.specimen


class SpecimensClientTest(unittest.TestCase):
    def setUp(self):
        self.process_ids1 = ["BOM1525-10", "BOM1528-10", "GBIR5337-13", "KKBNA817-05", "KKBNA820-05", "KKBNA821-05", 
                            "KKBNA824-05", "KKBNA827-05", "KKBNA830-05", "KKBNA831-05", "KKBNA834-05", "KKBNA840-05", 
                            "KKBNA851-05", "KKBNA852-05", "KKBNA854-05", "KKBNA855-05", "KKBNA856-05", "KKBNA857-05", 
                            "KKBNA859-05", "KKBNA861-05", "KKBNA866-05", "KKBNA867-05", "KKBNA872-05", "KKBNA875-05", 
                            "KKBNA986-06", "BOTW348-05", "GBIR3737-13", "GBIR4109-13", "GBIR5335-13", "KBNA876-04", 
                            "KBNA877-04", "KKBNA825-05", "KKBNA826-05", "KKBNA828-05", "KKBNA833-05", "KKBNA836-05", 
                            "KKBNA839-05", "KKBNA842-05", "KKBNA843-05", "KKBNA846-05", "KKBNA847-05", "KKBNA849-05", 
                            "KKBNA850-05", "KKBNA858-05", "KKBNA863-05", "KKBNA865-05", "KKBNA868-05", "KKBNA871-05", 
                            "KKBNA983-06", "KKBNA990-06", "BOTW347-05", "CDAMH045-05", "CDAMH062-05", "GBIR1167-08"]
        self.process_ids2 = self.process_ids1[0:int(len(self.process_ids1)/2)]
        self.specimen_client1 = pybold.specimen.SpecimensClient()
        self.specimen_client2 = pybold.specimen.SpecimensClient()

    def tearDown(self):
        pass


    def test_isinstance(self):
        msg = '{} should be an instance of SpecimensClient'
        self.assertIsInstance(self.specimen_client1, pybold.specimen.SpecimensClient, msg.format('specimen_client1'))
        self.assertIsInstance(self.specimen_client2, pybold.specimen.SpecimensClient, msg.format('specimen_client2'))
        
    def test_specimen_list_attribute(self):
        msg = 'SpecimensClient does not have a specimen_list attribute'
        self.assertTrue(hasattr(self.specimen_client1, 'specimen_list'), msg)
        self.assertTrue(hasattr(self.specimen_client2, 'specimen_list'), msg)
        msg = 'SpecimensClient.specimen_list must be of type list'
        self.assertIsInstance(self.specimen_client1.specimen_list, list, msg)
        self.assertIsInstance(self.specimen_client2.specimen_list, list, msg)
    
    def test_get(self):
        specimen_list = self.specimen_client1.get(ids='|'.join(self.process_ids1))
        msg = 'SpecimensClient.get() should return a list of specimens'
        self.assertIsNotNone(specimen_list, msg)
        self.assertIsInstance(specimen_list, list, msg )
        self.assertIsInstance(specimen_list[0], pybold.specimen.Specimen, msg)
        self.assertListEqual(specimen_list, self.specimen_client1.specimen_list, 'SpecimensClient.get() and SpecimensClient.specimen_list should be equal')
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()