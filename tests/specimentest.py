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
        self.process_ids = ["BOM1525-10", "BOM1528-10", "GBIR5337-13", "KKBNA817-05", "KKBNA820-05", "KKBNA821-05", 
                            "KKBNA824-05", "KKBNA827-05", "KKBNA830-05", "KKBNA831-05", "KKBNA834-05", "KKBNA840-05", 
                            "KKBNA851-05", "KKBNA852-05", "KKBNA854-05", "KKBNA855-05", "KKBNA856-05", "KKBNA857-05", 
                            "KKBNA859-05", "KKBNA861-05", "KKBNA866-05", "KKBNA867-05", "KKBNA872-05", "KKBNA875-05", 
                            "KKBNA986-06", "BOTW348-05", "GBIR3737-13", "GBIR4109-13", "GBIR5335-13", "KBNA876-04", 
                            "KBNA877-04", "KKBNA825-05", "KKBNA826-05", "KKBNA828-05", "KKBNA833-05", "KKBNA836-05", 
                            "KKBNA839-05", "KKBNA842-05", "KKBNA843-05", "KKBNA846-05", "KKBNA847-05", "KKBNA849-05", 
                            "KKBNA850-05", "KKBNA858-05", "KKBNA863-05", "KKBNA865-05", "KKBNA868-05", "KKBNA871-05", 
                            "KKBNA983-06", "KKBNA990-06", "BOTW347-05", "CDAMH045-05", "CDAMH062-05", "GBIR1167-08"]
        self.taxons = ["Archaeorhizomycetes", "Arthoniomycetes"]
        self.geographies = ["Canada", "United States"]
        self.bins = ["BOLD:AAA5125", "BOLD:AAA5126"]
        #self.containers = ["ACRJP", "ACRJI"]
        self.institutions = ["Biodiversity Institute of Ontario", "Agriculture and Agri-Food Canada"]
        self.researchers = ["Rodolphe Rougerie|Hai D. T. Nguyen"]

        self.specimen_client = pybold.specimen.SpecimensClient()

    def tearDown(self):
        self.process_ids = None
        self.taxons = None
        self.bins = None
        self.institutions = None
        self.researchers = None
        self.specimen_client = None


    def test_isinstance(self):
        self.assertIsInstance(self.specimen_client, pybold.specimen.SpecimensClient, 'specimen_client failed to instantiate as a SpecimensClient object.')
        
    def test_specimen_list_attribute(self):
        self.assertTrue(hasattr(self.specimen_client, 'specimen_list'), 'SpecimensClient does not have a specimen_list attribute')
        self.assertIsInstance(self.specimen_client.specimen_list, list, 'SpecimensClient.specimen_list must be of type list')
       
    def _test_get(self,taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, timeout=5):
        specimen_list = self.specimen_client.get(taxon=taxon, 
                                                 ids=ids, 
                                                 bins=bins, 
                                                 containers=containers, 
                                                 institutions=institutions, 
                                                 researchers=researchers, 
                                                 geo=geo,
                                                 timeout=timeout)
        msg = 'SpecimensClient.get() should return a list of specimens'
        self.assertIsNotNone(specimen_list, msg)
        self.assertIsInstance(specimen_list, list, msg )
        self.assertIsInstance(specimen_list[0], pybold.specimen.Specimen, msg)
        self.assertListEqual(specimen_list, self.specimen_client.specimen_list, 'SpecimensClient.get() and SpecimensClient.specimen_list should be equal')
        return specimen_list

    def test_get_by_processids(self):
        specimen_list = self._test_get(ids='|'.join(self.process_ids))
        for specimen in specimen_list:
            self.assertIn(specimen.process_id, self.process_ids, 'Returned Specimen was not part of original search criteria')
    
    def test_get_by_taxon_and_geo(self):
        specimen_list = self._test_get(taxon='|'.join(self.taxons), geo='|'.join(self.geographies))
        for specimen in specimen_list:
            found = False
            specimen_taxonomy = ' '.join(specimen.taxonomy.values())
            for taxon in self.taxons:
                if taxon in specimen_taxonomy:
                    found = True
                    break
            
            self.assertTrue(found, 'Specimen {} does not match the taxon search critera.'.format(specimen.process_id))
            
            found = False
            specimen_geography = [ specimen.geography.country, specimen.geography.province, specimen.geography.region ]
            for geo in self.geographies:
                if geo in specimen_geography:
                    found = True
                    break
            
            self.assertTrue(found, 'Specimen {} does not match the geography search criteria.'.format(specimen.process_id))
    
    def test_get_by_institution_and_researcher(self):
        specimen_list = self._test_get(institutions='|'.join(self.institutions), researchers='|'.join(self.researchers), timeout=10)
        for specimen in specimen_list:
            found = False
            specimen_researchers = ' '.join([str(specimen.record.taxonomy.identification_provided_by), str(specimen.record.collection_event.collectors)])
            for taxon in self.researchers:
                if taxon in specimen_researchers:
                    found = True
                    break
            
            self.assertTrue(found, 'Specimen {} does not match the researcher search critera.'.format(specimen.process_id))
            
            found = False
            specimen_institution = [ specimen.record.specimen_indentifiers.institution_storing ]
            for geo in self.institutions:
                if geo in specimen_institution:
                    found = True
                    break
            
            self.assertTrue(found, 'Specimen {} does not match the institution search criteria.'.format(specimen.process_id))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()