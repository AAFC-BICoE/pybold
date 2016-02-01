'''
Created on 2015-12-12

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca>
'''
import unittest

import pybold.sequence
import pybold.specimen
import pybold.tracefile


class Test(unittest.TestCase):


    def setUp(self):
        self.process_ids = ["BOM1525-10", "BOM1528-10", "GBIR5337-13", "KKBNA817-05", "KKBNA820-05", "KKBNA821-05", 
                            "KKBNA824-05", "KKBNA827-05", "KKBNA830-05", "KKBNA831-05", "KKBNA834-05", "KKBNA840-05", 
                            "KKBNA851-05", "KKBNA852-05", "KKBNA854-05", "KKBNA855-05", "KKBNA856-05", "KKBNA857-05"]
        self.record_ids = ["6283631", "6283648", "6283659", "6283679", "5064977", "5064989", "5065007", "5065009", 
                           "5065014", "6039995", "6039996", "6040002", "6040009", "6040013", "6040016", "6040035", 
                           "6040047", "6040055", "6040056", "6040057", "6040061", "6040062", "6040077", "6040078"]
        self.taxons = ["Archaeorhizomycetes", "Arthoniomycetes"]
        self.geographies = ["Canada", "United States"]
        self.bins = ["BOLD:AAA5125", "BOLD:AAA5126"]
        #self.containers = ["ACRJP", "ACRJI"]
        self.institutions = ["Biodiversity Institute of Ontario", "Agriculture and Agri-Food Canada"]
        self.researchers = ["Rodolphe Rougerie","Hai D. T. Nguyen"]
        self.markers = ["ITS", "COI"]
        # Instantiate the SequencesClient class to be used throughout all tests
        self.tracefiles_client = pybold.tracefile.TracefilesClient()

    def tearDown(self):
        self.process_ids = None
        self.record_ids = None
        self.taxons = None
        self.bins = None
        self.institutions = None
        self.researchers = None
        self.sequences_client = None
    
    def test_isinstance(self):
        self.assertIsInstance(self.tracefiles_client, pybold.tracefile.TracefilesClient, "tracefiles_client failed to instantiate as a TracefilesClient")
    
    def test_tracefile_list_attribute(self):
        self.assertTrue(hasattr(self.tracefiles_client, 'tracefile_list'), 'TracefilesClient does not have a tracefile_list attribute')
        self.assertIsInstance(self.tracefiles_client.tracefile_list, list, 'TracefilesClient.tracefile_list must be of type list')
    
    def test_get_process_ids(self):
        self.tracefiles_client.get(ids='|'.join(self.process_ids))
        process_ids = self.tracefiles_client.get_process_ids()
        for process_id in process_ids:
            self.assertIn(process_id, self.process_ids, 'Returned tracefile process_id ({}) was not part of original search criteria'.format(process_id))
    
    def test_get_sequences(self):
        self.tracefiles_client.get(ids='|'.join(self.process_ids))
        sequences = self.tracefiles_client.get_sequences()
        self.assertIsInstance(sequences, list, "TracefilesClient.get_sequences() should return a list of Sequence.")
        for sequence in sequences:
            self.assertIsInstance(sequence, pybold.sequence.Sequence, "TracefilesClient.get_sequences() should return a list of Sequence.")
    
    def test_get_specimens(self):
        self.tracefiles_client.get(ids='|'.join(self.process_ids))
        specimen_list = self.tracefiles_client.get_specimens()
        self.assertIsInstance(specimen_list, list, "TracefilesClient.get_specimens() should return a list of Specimen.")
        for specimen in specimen_list:
            self.assertIsInstance(specimen, pybold.specimen.Specimen, "TracefilesClient.get_specimens() should return a list of Specimen.")

    def _test_get(self,taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, marker=None, timeout=5):
        tracefile_list = self.tracefiles_client.get(taxon=taxon, 
                                                 ids=ids, 
                                                 bins=bins, 
                                                 containers=containers, 
                                                 institutions=institutions, 
                                                 researchers=researchers, 
                                                 geo=geo,
                                                 marker=marker,
                                                 timeout=timeout)
        msg = 'SequencesClient.get() should return a list of Tracefile objects'
        self.assertIsNotNone(tracefile_list, msg)
        self.assertIsInstance(tracefile_list, list, msg )
        self.assertIsInstance(tracefile_list[0], pybold.tracefile.Tracefile, msg)
        self.assertListEqual(tracefile_list, self.tracefiles_client.tracefile_list, 'TracefilesClient.get() and TracefilesClient.sequence_list should be equal')
        return tracefile_list
    
    def test_get_by_processids(self):
        tracefile_list = self._test_get(ids='|'.join(self.process_ids))
        for tracefile in tracefile_list:
            self.assertIn(tracefile.process_id, self.process_ids, 'Returned Tracefile process_id ({}) does not match original search criteria'.format(tracefile.process_id))

    def test_get_by_taxon_and_marker(self):
        tracefile_list = self._test_get(taxon='|'.join(self.taxons), marker='|'.join(self.markers))
        for tracefile in tracefile_list:
            found = False
            for marker in self.markers:
                if marker in tracefile.marker:
                    found = True
            self.assertTrue(found, "Returned Sequence was not part of original search criteria for markers {}.".format(', '.join(self.markers)))
             
        #TODO validate taxons

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()