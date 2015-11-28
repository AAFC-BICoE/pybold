'''
Created on 2015-11-16

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca
'''
from __builtin__ import super, str
import __builtin__
from lxml import objectify
import pickle
import threading 

from pybold import PUBLIC_API_URL, Endpoint
import pybold.sequence


class Specimen(object):
    '''
    classdocs
    '''

    def __init__(self, objectified_element=None, specimen_record_xml=None):
        if isinstance(specimen_record_xml, __builtin__.str):
            obj = objectify.fromstring(specimen_record_xml)
        elif isinstance(objectified_element, objectify.ObjectifiedElement):
            obj = objectified_element
        
        self.record = obj
        self.sequence = None
        
        super(Specimen, self).__init__()

    def get_taxonomy(self):
        if not hasattr(self.record, 'taxonomy'):
            return None
        
        ranks = ('phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species')
        taxonomy = {}
        for rank in ranks:
            name = None
            if hasattr(self.record.taxonomy, rank):
                name = getattr(self.record.taxonomy, rank).taxon.name
            
            taxonomy[rank] = name or ''

        return taxonomy
    
    def get_record_id(self):
        return int(self.record.record_id)
    
    def get_process_id(self):
        return str(self.record.processid)
    
    def get_tracefiles(self):
        if not hasattr(self.record, 'tracefiles'):
            return None
        
        return ( self.record.tracefiles.read[0], self.record.tracefiles.read[1] ) 
    
    def get_sequence(self):
        if self.sequence is None:
            self.sequence = pybold.sequence.SequencesClient().get(ids=self.get_process_id()).pop()
        
        return self.sequence
    
class SpecimensClient(Endpoint):
    '''
    classDocstring
    '''
    specimen_list = []
    endpoint_name = 'specimen'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(SpecimensClient, self).__init__()

    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None):
        result = super(SpecimensClient, self).get({'taxon': taxon, 
                                    'ids': ids, 
                                    'bin': bins, 
                                    'container': containers, 
                                    'institutions': institutions, 
                                    'researchers': researchers, 
                                    'geo': geo, 
                                    'format': 'xml'})
    
        bold_specimens = objectify.fromstring(result)
        #print bold_specimens.getchildren()
        for record in bold_specimens.record:
            self.specimen_list.append(Specimen(record))
    
        return self.specimen_list
    
    def get_taxonomies(self):
        taxonomies = {}
        for specimen in self.specimen_list:
            taxonomies[specimen.get_record_id()] = specimen.get_taxonomy()
        return taxonomies 
    
    def get_record_ids(self):
        ids = []
        for specimen in self.specimen_list:
            ids.append(specimen.get_record_id())    
        return ids
    
    def get_process_ids(self):
        ids = []
        for specimen in self.specimen_list:
            ids.append(specimen.get_process_id())
        
        return ids
    
    def get_sequences(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.sequence.SequencesClient(base_url = self.base_url).get(ids=ids_query)
    
if __name__ == "__main__":
    test = SpecimensClient()
    print test.url
    print test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.get_taxonomies()
    print test.get_record_ids()
    print test.get_sequences()
    print test.specimen_list.pop().get_tracefiles()