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
import pybold.tracefile


class Specimen(object):
    '''
    classdocs
    '''

    def __init__(self, objectified_element=None, specimen_record_xml=None):
        if isinstance(specimen_record_xml, __builtin__.str):
            obj = objectify.fromstring(specimen_record_xml)
        elif isinstance(objectified_element, objectify.ObjectifiedElement):
            obj = objectified_element
        elif objectified_element is None and specimen_record_xml is None:
            raise ValueError()
        
        
        self.record = obj
        self.sequence = None
        self.tracefiles = None
        
        super(Specimen, self).__init__()
    
    @property
    def taxonomy(self):
        if not hasattr(self.record, 'taxonomy'):
            return None
        
        ranks = ('phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species')
        taxonomy = {}
        for rank in ranks:
            name = None
            if hasattr(self.record.taxonomy, rank):
                name = str(getattr(self.record.taxonomy, rank).taxon.name)
            
            taxonomy[rank] = name or ''

        return taxonomy
    
    @property
    def record_id(self):
        return int(self.record.record_id)
    
    @property
    def process_id(self):
        return str(self.record.processid)
    
    @property
    def geography(self):
        class Geography():
            def __init__(self):
                self.country = None
                self.province = None
                self.region = None
                self.coordinates = (None, None)
        
        geo = Geography()
        try:
            geo.country = self.record.collection_event.country
            geo.province = self.record.collection_event.province
            geo.region = self.record.collection_event.region
            geo.coordinates = (self.record.collection_event.coordinates.lat, self.record.collection_event.coordinates.long)
        except (KeyError, AttributeError):
            pass
        
        return geo
    
    @property
    def sequence(self):
        # If property is not set, then call the setter
        if self.__sequence is None:
            # Provide the setter with a Sequence object using an API query that fetches the sequence
            self.sequence = pybold.sequence.SequencesClient().get(ids=self.process_id).pop()
        
        return self.__sequence
    
    @sequence.setter
    def sequence(self, sequence_obj):        
        self.__sequence =  sequence_obj
    
    @property
    def tracefiles(self):
        if self._tracefiles is None:
            self.tracefiles = pybold.tracefile.TracefilesClient().get(ids=self.process_id)
            
        return self._tracefiles
    
    @tracefiles.setter
    def tracefiles(self, traces_list):
        self._tracefiles = traces_list
    
class SpecimensClient(Endpoint):
    '''
    classDocstring
    '''
    ENDPOINT_NAME = 'specimen'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(SpecimensClient, self).__init__()
        
        self.specimen_list = []


    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, timeout=5):
        result = super(SpecimensClient, self).get({'taxon': taxon, 
                                    'ids': ids, 
                                    'bin': bins, 
                                    'container': containers, 
                                    'institutions': institutions, 
                                    'researchers': researchers, 
                                    'geo': geo, 
                                    'format': 'xml'}, timeout=timeout)
    
        bold_specimens = objectify.fromstring(result)
        #print bold_specimens.getchildren()
        for record in bold_specimens.record:
            self.specimen_list.append(Specimen(record))
    
        return self.specimen_list
    
    def get_taxonomies(self):
        taxonomies = {}
        for specimen in self.specimen_list:
            taxonomies[specimen.process_id] = specimen.taxonomy
        return taxonomies 
    
    def get_record_ids(self):
        ids = []
        for specimen in self.specimen_list:
            ids.append(specimen.record_id)    
        return ids
    
    def get_process_ids(self):
        ids = []
        for specimen in self.specimen_list:
            ids.append(specimen.process_id)
        
        return ids
    
    def get_sequences(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.sequence.SequencesClient(self.base_url).get(ids=ids_query)
    
    def get_tracefiles(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.tracefile.TracefilesClient(self.base_url).get(ids=ids_query)
    
if __name__ == "__main__":
    test = SpecimensClient()
    print test.url
    print test.get(taxon='ACRJP618-11|ACRJP619-11')
    print test.get_taxonomies()
    print test.get_record_ids()
    print test.get_sequences()
    print test.get_sequences().pop().seq
    print test.specimen_list[0].taxonomy
    print test.specimen_list[1].taxonomy
    print test.specimen_list[1].tracefiles
    #print test.specimen_list.pop().tracefiles()