'''
Created on 2015-11-22

@author: kandalafti
'''


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import StringIO
from exceptions import TypeError
import json
from lxml import objectify
import pickle

from pybold import Endpoint, PUBLIC_API_URL
import pybold.specimen


class Sequence(object):
    '''
    classdocs
    '''
    def __init__(self, sequence_record):
        '''
        Constructor
        '''
        if not isinstance(sequence_record, SeqRecord):
            raise TypeError()
        
        self.sequence_record = sequence_record
        self.specimen = None
        super(Sequence, self).__init__()
    
    @property
    def id(self):
        return self.sequence_record.id
    
    @property
    def seq(self):
        return self.sequence_record.seq
    
    @property
    def process_id(self):
        return self.id.split('|')[0]
    
    @property
    def identification(self):
        return self.id.split('|')[1]
    
    @property
    def marker(self):
        return self.id.split('|')[2]
    
    @property
    def accession(self):
        return self.id.split('|')[3]
    
    @property
    def specimen(self):
        if self.__specimen is None:
            self.specimen = pybold.specimen.SpecimensClient().get(ids=self.process_id).pop()
             
        return self.__specimen
    
    @specimen.setter
    def specimen(self, specimen_obj):
        self.__specimen = specimen_obj
    
class SequencesClient(Endpoint):
    '''
    Classdocs
    '''
    sequence_list = []
    endpoint_name = 'sequence'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(SequencesClient, self).__init__()
    
    
    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, marker=None):
        result = super(SequencesClient, self).get({'taxon': taxon, 
                                    'ids': ids, 
                                    'bin': bins, 
                                    'container': containers, 
                                    'institutions': institutions, 
                                    'researchers': researchers, 
                                    'geo': geo,
                                    'marker': marker })
    
        sequences_handle = StringIO.StringIO(result)
        
        for record in SeqIO.parse(sequences_handle, "fasta"):
            self.sequence_list.append(Sequence(record))
    
        return self.sequence_list
    
    def get_process_ids(self):
        ids = []
        for sequence in self.sequence_list:
            ids.append(sequence.process_id)
        
        return ids
    
    def get_specimens(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.specimen.SpecimensClient(self.base_url).get(ids=ids_query)
    
if __name__ == "__main__":
    test = SequencesClient()
    #print test.url
    test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.sequence_list[0].id
    print test.sequence_list[0].seq
    print test.sequence_list[0].specimen.taxonomy
    print test.sequence_list[1].specimen.taxonomy
    