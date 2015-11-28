'''
Created on 2015-11-22

@author: kandalafti
'''


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import StringIO
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
        self.sequence_record = sequence_record
        self.specimen = None
        super(Sequence, self).__init__()
        
    def get_id(self):
        return self.sequence_record.id
    
    def get_seq(self):
        return self.sequence_record.seq
    
    def get_process_id(self):
        return self.get_id().split('|')[0]
    
    def get_identification(self):
        return self.get_id().split('|')[1]
    
    def get_marker(self):
        return self.get_id().split('|')[2]
    
    def get_accession(self):
        return self.get_id().split('|')[3]
    
    def get_specimen(self):
        if self.specimen is None:
            self.specimen = pybold.specimen.SpecimensClient().get(ids=self.get_process_id()).pop()
             
        return self.specimen
    
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
            ids.append(sequence.get_process_id())
        
        return ids
    
    def get_specimens(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.specimen.SpecimensClient(self.base_url).get(ids=ids_query)
    
if __name__ == "__main__":
    test = SequencesClient()
    #print test.url
    test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.sequence_list[0].get_id()
    print test.sequence_list.pop().get_seq()
    print test.sequence_list.pop().get_specimen().record
    