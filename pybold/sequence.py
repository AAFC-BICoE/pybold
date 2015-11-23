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


class Sequence(object):
    '''
    classdocs
    '''
    def __init__(self, sequence_record):
        '''
        Constructor
        '''
        self.sequence_record = sequence_record
        super(Sequence, self).__init__()
        
    def get_record_id(self):
        return self.sequence_record.id
    
    def get_sequence(self):
        return self.sequence_record.seq
    
    
class Sequences(Endpoint):
    '''
    Classdocs
    '''
    sequence_list = []
    endpoint_name = 'sequence'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(Sequences, self).__init__()
    
    
    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, marker=None):
        result = super(Sequences, self).get({'taxon': taxon, 
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
    
if __name__ == "__main__":
    test = Sequences()
    #print test.url
    test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.sequence_list.pop().get_sequence()
    #print test.sequence_list.pop().get_record_id()
    