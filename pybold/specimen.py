'''
Created on 2015-11-16

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca
'''
from __builtin__ import super
from booby import Model, fields
import json
from lxml import objectify, etree
import pickle
from xml import sax
from  xml.sax.handler import ContentHandler

from pybold import PUBLIC_API_URL, Endpoint


class Specimen(object):
    '''
    classdocs
    '''
    record_id = fields.Integer()
    processid = fields.String()

    def __init__(self, obj):
        '''
        Constructor
        '''
        super(Specimen, self).__init__()
        self = obj
    
    def get_taxonomy(self):
        print pickle.dumps(self)
#          if self.taxonomy.phylum is not None:
#              return self.taxonomy.phylum.taxon.name
         
class Specimens(Endpoint):
    '''
    classDocstring
    '''
    endpoint_name = 'specimen'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(Specimens, self).__init__()

    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None):
        result = super(Specimens, self).get({'taxon': taxon, 
                                    'ids': ids, 
                                    'bin': bins, 
                                    'container': containers, 
                                    'institutions': institutions, 
                                    'researchers': researchers, 
                                    'geo': geo, 
                                    'format': 'xml'})
    
        records = objectify.fromstring(result)
        for record in records.record:
            self.append(Specimen(record))
        
if __name__ == "__main__":
    test = Specimens()
    print test.url
    test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.pop().get_taxonomy()
    print 'here'