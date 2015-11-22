'''
Created on 2015-11-16

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca
'''
from __builtin__ import super, str
import __builtin__
import json
from lxml import objectify
import pickle
import xmltodict

from pybold import PUBLIC_API_URL, Endpoint


class Specimen(objectify.ObjectifiedElement):
    '''
    classdocs
    '''
    
    def __new__(cls, objectified_element=None, specimen_record_xml=None):
        if isinstance(specimen_record_xml, __builtin__.str):
            return super(Specimen, cls).__new__(cls, specimen_record_xml = objectify.fromstring(specimen_record_xml))
        elif isinstance(objectified_element, objectify.ObjectifiedElement):
            return super(Specimen, cls).__new__(cls, objectified_element = objectified_element)
    

    def __init__(self, objectified_element=None, specimen_record_xml=None):
        if isinstance(specimen_record_xml, __builtin__.str):
            obj = objectify.fromstring(specimen_record_xml)
        elif isinstance(objectified_element, objectify.ObjectifiedElement):
            obj = objectified_element
        
        #print obj.getchildren()
        super(Specimen, self).__init__(obj)

    def get_taxonomy(self):
        #if self.taxonomy.phylum is not None:
         #   return self.taxonomy.phylum.taxon.name
         print "taxonomy"
     
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
    
        bold_specimens = objectify.fromstring(result)#.bold_records
        #print bold_specimens.getchildren()
        for specimen in bold_specimens.record:
            self.append(Specimen(specimen))
        
if __name__ == "__main__":
    test = Specimens()
    print test.url
    test.get(ids='ACRJP618-11|ACRJP619-11')

    #test.get(ids='ACRJP618-11')
    print test.pop().getchildren()
    #test.pop().get_taxonomy()
    #print 'here'