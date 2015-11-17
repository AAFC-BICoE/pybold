'''
Created on 2015-11-16

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca
'''
from __builtin__ import super
from booby import Model, fields
from lxml import objectify

from pybold import PUBLIC_API_URL, Endpoint


class Specimen(object):
    '''
    classdocs
    '''
    record_id = fields.Integer()
    processid = fields.String()

    def __init__(self, xml):
        '''
        Constructor
        '''
        super(Specimen, self).__init__()
        self = objectify.parse(xml)
    
    
         
class Specimens(Endpoint):
    '''
    classDocstring
    '''
    endpoint_name = 'specimen'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        super(Specimens, self).__init__()

        
if __name__ == "__main__":
    test = Specimens()
    print test.url