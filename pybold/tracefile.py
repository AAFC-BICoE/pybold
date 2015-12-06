'''
Created on 2015-12-02

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca>
@organization: Agriculture and Agri-Foods Canada
@group: Microbial Biodiversity Bioinformatics
@contact: mbb@agr.gc.ca 
'''
from StringIO import StringIO
import csv
import tarfile

from pybold import Endpoint, PUBLIC_API_URL


class Tracefile(object):
    '''
    classdocs
    '''


    def __init__(self, process_id, fileobj, marker, taxon, genbank_accession, filename):
        '''
        Constructor
        '''
        
        self.fileobj = fileobj
        self.process_id = process_id
        self.marker = marker
        self.taxon = taxon
        self.genbank_accession = genbank_accession
        self.filename = filename

        super(Tracefile, self).__init__()



class TracefilesClient(Endpoint):
    sequence_list = []
    ENDPOINT_NAME = 'trace'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        self.tracefile_list = []
        super(TracefilesClient, self).__init__()
    
    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, marker=None):
        '''
        @raise tarfile.ReadError: If the BOLD API return an invalid tarfile
        '''
         
        result = super(TracefilesClient, self).get({
                                    'taxon': taxon, 
                                    'ids': ids, 
                                    'bin': bins, 
                                    'container': containers, 
                                    'institutions': institutions, 
                                    'researchers': researchers, 
                                    'geo': geo,
                                    'marker': marker })
    
        chromat_tar = tarfile.open(mode='r:', fileobj=StringIO(result))
        
        # Build a dictionary containing the process_id, taxon, marker, genbank_accession, and filename 
        # based on the content of the TRACE_FILE_INFO.txt within the tar archive
        fileobj = chromat_tar.extractfile('TRACE_FILE_INFO.txt')
        tracefiles_d = {}
        for row in fileobj.readlines():
            row_traceinfo = row.rstrip().split("\t")
            if len(row_traceinfo) != 5:
                continue
            tracefiles_d[row_traceinfo[0]] = {'process_id': row_traceinfo[0],
                                                 'taxon': row_traceinfo[1],
                                                 'marker': row_traceinfo[2],
                                                 'genbank_accession': row_traceinfo[3],
                                                 'tracefile': row_traceinfo[4] 
                                                 }
        
        for member in chromat_tar.getmembers():
            if not member.isfile() or member.name == 'TRACE_FILE_INFO.txt':
                continue
            
            fileobj = chromat_tar.extractfile(member)
            
            for process_id in tracefiles_d.keys():
                if process_id not in member.name:
                    continue
                
                # Exit the loop because we found the process id in the tracefile name
                break
            
            if not process_id:
                raise ValueError()
            
            self.tracefile_list.append(Tracefile(process_id,
                                                  fileobj=fileobj, 
                                                  tracefiles_d[process_id]['marker'], 
                                                  tracefiles_d[process_id]['taxon'], 
                                                  tracefiles_d[process_id]['genbank_accession'],
                                                  member.name)
                                       )
            
            
        chromat_tar.close()
            
            
            
        
        #chromat_tar.list()
        
        
if __name__ == "__main__":
    test = TracefilesClient()
    test.get(ids='ACRJP618-11|ACRJP619-11')
    