'''
Created on 2015-12-02

@author: Iyad Kandalaft <iyad.kandalaft@agr.gc.ca>
@organization: Agriculture and Agri-Foods Canada
@group: Microbial Biodiversity Bioinformatics
@contact: mbb@agr.gc.ca 
'''
import StringIO
import os.path
import tarfile

from pybold import Endpoint, PUBLIC_API_URL
import pybold.sequence
import pybold.specimen


class Tracefile(object):
    '''
    classdocs
    '''

    def __init__(self, process_id, fileobj, marker, taxon, genbank_accession, filename):
        '''
        Constructor
        '''
        self.specimen = None
        self.sequence = None
        
        self.fileobj = fileobj
        self.process_id = process_id
        self.marker = marker
        self.taxon = taxon
        self.genbank_accession = genbank_accession
        self.filename = filename

        super(Tracefile, self).__init__()

    @property
    def format(self):
        return os.path.splitext(self.filename)[1].lstrip('.')

    @property
    def specimen(self):
        if self.__specimen is None:
            self.specimen = pybold.specimen.SpecimensClient().get(ids=self.process_id).pop()
             
        return self.__specimen
    
    @specimen.setter
    def specimen(self, specimen_obj):
        self.__specimen = specimen_obj

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

    def to_file(self, dir_path=None, filename=None):
        if dir_path is None:
            dir_path = os.path.curdir
        else:
            os.makedirs(dir_path, 0755)
            
        if filename is None:
            filename = self.filename

        original_pos = self.fileobj.tell()
        self.fileobj.seek(0)
        
        full_path = os.path.join(dir_path, filename)
        with open(full_path, 'w+') as fh:
            fh.write(self.fileobj.read())
        
        self.fileobj.seek(original_pos)
        
        return full_path
     
class TracefilesClient(Endpoint):
    ENDPOINT_NAME = 'trace'
    
    def __init__(self, base_url=PUBLIC_API_URL):
        self.base_url = base_url
        self.tracefile_list = []
        super(TracefilesClient, self).__init__()
    
    def get(self, taxon=None, ids=None, bins=None, containers=None, institutions=None, researchers=None, geo=None, marker=None, timeout=5):
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
                                    'marker': marker }, timeout=timeout)

        self._parse_tracefiles(result)
        
        return self.tracefile_list
        
    def _parse_tracefiles(self, response):
        chromat_tar = tarfile.open(mode='r:', fileobj=StringIO.StringIO(response))
        
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
                                                fileobj, 
                                                tracefiles_d[process_id]['marker'], 
                                                tracefiles_d[process_id]['taxon'], 
                                                tracefiles_d[process_id]['genbank_accession'],
                                                member.name
                                                )
                                       )
            
        chromat_tar.close()

    def get_process_ids(self):
        ids = []
        for tracefile in self.tracefile_list:
            ids.append(tracefile.process_id)
        
        return ids

    def get_sequences(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.sequence.SequencesClient(self.base_url).get(ids=ids_query)
        
    def get_specimens(self):
        ids_query = '|'.join(self.get_process_ids())
        return pybold.specimen.SpecimensClient(self.base_url).get(ids=ids_query)
        

if __name__ == "__main__":
    test = TracefilesClient()
    test.get(ids='ACRJP618-11|ACRJP619-11')
    print test.tracefile_list[0].fileobj.readline()
    print test.tracefile_list[0].format
    print test.tracefile_list[0].process_id
    print test.tracefile_list[0].marker
    print test.tracefile_list[0].taxon
    print test.tracefile_list[0].genbank_accession
    print test.tracefile_list[0].filename
    test.tracefile_list[0].to_file()
    print test.tracefile_list[0].sequence
    print test.tracefile_list[0].specimen
    print test.get_sequences()
    print test.get_specimens()
    