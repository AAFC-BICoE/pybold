PUBLIC_API_URL = 'http://www.boldsystems.org/index.php/API_Public/'
TAXON_API_URL = 'http://www.boldsystems.org/index.php/API_Tax/'

from exceptions import ValueError
import requests 
from urlparse import urljoin


class Endpoint(list):   
    endpoint_name = None
    base_url = None
    
    def __init__(self):
        super(Endpoint,self)

        self._set_base_url()

        
    def _set_base_url(self):
        if self.base_url is None or self.base_url == '':
            raise ValueError("BOLD API's base URL cannot be empty.")
        
        if not self.base_url.endswith('/'): self.base_url += '/'
        
        self.url = urljoin(self.base_url, self.endpoint_name)