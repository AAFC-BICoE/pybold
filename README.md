# pybold

Python bindings for Barcode of Life's Restful services (http://boldsystems.org/index.php/resources/api).

These python modules provide object oriented methods to query various endpoints from the Barcode of Life web services.

# Install using pip

The latest release can be installed using pip as follows:
  pip install pybold

Alternatively, the latest modifications in the source can be installed from github using pip
  pip install -e git+https://github.com/aafc-mbb/pybold#egg=pybold

# Getting Started

Generally, you will start with one of the client classes to query the web services and return a set of objects that match the search criteria.

'''
import pybold.specimen

taxons_list = ['Achaeorhizomycetes', 'Arthoniomycetes']

sc = pybold.specimen.SpecimensClient()
specimen_list = sc.get(taxons='|'.join(taxons_list))

for specimen in specimen_list:
  print specimen.process_id
'''
