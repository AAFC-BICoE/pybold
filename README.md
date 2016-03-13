# pybold

Python bindings for [Barcode of Life's Restful WebServices](http://boldsystems.org/index.php/resources/api).

These python modules provide object oriented methods to query various endpoints from the Barcode of Life web services.

# Install using pip

The latest release can be installed using pip as follows:
> pip install pybold

Alternatively, the latest modifications in the source can be installed from github using pip
> pip install -e git+https://github.com/aafc-mbb/pybold#egg=pybold

# Getting Started

Generally, you will start with one of the client classes to query the web services and return a set of objects that match the search criteria.

```
import pybold.specimen

taxons_list = ['Achaeorhizomycetes', 'Arthoniomycetes']

sc = pybold.specimen.SpecimensClient()
specimen_list = sc.get(taxons='|'.join(taxons_list))

for specimen in specimen_list:
  print specimen.process_id
```

Currently, there are 3 clients: SpecimensClient, SequencesClient, and TracefilesClient.  Each client has a _get_ method that will return a list of objects (Specimen, Sequence, Tracefile) corresponding to the search criteria.

As a convenience, the Specimen, Sequence, and Tracefile objects permit looking up the corresponding Specimen, Sequence, or Tracefile objects.

```
import pybold.tracefile

taxons_list = ['Achaeorhizomycetes', 'Arthoniomycetes']

tc = pybold.tracefile.TracefilesClient()
tracefile_list = tc.get(taxons='|'.join(taxons_list))

for tracefile in tracefile_list:
  print tracefile.sequence.seq
```
