from io import BytesIO
import zipfile
import tarfile
from urllib.request import urlopen
from collections import namedtuple
import os

tg = os.path.join(os.path.realpath(__file__))

Dataset = namedtuple('Dataset', ['name', 'instances', 'solutions'])

datasets = \
    [
        Dataset(name="qaplib", 
                instances="http://anjos.mgi.polymtl.ca/qaplib//data.d/qapdata.tar.gz", 
                solutions="http://anjos.mgi.polymtl.ca/qaplib//soln.d/qapsoln.tar.gz")
    ]

def process(url, to):
    resp = urlopen(url)
    bt = BytesIO(resp.read())
    f = None
    if ".tar." in dataset.instances or \
        dataset.instances.endswith(".tar"):
        # Data specified is a tarfile.
        f = tarfile.open(fileobj=bt, mode='r:*')
        f.extractall(to)
    elif dataset.instances.endswith(".zip"):
        # Data specified is a zipfile
        f = zipfile.ZipFile(bt)
        f.extractall(to)
    else:
        # Assume plain file.
        with open(to, 'w') as f:
            f.write(bt)

for dataset in datasets:
    c = os.path.dirname(os.path.realpath(__file__))
    path_instances = os.path.join(c, dataset.name, "instances")
    path_solutions = os.path.join(c, dataset.name, "solutions")
    process(dataset.instances, path_instances)
    process(dataset.solutions, path_solutions)
    
