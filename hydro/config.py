

from socket import gethostname

host_marble = {'cflags': ['-Wall', '-O3'], 'mpi4py_path': ''}
host_ria    = {'cflags': ['-fPIC', '-O3'], 'mpi4py_path':
                   '/home/jjz237/software/mpi4py/lib/python2.5/site-packages'}

try:
    host_config = {'marble.local': host_marble,
                   'PD161.PHYSICS.NYU.EDU': host_marble,
                   'ria.cosmo.fas.nyu.edu': host_ria
                   }[gethostname()]
except KeyError:
    print "Please add your machine, %s to the list of hosts in "\
        "hydro.config" % gethostname()
    exit()
