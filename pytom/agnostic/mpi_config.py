#!/usr/bin/env pytom

from pytom.agnostic.mpi import MPI
global mpi
mpi = None

def initialise_MPI():
    has_begun = False
    try:
        has_begun = mpi._begun
    except:
        pass

    if not has_begun:
        mpi = MPI()
        mpi.begin()