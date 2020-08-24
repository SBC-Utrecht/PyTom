#!/usr/bin/env pytom

from pytom.tompy.mpi import MPI

mpi = None

def initialise_MPI():
    has_begun = False
    try:
        has_begun = mpi._begun
    except:
        pass

    if not has_begun:
        global mpi
        mpi = MPI()
        mpi.begin()