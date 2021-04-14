##!/bin/bash

## This script merges desired classes from one or multiple particleLists of individual tomograms. It utilizes the pytom function pl.splitByClass
## Utrecht University, 19.11.2019, RE

#module load openmpi/2.1.1 pytom/0.971

from pytom.basic.structures import ParticleList
import os

## The variable nice_classes_list contains arrays specifying the tomogram and classes. Each array should contain the index of the desired tomograms in the first position, followed by the indices of the desired classes from that tomogram.
## In the example below, the first array [8,0,1,3,4] specifies the classes 0, 1, 3 and 4 within tomogram 8.

nice_classes_list = [[8,0,1,3,4], [9,0], [13,2], [14,0,1,2], [15,0,1,3], [16,0,1], [20,1], [21,3], [24,1], [27,2], [29,1], [30,1], [31,0], [32,1], [33,1], [34,1], [35,1], [36,1,4], [38,1], [43,3], [45,1], [46,1], [52,4], [53,1,3], [54,1], [79,3]]

pl=ParticleList()

## Below, adjust the path inside pl.fromXMLFile to point to the individual particle lists.

for nice_classes in nice_classes_list:
  pl.fromXMLFile('tomogram' + str(nice_classes_list[nice_classes][0]) + '/subtomograms/subtomos_bin6/AC3D/classified_pl_iter9.xml')
  lists=pl.splitByClass()


  for i in range(len(nice_classes_list[nice_classes])-1)
    ribos=ribos.append(lists[nice_classes_list[nice_classes][i + 1])

## specify the desired output path below
ribos.toXMLFile('AC3d_pickedClasses.xml')
