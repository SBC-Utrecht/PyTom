#!/usr/bin/env python 

# Convert PyTOM particle list from .xml to motl format.
# Run in Pytom environment.
# Tested in pytom/0.993
# Last edited 06-07-2021 by Sofie van Dorst

from pytom.basic.structures import ParticleList
import sys

print ('Usage is <programme.py> <xml_filename>.')

xml_file = sys.argv[1]
list_name = "".join(xml_file.split('.xml')[:-1]) 
motl_name = f"{list_name}.em"

pl = ParticleList()
pl.fromXMLFile(xml_file)
pl.toMOTL(motl_name)

print ("Your motl file has been created. Happy motling!")