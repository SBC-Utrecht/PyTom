#!/usr/bin/env python

'''
Created on Apr 4, 2011

@author: yuxiangchen
'''

if __name__ == '__main__':
    usage = "Usage: ./scriptname -l weighted_projections_directory -b binning_factor -p particleList [-s cube_size in pixel] -f [center_offset_in_each_dimension]"
    
    import sys, getopt
    
    if len(sys.argv) ==1:
        print "No argument is given!"
        print usage
        sys.exit() 
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hl:p:o:s:b:f:", ["help"])
    except getopt.GetoptError:
        print 'Command not right. Exit!'
        sys.exit()
    
    for o,a in opts:
        if o in ("-h", "--help"):
            print usage
            sys.exit()
        if o in ("-l"):
            projDir = a
        if o in ("-p"):
            plFilename = a
        if o in ("-b"):
            binning = int(a)
        if o in ("-s"):
            size = int(a)
        if o in ("-f"):
            offset = [int(i) for i in a.split(",")]
    
    from pytom.basic.structures import ParticleList
    pl = ParticleList('.')
    pl.fromXMLFile(plFilename)
    
    pl.setCenterInPickVolume(-offset[0], -offset[1], -offset[2]) # thomas has changed the sign of it!!!
    
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    projs = ProjectionList()
    projs.loadDirectory(projDir)
    projs.generateVolumes(pl, size, binning, False, True, False)
