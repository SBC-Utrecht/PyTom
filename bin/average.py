#!/usr/bin/env pytom

"""
Created on Aug 24, 2011

@author: hrabe
"""



if __name__ == '__main__':
# parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList
    from pytom.alignment.alignmentFunctions import distributeAverage
    from pytom_mpi import isInitialised
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Averge a particle. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/average.html',
                          authors='Thomas Hrabe',
                          options= [ScriptOption(['-p','--particleList'], 'XML particle list to be averaged.', arg=True, optional=False),
                                    ScriptOption(['-a','--average'], 'Resulting average filename.', arg=True, optional=False),
                                    ScriptOption(['--startIndex'], 'Average a subset of particles: start particle index.', arg=True, optional=True),
                                    ScriptOption(['--endIndex'], 'Average a subset of particles: last particle index.', arg=True, optional=True),
                                    ScriptOption(['--minimumScore'], 'Average a subset of particles: minimum score of particle. Does not work when start & end index are specified!', arg=True, optional=True),
                                    ScriptOption(['--maximumScore'], 'Average a subset of particles: maximum score of particle. Does not work when start & end index are specified!', arg=True, optional=True),
                                    ScriptOption(['--infoVolumes'], 'Generate info volumes like wedge volume.', arg=False, optional=True),
                                    ScriptOption(['--progressbarOff'], 'Display a progressbar. On by default', arg=False, optional=True),
                                    ScriptOption(['--fromAlignmentList'], 'Average from alignment list XML instead from particleListXML. Optional, off by default.', arg=True, optional=True),
                                    ScriptOption(['--subregion'], 'Average of particles whose coordinates are in [startX,startY,startZ,endX,endY,endZ].', arg=True, optional=True),
                                    ScriptOption(['--help'], 'Print this help.', arg=False, optional=True)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    
    try:
        particleListName, averageName, startIndex, endIndex, minimum, maximum, infoVolumes, progressbarOff, alignmentFileName ,subRegion, help= parse_script_options(sys.argv[1:], helper)
        
    except Exception as e:
        print e
        sys.exit()
        
    if help is True:
        print helper
        sys.exit()
        
    pl = None
    
    if particleListName:
        pl = ParticleList('/')
        pl.fromXMLFile(particleListName)
    
    if startIndex is not None and endIndex is not None:
        try:
            pl = pl[int(startIndex):int(endIndex)]
        except ValueError:
            print 'Your start or end index is not an integer: ', startIndex, endIndex
            print 'Aborting...'
            sys.exit()
    try:
        if subRegion and '[' in subRegion and ']' in subRegion:
            
            subRegion = [int(x) for x in subRegion.replace('[','').replace(']','').split(',')]
            
            newPL = ParticleList()
            
            for aParticle in pl:
                position = aParticle.getPickPosition()
                
                xOK = subRegion[0] <= position.getX() <= subRegion[3]
                yOK = subRegion[1] <= position.getY() <= subRegion[4]
                zOK = subRegion[2] <= position.getZ() <= subRegion[5]
                
                if xOK and yOK and zOK:
                    newPl.append(aParticle)
            
            if len(newPl) > 0:
                pl = newPl
    except:
        print 'Your subregion parameters seem to be incorrect. Use --help to look at their specification.'
        sys.exit()
        
    if startIndex is None and endIndex is None and minimum is not None and maximum is not None:
        newPl = ParticleList() 
    
        for aParticle in pl:
            if float(minimum) <= aParticle.getScore().getValue() <= float(maximum):
                newPl.append(aParticle)
                
        if len(newPl) > 0:
            pl = newPl
    
    if progressbarOff:
        progressbarOff = False
    else:
        progressbarOff = True
        
    if not particleListName and alignmentFileName:
        pl = ParticleList('/')
        pl.fromAlignmentList(alignmentFileName)
        
    if infoVolumes:
        infoVolumes = True
    else:
        infoVolumes = False
        
    distributeAverage(pl, averageName, progressbarOff, False, infoVolumes,sendEndMessage = True)
    
    