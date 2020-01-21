#!/usr/bin/env pytom
'''
Created on Jul 19, 2011

@author: hrabe
'''

    

if __name__ == '__main__':
    
    #this script is significantly linked to 
    #pytom/frontend/serverpages/createReconstructionJob.py
    #any changes like parameter names must be changed in the other script, too!

    import sys,getopt
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.basic.files import read_em_header
    from pytom.bin.createXMLsClosestMarker import *
    from pytom.gui.guiFunctions import datatypeAR

    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                      description='Reconstruct particles in a particle list. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/resonstructTomograms.html or\n\
                          http://www.pytom.org/doc/pytom/resonstructSubtomograms.html',
                      authors='Thomas Hrabe, FF',

                      options= [ScriptOption(['-t','--tomogram'], 'Reconstruct a tomogram. Specify name of tomogam here. You do not need a particle list for that!', arg=True, optional=True),
                                ScriptOption(['-p','--particleList'], 'XML particle list.', arg=True, optional=True),
                                ScriptOption(['-r','--prefixFileName'], 'Prefix of files that are selected.', arg=True, optional=True),
                                ScriptOption(['-l', '--logfileReconstruction'], 'Log files reconstruction, containing the positions of the markers.', arg=True, optional=True),
                                ScriptOption(['--projectionList'], 'XML projection list.', arg=True, optional=True),
                                ScriptOption(['--projectionDirectory'], 'Directory containing the projections.', arg=True, optional=True),
                                ScriptOption(['-w','--applyWeighting'], 'If projections are not weighted, apply weighting before. If omited, no weighting.', arg=True, optional=True),
                                ScriptOption(['-s','--size'], 'Size of particle cube / tomogram.', arg=True, optional=False),
                                ScriptOption(['-b','--coordinateBinning'], 'Binning factor of coordinates. If particle coordinates are determined in binned volume (with respect to projections) this binning factor needs to be specified.', arg=True, optional=True),
                                ScriptOption(['-o','--recOffset'], 'Cropping offset of the binned tomogram.', arg=True, optional=False),
                                ScriptOption(['--projBinning'], 'Bin projections BEFORE reconstruction. 1 is no binning, 2 will merge two voxels to one, 3 -> 1, 4 ->1 ...', arg=True, optional=True),
                                ScriptOption(['-m', '--metafile'], 'Supply a metafile to get tiltangles.', arg=True, optional=True),
                                ScriptOption(['-n', '--numProcesses'], 'number of parallel processes.', arg=True, optional=True),
                                ScriptOption(['--help'], 'Print this help.', arg=False, optional=True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
        
    particleList = None
    projectionDirectory = None
    aw = False
    
    try:
        tomogram, particleListXMLPath, prefix, logfile, projectionList, projectionDirectoryTemplate, aw, size, coordinateBinning, \
        recOffset, projBinning, metafile, numProcesses, help= parse_script_options(sys.argv[1:], helper)
    
    except Exception as e:
        print(e)
        sys.exit()
   
    if help:
        print(helper)
        sys.exit()



    size = [int(i) for i in size.split(",")]
    if len(size) == 1:
        tmp = size[0]
        size.append(tmp)
        size.append(tmp)

    if projBinning:
        projBinning = int(projBinning)
    else:
        projBinning = 1
        
    if coordinateBinning:
        coordinateBinning = float(coordinateBinning)
    else:
        coordinateBinning = 1
    if aw is None:
        aw = 0
    else:
        aw = int(aw)

    if recOffset:
        recOffset = [int(i) for i in recOffset.split(",")]
    else:
        recOffset = [0.,0.,0.]
    
    try:
        numProcesses = int(numProcesses)
    except:
        numProcesses = 0

    print(particleListXMLPath)
    xmlsBasedOnClosestMarker = extractParticleListsClosestToRefMarker(particleListXMLPath, logfile,
                                                                      binning_factor=coordinateBinning,
                                                                      projDirTemplate=projectionDirectoryTemplate)

    if not prefix:
        prefix = 'sorted_aligned'


    if tomogram:
        vol = projections.reconstructVolume( dims=size, reconstructionPosition=recOffset,
            binning=projBinning, applyWeighting=aw)
        vol.write(tomogram)
        
    else:
        print(xmlsBasedOnClosestMarker)
        for markerIndex, particleListXMLPath in xmlsBasedOnClosestMarker:
            markerIndex = int(markerIndex)
            projections = ProjectionList()
            projectionDirectory = projectionDirectoryTemplate.replace('_CLOSEST_', '_{:04d}_'.format(markerIndex))
            alignResultFile = os.path.join(projectionDirectory, 'alignmentResults.txt')

            if not os.path.exists(alignResultFile):
                alignResultFile = ''
            else:
                alignmentResults = loadstar(alignResultFile, dtype=datatypeAR)
                projectionsFileNames = alignmentResults['FileName']
                projectionDirectory = os.path.dirname(projectionsFileNames[0])
                print(projectionDirectory,checkDirExists(projectionDirectory) )
                prefix = os.path.basename(projectionsFileNames[0]).split('_')[0] + '_'
            if checkFileExists(projectionList):
                projections.fromXMLFile(projectionList)
            elif checkDirExists(projectionDirectory):
                projections.loadDirectory(projectionDirectory, metafile=metafile, prefix=prefix)
            else:
                raise RuntimeError('Neither projectionList existed nor the projectionDirectory you specified! Abort')

            # transform the cropping offset
            if len(projections) == 0:
                print(markerIndex, projectionDirectory, metafile, prefix)
                continue

            tmp = projections[0]
            sx = tmp.getXSize() # here should be the size of original projection!
            sy = tmp.getYSize()
            recOffset2 = [0,0,0]
            recOffset2[0] = -sx/2 + recOffset[0]*coordinateBinning
            recOffset2[1] = -sy/2 + recOffset[1]*coordinateBinning
            recOffset2[2] = -sx/2 + recOffset[2]*coordinateBinning

            # set particle list in order to reconstruct subtomograms
            particleList = ParticleList()

            try:
                particleList.fromXMLFile(particleListXMLPath)
            except RuntimeError:
                print('Error reading particleList XML file! Abort')
                sys.exit()

            from pytom.basic.structures import PickPosition
            for particle in particleList:
                pickPosition = particle.getPickPosition()
                x = (pickPosition.getX() * coordinateBinning + recOffset2[0])
                y = (pickPosition.getY() * coordinateBinning + recOffset2[1])
                z = (pickPosition.getZ() * coordinateBinning + recOffset2[2])

                particle.setPickPosition( PickPosition(x=x, y=y, z=z))

            projections.reconstructVolumes(particles=particleList, cubeSize=int(size[0]), \
                                           binning=projBinning, applyWeighting = aw, \
                                           showProgressBar = True,verbose=False, \
                                           preScale=projBinning,postScale=1, num_procs=numProcesses, alignResultFile=alignResultFile)



            del projections
