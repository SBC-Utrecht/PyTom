
"""
script for cropping particles
"""

def writeCroppedParticles(particleListName, output, center, cubesize):
    """
    @param particleListName: name of particle list
    @type particleListName: str
    @param output: Name of output particles
    @type output: str
    @param center: center of output particles in template orientation
    @type center: list
    @param cubesize: Size of output particles in pixel
    @type cubesize: int

    """
    from pytom.basic.structures import ParticleList, Particle, Shift
    from pytom_volume import transformSpline as transform
    from pytom_volume import subvolume, vol


    pl = ParticleList()
    pl.fromXMLFile(filename=particleListName)
    #copy particle list for list of cropped particles
    pl_new = pl.copy()
    pvol = pl[0].getVolume()
    sizeX = pvol.sizeX() 
    sizeY = pvol.sizeY()
    sizeZ = pvol.sizeZ() 
    pvol_ali = vol(sizeX, sizeY, sizeZ) 
    subV = vol(cubesize, cubesize, cubesize)

    sub_startX = center[0]-cubesize/2
    sub_startY = center[1]-cubesize/2
    sub_startZ = center[2]-cubesize/2
    if (sub_startX < 0) or (sub_startY < 0) or (sub_startZ < 0):
        raise ValueError('cubesize too large :(')

    for (ipart, part) in enumerate(pl):
        pvol_ali.setAll(0) 
        subV.setAll(0)
        pvol = part.getVolume()
        rot = part.getRotation()
        rotinvert = rot.invert()
        shiftV = part.getShift() 
        transform(pvol, pvol_ali, rotinvert[0], rotinvert[1], rotinvert[2], 
                  sizeX/2, sizeY/2, sizeZ/2, -shiftV[0], -shiftV[1], -shiftV[2], 0, 0, 0) 
        # box out subvolume
        subV = subvolume(pvol_ali,  sub_startX, sub_startY, sub_startZ, cubesize, cubesize, cubesize)
        transform(subV, subV, rot[0], rot[1], rot[2], cubesize/2, cubesize/2, cubesize/2, 0, 0, 0, 0, 0, 0)
        fname = part.getFilename()
        idx = fname.split('_')[-1].split('.')[0] 
        nfname = output+'_'+idx+'.em'
        print "write file " + nfname
        subV.write(nfname)
        pl_new[ipart].setFilename(newFilename=nfname)
        pl_new[ipart].setShift(shift=Shift(0,0,0))
    return pl_new

                  
                  
if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options=[ScriptOption(['--particleList'], 'Name of particle list', arg=True, optional=False),
             ScriptOption(['--output'], 'Name of output particles (<output>_<index>.em)', arg=True, optional=False),
             ScriptOption(['--center'], 'Center of output particles in template orientation (3dim vec:x,y,z) starting at 0', 
                          arg=True, optional=True),
             ScriptOption(['--cubesize'], 'Size of output particles in pixel', arg=True, optional=False),
             ScriptOption(['--outParticleList'], 'Name of output particle list', arg=True, optional=True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Script for cropping particles',
                          authors='Friedrich Foerster, Stefan Pfeffer',
                          options=options)
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    particleListName, output, center, cubesize, outParticleListName = parse_script_options(sys.argv[1:], helper)

    if not center:
        center = [0,0,0]
    else:
        center = [int(center.split(',')[0]), int(center.split(',')[1]), int(center.split(',')[2])]
    if not output:
        raise ValueError("provide output filename")
    if not cubesize:
        raise ValueError("provide cubesize of new volumes")
    else:
        cubesize = int(cubesize)
    if not outParticleListName:
        outParticleListName='particleListCropped.xml'

    pl_new = writeCroppedParticles(particleListName=particleListName, 
                                   output=output, center=center, 
                                   cubesize=cubesize)
    pl_new.toXMLFile(outParticleListName)



