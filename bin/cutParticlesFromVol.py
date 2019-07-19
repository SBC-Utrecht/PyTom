#!/usr/bin/env pytom

'''
Created on Apr 4, 2011

@author: yuxiangchen
'''

if __name__ == '__main__':
    import sys
    
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom_volume import read, subvolume
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.basic.structures import ParticleList
    from pytom.tools.files import checkDirExists
    import os
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Cut particles out from a volume, given the particle list. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/cutParticles.html',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-v','--volume'], 'Volume.', arg=True, optional=False),
                                    ScriptOption(['-p','--particleFile'], 'Particle list.', arg=True, optional=False),
                                    ScriptOption(['-c','--cubeSize'], 'Cube size along each dimension.', arg=True, optional=False),
                                    ScriptOption(['--help'], 'Print this help.', arg=False,optional=True)])

    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        volFilename, plFilename, cubeSize, help = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if help is True:
        print helper
        sys.exit()

    cubeSize = int(cubeSize)

    particleList = ParticleList()
    try:
        particleList.fromXMLFile(plFilename)
    except:
        from pytom.localization.structures import readParticleFile
        particles = readParticleFile(plFilename)
        
        particleList = ParticleList()
        for particle in particles:
            particleList.append(particle.toParticle())
            
    particlePath = particleList[0].getFilename()
    particleFolder = particlePath[0:particlePath.rfind('/')]
    
    if not checkDirExists(particleFolder):
        os.makedirs(particleFolder)
    
    prog = FixedProgBar(0, len(particleList)-1, '')

    newParticleList = ParticleList()
    
    vol = read(volFilename)
    volX = vol.sizeX()
    volY = vol.sizeY()
    volZ = vol.sizeZ()

    i = 0
    for particle in particleList:
        i = i+1 # even if some particles are skipped, the progress shall be updated
        prog.update(i)
        
        pi = particle.getPickPosition()
        try:
            x = int(pi.getX()-cubeSize/2)
            y = int(pi.getY()-cubeSize/2)
            z = int(pi.getZ()-cubeSize/2)
            
            if x < cubeSize or y < cubeSize or z < cubeSize or\
               x+cubeSize > volX or y+cubeSize > volY or z+cubeSize > volZ:
                print 'Coordinate out of bounds (',x,y,z,') for '
                print particle
                print 'Particle could not be cut out from origin volume!'
                print ''
                continue
            
            v = subvolume(vol, x, y, z, cubeSize, cubeSize, cubeSize) # faster this way
            
            newParticleList.append(particle) # this part should be inside the try block
            v.write(particle.getFilename())
        except:
            print 'Error for'
            print particle
            print 'Particle could not be cut out from origin volume!'
            print ''
            continue


    if len(particleList) != len(newParticleList):
        new_plFilename = None
        if '/' in plFilename:
            new_plFilename = plFilename[:-4] + 'New.xml'
        else:
            new_plFilename = 'new_'+plFilename
        
        print 'Length of particle list has been changed. The new particle list is saved as', new_plFilename
        newParticleList.toXMLFile(new_plFilename)

