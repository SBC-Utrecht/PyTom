

def create_RandomParticleList( reffile, pl_filename='pl.xml', pdir='./testparticles', nparticles=10):
    """
    @param reffile: reference file
    @type reffile: C{str}
    @param nparticles: number of particles (default: 10)
    @type nparticles: C{int}
    @param pl_filename: particle list filename
    @type pl_filename: C{str}
    @param pdir: particle directory
    @type pdir: C{str}
    @return: particleList
    @rtype: L{pytom.basic.ParticleList}
    
    """
    from pytom.basic.structures import Particle, ParticleList, Rotation, Shift
    from pytom_volume import vol, rotate, shift, read
    from pytom.basic.transformations import general_transform_crop
    from pytom.simulation.whiteNoise import add as addNoise
    import random
    from os import mkdir
    from pytom.score.score import FLCFScore as score

    try:
        mkdir(pdir)
    except FileExistsError:
        print('directory '+pdir+' existed already - using this one')
    random.seed(0)

    pl = ParticleList( directory='./')
    ref = read(reffile)
    for ii in range(0, nparticles):
        rot = Rotation(random.uniform(0,360), random.uniform(0,360), 
                       random.uniform(0,180))
        shift = Shift( x=random.uniform(-5,5), y=random.uniform(-5,5), 
                       z=random.uniform(-5,5))
        rotvol = general_transform_crop(v=ref, rot=rot, shift=shift, 
                 scale=None, order=[0,1,2])
        # add some noise
        noisy = addNoise(volume=rotvol,SNR=1)
        fname = pdir + '/particle_' + str(ii) + '.em'
        noisy.write( fname)
        p = Particle(filename=fname, rotation=rot, shift=shift, wedge=None, 
                     className=0, pickPosition=None, score=score, sourceInfo=None)
        p.setScoreValue(0.0)
        pl.append( particle=p)

    pl.setFileName( filename=pl_filename) 
    pl.toXMLFile(filename=pl_filename)
    return pl

def cleanUp_RandomParticleList( pl_filename='pl.xml', pdir='./testparticles'):
    """
    remove directories 
    """
    from os import remove, rmdir
    from pytom.basic.structures import ParticleList

    pl = ParticleList()
    pl.fromXMLFile(filename=pl_filename)
    for part in pl:
        remove(part.getFilename())
    rmdir(pdir)
    remove(pl_filename)

def create_TiltSeries(data, tiltAngles, outputfolder='./'):
    from pytom.tompy.io import read, write
    from pytom.tompy.transform import rotate_axis
    import os

    if data.__class__ == str:
        data = read(data)

    if not os.path.exists(outputfolder): os.mkdir(outputfolder)

    for n, tiltAngle in enumerate(tiltAngles):
        outname = os.path.join(outputfolder, 'sorted_{:03d}.mrc'.format(n))
        write(outname, rotate_axis(data, tiltAngle, axis='y').sum(axis=2),tilt_angle=tiltAngle)