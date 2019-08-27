#!/usr/bin/env pytom

'''
Created on Jun 9, 2010

@author: chen
'''

def usage():
    print("Usage: ./scriptname -j jobFilename -r resultFilename -o orientationFilename -n maxNumOfParticles [-v minScoreValue] -s sizeOfParticleInRadius -p particleList -m  [-w (write to disk, length along each dimension)] [-g margin]")


def extractCandidates(jobFilename='', resultFilename='', orientFilename='', sizeParticle=None, maxNumParticle=0,
                      minScore=-1, write2disk=0, margin=None,mask=None):
    # construct the original job from the xml file
    if jobFilename=='':
        jobFilename='JobInfo.xml'
        
    from pytom.localization.peak_job import PeakJob
    job = PeakJob()
    job.fromXMLFile(jobFilename)
    
    from pytom.localization.extraction_peak_result import ExPeakResult
    res = ExPeakResult()
    res.volFilename = job.volume.getFilename()


    from pytom.basic.files import read
    from pytom_numpy import vol2npy
    from copy import deepcopy


    if resultFilename=='':
        res.resultFilename='scores.em'
    else:
        res.resultFilename = resultFilename
    if orientFilename=='':
        res.orientFilename='angles.em'
    else:
        res.orientFilename = orientFilename

    if mask:
        resultfile, maskfile = read(res.resultFilename), read(mask)
        resultdata, maskdata = deepcopy(vol2npy(resultfile)), deepcopy(vol2npy(maskfile))

        import mrcfile
        masked_data = (resultdata*maskdata)
        #masked_data[masked_data < 0.00001] = resultdata.median()

        mrcfile.new(res.resultFilename.replace('.em', '_masked.em'),masked_data.T, overwrite=True)
        res.resultFilename = res.resultFilename.replace('.em', '_masked.em')


    res.angleList = job.rotations[:]
    res.score = job.score.__class__
    
    if maxNumParticle <= 0:
        return None
    
    res.readAll()
    
    if sizeParticle==None:
        ref = job.reference.getVolume()
        sizeParticle = [ref.sizeX(),ref.sizeY(),ref.sizeZ()]
    
    particleList = res.findParticles(sizeParticle,maxNumParticle,minScore,write2disk,margin, offset=job.subregion[:3])
    
    return particleList

if __name__ == '__main__':
    import sys, os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
     
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Extract candidate molecules from localization result.',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-j','--jobFile'], 'Localization job XML file.', arg=True, optional=False),
                                    ScriptOption(['-q', '--mask'], 'Mask file (mask.em).',arg=True, optional=True),
                                    ScriptOption(['-r', '--result'], 'File with score coefficients (score.em).',arg=True, optional=False),
                                    ScriptOption(['-o','--orientation'], 'File with orientation indices (angles.em).', arg=True, optional=False),
                                    ScriptOption(['-n','--numberCandidates'], 'Number of candidates to extract.', arg=True, optional=False),
                                    ScriptOption(['-s','--size'], 'Radius around potential candidate that will be ignored during further processing.', arg=True, optional=False),
                                    ScriptOption(['-p','--particleList'], 'Name of particle list XML file.', arg=True, optional=False),
                                    ScriptOption(['-t','--particlePath'], 'Path prepended to each particle.', arg=True, optional=True),
                                    ScriptOption(['-v','--minimalScoreValue'], 'Minimal score value to which to extract.', arg=True, optional=True),
                                    ScriptOption(['-m','--motlList'], 'Write a MOTL file with candidates. The name of the file will be an extension of the particle list with .em.', arg=True,optional=True),
                                    ScriptOption(['-g','--margin'], 'Size of outer margin that will be ignored for potential candidates.', arg=True,optional=True),
                                    ScriptOption(['-w','--sizeCubes'], 'If specified, it will cut out candidates from the original tomogram with the specified size.', arg=True,optional=True),
                                    ScriptOption(['--scale'], 'Scale coordinates by a factor. Set > 1 to adjust to larger volumes. Use 2 if the localization tomo was 1x binned.', arg=True, optional=True),
                                    ScriptOption(['--help'], 'Print this help.', arg=False,optional= True)])
    
    if len(sys.argv) ==1:
        print(helper)
        sys.exit()
    
    
    try:
        jobFilename, maskFile, resultFilename, orientFilename, maxNumParticle, sizeParticle, plFilename, particlePath, minScore, motlFilename, margin, write2disk, scale, help = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
        
    if minScore == None:
        minScore = -1
    else:
        minScore = float(minScore)
    if write2disk == None:
        write2disk = 0
    if margin.__class__ == str:
        margin = int(margin)
    if particlePath is None:
        particlePath = './'
    if not scale is None:
        scale = float(scale)
    else:
        scale = 1.0


    res=extractCandidates(jobFilename,resultFilename,orientFilename,int(sizeParticle),int(maxNumParticle),minScore,int(write2disk),margin, mask = maskFile)
    if not plFilename and not motlFilename:
        raise RuntimeError('You must specify at least a particle list or a motl file as result of this script!')
    
    if plFilename:
        from pytom.basic.structures import ParticleList
        pl = ParticleList()
        
        from pytom.localization.peak_job import PeakJob
        job = PeakJob()
        job.fromXMLFile(jobFilename)

        wedge = job.wedge

        if particlePath[-1] != os.sep:
            particlePath += os.sep
        
        for particle in res:
            newParticle = particle.toParticle()
            newParticle.setWedge(wedge)
            newParticle.setFilename(particlePath + newParticle.getFilename())

            if scale != 1.0:
                pi = newParticle.getPickPosition()
                pi.setX(pi.getX() * scale)
                pi.setY(pi.getY() * scale)
                pi.setZ(pi.getZ() * scale)
                newParticle.setPickPosition(pi)

            pl.append(newParticle)
        
        pl.toXMLFile(plFilename)
        
    
    if motlFilename:
        from pytom.basic.structures import ParticleList
        
        pl = ParticleList()
        for newParticle in res:

            if scale != 1.0:
                pi = newParticle.getPickPosition()
                pi.setX(pi.getX() * scale)
                pi.setY(pi.getY() * scale)
                pi.setZ(pi.getZ() * scale)
                newParticle.setPickPosition(pi)

            pl.append(newParticle.toParticle())
            
        pl.toMOTL(motlFilename)
        
