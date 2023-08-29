#!/usr/bin/env pytom
'''
Created on Jan 28, 2013

@author: thrabe
'''


if __name__ == '__main__':
    # parse command line arguments

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.agnostic.io import read_size
    helper = ScriptHelper(sys.argv[0].split('/')[-1], 
                          description='Create a localization job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/localization.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-v','--volume'], 'Volume : the big volume', arg=True, optional=False),
                                   ScriptOption(['-r','--reference'], 'Reference : the molecule searched', arg=True, optional=False),
                                   ScriptOption(['-m','--mask'], 'Mask : a mask ', arg=True, optional=False),
                                   ScriptOption(['--spherical_mask'], 'Set this option if the mask is spherical to '
                                                                      'speed up the calculation.',
                                                arg=False, optional=True),
                                   ScriptOption(['--wedge1'], 'Wedge : first tilt angle. Must be 90-tilt!', arg=True, optional=False),
                                   ScriptOption(['--wedge2'], 'Wedge : second tilt angle.  Must be 90-tilt!', arg=True, optional=False),
                                   ScriptOption(['-a','--angles'], '''Angles : name of angle list. Either : 
                                    angles_50_100.em
                                    angles_38.53_256.em
                                    angles_35.76_320.em
                                    angles_25.25_980.em
                                    angles_19.95_1944.em
                                    angles_18_3040.em    
                                    angles_12.85_7112.em    
                                    angles_11_15192.em    
                                    angles_07_45123.em
                                    angles_3_553680.em
                                    ''', arg=True, optional=False),
                                   ScriptOption(['-d','--destination'], 'Destination : destination directory', arg=True, optional=False),
                                   ScriptOption(['--zstart'], 'Start of z height subregion to search.', arg=True,
                                                optional=True),
                                   ScriptOption(['--zend'], 'Start of z height subregion to search.', arg=True,
                                                optional=True),
                                   ScriptOption(['-b','--band'], 'Lowpass filter : band - in pixels', arg=True,
                                                optional=True),
                                   ScriptOption(['-j','--jobName'], 'Specify job.xml filename', arg=True, optional=False),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])

    if len(sys.argv) <= 2:
        print(helper)
        sys.exit()
    try:
        volume, reference, mask, mask_is_spherical, wedge1,wedge2,angles,destination,zstart,zend,band,jobName,\
        help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
        
    if help is True:
        print(helper)
        sys.exit()

    if mask_is_spherical is None:
        mask_is_spherical = False
    
    if not checkFileExists(volume):
        raise RuntimeError('Volume file ' + volume + ' does not exist!')
    
    if not checkFileExists(reference):
        raise RuntimeError('Reference file ' + reference + ' does not exist!')
    
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')
    
    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')

    subregion = [0, 0, 0, 0, 0, 0]
    try:
        zstart, zend = int(zstart), int(zend)
        x, y, z = map(int, read_size(volume))
        if z >= zend > zstart:
            subregion = [0, 0, zstart, x, y, zend - zstart]
    except Exception as e:
        print('Non or invalid z search range provided.')
    
    from pytom.basic.structures import Mask,Reference,Wedge,BandPassFilter
    from pytom.localization.structures import Volume
    from pytom.angles.globalSampling import GlobalSampling
    from pytom.basic.score import FLCFScore
    from pytom.localization.peak_job import PeakJob
    
    v = Volume(volume, subregion=subregion)
    r = Reference(reference)
    m = Mask(mask)
    w = Wedge([float(wedge1),float(wedge2)])
    a = GlobalSampling(angles)
    bp = BandPassFilter(0,float(band),0) if band is not None else None
    
    job = PeakJob(volume=v, reference=r, mask=m, wedge=w, rotations=a, score=FLCFScore(), jobID=0, members=1,
                  dstDir=destination, bandpass=bp)
    
    job.toXMLFile(jobName)
