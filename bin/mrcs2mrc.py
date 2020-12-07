#!/usr/bin/env pytom

import mrcfile
import sys, os
from pytom.tompy.io import read, write
from pytom.tompy.mpi import MPI

#mpi = MPI()

def extract_single_image(dataSlice, sliceId, out_name, tiltangle, origdir, outdir, prefix):
    if origdir:
        outname = os.path.join(outdir, out_name.replace('sorted_', prefix))
    else:
        outname = os.path.join(outdir, '{}{:02d}.mrc'.format(prefix, sliceId))
    print(f'extracted {os.path.basename(outname)} into {outdir}')
    # write(outname, data[:,:,sliceId], tilt_angle=tiltangle)
    mrcfile.new(outname, dataSlice.astype('float32'), overwrite=True)


if __name__=='__main__':
    
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList
    from numpy import arange

    options = [ScriptOption(['-f', '--fileName'], 'Filename of mrc stack.', True, False),
               ScriptOption(['-t', '--targetDir'], 'Folder in which the output files are saved', True, False),
               ScriptOption(['-p', '--prefix'], 'Prefix to filename. Default name of file.', True, False),
               ScriptOption(['-o', '--origDir'], 'Directory from which images in the stack originate. '
                                                'It will use the prefix "sorted" to select the output file names.',
                            True, True),
               ScriptOption(['-i', '--angularIncrement'], 'Angular Increment between tilts (typically 2 or 3 degrees). Default 2 degrees.', True, True),
               ScriptOption(['-s', '--startingAngle'], 'Lowest angle in tiltseries. Default -60 degrees.', True, True),
               ScriptOption(['-e', '--endingAngle'], 'Highest angle in tiltseries. Default 60 degrees .', True, True),
               ScriptOption(['--metaFile'], 'Pytom-specific meta file containing the tiltangles of the images.', True, True),
               ScriptOption(['-m', '--mdoc'], 'Create truncated mdoc files.', False, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name             
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)
    # mpi.begin()

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, outdir, prefix, origdir, increment, startAngle, endAngle, metafile, mdoc, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)

        print(helper)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()


    if mdoc is None: mdoc =False
    else: mdoc = True

    if os.path.exists(filename):
        data = mrcfile.open(filename, permissive=True).data.copy()
    else:
        print(helper)
        sys.exit()

    if prefix is None:
        prefix = os.path.splitext(os.path.basename(filename))[0]
    if not prefix[-1] == '_': prefix += '_'


    increment = 2 if increment is None else float(increment)
    startAngle = -60 if startAngle is None else float(startAngle)
    endAngle = 60 if endAngle is None else float(endAngle)
    angles = arange(startAngle, endAngle + increment/2, increment)


    if not metafile is None and os.path.exists(metafile):
        from pytom.basic.datatypes import DATATYPE_METAFILE
        from pytom.gui.guiFunctions import loadstar
        metadata = loadstar(metafile, dtype=DATATYPE_METAFILE)
        angles = metadata['TiltAngle']

    x,y,z = data.shape
    sliceNR = min(x,y,z)

    if len(angles) != sliceNR:
        raise Exception('Angles and number of images do not coincide.')

    if mdoc:
        mdocname = os.path.join(outdir, os.path.basename(filename).replace('.mrc', '.mdoc').replace('.ali','.mdoc').replace('.st','.mdoc'))
        print(mdocname)
        mdocfile = open(mdocname, 'w')

    d = '''[ZValue = {}]
TiltAngle = {}
SubFramePath = X:\{}

'''


    if origdir:
        out_names = sorted([fname for fname in os.listdir(origdir) if fname.startswith('sorted')])
        if len(out_names) != sliceNR:
            print(out_names)
            raise Exception('Number of files in orig dir is different from the number of defined images in stack.')

    else:
        out_names = [os.path.join(outdir, '{}{:02d}.mrc'.format(prefix, sliceId)) for sliceId in range(sliceNR)]

    out = []
    for sliceId in range(sliceNR):
        out.append((data[sliceId,:,:], sliceId, out_names[sliceId], angles[sliceId], origdir, outdir, prefix))
        #tiltangle = angles[sliceId]
        #if origdir:
        #    outname = os.path.join(outdir, out_names[sliceId].replace('sorted_', prefix))
        #else:
        #    outname = os.path.join(outdir, '{}{:02d}.mrc'.format(prefix, sliceId))
        #print(f'extracted {os.path.basename(outname)} into {outdir}')
        #write(outname, data[:,:,sliceId], tilt_angle=tiltangle)
        #mrcfile.new(outname, data[sliceId, :,:].astype('float32'), overwrite=True)
        if mdoc:
            mdocfile.write(d.format(sliceId, angles[sliceId], os.path.basename(out_names[sliceId])))

        extract_single_image(*out[-1])

    #mpi.parfor(extract_single_image, out)
    #mpi.end()

    if mdoc: mdocfile.close()
    
