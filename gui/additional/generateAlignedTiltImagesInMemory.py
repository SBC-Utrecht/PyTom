from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR
import glob
import mrcfile
import os
import numpy

def determine_tiltangles(image_list, metafile):
    basename_images = [os.path.basename(image) for image in image_list]
    identifier = [f.split('.')[-2].split('_')[-1] for f in basename_images]
    positions = map(int, identifier)
    metadata = numpy.loadtxt(metafile, dtype=datatype)
    tilt_angles = metadata['TiltAngle'][positions]
    return tilt_angles

def toProjectionStackFromAlignmentResultsFile( alignmentResultsFile, weighting=None, lowpassFilter=0.9, binning=None, num_procs=1):
    """read image and create aligned projection stack, based on the results described in the alignmentResultFile.

       @param alignmentResultsFile: result file generate by the alignment script.
       @type datatypeAR: gui.guiFunction.datatypeAR
       @param weighting: weighting (<0: analytical weighting, >1 exact weighting (value corresponds to object diameter in pixel AFTER binning)
       @type weighting: float
       @param lowpassFilter: lowpass filter (in Nyquist)
       @type lowpassFilter: float
       @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.

       @author: GS
    """
    import numpy
    from pytom_numpy import vol2npy
    from pytom.basic.files import read_em, write_em
    from pytom.basic.functions import taper_edges
    from pytom.basic.transformations import general_transform2d
    from pytom.basic.fourier import ifft, fft
    from pytom.basic.filter import filter as filterFunction, bandpassFilter
    from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift
    from pytom_volume import complexRealMult, vol, paste
    import pytom_freqweight
    from pytom.basic.transformations import resize
    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR
    from pytom.gui.additional.reconstructionStructures import Projection, ProjectionList
    from pytom_numpy import vol2npy

    alignmentResults = numpy.loadtxt(alignmentResultsFile, dtype=datatypeAR)
    imageList = alignmentResults['FileName']
    tilt_angles = alignmentResults['TiltAngle']
    print(tilt_angles)
    a = mrcfile.open(imageList[0],permissive=True)
    imdim = a.data.T.shape[0]

    if binning:
        imdim = int(float(imdim) / float(binning) + .5)
    else:
        imdim = imdim

    print('imdim', imdim)
    sliceWidth = imdim

    # pre-determine analytical weighting function and lowpass for speedup
    if (weighting != None) and (weighting < -0.001):
        w_func = fourierFilterShift(rampFilter(imdim, imdim))

    # design lowpass filter
    if lowpassFilter:
        if lowpassFilter > 1.:
            lowpassFilter = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
        # weighting filter: arguments: (angle, cutoff radius, dimx, dimy,
        lpf = pytom_freqweight.weight(0.0, lowpassFilter * imdim // 2, imdim, imdim // 2 + 1, 1,
                                      lowpassFilter / 5. * imdim)
        # lpf = bandpassFilter(volume=vol(imdim, imdim,1),lowestFrequency=0,highestFrequency=int(lowpassFilter*imdim/2),
        #                     bpf=None,smooth=lowpassFilter/5.*imdim,fourierOnly=False)[1]
        print(lpf)

    projectionList = ProjectionList()
    for n, image in enumerate(imageList):
        atx = alignmentResults['AlignmentTransX'][n]
        aty = alignmentResults['AlignmentTransY'][n]
        rot = alignmentResults['InPlaneRotation'][n]
        mag = alignmentResults['Magnification'][n]

        projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx,alignmentTransY=aty,
                                alignmentRotation=rot,alignmentMagnification=mag)
        projectionList.append(projection)

    stack = vol(imdim, imdim, len(imageList))
    stack.setAll(0.0)

    phiStack = vol(1, 1, len(imageList))
    phiStack.setAll(0.0)

    thetaStack = vol(1, 1, len(imageList))
    thetaStack.setAll(0.0)

    offsetStack = vol(1, 2, len(imageList))
    offsetStack.setAll(0.0)

    for (ii, projection) in enumerate(projectionList):
        if projection._filename.split('.')[-1] == 'st':
            from pytom.basic.files import EMHeader, read
            idx = projection._index
            image = read(file=projection._filename,
                         subregion=[0, 0, idx - 1, imdim, imdim, 1],
                         sampling=[0, 0, 0], binning=[0, 0, 0])
            if not (binning == 1) or (binning == None):
                image = resize(volume=image, factor=1 / float(binning))[0]
        else:
            # read projection files
            from pytom.basic.files import EMHeader, read, read_em_header
            print(type(str(projection._filename)))
            image = read(str(projection._filename))
            image = resize(volume=image, factor=1 / float(binning))[0]

        if lowpassFilter:
            filtered = filterFunction(volume=image, filterObject=lpf, fourierOnly=False)
            image = filtered[0]

        tiltAngle = projection._tiltAngle

        # normalize to contrast - subtract mean and norm to mean
        immean = vol2npy(image).mean()
        image = (image - immean) / immean

        # smoothen borders to prevent high contrast oscillations
        image = taper_edges(image, imdim // 30)[0]

        # transform projection according to tilt alignment
        transX = projection._alignmentTransX / binning
        transY = projection._alignmentTransY / binning
        rot = float(projection._alignmentRotation)
        mag = float(projection._alignmentMagnification)

        print(transX,transY,rot,mag)
        image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=[2, 1, 0], crop=True)

        # smoothen once more to avoid edges
        image = taper_edges(image, imdim // 30)[0]

        # analytical weighting
        if (weighting != None) and (weighting < 0):
            image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))

        elif (weighting != None) and (weighting > 0):
            w_func = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
            image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))

        thetaStack(int(round(projection.getTiltAngle())), 0, 0, ii)
        offsetStack(int(round(projection.getOffsetX())), 0, 0, ii)
        offsetStack(int(round(projection.getOffsetY())), 0, 1, ii)
        paste(image, stack, 0, 0, ii)
        fname = '/Users/gijs/Documents/PostDocUtrecht/Data/Juliette/sorted_aligned_novel_{:02d}.mrc'.format(ii)
        write_em(fname.replace('mrc','em'), image)
        mrcfile.new(fname, vol2npy(image).copy().astype('float32').T, overwrite=True)

    return [stack, phiStack, thetaStack, offsetStack]




if __name__ == '__main__':
    import sys
    alignmentResultFile = sys.argv[1]

    toProjectionStackFromAlignmentResultsFile(alignmentResultFile, binning=1, weighting=0, lowpassFilter=0.9)
