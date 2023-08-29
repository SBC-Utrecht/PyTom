import mrcfile
import copy


def writeAlignedProjections(TiltSeries_, weighting=None,
                            lowpassFilter=None, binning=None,verbose=False, write_images=True, order=(2,1,0)):
    """write weighted and aligned projections to disk1

       @param TiltSeries_: Tilt Series
       @type TiltSeries_: reconstruction.TiltSeries
       @param weighting: weighting (<0: analytical weighting, >1 exact weighting (value corresponds to object diameter in pixel AFTER binning)
       @type weighting: float
       @param lowpassFilter: lowpass filter (in Nyquist)
       @type lowpassFilter: float
       @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.

       @author: FF
    """
    import numpy as np
    from pytom.lib.pytom_numpy import vol2npy
    from pytom.basic.files import read_em, write_em
    from pytom.basic.functions import taper_edges
    from pytom.basic.transformations import general_transform2d
    from pytom.basic.fourier import ifft, fft
    from pytom.basic.filter import filter as filterFunction, bandpassFilter
    from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift, rotateFilter
    from pytom.lib.pytom_volume import complexRealMult, vol, pasteCenter, sum
    import pytom.lib.pytom_freqweight as pytom_freqweight
    from pytom.basic.transformations import resize
    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatypeAR
    import os
    from pytom.basic.files import EMHeader, read, read_em_header
    from pytom.agnostic.io import read_size
    if binning:
        imdim = int(float(TiltSeries_._imdim)/float(binning)+.5)
    else:
        imdim = TiltSeries_._imdim

    imdimX = int(np.around(read_size(TiltSeries_._ProjectionList[0]._filename, 'x')/binning,0))
    imdimY = int(np.around(read_size(TiltSeries_._ProjectionList[0]._filename, 'y')/binning,0))

    imdim = int(max(imdimX, imdimY))

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
        lpf = pytom_freqweight.weight(0.0,lowpassFilter*imdim/2, imdim, imdim//2+1,1, lowpassFilter/5.*imdim)
        # lpf2 = pytom_freqweight.weight(0.0, lowpassFilter * imdimY / 2, imdimX, imdimY // 2 + 1, 1,
        #                               lowpassFilter / 5. * imdimY)
        #lpf = bandpassFilter(volume=vol(imdim, imdim,1),lowestFrequency=0,highestFrequency=int(lowpassFilter*imdim/2),
        #                     bpf=None,smooth=lowpassFilter/5.*imdim,fourierOnly=False)[1]

    tilt_angles = []

    for projection in TiltSeries_._ProjectionList:
        tilt_angles.append( projection._tiltAngle )
    tilt_angles = sorted(tilt_angles)

    #q = np.matrix(abs(np.arange(-imdim//2, imdim//2)))

    alignmentResults = np.zeros((len(TiltSeries_._ProjectionList)),dtype=datatypeAR)
    alignmentResults['TiltAngle'] = tilt_angles

    for (ii,projection) in enumerate(TiltSeries_._ProjectionList):

        alignmentResults['FileName'][ii] = os.path.join(os.getcwd(), projection._filename)
        transX = -projection._alignmentTransX / binning
        transY = -projection._alignmentTransY / binning
        rot = -(projection._alignmentRotation + 90.)
        mag = projection._alignmentMagnification

        alignmentResults['AlignmentTransX'][ii] = transX
        alignmentResults['AlignmentTransY'][ii] = transY
        alignmentResults['InPlaneRotation'][ii] = rot
        alignmentResults['Magnification'][ii] = mag

        if write_images:
            if projection._filename.split('.')[-1] == 'st':
                from pytom.basic.files import EMHeader, read
                header = EMHeader()
                header.set_dim(x=imdim, y=imdim, z=1)
                idx = projection._index
                if verbose:
                    print("reading in projection %d" % idx)
                image = read(file=projection._filename, subregion=[0,0,idx-1,TiltSeries_._imdim,TiltSeries_._imdim,1],
                             sampling=[0,0,0], binning=[0,0,0])
                if not (binning == 1) or (binning == None):
                    image = resize(volume=image, factor=1/float(binning))[0]
            else:
                # read projection files
                from pytom.basic.files import EMHeader, read, read_em_header
                image = read(projection._filename)
                image = resize(volume=image, factor=1 / float(binning))[0]
                if projection._filename[-3:] == '.em':
                    header = read_em_header(projection._filename)
                else:
                    header = EMHeader()
                    header.set_dim(x=imdim, y=imdim, z=1)


            tiltAngle = projection._tiltAngle
            header.set_tiltangle(tiltAngle)

            # # 5 -- Optional Low Pass Filter
            # if lowpassFilter:
            #     filtered = filterFunction( volume=image, filterObject=lpf2, fourierOnly=False)
            #     image = filtered[0]

            # 1 -- normalize to contrast - subtract mean and norm to mean
            immean = sum(image) / image.numelem()
            image = (image - immean)/immean


            # 2 -- smoothen borders to prevent high contrast oscillations
            image = taper_edges(image, imdim//30)[0]
            # transform projection according to tilt alignment


            if projection._filename.split('.')[-1] == 'st':
                newFilename = (TiltSeries_._alignedTiltSeriesName+"_"+str(projection.getIndex())+'.em')
            else:
                newFilename = (TiltSeries_._alignedTiltSeriesName+"_"+str(projection.getIndex())
                               +'.'+TiltSeries_._tiltSeriesFormat)
            if verbose:
                tline = ("%30s" %newFilename)
                tline = tline + (" (tiltAngle=%6.2f)" % tiltAngle)
                tline = tline + (": transX=%6.1f" %transX)
                tline = tline + (", transY=%6.1f" %transY)
                tline = tline + (", rot=%6.2f" %rot)
                tline = tline + (", mag=%5.4f" %mag)
                print(tline)

            # 3 -- square if needed
            if imdimY != imdimX:
                print('Squared image to larger dimension')
                newImage = vol(imdim, imdim,1)
                newImage.setAll(0)
                pasteCenter(image, newImage)
                image = newImage

            # 4 -- Rotate

            image = general_transform2d(v=image, rot=rot, shift=[transX,transY], scale=mag, order=order, crop=True)

            # 5 -- Optional Low Pass Filter
            if lowpassFilter:
                filtered = filterFunction( volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            # 6 -- smoothen once more to avoid edges
            image = taper_edges(image, imdim//30)[0]

            # 7 -- weighting
            if (weighting != None) and (weighting < 0):
                image = (ifft( complexRealMult( fft( image), w_func) ) / (image.size_x()*image.size_y()*image.size_z()) )
            elif (weighting != None) and (weighting > 0):
                if weighting > 1.5:
                    w_func = fourierFilterShift(rotateFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
                else:
                    w_func = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))

                image = (ifft( complexRealMult( fft( image), w_func) )/
                      (image.size_x()*image.size_y()*image.size_z()) )
            header.set_tiltangle(tilt_angles[ii])

            if newFilename.endswith('.mrc'):
                data = copy.deepcopy(vol2npy(image))
                mrcfile.new(newFilename,data.T.astype(np.float32), overwrite=True)
            else:
                write_em(filename=newFilename, data=image, header=header)

            if verbose:
                tline = ("%30s written ..." %newFilename)

    outname = os.path.join(os.path.dirname(TiltSeries_._alignedTiltSeriesName), 'alignmentResults.txt')
    np.savetxt(outname, alignmentResults, fmt=fmtAR, header=headerAlignmentResults)
    print('Alignment successful. See {} for results.'.format(outname))
