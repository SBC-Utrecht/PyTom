'''
started on Feb 02, 2013

@author: foerster
'''


from pytom.basic.structures import PyTomClass
from pytom_volume import vol

class Image(PyTomClass):
    """
    2D Image 
    """
    def __init__(self, filename=None, boxCoords=[0,0,0], dims=[32,32,1],
                 shiftX=0, shiftY=0, appliedShiftX=0, appliedShiftY=0, rotation=0,
                 index = 0, verbose=False):
        """
        @param filename:
        @type filename: L{str}
        @param boxCoords: coordinates in file (lower bounds)
        @type boxCoords: L{list} (2-dim or 3-dim)
        @param dims: dimensions of subframe
        @type dims: 2-dim or 3-dim list
        @param shiftX: shift in X (if determined previously)
        @type shiftX: float
        @param shiftY: shift in Y (if determined previously)
        @type shiftY: L{float}
        @param rotation: rotation in degrees (if determined previously)
        @type rotation: L{float}
        @param appliedShiftX: applied shift in X (if any shift applied to image)
        @type appliedShiftX: L{float}
        @param appliedShiftY: applied shift in Y (if any shift applied to image)
        @type appliedShiftY: L{float}
        @param index: index of image for other applications
        @type index: L{int}
        """
        self.verbose = verbose
        from pytom_volume import read, vol
        if filename:
            self.read(filename=filename, boxCoords=boxCoords, dims=dims)
        else:
            self.data=vol(dims[0],dims[1],dims[2])
        self.shiftX = shiftX
        self.shiftY = shiftY
        self.appliedShiftX = appliedShiftX
        self.appliedShiftY = appliedShiftY
        self.rotation = rotation
        self.boxCoords = boxCoords
        self.dims = dims
        self.index = 0

    def read(self, filename, boxCoords=[0, 0, 0], dims=[32, 32, 1]):
        """
        read 2D image or part of it from image or volume file

        @param filename: filename of image / vol
        @type filename: L{str}
        @param boxCoords: coordinates in file (lower bounds)
        @type boxCoords: 2-dim or 3-dim list
        @param dims: dimensions of subframe
        @type dims: 2-dim or 3-dim list
        """
        from pytom_volume import read

        self.data = read(filename, int(boxCoords[0]), int(boxCoords[1]), int(boxCoords[2]),
                         int(dims[0]), int(dims[1]), int(dims[2]), 0, 0, 0, 1, 1, 1)
        if self.verbose:
            print(
            "read file <" + filename + "> at position " + str(boxCoords[0]) + ", " + str(boxCoords[1]) + ", " + str(
                boxCoords[2]))
        self.boxCoords = boxCoords
        self.dims = dims

    def copy(self):
        """
        make a copy of an image
        @return: image copy
        @rtype: L{Image}
        """
        imageCopy = Image(filename=None, boxCoords=self.boxCoords, dims=self.dims,
                          shiftX=self.shiftX, shiftY=self.shiftY,
                          appliedShiftX=self.appliedShiftX,
                          appliedShiftY=self.appliedShiftY, rotation=self.rotation)
        imageCopy.data.copyVolume(self.data)

        return imageCopy

    def normalize(self, normtype="StdMeanInMask", mask=None, p=None):
        """
        normalize image

        @param normtype: normalization type
        @type normtype: str ("StdMeanInMask", "StdMean")
        @param mask: mask volume
        @type mask: L{ptom_volume.vol}
        @param p: sum of gray values in mask (if pre-computed)
        @type p: L{int} or L{float}
        """
        if normtype == "StdMeanInMask":
            from pytom.basic.normalise import normaliseUnderMask

            if mask == None:
                raise ValueError("StdMeanInMask normalization requires mask!")
            # spherical mask
            if (type(mask) == float) or (isinstance(mask, (int, long))):
                if mask <= 0:
                    raise ValueError("Value for mask radius must be > 0!")
                from pytom.basic.functions import initSphere

                mask = initSphere(sizeX=self.data.sizeX(), sizeY=self.data.sizeY(),
                                  sizeZ=self.data.sizeZ(), radius=mask,
                                  smooth=mask/10., maxradius=0, cent=None)
            # user-specified generic mask
            else:
                if type(mask) != vol:
                    raise TypeError("Mask must be pytom_volume.vol")
                if ((mask.sizeX() != self.data.sizeX()) or
                        (mask.sizeY() != self.data.sizeY()) or
                        (mask.sizeZ() != self.data.sizeZ())):
                    raise ValueError("Mask have same dimension as image")

            normvol, p = normaliseUnderMask(volume=self.data, mask=mask, p=p)
            self.data = normvol
            #return number of voxels in volume for speed-up
            return p

        if normtype == "StdMean":
            from pytom.basic.normalise import mean0std1
            mean0std1(self.data)
            # return p for consistency
            return None

    def taper_edges(self, width, taper_mask=None):
        """
        taper edges of image

        @param width: width of tapered edge
        @type width: L{float} or L{int}
        @param taper_mask: mask for tapering - if None it will be generated
        @type taper_mask: L{pytom_volume.vol}
        @return: taper_mask
        @rtype: L{pytom_volume.vol}
        """
        from pytom.basic.functions import taper_edges

        (self.data, taper_mask) = taper_edges(image=self.data, width=width, taper_mask=taper_mask)

        return taper_mask


    def bandpass(self, lowfreq, hifreq, smooth=0., bpf=None):
        """
        bandpass filter image

        @param lowfreq: lowest frequency of filter (=hi-pass)
        @type lowfreq: L{float}
        @param hifreq: highest frequency of filter (=low-pass)
        @type hifreq: L{float}
        @param smooth: smoothing
        @type smooth: L{float}
        @param bpf: bandpass filter object

        @return: bandpass filter object
        @rtype: L{pytom.basic.structures.BandPassFilter}
        """
        from pytom.basic.filter import bandpassFilter

        res = bandpassFilter(volume=self.data, lowestFrequency=lowfreq,
                             highestFrequency=hifreq, bpf=None, smooth=smooth, fourierOnly=False)
        self.data = res[0]
        bpf = res[1]
        return bpf

    def bandpass_in_nyquist(self, lowfreq, hifreq, smooth=0., bpf=None):
        """
        bandpass filter image with frequencies specified in Nyquist

        @param lowfreq: lowest frequency of filter in Nyquist (=hi-pass)
        @param hifreq: highest frequency of filter in Nyquist (=low-pass)
        @param smooth: smoothing in Nyquist
        @param bpf: bandpass filter object

        @return: bandpass filter object
        @rtype: L{pytom.basic.structures.BandPassFilter}
        """
        from pytom.basic.filter import bandpassFilter

        res = bandpassFilter(volume=self.data,
                             lowestFrequency=lowfreq * self.dims[0] / 2.,
                             highestFrequency=hifreq * self.dims[0] / 2., bpf=None,
                             smooth=smooth * self.dims[0] / 2., fourierOnly=False)
        self.data = res[0]
        bpf = res[1]
        return bpf


class ImageStack(PyTomClass):
    """
    Stack of 2D images
    """

    def __init__(self, verbose=False):
        """
        """
        self.images = []
        self.imageCopies = []
        # store average here
        self.averageData = None
        # number of voxels in mask (if used for processing)
        self.p = None
        self.dimX = 0
        self.dimY = 0
        self.verbose = verbose


    def addImageFromFile(self, filename, boxCoords=[0, 0, 0], dims=[32, 32, 1],
                         shiftX=0, shiftY=0, rotation=0, appliedShiftX=0, appliedShiftY=0,
                         index=0):
        """
        add image from file

        @param filename: filename of image / vol
        @type filename: L{str}
        @param boxCoords: coordinates of upper left corner of boxed frame
        @type boxCoords: 3-dim list
        @param dims: dimensions of frame
        @type dims: L{list} (3-dim)
        @param shiftX: shift in X (if determined previously)
        @type shiftX: L{float}
        @param shiftY: shift in Y (if determined previously)
        @type shiftY: L{float}
        @param rotation: rotation in degrees (if determined previously)
        @type rotation: L{float}
        @param appliedShiftX: applied shift in X (if any shift applied to image)
        @type appliedShiftX: L{float}
        @param appliedShiftY: applied shift in Y (if any shift applied to image)
        @type appliedShiftY: L{float}
        @param index: index of image
        @type index: L{int}
        """
        image = Image(filename=filename, boxCoords=boxCoords, dims=dims,
                      shiftX=shiftX, shiftY=shiftY, appliedShiftX=appliedShiftX,
                      appliedShiftY=appliedShiftY, rotation=rotation,
                      index=index, verbose=self.verbose)
        self.append(image)

    def append(self, image):
        """
        append 2D image to ImageStack

        @param image: 2D image
        @type image: L{ptom_volume.vol}
        """
        if (type(image) != Image):
            print type(image)
            raise TypeError("Input must be Image")
        else:
            self.images.append(image)
            # copy working copy
            self.imageCopies.append(image.copy())
        if self.dimX == 0:
            self.dimX = image.data.sizeX()
            self.dimY = image.data.sizeY()

    def write(self, filename):
        """
        write ImageStack into file

        @param filename: name of stack
        @type filename: L{str}
        """
        from pytom_volume import vol

        stack = vol(self.dimX, self.dimY, len(self.images))
        for ii in range(0, len(self.images)):
            for ix in range(0, self.dimX):
                for iy in range(0, self.dimY):
                    v = self.images[ii].data.getV(ix, iy, 0)
                    stack.setV(v, ix, iy, ii)
        stack.write(filename)

    def writeWorkingCopies(self, filename):
        """
        write working images (imagesCopies) in ImageStack into file

        @param filename: name of stack
        @type filename: L{str}
        """
        from pytom_volume import vol

        stack = vol(self.dimX, self.dimY, len(self.images))
        for ii in range(0, len(self.images)):
            for ix in range(0, self.dimX):
                for iy in range(0, self.dimY):
                    v = self.imageCopies[ii].data.getV(ix, iy, 0)
                    stack.setV(v, ix, iy, ii)
        stack.write(filename)

    def average(self, mask=None):
        """
        average Image Stack

        @param mask: mask is multiplied with average if specified
        @type mask: L{pytom_volume.vol}
        @return: average
        @rtype: L{ptom_volume.vol}
        """
        from pytom.basic.transformations import general_transform2d

        if self.averageData == None:
            self.averageData = vol(self.dimX, self.dimY, 1)
        self.averageData.setAll(0.)
        for ii in range(0, len(self.images)):
            self.imageCopies[ii].data = general_transform2d(
                v=self.images[ii].data, rot=-self.images[ii].rotation,
                shift=[-self.images[ii].shiftX, -self.images[ii].shiftY],
                scale=1., order=[2, 1, 0], crop=True)
            self.averageData = self.averageData + self.imageCopies[ii].data
        if mask:
            self.averageData = self.averageData * mask
        return self.averageData

    def writeAverage(self, filename):
        """
        write average to file
        @param filename: name of file
        @type filename: L{str}
        """
        self.averageData.write(filename)


    def subtractImageFromAverage(self, ii, mask=None):
        """
        remove a stack image from current average - working copy is used!

        @param ii: index of image to be subtracted
        @type ii: L{int}
        @param mask: mask is multiplied with average if specified
        @type mask: L{pytom_volume.vol}
        @return: average of remaining images
        @rtype: L{pytom_volume.vol}
        """
        av = vol(self.dimX, self.dimY, 1)
        av.copyVolume(self.averageData)
        if mask:
            av.__sub__(self.imageCopies[ii].data * mask)
        else:
            av.__sub__(self.imageCopies[ii].data)
        return av

    def exMaxAlign(self, niter=10, mask=None):
        """
        expectation maximization alignment of imageStack

        @param niter: number of iterations
        @type niter: L{int}
        @param mask: mask applied to average before alignment
        @type mask: L{pytom_volume.vol}
        @author: FF
        """
        from pytom.basic.correlation import nXcf, subPixelPeak
        from pytom_volume import peak

        for iexMax in range(0, niter):
            if self.verbose:
                print "Starting Ex-Max iteration " + str(iexMax + 1) + " of " + str(niter)
            meanX = 0.
            meanY = 0.
            #Max - step
            self.average(mask=mask)
            # Ex-step
            for ii in range(0, len(self.images)):
                ccf = nXcf(volume=self.images[ii].data,
                           template=self.subtractImageFromAverage(ii=ii, mask=mask))
                pos = peak(ccf)
                peakinfo = subPixelPeak(scoreVolume=ccf, coordinates=pos,
                                        cubeLength=8, verbose=self.verbose)
                peakval = peakinfo[0]
                pos = peakinfo[1]
                self.images[ii].shiftX = float(pos[0] - self.dimX)
                self.images[ii].shiftY = float(pos[1] - self.dimY)
                meanX = meanX + self.images[ii].shiftX
                meanY = meanY + self.images[ii].shiftY
            # set mean shifts to 0
            meanX = meanX / len(self.images)
            meanY = meanY / len(self.images)
            for ii in range(0, len(self.images)):
                self.images[ii].shiftX = self.images[ii].shiftX - meanX
                self.images[ii].shiftY = self.images[ii].shiftY - meanY

    def normalize(self, normtype="StdMeanInMask", mask=None):
        """
        normalize each image in ImageStack

        @param normtype: normalization type
        @type normtype: str ("StdMeanInMask", )
        @param mask: mask used for normalization
        @type mask: L{pytom_volume.vol}
        """
        for ii in range(0, len(self.images)):
            self.p = self.images[ii].normalize(normtype=normtype, mask=mask,
                                               p=self.p)
            if self.verbose:
                print "normalized image " + str(ii + 1) + " of " + str(len(self.images))

    def taper_edges(self, width, taper_mask=None):
        """
        taper edges of image stack
        @param width: width of tapered edge
        @type width: L{float} or L{int}
        @param taper_mask: mask for tapering - if None it will be generated
        @type taper_mask: L{pytom_volume.vol}

        @return taper_mask
        @rtype: L{pytom_volume.vol}
        """
        for ii in range(0, len(self.images)):
            taper_mask = self.images[ii].taper_edges(width=width, taper_mask=taper_mask)
            if self.verbose:
                print "edge tapered for image " + str(ii + 1) + " of " + str(len(self.images))
        return taper_mask

    def bandpass(self, lowfreq, hifreq, smooth=0., bpf=None):
        """
        bandpass image stack

        @param lowfreq: lowest frequency of filter in pixel (=hi-pass)
        @param hifreq: highest frequency of filter in pixel (=low-pass)
        @param smooth: smoothing in pixel
        @param bpf: bandpass filter object

        @return: bandpass filter object
        @rtype: L{pytom.basic.structures.BandPassFilter}
        """
        for ii in range(0, len(self.images)):
            bpf = self.images[ii].bandpass(lowfreq=lowfreq, hifreq=hifreq,
                                           smooth=smooth, bpf=bpf)
            if self.verbose:
                print "bandpass filtered image " + str(ii + 1) + " of " + str(len(self.images))

        return bpf

    def bandpass_in_nyquist(self, lowfreq, hifreq, smooth=0., bpf=None):
        """
        bandpass image stack with frequencies specified in Nyquist

        @param lowfreq: lowest frequency of filter in Nyquist (=hi-pass)
        @param hifreq: highest frequency of filter in Nyquist (=low-pass)
        @param smooth: smoothing in Nyquist
        @param bpf: bandpass filter object

        @return: bandpass filter object
        @rtype: L{pytom.basic.structures.}
        """
        for ii in range(0, len(self.images)):
            bpf = self.images[ii].bandpass_in_nyquist(lowfreq=lowfreq, hifreq=hifreq,
                                                      smooth=smooth, bpf=bpf)
        return bpf

