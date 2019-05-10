
from pytom.gui.guiFunctions import datatype
import os
import numpy as np
from tompy.mpi import MPI
import mrcfile


mpi = MPI()


def CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binning_factor):
    """
    This function corrects the astigmatic defocus gradient of a tiltseries.

    INPUT
    dz1 --- defocus of tiltseries first principal component
    dz2 --- defocus of tiltseries second principal component
    alpha --- astigmatism angle
    alg --- alignment parameter of the tiltseries
    PathTS --- path of tiltseries
    NameTS --- name of projections of tiltseries
    PathTScor --- path of the corrected tiltseries
    NameTScor --- name of the projections of the corrected tiltseries
    gs --- grid spacing /in pixel
    fs --- fieldsize /in pixel
    Objectpixelsize --- /in nm
    Voltage --- /in kV
    Cs --- /in mm
    Amplitude/Phase ratio

    OUTPUT
    CTF corrected tiltseries
    """
    
    # NumOfProj
    import glob

    fnames = [line for line in sorted(glob.glob(uprefix+"*"))]
    
    NumOfProj = len(glob.glob(uprefix+"*"))
    args = []

    for p, fname in enumerate(fnames):
        # Corrected projection name
        new_fname = cprefix + os.path.basename(fname).split("_")[-1]

        args.append((fname, new_fname, p, metafile, gs, fs, binning_factor))
        #CorrectProjection_proxy(fname, new_fname, p, metafile, gs, fs, binning_factor)
    # Parallelization
    res = mpi.parfor(CorrectProjection_proxy, args)

def CalculateDefocusModel(dz0,Objectpixelsize,Imdim,Tiltangle,Tiltaxis,tx=0,ty=0,sfx=1,sfy=1):
    """
This program models a defocus plane.
The defocus plane rely on the assumption that the defocus is constant parallel 
to the y-axis.

INPUT
dz0 - mean defocus of the plane /in mu m
Objectpixelsize --- Objectpixelsize of the projection /in nm
Imdim --- image dimension of the projection and of the defocus plane
Tiltangle --- tiltangle of the projection /in deg
Tiltaxis --- tiltaxis of the projection /in deg
tx --- translation of the projection in x-direction /in pixel
ty --- translation of the projection in y-direction /in pixel
sfx --- scaling of the projection in x-direction
sfy --- scaling of the projection in y-direction

OUTPUT
dzp --- defocus plane /in mu m
dzm --- central line of defocus plane parallel to the x-axis /in mu m
    """
    # Calculate defocus plane

    dzp = dz0 + (1./1000) * np.arange(-Imdim, Imdim)*Objectpixelsize * np.tan(Tiltangle*np.pi/180)
    dzp = np.transpose(np.tile(dzp, [2*Imdim,1]))

    # Inverse transformation
    if Tiltaxis:
        dzp = AlignProjectionINV(dzp,Tiltaxis,tx,ty,sfx,sfy)

    # Cut defocus plane
    dzp = dzp[Imdim-Imdim//2:2*Imdim-Imdim//2,Imdim-Imdim//2:2*Imdim-Imdim//2]

    # Cut central line
    dzm = dzp[:,Imdim//2]

    return dzp

def AlignProjectionINV(proj,alpha,tx,ty,sfx,sfy,SmoothFlag=1):
    """
    This function aligns a projection in backward direction from the aligned
    version to the original micrograph.
    The tranformation takes into account the following parameter:
    translation, rotation and scaling in x and y-direction.
    Interpolation method is bicubic.
    Additionally the image borders can be smoothed, to prevent high contrast
    oscillations if the projection is fourier filtered and for the reconstruction that
    the outer area looks smooth without edges.
    Image areas without information are filled with the mean value of the projection.

    INPUT
    proj - projection
    alpha - in-plane rotation, tilt axis is rotated parallel to the y-axis
    tx - translation in x-direction
    ty - translation in y-direction
    sfx - scaling in x-direction
    sfy - scaling in y-direction
    SmoothFlag - 1 smoothing is applied
               - 0 no smoothing

    OUTPUT
    projt - transformed projection
    """
    if SmoothFlag:
        proj = SmoothEdges(proj,int(np.round(proj.shape[0]/16.)),int(np.round(proj.shape[0]/16.)),np.mean(proj))

    # Build tform
    # (1) Move
    m = np.matrix([[1,0,0], [0,1,0], [-ty,-tx,1]])

    # (2) Rotate
    alpha = alpha*np.pi/180
    r = np.matrix([[np.cos(-alpha),-np.sin(-alpha),0], [np.sin(-alpha),np.cos(-alpha),0], [0,0,1]])

    # (3) Scale
    s = np.matrix([[1./sfy,0,0], [0,1./sfx,0], [0,0,1]])

    T = np.linalg.inv(m*r*s)
    T[0,2]=0
    T[1,2]=0
    T[2,2]=1

    # Tranform projection
    projt = imtransform(proj,T)

    # Smooth edges
    if SmoothFlag:
        projt = SmoothEdges(projt,int(np.round(projt.shape[0]/16.)),int(np.round(projt.shape[0]/16.)),np.mean(projt))

    return projt

def SmoothEdges(image2,border=None,sigma=None,imdev=None):
    """
    """
    image = np.array(image2, copy=True)
    if border is None:
        border = int(np.round(image.shape[0]/16))

    if sigma is None:
        sigma = border

    ix_up, iy_up = image.shape

    if imdev is None:
        imdev = (np.sum(image[0,:]) + np.sum(image[ix_up-1,:]) + np.sum(image[1:ix_up-1,0]) + np.sum(image[1:ix_up-1,iy_up-1]))/(2*(ix_up+iy_up)-2)

    xx = np.arange(1, border+1, dtype='float32')
    fact = (np.exp(-2*(xx**2)/(sigma**2))-np.exp(-2))/(1-np.exp(-2))
    if sigma < border:
        fact[int(np.floor(sigma))-1:int(border)] = 0

    for ix in range(border):
        if iy_up == 1:
            image[ix,0] = imdev + (image[ix,0]-imdev)*fact[border-ix-1]
            image[ix_up-1,0] = imdev + (image[ix_up-1,0]-imdev)*fact[border-ix-1]
        else:
            image[ix,ix:iy_up] = imdev + (image[ix,ix:iy_up]-imdev)*fact[border-ix-1]
            image[ix_up-1,ix:iy_up] = imdev + (image[ix_up-1,ix:iy_up]-imdev)*fact[border-ix-1]
            image[ix+1:ix_up,ix] = imdev + (image[ix+1:ix_up,ix]-imdev)*fact[border-ix-1]
            image[ix+1:ix_up,iy_up-1] = imdev + (image[ix+1:ix_up,iy_up-1]-imdev)*fact[border-ix-1]
            iy_up = iy_up-1

        ix_up = ix_up-1


    return image

def imtransform(proj, T):
    from scipy import mgrid
    cx = proj.shape[0] / 2
    cy = proj.shape[1] / 2
    grid = mgrid[-float(cx):proj.shape[0]-cx, -float(cy):proj.shape[1]-cy]
    temp = grid.reshape((2, grid.size / 2))

    x, y = temp
    # temp2 = np.array([y, x, np.ones(x.shape)])
    temp2 = np.array([x, y, np.ones(x.shape)])

    temp3 = np.dot(temp2.transpose(), T)
    # yy, xx, zz = temp3.transpose()
    xx, yy, zz = temp3.transpose()
    grid = np.reshape(np.array([xx, yy]), grid.shape)
    grid[0] += cx
    grid[1] += cy

    from scipy.ndimage import map_coordinates
    projt = map_coordinates(proj, grid, order=2, cval=np.mean(proj, dtype='double'))

    return projt

def CorrectProjection_proxy(fname, new_fname, p, metafile, gs, fs, binning_factor):
    """
    """
    print('Correct projection:', fname)

    # load the metadata
    from numpy import loadtxt
    metadata = loadtxt(metafile, dtype=datatype)
    
    # Alignment parameter
    Tiltangles = metadata['TiltAngle']
    Tiltaxis = metadata['InPlaneRotation']

    dz1 = metadata['DefocusU']
    dz2 = metadata['DefocusV']
    alpha = metadata['DefocusAngle']

    Voltage = metadata['Voltage'][p]
    Cs = metadata['SphericalAberration'][p]
    A = metadata['AmplitudeContrast'][p]
    Imdim = metadata['ImageSize'][p]

    if Imdim == 0: Imdim = 3710

    Objectpixelsize = metadata['PixelSpacing'][p] * 10. * metadata['Magnification'][p] * binning_factor
                    
    from tompy.io import read, write

    # Load projection
    proj = read(fname)
    proj = np.squeeze(proj) # squeeze it to 2D

    # Create defocus plane dz1
    dz1p = CalculateDefocusModel(dz1[p],Objectpixelsize,Imdim,Tiltangles[p],Tiltaxis[p])
    
    # Create defocus plane dz2
    dz2p = CalculateDefocusModel(dz2[p],Objectpixelsize,Imdim,Tiltangles[p],Tiltaxis[p])
    
    # Create astigmatism angle plane
    alphap = (alpha[p]+Tiltaxis[p]*0)*np.ones((Imdim,Imdim))
    
    # !!! ADDED Tiltaxis(1,p) to astigmatsm angle to correct for Tiltaxis rotation during CTF determination !!! -- SP 7.7.16
    # originally: alphap = alpha[p]*np.ones((Imdim,Imdim))

    projc = CorrectProjection(proj,dz1p,dz2p,alphap,gs,fs,Objectpixelsize,Voltage,Cs,A)

    # Save projection
    mrcfile.new(new_fname, projc.T,overwrite=True)
    #write(new_fname, projc, Tiltangles[p])

    return True

def CorrectProjection(proj, dzp1, dzp2, alphap, gs, fs, Objectpixelsize, Voltage, Cs, A):
    """
    This function corrects a defocus gradient on a projection with phase-flipping.
    A patch of size fs is cutted out and corrected, then a region of size gs
    is cutted out from the patch and builds up the corrected projection.

    INPUT
    proj    --- projection with defocus gradient
    dz1p    --- defocus plane which describes the defocus gradient
                in the direction of the first principal axis
    dz2p    --- defocus plane which describes the defocus gradient
                in the direction of the second principal axis,
                if empty dzp2 equals dzp1
    alpha   --- astigmatism angle plane,
                if empty astigmatism angle is zero
    gs      --- grid spacing [2,4,6,...] is the size of the area which is
                corrected with a constant ctf. The ctf value is taken
                from the center of this area. This area builds up the
                corrected projection.
    fs      --- fieldsize [2,4,6,...] & (fs>=gs) is the size of the area which
                is extracted from the projection and corrected with a
                constant ctf. Fieldsize is also the size of the modelled ctf.
    Objectpixelsize --- /in nm
    Voltage --- /in kV
    Cs --- /in mm
    A ---

    OUTPUT
    projc --- corrected projection
    """
    # Initialize corrected projection
    projc = np.copy(proj)

    # Prepare coordinate list
    x = np.arange(gs//2-1, proj.shape[0]-gs//2, gs)
    ind = ~np.logical_or(x<fs//2-1, x>proj.shape[0]-fs//2-1)
    x = x[ind]
    y = np.arange(gs//2-1, proj.shape[1]-gs//2, gs)
    ind = ~np.logical_or(y<fs//2-1, y>proj.shape[1]-fs//2-1)
    y = y[ind]

    for xx in range(len(x)):
        for yy in range(len(y)):
            # Model CTF
            ctf = ModelCTF(dzp1[x[xx],y[yy]],dzp2[x[xx],y[yy]],alphap[x[xx],y[yy]],A,Objectpixelsize,Voltage,Cs,fs)
            
            # Extract correction patch
            # smooth edges
            # and do the phase flipping
            tmp = np.array(proj[x[xx]-fs//2+1:x[xx]+fs//2+1,y[yy]-fs//2+1:y[yy]+fs//2+1], copy=True)
            tmp2 = SmoothEdges(tmp)
            cpatch = PhaseDeconv(tmp2,ctf)
            
            # Cut correction patch
            cpatchc = cpatch[fs//2-gs//2:fs//2+gs//2,fs//2-gs//2:fs//2+gs//2]
            
            # Paste correction patch into corrected projection
            projc[x[xx]-gs//2+1:x[xx]+gs//2+1,y[yy]-gs//2+1:y[yy]+gs//2+1] = cpatchc


    return projc

def ModelCTF(Dz1,Dz2,alpha0,A,Objectpixelsize,Voltage,Cs,Imdim):
    """
    Models a CTF.
    """
    # Units -> SI units
    Dz1 = Dz1*10**(-6)
    Dz2 = Dz2*10**(-6)
    alpha0 = alpha0*(np.pi/180)
    Objectpixelsize = Objectpixelsize*10**(-9)
    Voltage = Voltage*1000
    Cs = Cs*10**(-3)

    # Electron wavelength
    lmbd = np.sqrt(150.4/((Voltage*(1+Voltage/1022000.))))*10**(-10)

    # Interpolation field
    f = 1/(2.*Objectpixelsize)
    x, y = np.meshgrid(np.arange(-f, f, 2.*f/Imdim), np.arange(-f, f, 2.*f/Imdim))

    # Wave vector and direction
    k = np.sqrt(x**2+y**2)
    alpha = np.arctan2(x,y)+np.pi

    # Defocus
    deltaDz = 0.5*(Dz1 + Dz2 + (Dz1 - Dz2)*np.cos(2*(alpha-alpha0)))

    # Phase shift due to lens aberrations and defocus 
    chi = np.pi*lmbd*k**2*(deltaDz - 0.5*lmbd**2*k**2*Cs)

    # Contrast transfer function
    ctf = -np.sqrt(1-A**2)*np.sin(chi)+A*np.cos(chi)
    
    # CHANGED sign for amplitude contrast to fix phase flip for low frequencies -- SP 7.9.16
    # originally: ctf = -np.sqrt(1-A**2)*np.sin(chi)-A*np.cos(chi)

    return ctf

def PhaseDeconv(proj,otf):
    """
    """
    # FT
    projfft = np.fft.fftshift(np.fft.fft2(proj))

    # Prepare OTF
    otf = np.sign(otf)

    # Set sign of center pixel to (+)
    otf[otf.shape[0]//2,otf.shape[1]//2] = 1

    # Phase deconvolution
    projfft = projfft*otf

    # IFT
    projdeconv = np.real(np.fft.ifft2(np.fft.fftshift(projfft)))

    return projdeconv


if __name__ == '__main__':
    # parse args
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-u", dest="uprefix",
                      help="Name prefix of uncorrected projections")
    parser.add_option("-c", dest="cprefix",
                      help="Name prefix of corrected projections")
    parser.add_option('--gridSpacing', dest='gridSpacing', help='Grid Spacing')
    parser.add_option('--fieldSize', dest='fieldSize', help='Field Size')
    parser.add_option("--metafile", dest="metafile",
                      help="Filename of metadata file", metavar="FILE")
    parser.add_option("--binningFactor", dest="binningFactor",
                      help="Object pixelsize")

    (options, args) = parser.parse_args()

    metafile = options.metafile
    uprefix  = options.uprefix
    cprefix  = options.cprefix

    gs = int(options.gridSpacing)
    fs = int(options.fieldSize)
    binningFactor = int(options.binningFactor)

    # start the clustering
    mpi.begin()

    CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binningFactor)

    mpi.end()


