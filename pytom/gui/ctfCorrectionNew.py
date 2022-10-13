from cupyx.scipy.fftpack.fft import get_fft_plan
from cupyx.scipy.fftpack.fft import fftn as fftnP
from cupyx.scipy.fftpack.fft import ifftn as ifftnP
from cupyx.scipy.ndimage import map_coordinates
from pytom.gui.guiFunctions import datatype
import os
import numpy as np
#from agnostic.mpi import MPI
import mrcfile
import time

#mpi = MPI()


def CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binning_factor, rotation_angle, def_grad_strip=2.5):
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
    from pytom.gui.guiFunctions import loadstar

    fnames = [line for line in sorted(glob.glob(uprefix+"*"))]
    
    NumOfProj = len(glob.glob(uprefix+"*"))
    args = []
    metadata = loadstar(metafile, dtype=datatype)

    for p, fname in enumerate(fnames):
        # Corrected projection name
        new_fname = cprefix + os.path.basename(fname).split("_")[-1]

        #args.append((fname, new_fname, p, metafile, gs, fs, binning_factor, rotation_angle))
        CorrectProjection_proxy(fname, new_fname, p, metadata, gs, fs, binning_factor, rotation_angle, def_grad_strip)

    # Parallelization
    #res = mpi.parfor(CorrectProjection_proxy, args)

def SmoothEdgesMaskBased(image2, mask=None):
    #import time
    s = time.time()
    image = np.array(image2, copy=True)
    border = int(np.around(image.shape[0]/16))
    sigma = border

    ix_up, iy_up = image.shape
    imdev = (np.sum(image[0,:]) + np.sum(image[ix_up-1,:]) + np.sum(image[1:ix_up-1,0]) + np.sum(image[1:ix_up-1,iy_up-1]))/(2*(ix_up+iy_up)-2)

    #xx = np.arange(1, border+1, dtype=np.float32)
    #fact = (np.exp(-2*(xx**2)/(sigma**2))-np.exp(-2))/(1-np.exp(-2))

    #if mask is None:
    #mask = np.zeros(image2.shape)
    #x,y = image2.shape
    #for i in range(border+1):
    #    mask[i:x-i,i:y-i] = fact[border-i-1]
    #    if border == i:
    #        mask[i:x-i,i:y-i] = 1

    image = (image-imdev)* mask + imdev
    #print(time.time()-s)
    return image

def CorrectProjection_proxy(fname, new_fname, p, metadata, gs, fs, binning_factor, rotation_angle, def_grad_strip=2.5):
    """
    """
    print('Correct projection:', fname)

    # load the metadata
    #from numpy import loadtxt

    # Alignment parameter
    Tiltaxis        = rotation_angle
    Tiltangle       = metadata['TiltAngle'][p]
    dz1             = -1*metadata['DefocusU'][p]
    dz2             = -1*metadata['DefocusV'][p]
    alpha           = metadata['DefocusAngle'][p]
    Voltage         = metadata['Voltage'][p]
    Cs              = metadata['SphericalAberration'][p]
    A               = metadata['AmplitudeContrast'][p]
    Imdim           = metadata['ImageSize'][p]
    Objectpixelsize = metadata['PixelSpacing'][p] * 0.1 * binning_factor

    from pytom.agnostic.io import read, write
    from pytom.agnostic.tools import paste_in_center
    # Load projection


    # Fast options are 3840, 3888, 4096
    size = 3888

    projInit = np.array(read(fname))
    projInit = np.squeeze(projInit) # squeeze it to 2D
    proj = np.zeros((size,size))
    proj = paste_in_center(projInit, proj)

    B = (proj.shape[0]-projInit.shape[0])//2

    sx, sy = proj.shape[0] // 2, proj.shape[1] // 2
    x, y = np.meshgrid(np.arange(-sx, sx), np.arange(-sx, sx))

    import time
    s = time.time()
    projc = CorrectProjection2(proj, dz1, dz2, alpha + Tiltaxis, Tiltangle, Tiltaxis, dfgrad=def_grad_strip, Border=B,
                               ObjectPixelSize=Objectpixelsize,Voltage=Voltage,CS=Cs,AmpContrast=A, x=x, y=y)
    print(time.time()-s)
    # Save projection
    try: mrcfile.new(new_fname, projc.T.get().astype('float32'),overwrite=True)
    except: mrcfile.new(new_fname, projc.T.astype('float32'), overwrite=True)

    #write(new_fname, projc, Tiltangles[p])

    return True

def GenerateMask(image,border=None):
    mask = np.zeros((image.shape[1],image.shape[1]))
    x2, y2 = mask.shape
    if border is None: border = int(np.around(mask.shape[1]/16))
    sigma=border
    xxx = np.arange(1, border + 1, dtype=np.float32)
    fact = (np.exp(-2 * (xxx ** 2) / (sigma ** 2)) - np.exp(-2)) / (1 - np.exp(-2))

    for i in range(border + 1):
        mask[i:x2 - i, i:y2 - i] = fact[border - i - 1]
        if border == i:
            mask[i:x2 - i, i:y2 - i] = 1
    return mask

def GenerateFixedArraysCTF(size, Objectpixelsize, alpha0):
    f = 1 / (2. * Objectpixelsize * 10 ** (-9))
    xz, yz = np.meshgrid(np.arange(-f, f, 2. * f / size)[:size][size // 2 - 1:], np.arange(-f, f, 2. * f / size)[:size])
    # Wave vector and direction
    k = np.sqrt(xz ** 2 + yz ** 2)
    k2 = k ** 2
    k4 = k2 ** 2
    alpha = np.arctan2(xz, yz) + np.pi
    cosalpha = np.cos(2 * (alpha - alpha0))
    return k2, k4, alpha, cosalpha

def CorrectProjection2(proj, dzp1, dzp2, astigangle, tiltangle, tiltaxis, Border=145, ObjectPixelSize=0.175, Voltage=200, CS=2.7,
                      AmpContrast=0.08, dfgrad=2.5, overlap=0.15, x=None, y=None):

    import time

    s = time.time()

    astigangle = np.deg2rad(astigangle) + np.pi/2
    tiltangle = np.deg2rad(tiltangle)
    tiltaxisrad = np.deg2rad(tiltaxis)


    sx, sy = proj.shape[0]//2, proj.shape[1]//2

    border = Border
    weights = np.zeros_like(proj, dtype=np.float32)
    corr_proj = np.zeros_like(proj, dtype=np.float32)


    #x, y = np.meshgrid(np.arange(-sx, sx), np.arange(-sx, sx))
    x = x.astype(np.float32) * np.cos(tiltaxisrad) + y.astype(np.float32) * np.sin(tiltaxis)
    del y

    # Calculate largest perpendicular distance to tiltaxis in image (excluding border)
    rcent = np.sqrt((sy - border) ** 2 + (sx - border) ** 2)
    max_distance = np.cos((45 - (tiltaxis % 90)) * np.pi / 180) * rcent * -1
    sang = np.sin(np.abs(tiltangle)) * ObjectPixelSize

    # Total gradient of defocus across image
    df = 2 * np.abs(max_distance) * sang
    steps = int(max(1, np.ceil(np.abs(df) / dfgrad)))
    stepsize = 2 * np.abs(max_distance) / steps
    overlap = (stepsize / 2 * overlap) // 1

    cur = max_distance
    dzp1 -= df*1E-3 / 2
    dzp2 -= df*1E-3 / 2
    np.cuda.stream.get_current_stream().synchronize()

    mask = GenerateMask(proj, border=border+128)
    tmp2 = SmoothEdgesMaskBased(proj, mask)
    projfft = np.fft.rfft2(tmp2)


    k2, k4, alphar, cosalpha = GenerateFixedArraysCTF(projfft.shape[0], ObjectPixelSize, astigangle)

    import time

    stepsizeDF = stepsize*sang*1E-3

    for i in range(steps):
        #print(dzp1+stepsizeDF/2, dzp2+stepsizeDF/2, astigangle, AmpContrast, ObjectPixelSize, Voltage, CS)
        ctf = ModelCTF(dzp1+stepsizeDF/2, dzp2+stepsizeDF/2, astigangle, AmpContrast, ObjectPixelSize, Voltage, CS, projfft.shape[1], projfft.shape[0],
                       k4, k2, alphar, cosalpha)
        cpatch = PhaseDeconv(projfft, ctf)

        MIN = x >= cur - overlap
        MAX = x < cur + fs + overlap
        strip = (1 * MIN * MAX).astype(np.float32)
        corr_proj += cpatch * strip
        weights += strip

        cur += stepsize
        dzp1 += stepsizeDF
        dzp2 += stepsizeDF

    weights[weights < 1] = 1
    return (corr_proj / weights)[border:-border,border:-border]

def ModelCTF(Dz1,Dz2,alpha0,A,Objectpixelsize,Voltage,Cs,Imdim, Imdim2, k4, k2, alpha, cc):
    """
    Models a CTF.
    """
    print(Dz1, Dz2)

    # Units -> SI units
    Dz1 = Dz1*10**(-6)
    Dz2 = Dz2*10**(-6)
    alpha0 = alpha0 #*(np.pi/180)
    Objectpixelsize = Objectpixelsize*10**(-9)
    Voltage = Voltage*1000
    Cs = Cs*10**(-3)

    # Electron wavelength
    lmbd = np.sqrt(150.4/((Voltage*(1+Voltage/1022000.))))*10**(-10)

    # Interpolation field

    #f = 1/(2.*Objectpixelsize)
    #x, y = np.meshgrid(np.arange(-f, f, 2.*f/Imdim), np.arange(-f, f, 2.*f/Imdim2))

    # Wave vector and direction
    #k = np.sqrt(x**2+y**2)
    #alpha = np.arctan2(x,y)+np.pi

    # Defocus
    deltaDz = 0.5*(Dz1 + Dz2 + (Dz1 - Dz2)*cc)#np.cos(2*(alpha-alpha0)))

    # Phase shift due to lens aberrations and defocus 
    chi = np.pi*lmbd*k2*deltaDz - 0.5*np.pi*lmbd**3*Cs*k4

    # Contrast transfer function
    ctf = -np.sqrt(1-A**2)*np.sin(chi)+A*np.cos(chi)

    #mrcfile.new('rockit.mrc', ctf.get().astype('float32'),overwrite=True)
    # CHANGED sign for amplitude contrast to fix phase flip for low frequencies -- SP 7.9.16
    # originally: ctf = -np.sqrt(1-A**2)*np.sin(chi)-A*np.cos(chi)

    return ctf

def PhaseDeconv(projfft,ctf):
    """
    """
    # FT

    ctf2 = np.fft.fftshift(ctf,axes=0) #[:,:ctf.shape[0]//2+1]
    #mrcfile.new('jumbo.mrc', (ctf2).get().astype('float32'), overwrite=True)
    otf = np.sign(ctf2)


    # Set sign of center pixel to (+)
    otf[otf.shape[0]//2,-1] = 1
    print(ctf2[0,0])
    # IFT of Deconvolved Phase
    #projdeconv = (np.fft.irfft2(projfft*otf))
    import time
    projdeconv = np.fft.irfft2(projfft*ctf2)
    return projdeconv


if __name__ == '__main__':
    # parse args
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-u", dest="uprefix",
                      help="Name prefix of uncorrected projections")
    parser.add_option("-c", dest="cprefix",
                      help="Name prefix of corrected projections")
    parser.add_option('--rotationAngle', dest='rotationAngle', help='In-plane Rotation Angle of the tiltaxis. Please note that Alignment corrects for 180')
    parser.add_option('--gridSpacing', dest='gridSpacing', help='Grid Spacing')
    parser.add_option('--fieldSize', dest='fieldSize', help='Field Size')
    parser.add_option("--metafile", dest="metafile",
                      help="Filename of metadata file", metavar="FILE")
    parser.add_option("--binningFactor", dest="binningFactor",
                      help="Object pixelsize")
    parser.add_option("--defocusGradient", dest="defocusGradient",
                      help="Defocus gradient across strip (nm)")

    (options, args) = parser.parse_args()

    metafile = options.metafile
    uprefix  = options.uprefix
    cprefix  = options.cprefix
    rotangle = int(options.rotationAngle)
    gs = int(options.gridSpacing)
    fs = int(options.fieldSize)
    binningFactor = int(options.binningFactor)
    defocusGradient = float(options.defocusGradient) if not options.defocusGradient is None else 5.
    global np

    if 1:
        import cupy as np

    else:
        import numpy as np
        global map_coordinates
        from scipy.ndimage import map_coordinates

    #np.cuda.Device(1).use()
    # start the clustering
    #mpi.begin()
    #print(rotangle)
    s = time.time()
    CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binningFactor, rotangle, defocusGradient)
    print(f'total time: {time.time()-s}')
    #mpi.end()


