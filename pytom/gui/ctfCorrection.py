from cupyx.scipy.fftpack.fft import get_fft_plan
from cupyx.scipy.fftpack.fft import fftn as fftnP
from cupyx.scipy.fftpack.fft import ifftn as ifftnP
from cupyx.scipy.ndimage import map_coordinates
from pytom.gui.guiFunctions import datatype
import os
import numpy as np
from pytom.agnostic.mpi import MPI
import mrcfile
import time

#mpi = MPI()



def CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binning_factor,rotation_angle):
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
        #args.append((fname, new_fname, p, metafile, gs, fs, binning_factor, rotation_angle))
        CorrectProjection_proxy(fname, new_fname, p, metafile, gs, fs, binning_factor, rotation_angle)
        break
    # Parallelization
    #res = mpi.parfor(CorrectProjection_proxy, args)

def CalculateDefocusModel(dz0, Objectpixelsize, Imdim, Tiltangle, Tiltaxis, tx=0,ty=0,sfx=1,sfy=1):
    """
    This program models a defocus plane.
    The defocus plane rely on the assumption that the defocus is constant parallel 
    to the y-axis.

    INPUT
    @param dz0: mean defocus of the plane in mu m
    @type dz0: C{float}
    @param Objectpixelsize: Objectpixelsize of the projection /in nm
    @type Objectpixelsize: C{float}
    @param Imdim: image dimension of the projection and of the defocus plane
    @type Imdim: C{int}
    @param Tiltangle: tiltangle of the projection /in deg
    @type Tiltangle: C{float}
    @param Tiltaxis: tiltaxis of the projection /in deg
    @type Tiltaxis: C{float}
    @param tx: translation of the projection in x-direction /in pixel
    @type tx: C{float}
    @param ty: translation of the projection in y-direction /in pixel
    @type ty: C{float}
    @param sfx: scaling of the projection in x-direction
    @type sfx: C{float}
    @param sfy: scaling of the projection in y-direction
    @type sfy: C{float}

    OUTPUT
    @return dzp: defocus plane /in mu m
    @rtype: C{list}
    """
    # Calculate defocus plane
    import numpy as np
    dzp = dz0 + (1./1000) * np.arange(-Imdim, Imdim)*Objectpixelsize * np.tan( Tiltangle * np.pi / 180)
    dzp = np.transpose(np.tile(dzp, [2*Imdim,1]))

    # Inverse transformation
    if Tiltaxis:
        print('inv transform')
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
    @param proj: projection
    @type proj: 
    @param alpha: in-plane rotation, tilt axis is rotated parallel to the y-axis
    @param tx: translation in x-direction
    @param ty: translation in y-direction
    @param sfx: scaling in x-direction
    @param sfy: scaling in y-direction
    @param SmoothFlag: 1 smoothing is applied - 0 no smoothing


    OUTPUT
    projt - transformed projection
    """

    import numpy as np
    if SmoothFlag:
        proj = SmoothEdges(proj,int(np.around(proj.shape[0]/16.)),int(np.around(proj.shape[0]/16.)),np.mean(proj))


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
        projt = SmoothEdges(projt,int(np.around(projt.shape[0]/16.)),int(np.around(projt.shape[0]/16.)),np.mean(projt))

    return projt

def SmoothEdges(image2,border=None,sigma=None,imdev=None):
    """
    @param image2: input image
    @type image2: numpy array
    @param border: border in pixel (default: 1/16 of image dimensions)
    @type border: C{int}
    @param sigma: smoothing decay at edge in pixel (default == border)
    @type sigma C{int}
    @param imdev: mean value at edge of image
    @type imdev: C{float}
    @return: image
    @rtype image: numpy array
    """
    image = np.array(image2, copy=True)
    #if border is None:
    border = int(np.around(image.shape[1]/16))

    #if sigma is None:
    sigma = border

    ix_up, iy_up = image.shape

    #if imdev is None:
    imdev = (np.sum(image[0,:]) + np.sum(image[ix_up-1,:]) + np.sum(image[1:ix_up-1,0]) + np.sum(image[1:ix_up-1,iy_up-1]))/(2*(ix_up+iy_up)-2)


    xx = np.arange(1, border+1, dtype='float32')
    fact = (np.exp(-2*(xx**2)/(sigma**2))-np.exp(-2))/(1-np.exp(-2))
    #print(fact.shape, image.shape)
    #if sigma < border:
    #    fact[int(np.floor(sigma))-1:int(border)] = 0
    image2 = image -imdev
    for ix in range(border):
        ff = fact[border-ix-1]
        #if iy_up == 1:
        #    print('hit')
        #    image[ix,0]               = imdev + (image[ix,0]-imdev)*ff
        #    image[ix_up-1,0]          = imdev + (image[ix_up-1,0]-imdev)*ff
        #else:
        image[ix,ix:iy_up]        = imdev + (image2[ix,ix:iy_up])*ff
        image[ix_up-1,ix:iy_up]   = imdev + (image2[ix_up-1,ix:iy_up])*ff
        image[ix+1:ix_up,ix]      = imdev + (image2[ix+1:ix_up,ix])*ff
        image[ix+1:ix_up,iy_up-1] = imdev + (image2[ix+1:ix_up,iy_up-1])*ff
        iy_up = iy_up-1
        ix_up = ix_up-1

    return image


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


def imtransform(proj, T):
    """
    transform image 
    @param proj: input image
    @type proj: numpy array
    @param T: transformation matrix (?)
    @type T:
    @return: transformed image
    @rtype: numpy array
    """
    from scipy import mgrid
    import numpy

    cx = proj.shape[0] // 2
    cy = proj.shape[1] // 2
    grid = np.mgrid[-float(cx):proj.shape[0]-cx, -float(cy):proj.shape[1]-cy]
    temp = grid.reshape((2, grid.size // 2))
    x, y = temp
    # temp2 = np.array([y, x, np.ones(x.shape)])
    temp2 = np.array([x, y, np.ones(x.shape)])
    temp3 = np.dot(temp2.transpose(), np.array(numpy.array(T)))

    # yy, xx, zz = temp3.transpose()
    xx, yy, zz = temp3.transpose()
    grid = np.reshape(np.array([xx, yy]), grid.shape)
    grid[0] += cx
    grid[1] += cy

    del temp, temp2, temp3,x,y

    #from scipy.ndimage import map_coordinates
    projt = map_coordinates(np.array(proj), grid.reshape(len(grid),-1), cval=np.mean(proj, dtype='double')).reshape(grid.shape[1:])

    return np.array(projt)

def CorrectProjection_proxy(fname, new_fname, p, metafile, gs, fs, binning_factor, rotation_angle):
    """
    @param fname: filename
    @type fname: C{str}
    @param new_fname:
    @type new_fname: C{str}
    @param p:
    @param metafile: star-file with metadata (contains defoci (long-/short), astig agnle alpha, 
                     voltage, Cs, Amplitude-contrast, Imdim, PixelSpacing, Magnification
    @type metafile: C{str}
    @param gs: grid spacing [2,4,6,...] is the size of the area which is
               corrected with a constant ctf. The ctf value is taken
               from the center of this area. This area builds up the
               corrected projection.
    @param fs: fieldsize [2,4,6,...] & (fs>=gs) is the size of the area which
               is extracted from the projection and corrected with a
               constant ctf. Fieldsize is also the size of the modelled ctf.
    @param binning_factor: de-magnfication factor
    @type binning_factor: C{float}
    @param rotation_angle:
    @return:

    """
    print('Correct projection:', fname)

    # load the metadata
    #from numpy import loadtxt
    from pytom.gui.guiFunctions import loadstar
    metadata = loadstar(metafile, dtype=datatype)
    
    # Alignment parameter
    Tiltangles = metadata['TiltAngle']
    Tiltaxis = rotation_angle
    directions = {0: 'horizontal', 1: 'vertical'}
    direction = directions[ int(np.around(Tiltaxis/90)) % 2 ]

    dz1 = metadata['DefocusU']
    dz2 = metadata['DefocusV']
    alpha = metadata['DefocusAngle']

    Voltage = metadata['Voltage'][p]
    Cs = metadata['SphericalAberration'][p]
    A = metadata['AmplitudeContrast'][p]
    Imdim = metadata['ImageSize'][p]

    if Imdim == 0: Imdim = 3710
    
    Objectpixelsize = metadata['PixelSpacing'][p] * 0.1 * binning_factor
                    
    from pytom.agnostic.io import read, write

    # Load projection
    proj = np.array(read(fname))
    proj = np.squeeze(proj) # squeeze it to 2D

    # Create defocus plane dz1
    dz1p = CalculateDefocusModel(dz1[p],Objectpixelsize,Imdim,Tiltangles[p],Tiltaxis)
    
    # Create defocus plane dz2
    dz2p = CalculateDefocusModel(dz2[p],Objectpixelsize,Imdim,Tiltangles[p],Tiltaxis)
  
    # Create astigmatism angle plane
    alphap = (alpha[p]+Tiltaxis)*np.ones((Imdim,Imdim))
    
    # !!! ADDED Tiltaxis(1,p) to astigmatsm angle to correct for Tiltaxis rotation during CTF determination !!! -- SP 7.7.16
    # originally: alphap = alpha[p]*np.ones((Imdim,Imdim))

    projc = CorrectProjection(proj,dz1p,dz2p,alphap,gs,fs,Objectpixelsize,Voltage,Cs,A, direction=direction)

    # Save projection
    try: mrcfile.new(new_fname, projc.T.get(),overwrite=True)
    except: mrcfile.new(new_fname, projc.T, overwrite=True)

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
    xz, yz = np.meshgrid(np.arange(-f, f, 2. * f / size)[size // 2 - 1:], np.arange(-f, f, 2. * f / size))
    # Wave vector and direction
    k = np.sqrt(xz ** 2 + yz ** 2)
    k2 = k ** 2
    k4 = k2 ** 2
    alpha = np.arctan2(xz, yz) + np.pi
    cosalpha = np.cos(2 * (alpha - alpha0))
    return k2, k4, alpha, cosalpha

def CorrectProjection(proj, dzp1, dzp2, alphap, gs, fs, Objectpixelsize, Voltage, Cs, A, direction='horizontal'):
    """
    This function corrects a defocus gradient on a projection with phase-flipping.
    A patch of size fs is cutted out and corrected, then a region of size gs
    is cutted out from the patch and builds up the corrected projection.

    INPUT
    @param proj: projection with defocus gradient
    @type proj: numpy array
    @param dz1p: defocus plane which describes the defocus gradient
                 in the direction of the first principal axis
    @type dz1p: C{float}
    @param dz2p: defocus plane which describes the defocus gradient
                in the direction of the second principal axis,
                if empty dzp2 equals dzp1
    @type dz2p: C{float}
    @param alpha: astigmatism angle plane,
                  if empty astigmatism angle is zero
    @type alpha: C{float}
    @param gs: grid spacing [2,4,6,...] is the size of the area which is
               corrected with a constant ctf. The ctf value is taken
               from the center of this area. This area builds up the
               corrected projection.
    @param fs: fieldsize [2,4,6,...] & (fs>=gs) is the size of the area which
               is extracted from the projection and corrected with a
               constant ctf. Fieldsize is also the size of the modelled ctf.
    @param Objectpixelsize: pixel size /in nm
    @type Objectpixelsize: C{float}
    @param Voltage: /in kV
    @type Voltage: C{float}
    @param Cs: /in mm
    @type Cs: C{float}
    @param A: Amplitude contrast fraction (e.g., 0.1)
    @type A: C{float}

    OUTPUT
    @return: corrected projection
    @rtype: numpy array
    """

    import numpy
    #proj2 = np.zeros((3712,3712), dtype=np.float32)
    #proj2[1:-1,1:-1] = proj[:,:]
    #proj = proj2
    # Initialize corrected projection
    projc = np.copy(proj)

    # Prepare coordinate list
    x = np.arange(gs//2-1, proj.shape[0]-gs//2, gs)
    ind = ~np.logical_or(x<fs//2-1, x>proj.shape[0]-fs//2-1)
    x = x[ind]
    y = np.arange(gs//2-1, proj.shape[1]-gs//2, gs)
    ind = ~np.logical_or(y<fs//2-1, y>proj.shape[1]-fs//2-1)
    y = y[ind]

    mask = GenerateMask(projc, border=fs//2)
    #tmp = projc#[:, y[yy] - fs // 2 + 1:y[yy] + fs // 2 + 1]  # , copy=True)
    tmp2 = SmoothEdgesMaskBased(projc,mask)
    projfft = np.fft.rfft2(tmp2)
    
    fftplan = get_fft_plan(projfft)

    k2, k4, alpha, cosalpha = GenerateFixedArraysCTF(projfft.shape[0], Objectpixelsize, alphap[x[0], y[0]])

    sx, sy = np.meshgrid(np.arange(0, projc.shape[0]), np.arange(0, projc.shape[1]))
    alfa = alphap[0][0] * np.pi / 180 + np.pi/2
    sx = sx.astype(np.float32) * np.cos(alfa) + sy.astype(np.float32) * np.sin(alfa)

    projc *= 0
    uno = np.ones_like(projc)
    w = np.zeros_like(uno)



    if direction == 'horizontal':
        for xx in range(len(x)):
            for yy in range(len(y)):
                b,a = int(x[xx]), int(y[yy])
                print(direction, dzp1[a,b], dzp2[a,b])
                ctf = ModelCTF(dzp1[a,b], dzp2[a,b]-2, alphap[a,b], A, Objectpixelsize, Voltage, Cs,
                               projfft.shape[1], projfft.shape[0], k4, k2, alpha, cosalpha)
                #mrcfile.new('fft_im.mrc', np.abs(projfft).get().astype('float32'), overwrite=True)
                #mrcfile.new('ctf.mrc', ctf.get().astype('float32'), overwrite=True)
                #import sys
                #sys.exit()

                cpatch = PhaseDeconv(projfft, ctf)
                #print(cpatch.sum(), (cpatch * ((sx >= b-gs//2+1) * (sx < b+gs//2+1) )).sum())
                #projc += cpatch * (sx >= b-gs//2+1) * (sx < b+gs//2+1)

                u = uno * (sx >= a - gs // 2 + 1 -2) * (sx < a + gs // 2 + 1 +2)
                w += u
                projc += cpatch * u
                print(u.sum(), w.sum(), b, a, projc.sum(), cpatch.sum())
                #print(cpatch.shape)
                #projc[b-gs//2+1:b+gs//2+1, :] = cpatch[b-gs//2+1:b+gs//2+1, :]
            break
    else:
        for xx in range(len(x)):
            for yy in range(len(y)):
                a, b = int(x[xx]), int(y[yy])
                print(direction, dzp1[a,b], dzp2[a,b])
                ctf = ModelCTF(dzp1[a,b], dzp2[a,b], alphap[a,b], A, Objectpixelsize, Voltage, Cs,
                               projfft.shape[1], projfft.shape[0], k4, k2, alpha, cosalpha)
                cpatch = PhaseDeconv(projfft, ctf)
                u = uno * (sx >= b-gs//2+1) * (sx < b+gs//2+1)
                w += u
                projc += cpatch * u
                print(cpatch.sum(), a, b, gs, sx.max(), sx.min(), (cpatch * u).sum())

                #projc[:,b-gs//2+1:b+gs//2+1] = cpatch[:,b-gs//2+1:b+gs//2+1]
                #print(cpatch.shape)
            break
    '''
    for xx in range(len(x)):
        for yy in range(len(y)):
            #s = time.time()
            #temp = np.zeros_like(projc)
            # Model CTF
            #try:
            #    DX, DY = numpy.around(float(dzp1[x[xx], y[yy]]), 8), numpy.around(float(dzp2[x[xx], y[yy]]), 8)
            #    ctf,otf = d[DX][DY]
            #except:
            print(f'defocus = {dzp1[x[xx],y[yy]]}')
            ctf = ModelCTF(dzp1[x[xx],y[yy]], dzp2[x[xx],y[yy]], alphap[x[xx],y[yy]], A, Objectpixelsize, Voltage, Cs,
                           projfft.shape[1], projfft.shape[0], k4, k2, alpha, cosalpha)
            #np.cuda.stream.get_current_stream().synchronize()
            #print('time modelctfgen:\t', time.time()-s)
            #if otf is None: print(otf)
            # Extract correction patch
            # smooth edges
            # and do the phase flipping

            #temp[x[xx]-fs//2+1:x[xx]+fs//2+1,y[yy]-fs//2+1:y[yy]+fs//2+1] = tmp2
            #ss = time.time()
            cpatch = PhaseDeconv(projfft, ctf)
            np.cuda.stream.get_current_stream().synchronize()


            #print('time deconvolve:\t', time.time() -ss )
            # Cut correction patch
            cpatchc = cpatch[fs//2-gs//2:fs//2+gs//2,fs//2-gs//2:fs//2+gs//2]

            # Paste correction patch into corrected projection
            projc[:,y[yy]-gs//2+1:y[yy]+gs//2+1] = cpatch[:,y[yy]-gs//2+1:y[yy]+gs//2+1]
            #projc[:,:16] = cpatch[:,:16]
            #np.cuda.stream.get_current_stream().synchronize()
            #print('time copy:\t\t', time.time() - sss, alphap[x[xx],y[yy]] )
            
            DX,DY = numpy.around(float(dzp1[x[xx],y[yy]]),8), numpy.around(float(dzp2[x[xx],y[yy]]),8)

            if not DX in d.keys():
                d[DX] = {}

            if not DY in d[DX].keys():
                d[DX][DY] = [ctf, otf]
            
            #print('time total:\t\t', time.time()-s, '\n')
        break
    '''
    #print(projc.shape)
    mrcfile.new('jj.mrc', w.get().T, overwrite=True)
    w[w<1] =1
    return (projc* 1 * (mask>1-1E-6) / w).astype(np.float32)

def ModelCTF(Dz1,Dz2,alpha0,A,Objectpixelsize,Voltage,Cs,Imdim, Imdim2, k4, k2, alpha, cc):
    """
    Models a CTF.
    INPUT
    @param Dz1: defocus plane which describes the defocus gradient
                in the direction of the first principal axis
    @type Dz1: C{float}
    @param Dz2: defocus plane which describes the defocus gradient
                in the direction of the second principal axis,
                if empty dzp2 equals dzp1
    @type Dz2: C{float}
    @param alpha0: astigmatism angle plane
    @type alpha0: C{float}
    @param A: Amplitude contrast fraction (e.g., 0.1)
    @type A: C{float}
    @param Objectpixelsize: pixel size /in nm
    @type Objectpixelsize: C{float}
    @param Voltage: /in kV
    @type Voltage: C{float}
    @param Cs: /in mm
    @type Cs: C{float}
    @param Imdim: image dimension
    @type Imdim: C{int}
    @return: CTF in Fourier space
    @rtype: L{numpy.ndarray}
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

    ctf = -np.sqrt(1.-A**2)*np.sin(chi)+A*np.cos(chi)
    
    # CHANGED sign for amplitude contrast to fix phase flip for low frequencies -- SP 7.9.16
    # originally: ctf = -np.sqrt(1-A**2)*np.sin(chi)-A*np.cos(chi)

    return ctf

def PhaseDeconv(projfft,ctf):
    """
    Phase deconvolution of image
    @param proj: input projection
    @type proj: L{numpy.ndarray}
    @param otf: object transfer function 
    @type otf: L{numpy.ndarray}
    @return: phase flipped image
    @rtype: L{numpy.ndarray}

    """
    # FT

    ctf2 = np.fft.fftshift(ctf,axes=0) #[:,:ctf.shape[0]//2+1]
    #mrcfile.new('jumbo.mrc', (ctf2).get().astype('float32'), overwrite=True)
    otf = np.sign(ctf2)
    # Set sign of center pixel to (+)
    otf[otf.shape[0]//2,-1] = 1


    # IFT of Deconvolved Phase
    #projdeconv = (np.fft.irfft2(projfft*otf))

    projdeconv = np.fft.irfft2(projfft*otf)


    return projdeconv


if __name__ == '__main__':
    # parse args
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-u", dest="uprefix",
                      help="Name prefix of uncorrected projections")
    parser.add_option("-c", dest="cprefix",
                      help="Name prefix of corrected projections")
    parser.add_option('--rotationAngle', dest='rotationAngle', 
                      help='In-plane Rotation Angle of the tiltaxis. Please note that Alignment corrects for 180')
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
    rotangle = int(options.rotationAngle)
    gs = int(options.gridSpacing)
    fs = int(options.fieldSize)
    binningFactor = int(options.binningFactor)

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
    CorrectTiltseries(metafile, uprefix, cprefix, gs, fs, binningFactor, rotangle)
    print(f'total time: {time.time()-s}')
    #mpi.end()


