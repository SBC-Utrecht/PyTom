"""
Gijs' tilt series simulation - used for Ilja's contest paper
"""
import numpy
import sys
# from chimera import runCommand
from pytom.gui.mrcOperations import *
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *
from pytom.basic.transformations import rotate
from pytom_volume import read
from pytom_numpy import vol2npy
import os
import mrcfile
from numpy import *
from numpy.fft import *
import os
from scipy.ndimage import gaussian_filter
import scipy
from nufft.reconstruction import fourier_2d1d_iter_reconstruct
import numpy
from pytom.tools.script_helper import ScriptHelper, ScriptOption
from pytom.tools.parse_script_options import parse_script_options
from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
from pytom.reconstruction.reconstructionStructures import *
from pytom.basic.files import *
import mrcfile
import numpy as xp

# from scipy.ndimage import rotate

# Dictionary containing physical constants required for calculation
phys_const_dict = {
    "c": 299792458, # m/s
    "el": 1.60217646e-19, # C
    "h": 6.62606896e-34, # J*S
    "h_ev": 4.13566733e-15, # eV*s
    "h_bar": 1.054571628e-34, # J*s
    "h_bar_ev": 6.58211899e-16, # eV*s

    "na": 6.02214179e23, # mol-1
    "re": 2.817940289458e-15, # m
    "rw": 2.976e-10, # m

    "me": 9.10938215e-31, # kg
    "me_ev": 0.510998910e6, # ev/c^2
    "kb": 1.3806503e-23, # m^2 kgs^-2 K^-1

    "eps0": 8.854187817620e-12 # F/m
}

# function calculates the electron wavelength given a certain accelaration voltage V
def wavelength_eV2m(V):
    h = phys_const_dict["h"]
    e = phys_const_dict["el"]
    m = phys_const_dict["me"]
    c = phys_const_dict["c"]

    Lambda = h/xp.sqrt(e*V*m*(e/m*V/c**2 + 2 ))

    return Lambda



def tom_bandpass(image, low, hi, smooth=0):
    try:
        s1, s2, s3 = image.shape
    except:
        s1, s2 = image.shape
        s3 = 1
    scf = 1 / (s1 * s2 * s3)

    if smooth < 0.001:
        [x, y, z] = meshgrid(arange(-floor(s1 / 2), -floor(s1 / 2) + s1 - 1),
                             arange(-floor(s2 / 2), -floor(s2 / 2) + s2 - 1),
                             arange(-floor(s3 / 2), -floor(s3 / 2) + s3 - 1))

        r = sqrt(x ** 2 + y ** 2 + z ** 2)
        lowp = (r <= hi)
        highp = (r >= low)
        image = fftshift(fftn(image))
        image = scf * real(ifftn(fftshift(highp * lowp * image)))

    else:
        image = fftshift(fftn(image));
        if low > 0:
            image = real(ifftn(fftshift(spheremask(image, hi, smooth) - spheremask(image, low, smooth))))
        else:
            image = real(ifftn(fftshift(spheremask(image, hi, smooth))))

    return image


def tom_error(a, m=0, v=1):
    if len(a.shape) == 1:
        s1 = a.shape[0]
        c = a + [sqrt(v) * randn(s1) + m]
    elif len(a.shape) == 2:
        s1, s2 = a.shape
        c = a + [sqrt(v) * randn(s1, s2) + m]
    elif len(a.shape) == 3:
        s1, s2, s3 = a.shape
        c = a + [sqrt(v) * randn(s1, s2, s3) + m]
    return c


def spheremask(vol, radius, smooth=0):
    shape = vol.shape

    mask = numpy.ones_like(vol).astype(numpy.float32)
    if len(shape) == 2:
        x, y = meshgrid(arange(-shape[0] / 2, shape[0] / 2, 1), arange(-shape[1] / 2, shape[1] / 2))
        r = sqrt(x ** 2 + y ** 2)
        # imshow(mask[:,:])
        # show()
    if len(shape) == 3:
        x, y, z = meshgrid(arange(-shape[1] / 2, shape[1] / 2, 1), arange(-shape[0] / 2, shape[0] / 2),
                           arange(-shape[2] / 2, shape[2] / 2, 1))
        r = sqrt(x ** 2 + y ** 2 + z ** 2)
        # imshow(mask[mask.shape[0]//2,:,:])
        # show()

    mask[r > radius] = 0

    if smooth > 0.01:
        mask = gaussian_filter(mask, smooth)
    return vol * mask


#def wavelength_eV2nm(ev):
#    # h / sqrt( 2 * me * el * ev * 1)
#    return 6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * ev * 1) ** 0.5


def create_ctf(Dz, vol, pix_size, voltage=200E3, Cs=2.7E-3, sigma=0.):
    '''
    %TOM_CREATE_CTF calculates 2D or 3D CTF (pure phase contrast)
    %
    %   Note: only tested for even dimensions!
    %
    %   [ctf_out amplitude] = tom_create_ctf(Dz, vol, pix_size, voltage, Cs, sigma)
    %
    %PARAMETERS
    %  INPUT
    %   Dz       : Defocus (<0 underfocus, >0 overfocus) (in \mu m);
    %   vol      : Volume (or image)
    %   pix_size : pixel size (in nm) (default: 0.72 nm)
    %   voltage  : accelerating Voltage (in eV) (default: 300 kV)
    %   Cs       : sperical aberration
    %   Sigma    : envelope of ctf (optional). If chosen decay
    %              ctf ~exp(-(freq/sigma)^2) is assumed. SIGMA is in units of Nyqvist.
    %              => ctf(sigma) = 1/e
    %
    %  OUTPUT
    %   ctf_out  : output containing the centrosymmetric ctf. It can be used as
    %              a fourier filter.
    %
    %EXAMPLE
    %   im = tom_emread('proteasome.em');
    %   ctf = tom_create_ctf(-4.4,im.Value,1.0240, 120,2);
    %   tom_imagesc(ctf)
    %
    %REFERENCES
    %
    %SEE ALSO
    %   TOM_CTF TOM_FOURIER TOM_IFOURIER
    %
    %
    %   created by FF 09/14/04
    %
    %   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
    %   Journal of Structural Biology, 149 (2005), 227-234.
    %
    %   Copyright (c) 2004-2007
    %   TOM toolbox for Electron Tomography
    %   Max-Planck-Institute of Biochemistry
    %   Dept. Molecular Structural Biology
    %   82152 Martinsried, Germany
    %   http://www.biochem.mpg.de/tom

    error(nargchk(2,6,nargin));

    if nargin<5
      Cs=2*10^(-3);
    else
       Cs=Cs*10^(-3);
    end;
    if nargin<4
        voltage=300000;
    else
        voltage = voltage * 1000;
    end;

    if nargin <3
        pix_size=0.72*10^(-9);
    else
        pix_size=pix_size*10^(-9);
    '''
    Dz = Dz * 10 ** (-6)

    Lambda = wavelength_eV2m(voltage)

    Ny = 1 / (2 * pix_size)

    if len(vol.shape) > 2:
        R, Y, Z = xp.meshgrid(xp.arange(-Ny, Ny, 2 * Ny / vol.shape[0]), xp.arange(-Ny, Ny, 2 * Ny / vol.shape[1]),
                           xp.arange(-Ny, Ny, 2 * Ny / vol.shape[2]))
        r = xp.sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[1])))
        r = xp.sqrt(R ** 2 + Y ** 2)

    vol = xp.sin(pi / 2 * (Cs * Lambda ** 3 * r ** 4 - 2 * Dz * Lambda * r ** 2))
    amplitude = xp.cos(pi / 2 * (Cs * Lambda ** 3 * r ** 4 - 2 * Dz * Lambda * r ** 2))

    if sigma:
        vol = vol * xp.exp(-(r / (sigma * Ny)) ** 2)
        amplitude = amplitude * xp.exp(-(r / (sigma * Ny)) ** 2)

    return vol, amplitude


def tom_dev(image):
    return image.mean(), image.max(), image.min(), image.std(), image.std() ** 2


def calcCTF(defocus, image, pixelsize, voltage, spherical_aberration, sigma_decay_ctf, amplitude_contrast):
    ctf1, amplitude = create_ctf(defocus, image, pixelsize, voltage, spherical_aberration, sigma_decay_ctf)

    ctf1 = -ctf1
    amplitude = -amplitude
    ctf = (1-amplitude_contrast) * ctf1 + amplitude_contrast * amplitude
    return ctf

def simutomo(dens, defocus, pixelsize, angles, SNR, the, psi, phi, globalsigma=0, gpu=False, pytom=False,
             outputFolder='./', modelID=0, voltage=200e3, amplitude_contrast=0.07, spherical_aberration=2.7E-3, sigma_decay_ctf=0.4):
    '''%
    %   simtomo = av3_simutomo(dens, tiltrange, tiltincr, SNR, the, psi, phi, globalsigma)
    %
    %   simulate a tomogram of a given volume realistically.
    %
    %   fixed CTF parametersacceleration voltage=300kV, C_s=2, b-factor=0.4
    %
    %   INPUT
    %   dens        density (3 dim volme)
    %   def         defocus (in mum)
    %   pixelsize   pixel size (in A)
    %   tiltrange   tilt range ([-60,60]) (in deg)
    %   tiltincr    tilt increment (in deg)
    %   SNR         signal-to-noise ratio (vari_signal/var_noise)
    %   the         rotation of dens - theta angle
    %   psi         psi angle
    %   phi         phi angle
    %   globalsigma approximate stdev of projs by mini-tilt series? - takes
    %               much longer ...
    %
    %
    %   Copyright (c) 2005-2010
    %   Max-Planck-Institute for Biochemistry
    %   Dept. Molecular Structural Biology
    %   82152 Martinsried, Germany
    %   http://www.biochem.mpg.de/foerster
    %'''

    if gpu:
        rotate = gpu_rotate
    elif pytom:
        from pytom.basic import rotate
    else:
        from scipy.ndimage import rotate

    # rotate volume according to specified euler
    # rotvol = rotate(dens,phi)

    rotvol = dens

    # calculate CTF

    ctf = calcCTF(defocus, (dens.sum(axis=0)), pixelsize, voltage, spherical_aberration, sigma_decay_ctf, amplitude_contrast)

    mask = spheremask(ones(dens.shape), dens.shape[0] / 2 - 3, (dens.shape[0] / 2 - 3) / 20)
    sigma = 0
    irun = 0

    # pre-calculate average stdv for a whole tilt series
    if globalsigma > 0.0001:
        for ipro in range(-90, 90, tiltincr):
            irun = irun + 1
            # calc projection
            proj = squeeze(sum(gpu_rotate(rotvol, [270, 90, ipro]), 3))
            # multiply by ctf
            proj = real(ifftn(ifftshift(ctf * fftshift(fftn(proj)))))
            [mv, mn, mx, stv, tmp] = tom_dev(proj)
            sigma = sigma + tmp

        sigma = sigma / irun
    else:
        # calc 0 deg projection
        proj = squeeze(rotvol.sum(axis=0))
        # multiply by ctf
        proj = real(ifftn(fftshift(ctf * fftshift(fftn(proj)))))
        [mv, mn, mx, stv, sigma] = tom_dev(proj)

    # take care of signal reduction due to ctf
    # ratio of white noise with and without mult by ctf
    [mv, mn, mx, stv, corrfac] = tom_dev(real(ifftn(fftshift(
        ctf * spheremask(fftshift(fftn(tom_error(zeros((ctf.shape[0], ctf.shape[1])), 0, 1.)[0, :, :])), 10, 0)))))
    [mv, mn, mx, stv, whitenoise] = tom_dev(real(ifftn(
        fftshift(spheremask(fftshift(fftn(tom_error(zeros((ctf.shape[0], ctf.shape[1])), 0, 1.)[0, :, :])), 10, 0)))))
    corrfac = whitenoise / corrfac
    sigma = sigma * corrfac

    # generate weighting function for reconstruction
    s1, s2, s3 = rotvol.shape
    [x, y] = meshgrid(arange(-s2 / 2, s2 / 2), arange(-s3 / 2, s3 / 2))
    r = sqrt(x ** 2 + y ** 2)
    mask = ones((s2, s3))
    mask[r >= (s2 // 2) - 1] = 0

    # imshow(mask)
    # show()
    # mask = (abs(x)*spheremask(ones_like(rotvol),s1/2-1,1))

    # simulate tomogram
    simtomo = zeros_like(rotvol)
    sx, sy, sz = dens.shape

    projections = []
    nfp = []
    for n, ipro in enumerate(arange(tiltrange[0], tiltrange[1] + 1, tiltincr)):
        proj = rot90(noisefree_projections[n])  # squeeze(rotate(rotvol,ipro).sum(axis=2))

        # ctf dependent contribution of noise
        noisy = tom_error(proj, 0, 0.5 / SNR * sigma)[0, :, :]

        nx, ny = noisy.shape
        # noisy = noisy[nx/2-sx/2:nx/2+sx/2,nx/2-sx/2:nx/2+sx/2]

        # ctf independent contribution of noise
        bb = tom_error(zeros(noisy.shape), 0, 0.5 / SNR * sigma)[0, :, :]
        bg = tom_bandpass(bb, 0, 1, 0.2 * noisy.shape[0])

        # from matplotlib.colors import LogNorm
        # imshow(abs(fftshift(fftn(bg))),norm=LogNorm())
        # show()
        # add both contributions (ratio 1:1) in Fourier space
        tmp = fftshift(fftn(noisy)) + fftshift(fftn(bg))
        tmp = ifftshift(tmp * mask)
        noisy = real(ifftn(tmp))
        # print noisy.shape
        # fig,ax = subplots(1,1,figsize=(10,10))
        # ax.imshow(noisy,cmap='binary')#, vmax=median(final))
        # ax.set_yticks([])
        # ax.set_xticks([])
        # savefig('/Users/gijs/Desktop/simulation2_SNR_{:4.3f}.png'.format(SNR))
        # show()
        # break
        # back projection
        # tom_backproj3dc(simtomo, single(noisy), single(0), single(-ipro), [single(0),single(0),single(0)])
        # convert_numpy_array3d_mrc(proj,'noisyProjections/simulated_proj_model{:02d}_{}.mrc'.format(modNR, n+1))
        # print "{:4d} {:8.2f} {:8.2f}".format( n, noisy.mean(), proj.mean() )

        if not os.path.exists(f'{outputFolder}/model_{modelID}/noisyProjections/'):
            os.mkdir(f'{outputFolder}/model_{modelID}/noisyProjections/')
        out = mrcfile.new(
            f'{outputFolder}/model_{modelID}/noisyProjections/simulated_proj_model{modelID}_{n+1}.mrc',
            noisy.astype(float32), overwrite=True)
        out.close()
        projections.append(noisy)
        nfp.append(proj * -1.)

    # sys.exit()

    # simtomo =  fourier_2d1d_iter_reconstruct(projections, list(range(tiltrange[0],tiltrange[1]+1,tiltincr)), iter)
    # simtomo =  fourier_2d1d_iter_reconstruct(nfp, list(range(tiltrange[0],tiltrange[1]+1,tiltincr)), iter)
    # normalize to stdv
    # [mnv, maxv, minv, stv, dummy] = tom_dev(simtomo)
    # simtomo = (simtomo-mnv)/stv
    # return simtomo
    return ''


def quat_mult(q1, q2):
    r"""
    Return the product of two quaternions

    Args:
       :q1 (array): Length-4 array [:math:`w_1`, :math:`x_1`, :math:`y_1`, :math:`z_1`] that represents the first quaternion

       :q2 (array): Length-4 array [:math:`w_2`, :math:`x_2`, :math:`y_2`, :math:`z_2`] that represents the second quaternion
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return numpy.array([w, x, y, z])


def quat_vec_mult(q, v):
    r"""
    Return the product of a quaternion and a vector

    Args:
       :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] that represents the quaternion

       :v (array): Length-3 array [:math:`v_x`, :math:`v_y`, :math:`v_z`] that represents the vector
    """
    q2 = numpy.array([0., v[0], v[1], v[2]])
    return quat_mult(quat_mult(q, q2), quat_conj(q))[1:]


def quat_conj(q):
    r"""
    Return the conjugate quaternion as a length-4 array [w,-ix,-jy,-kz]

    Args:
       :q (array): Numpy array :math:`[w,x,y,z]` that represents the quaternion
    """
    iq = q.copy()
    iq[1:] = -iq[1:]
    return iq


def rotate_quat(v, q):
    r"""
    Return rotated version of a given vector by a given quaternion

    Args:
       :v (array): Length-3 array :math:`[v_x,v_y,v_z]` that represents the vector

       :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] that represents the quaternion
    """
    return quat_vec_mult(q, v)


def dotproduct(v1, v2):
    return numpy.sum((a * b) for a, b in zip(v1, v2))


def rand_quat():
    r"""
    Obtain a uniform random rotation in quaternion representation ([Shoemake1992]_ pages 129f)
    """
    x0, x1, x2 = numpy.random.random(3)
    theta1 = 2. * numpy.pi * x1
    theta2 = 2. * numpy.pi * x2
    s1 = numpy.sin(theta1)
    s2 = numpy.sin(theta2)
    c1 = numpy.cos(theta1)
    c2 = numpy.cos(theta2)
    r1 = numpy.sqrt(1 - x0)
    r2 = numpy.sqrt(x0)
    q = numpy.array([s1 * r1, c1 * r1, s2 * r2, c2 * r2])
    return q


def euler_from_quat(q, rotation_axes="zxz"):
    if len(rotation_axes) != 3:
        print("Error: rotation_axes = %s is an invalid input." % rotation_axes)
        return
    for s in rotation_axes:
        if s not in "xyz":
            print("Error: rotation_axes = %s is an invalid input." % rotation_axes)
            return
    i1 = 0 if rotation_axes[0] == "x" else 1 if rotation_axes[0] == "y" else 2 if rotation_axes[0] == "z" else None
    i2 = 0 if rotation_axes[1] == "x" else 1 if rotation_axes[1] == "y" else 2 if rotation_axes[1] == "z" else None
    i3 = 0 if rotation_axes[2] == "x" else 1 if rotation_axes[2] == "y" else 2 if rotation_axes[2] == "z" else None
    v3 = numpy.array([0., 0., 0.])
    v3[i3] = 1.
    v3r = rotate_quat(v3, q)
    if ((i1 == 0) and (i2 == 2) and (i3 == 0)) or \
            ((i1 == 1) and (i2 == 0) and (i3 == 1)) or \
            ((i1 == 2) and (i2 == 1) and (i3 == 2)):
        e0 = numpy.arctan2(v3r[(i1 + 2) % 3], v3r[(i1 + 1) % 3])
        e1 = numpy.arccos(v3r[i1])
    elif ((i1 == 0) and (i2 == 2) and (i3 == 1)) or \
            ((i1 == 1) and (i2 == 0) and (i3 == 2)) or \
            ((i1 == 2) and (i2 == 1) and (i3 == 0)):
        e0 = numpy.arctan2(v3r[(i1 + 2) % 3], v3r[(i1 + 1) % 3])
        e1 = -numpy.arcsin(v3r[i1])
    elif ((i1 == 0) and (i2 == 1) and (i3 == 0)) or \
            ((i1 == 1) and (i2 == 2) and (i3 == 1)) or \
            ((i1 == 2) and (i2 == 0) and (i3 == 2)):
        e0 = numpy.arctan2(v3r[(i1 + 1) % 3], -v3r[(i1 + 2) % 3])
        e1 = numpy.arccos(v3r[i1])
    else:
        e0 = numpy.arctan2(-v3r[(i1 + 1) % 3], v3r[(i1 + 2) % 3])
        # The reference states this:
        # e1 = -numpy.arcsin(v3r[i1])
        # The tests only pass with the inverse sign, so I guess this is a typo.
        e1 = numpy.arcsin(v3r[i1])
    q1 = numpy.array([numpy.cos(e0 / 2.), 0., 0., 0.])
    q1[1 + i1] = numpy.sin(e0 / 2.)
    q2 = numpy.array([numpy.cos(e1 / 2.), 0., 0., 0.])
    q2[1 + i2] = numpy.sin(e1 / 2.)
    q12 = quat_mult(q1, q2)
    v3n = numpy.array([0., 0., 0.])
    v3n[(i3 + 1) % 3] = 1.
    v3n12 = quat_vec_mult(q12, v3n)
    v3nG = quat_vec_mult(q, v3n)
    e2_mag = numpy.arccos(dotproduct(v3n12, v3nG))
    vc = crossproduct(v3n12, v3nG)
    m = dotproduct(vc, v3r)
    e2 = numpy.sign(m) * e2_mag
    return numpy.array([e0, e1, e2])


def crossproduct(a, b):
    c = numpy.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])
    return c


def generate_map(pdbid, id):
    print(pdbid)

    runCommand('open {}'.format(pdbid))
    runCommand('addh')
    runCommand('molmap #{} 3 symmetry biomt'.format(id))
    runCommand('volume #{} save {}_model.mrc'.format(id, pdbid))
    return '{}_model.mrc'.format(pdbid)


def addNoise(noise_size, dim):
    noise_no_norm = abs(ifftn(fftshift(fftn(random((noise_size, noise_size, noise_size)))), [dim] * 3))
    noise = 1. + 0.1 * noise_no_norm / abs(noise_no_norm).max()
    return noise


# a = read_mrc('1bxn_chimera.mrc')

def addStructuralNoise(model, th=1.3):
    dx, dy, dz = model.shape
    data2 = model.copy()
    data2[data2 > th] = th

    data4 = (1 - data2 / th)

    data5 = numpy.ones((dx, dy, dz))
    data5[:, :] = data4
    a2 = numpy.zeros((dx, dy, dz))
    a2[:, :] = model
    # noise = addNoise(dx*9//10, dz)
    # noise = noise[:dx,:dy,:dz]
    noise = normal(0.94, 0.05, (dx, dy, dz))

    final = a2 + data5 * noise
    # final = numpy.ones((128,128,128))

    # final = data4

    # imshow(final[64,:,:],cmap='binary')
    # imshow(final.sum(axis=0))
    # show()
    print(median(noise[final > 1]))
    return final


def crop(data, factor):
    s = (data.shape[0] % factor) // 2
    d = (data.shape[0] % factor) % 2
    data = data[s:data.shape[0] - s - d, s:data.shape[0] - s - d, s:data.shape[0] - s - d]

    ds = data.shape[0] // factor
    image_size = data.shape[0]
    binned = data.reshape(ds, image_size // ds, ds, image_size // ds, ds, image_size // ds).mean(-1).mean(1).mean(-2)

    # fig,ax = subplots(1,2,figsize=(10,5))
    # ax[0].imshow(binned.sum(axis=0))
    ft = fftshift(fftn(data))
    x, y, z = numpy.array(ft.shape) // 2
    ff = factor
    ft = ft[int(x - x // ff):int(x + x // ff), int(y - y // ff):int(y + y // ff), int(z - z // ff):int(z + z // ff)]
    particle = abs(ifftn(fftshift(ft)))
    # ax[1].imshow(particle.sum(axis=0))
    # show()

    return binned

def createramp(sz,dim):
# Takes a list of dimensions sizes (sz) and a specfic dimension (dim) as input.
# Output returns a 2d matrix with a gradient of this dimensions towards the center
# in the specified direction. Adapted from dip_image matlab toolbox.
    if len(sz) >= dim and sz[dim] > 1:
        x = sz[dim] - 1

        # make a floating point gradient with true origin
        x = xp.arange(-x/2, x/2 + 1, 1)

        sz[dim] = 1
        xsz = xp.ones(len(sz),xp.int32) # np.ones() and np.int8
        xsz[dim] = len(x)
        x = xp.reshape(x, tuple(xsz)) # np.reshape()
    else:
        x = 0

    return xp.tile(x, sz)

def rr(x = 256, y = 256):
# Takes size for x and y dimension as input and returns a 2d matrix with a gradient of angular coordinates towards the
# center of the array. Function addapted from dip_image matlab toolbox.
    sz = [x, y]
    out = createramp(sz, 0)**2
    out = out + createramp(sz, 1)**2
    return out**0.5


def generate_model(outputFolder, modelID, listpdbs, SIZE, waterdensity=1., numberOfParticles=1000):
    if not os.path.exists(os.path.join(outputFolder, f'model_{modelID}')):
        os.mkdir(os.path.join(outputFolder, f'model_{modelID}'))

    # ['3cf3', '1s3x','1u6g','4b4t','1qvr','3h84','2cg9','3qm1','3gl1','3d2f','4d8q','1bxn']

    factor = 10
    skipped = 0
    models = []
    dims = []
    for n, pdbid in enumerate(listpdbs):
        if not os.path.exists('{}/{}_model.mrc'.format(outputFolder, pdbid)):
            fname = generate_map(pdbid, n + skipped)
        else:
            fname = '{}/{}_model.mrc'.format(outputFolder, pdbid)
            skipped += 1
        if 1 or not os.path.exists('{}/{}_model_binned.mrc'.format(outputFolder, pdbid)):
            m = read_mrc(fname)
            size = max(m.shape)
            m2 = numpy.zeros((size, size, size))
            dx, dy, dz = m.shape
            sx, ex = (size - dx) // 2, size - int(ceil((size - dx) / 2.))
            sy, ey = (size - dy) // 2, size - int(ceil((size - dy) / 2.))
            sz, ez = (size - dz) // 2, size - int(ceil((size - dz) / 2.))
            m2[sx:ex, sy:ey, sz:ez] = m
            m_r = crop(m2, factor)
            # print m_r.sum(), m2.sum(), m_r.max(), m2.max()
            # m_str = addStructuralNoise(m_r,m_r.shape[0])
            # print median(m_str[m_str > 1.]),median(m_r[m_r > 1.])
            convert_numpy_array3d_mrc(m_r, '{}/{}_model_binned.mrc'.format(outputFolder, pdbid))


        else:
            m_r = read_mrc('{}/{}_model_binned.mrc'.format(outputFolder, pdbid))

        models.append(m_r)
        dims.append(m_r.shape[0])

    X, Y, Z = 200, SIZE, SIZE
    cell = numpy.zeros((X, Y, Z))

    cell[:, :, :] = waterdensity

    mask = numpy.zeros_like(cell)

    #if addStructNoise:
    #    cell = addStructuralNoise(cell)
    volumes = []
    for i in range(len(models)):
        vol = read('{}_model_binned.mrc'.format(listpdbs[i]))
        volumes.append(vol)

    total = [0, ] * len(models)
    particleNr = 0
    outline = ''
    for i in range(numberOfParticles):
        if i % 200 == 0:
            print(i)
        a = numpy.random.randint(0, len(models))
        q = rand_quat()
        euler = euler_from_quat(q)
        euler *= 180. / pi

        mdl_vol = rotate(volumes[a], float(euler[0]), x=float(euler[1]), z2=float(euler[2]))
        mdl = vol2npy(mdl_vol)
        # mdl = rotate(models[a],randint(0,360),reshape=False)
        xx, yy, zz = mdl.shape

        for i in range(900):
            loc_x = numpy.random.randint(xx // 2 + 1, X - xx // 2 - 1)
            loc_y = numpy.random.randint(yy // 2 + 1, Y - yy // 2 - 1)
            loc_z = numpy.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            if mask[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
               loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
               loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2].sum() == 0:
                break
        if i > 898:
            continue
        # print ('place particle {} at: {},{},{}'.format(particleNr,loc_x,loc_y,loc_z))

        mask[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
        loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
        loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2] = 1

        cell[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
        loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
        loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2] = mdl

        if abs(512 - loc_y) < 256 and abs(512 - loc_z) < 256:
            outline += "{} {:4d} {:4d} {:4d} {:4.0f} {:4.0f} {:4.0f}\n".format(listpdbs[a], loc_x, loc_y - 256,
                                                                               loc_z - 256, euler[0], euler[1],
                                                                               euler[2])

        particleNr += 1

        total[a] += 1

    cell = addStructuralNoise(cell)

    grandcell = zeros((750, SIZE, SIZE))
    grandcell[275:475, :, :] = cell

    convert_numpy_array3d_mrc(grandcell, f'{outputFolder}/model_{modelID}/grandmodel_{modelID}.mrc')

    convert_numpy_array3d_mrc(cell[:, SIZE // 4:-SIZE // 4, SIZE // 4:-SIZE // 4],
                              f'{outputFolder}/model_{modelID}/grandmodel_{modelID}.mrc')

    outfile = open(f'{outputFolder}/model_{modelID}/particle_locations_model_{modelID}.txt', 'w')
    outfile.write(outline)

    return grandcell

def generate_projections(grandcell, angles, outputFolder='./', modelID=0, pixelSize=1e-9, voltage = 200e3,
                         sigmaDecayCtf=0.4, sphericalAberration=2.7E-3, multislice=None, msdz=5e-9, heightBox=1200,
                         amplitudeContrast=0.07, defocus=2, sigmaDecayCTF=0.4):

    from voltools import StaticVolume, Interpolations
    import cupy
    import numpy as xp

    SIZE = grandcell.shape[1]
    print('size of grandcell: ', SIZE)
    if grandcell.shape[0] >= heightBox:
        raise Exception('Your model is larger than than the box that we will rotate (heightBox parameter)')

    offset = (heightBox - grandcell.shape[0])//2

    volume = xp.zeros((heightBox, SIZE, SIZE), dtype=float32)
    volume[offset:-offset, :, :] = grandcell[:, :, :] # place the 750 long grandcell in a 1200 long volume to accommodate rotation

    gpu_volume = cupy.array(volume)
    d_vol = StaticVolume(gpu_volume, interpolation=Interpolations.FILT_BSPLINE) # turn into Voltools volume

    noisefree_projections = xp.zeros((len(angles), SIZE//2, SIZE//2), dtype=float64)

    #TODO allow for different defocus per tilt image. Now single defocus for all images
    ctf = calcCTF(defocus, xp.zeros((SIZE,SIZE)), pixelSize, voltage, sphericalAberration, sigmaDecayCTF, amplitudeContrast)

    for n, angle in enumerate(angles):

        print('simulating projection from tilt angle ',angle)

        rotated_volume = d_vol.transform(rotation=(0, angle, 0), rotation_order='szyx', rotation_units='deg') # if center = None: center = np.divide(volume.shape, 2, dtype=np.float32)

        if not multislice:
            projected_tilt_image = rotated_volume.sum(axis=0).get() # remove .get() if fully cupy
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*ctf))**2
            #TODO is abs()**2 the way to go for exit wave field?

        else:
            # raise Exception('Not implemented yet!!!')
            # Multislice wave propagation only needs to be done in a 512 by 512 area as that produces the final image.
            # calculate the number of slices based on the msdz

            zheight = heightBox * pixelSize # thickness of volume in nm???

            # Determine the number of slices
            if msdz >= zheight:
                n_slices = 1
                msdz = zheight
                print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
            elif msdz <= pixelSize:
                n_slices = heightBox
                msdz = pixelSize
                print('The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
            else:
                n_slices = int(xp.ceil(xp.around(zheight/msdz,3)))

            # Allocate space for multislice projection
            projected_potent_ms = xp.zeros((n_slices, SIZE, SIZE), dtype=complex)

            px_per_slice = int(xp.ceil(xp.around(heightBox / n_slices,3)))
            # PROJECTED POTENTIAL whithin slices (phase grating)
            for ii in range(n_slices):
                projected_potent_ms[ii, :, :] = rotated_volume[ii*px_per_slice : (ii+1)*px_per_slice - 1, :, :].mean(axis=0).get() # remove .get() if fully cupy


            # zheight of slices in nm for Fresnel propagator
            dzprop = px_per_slice * pixelSize
            # wavelength
            Lambda = wavelength_eV2m(voltage)
            # relative mass
            relmass = phys_const_dict["me"] + phys_const_dict["el"] * voltage / (phys_const_dict["c"]**2)
            # sigma_transfer
            sig_transfer = 2 * pi * relmass * phys_const_dict["el"] * Lambda / (phys_const_dict["h"]**2)
            # TRANSMISSON FUNCTION: the slice thickness is constant
            psi_t = xp.exp(1j * sig_transfer * projected_potent_ms * dzprop)


            # Determine Fresnel Propagator for slice

            # Get value of q in Fourier space because the Fresnel propagator needs them
            xwm = pixelSize * SIZE  # pixelsize for multislice * size sample
            q_true_pix_m = 1 / xwm
            q_m = rr(SIZE, SIZE) * q_true_pix_m  # frequencies in Fourier domain

            # FRESNEL PROPAGATOR
            P = xp.exp(-1j * pi * Lambda * (q_m**2) * dzprop)


            # PsiExit = multislice(psi_t, Nm, n, Lambda , q_m, dzprop); # function in matlab code
            # MULTISLICE
            psi_multislice = xp.zeros((SIZE, SIZE), dtype=complex) + 1 # should be complex datatype

            num_px_last_slice = heightBox % px_per_slice

            for ii in range(n_slices-min(1,num_px_last_slice)):
                #TODO ADD APPLICATION OF PYTOM.TOMPY
                psi_multislice = xp.fft.ifftn(xp.fft.fftn( psi_multislice * xp.squeeze(psi_t[ii, :, :]) )*P)

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                dzprop_end = num_px_last_slice * pixelSize
                psi_t[-1, :, :] = xp.exp(1j * sig_transfer * projected_potent_ms[-1, :, :] * dzprop_end)
                P_end = xp.exp(-1j * pi * Lambda * (q_m ** 2) * dzprop_end)
                psi_multislice = xp.fft.ifftn(xp.fft.fftn( psi_multislice * xp.squeeze(psi_t[-1, :, :]) )*P_end)

            # To make clear that we obtain the exit_wave here! Could be removed.
            psi_exit = psi_multislice

            # GET INTENSITIES IN IMAGE PLANE

            # intensity = modulo squared (exit wave)

            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(psi_exit) * ctf)) ** 2

            print('complex part of projection: ', xp.abs(psi_exit.imag).mean())

                # noisefree_projections[n] = gpu_proj[300:-300,300:-300]

        noisefree_projections[n] = projected_tilt_image[SIZE//4:-SIZE//4, SIZE//4:-SIZE//4]

        # =  create_projections(dens, range(-60, 61, 3))
        # noisefree_projections[n] /= noisefree_projections[n][:,140:160].mean()

    if not os.path.exists('{}/model_{}'.format(outputFolder, modelID)): os.mkdir(
        '{}/model_{}'.format(outputFolder, modelID))

    a = mrcfile.new('{}/model_{}/projections_grandmodel_{}.mrc'.format(outputFolder, modelID, modelID),
                    data=noisefree_projections.astype(xp.float32), overwrite=True) # if cupy noisefree_projections.get()
    a.close()
    del d_vol # + del volume?
    return noisefree_projections

def add_effects_microscope(outputFolder, modelID, grandcell, noisefree_projections, defocus, pixelSize, angles, SNR):

    for n in range(len(noisefree_projections)):
        noisefree_projections[n] = rot90(noisefree_projections[n], 3)

    
    diff = (grandcell.shape[-1] - noisefree_projections.shape[-1])

    simtomo = simutomo(grandcell[:, diff // 2:-diff // 2, diff // 2:-diff // 2], defocus, pixelSize,
                       angles, SNR, 0, 0, 0, outputFolder=outputFolder, modelID=modelID)

    if simtomo:
        convert_numpy_array3d_mrc(simtomo, f"{outputFolder}/simutomo_{modelID}.mrc")

def reconstruct_tomogram(prefix, suffix, start_idx, end_idx, volsize, angles, outputFolder, modelID, weighting=-1):
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()

    for i in range(start_idx, end_idx+1):
        p = Projection(prefix+str(i)+suffix,tiltAngle=angles[i-1])
        projections.append(p)

    outputname = os.path.join(outputFolder, f'tomogram_model_{modelID}.em') 

    vol = projections.reconstructVolume( dims=vol_size, reconstructionPosition=[0,0,0], binning=1, applyWeighting=weighting)
    vol.write(outputname)
    os.system('em2mrc.py -f {} -t {}'.format(outputname, os.path.dirname(outputname)) ) 
    os.system(f'rm {outputname}')

if __name__ == '__main__':

    import configparser
    from pytom.gui.guiFunctions import loadstar, datatype
    from pytom.gui.mrcOperations import read_mrc

    config = configparser.ConfigParser()
    try:
        config.read_file(open('simulation.conf'))
    except Exception as e:
        print(e)
        raise Exception('Could not open config file.')

    print(config.sections())

    try:
        outputFolder = config['General']['OutputFolder']
        modelID = int(config['General']['ModelID'])

        # meta file
        metadata = loadstar(config['General']['MetaFile'], dtype=datatype)
        angles = metadata['TiltAngle'] # specified in degrees
        defocus = metadata['DefocusU'][0] * 1E-6 # defocus in um
        voltage = metadata['Voltage'][0] * 1E3 # voltage in keV
        sphericalAberration = metadata['SphericalAberration'][0] * 1E-3 # spherical aberration in mm
        amplitudeContrast = metadata['AmplitudeContrast'][0] # fraction of amplitude contrast
        pixelSize = metadata['PixelSpacing'][0] * 1E-9 # pixel size in nm
    except Exception as e:
        print(e)
        raise Exception('No general parameters specified in the config file.')

    if 'GenerateModel' in config.sections():
        generateModel = True
        try:
            listpdbs = eval(config['GenerateModel']['Models'])
            SIZE = int(config['GenerateModel']['Size'])
            # instead: modelSize = int(config['GenerateModel']['Size'])
            # iceThickness = int(config['GenerateModel']['IceThickness'])
            waterdensity = float(config['GenerateModel']['WaterDensity'])
            numberOfParticles = int(config['GenerateModel']['NumberOfParticles'])
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters.')
    else:
        generateModel = False

    if 'GenerateProjections' in config.sections():
        generateProjections = True
        try:
            multislice = bool(config['GenerateProjections']['MultiSlice'])
            if multislice:
                msdz = float(config['GenerateProjections']['MultiSliceSize']) * 1E-9 # multislice step size in nm
            heightBox = int(config['GenerateProjections']['HeightBox'])
            sigmaDecayCTF = float(config['GenerateProjections']['SigmaDecayCTF'])
            if not(generateModel):
                grandmodelFile = config['GenerateProjections']['GrandModelFile']
        except Exception as e:
            print(e)
            raise Exception('Missing generate projection parameters.')
    else:
        generateProjections = False

    if 'AddEffectsMicroscope' in config.sections():
        addEffectsMicroscope = True
        try:
            SNR = float(config['AddEffectsMicroscope']['SNR'])
            if not(generateProjections):
                projectionsFile = config['AddeffectsMicroscope']['ProjectionsFile']
        except Exception as e:
            print(e)
            raise Exception('Missing add effects microscope parameters.')
    else:
        addEffectsMicroscope = False

    if 'ReconstructTomogram' in config.sections():
        reconstructTomogram = True
        try:
            prefix = config['ReconstructTomogram']['Prefix']
            suffix = config['ReconstructTomogram']['Suffix']
            start = int(config['ReconstructTomogram']['StartIdx'])
            end = int(config['ReconstructTomogram']['EndIdx'])
            weighting = int(config['ReconstructTomogram']['Weighting'])
            sizeRecon = int(config['ReconstructTomogram']['SizeRecon'])
        except Exception as e:
            print(e)
            raise Exception('Missing reconstruct tomogram parameters.')
    else:
        reconstructTomogram = False

    # Generate or read a grand model
    if generateModel:
        print('Generating model')
        grandcell = generate_model(outputFolder, modelID, listpdbs, SIZE, waterdensity, numberOfParticles)
    elif grandmodelFile:
        grandcell = read_mrc(grandmodelFile)
    else:
        grandcell = None

    # Generate or read noise-free projections
    if generateProjections and not (grandcell is None):
        print('Simulating projections')
        noisefree_projections = generate_projections(grandcell, angles, outputFolder=outputFolder, modelID=modelID,
                                                     pixelSize=pixelSize, voltage=voltage, multislice=multislice, msdz=msdz)
    elif projectionsFile:
        noisefree_projections = read_mrc(projectionsFile)
    else:
        noisefree_projections = None

    # Add effects of the microscope to the noise-free projections
    if addEffectsMicroscope and not (grandcell is None) and not (noisefree_projections is None):
        print('Adding effects of microscope')
        add_effects_microscope(outputFolder, modelID, grandcell, noisefree_projections, defocus, pixelSize, angles, SNR)

    # Reconstruct tomogram
    if reconstructTomogram:
        print('Reconstructing tomogram')
        # prefix  = os.path.join(outputFolder, f'model_{modelID}/noisyProjections/simulated_proj_model{modelID}_')
        # suffix  = '.mrc'
        vol_size = [sizeRecon,sizeRecon,sizeRecon]
        reconstruct_tomogram(prefix, suffix, start, end, vol_size, angles, outputFolder, modelID, weighting=weighting)
