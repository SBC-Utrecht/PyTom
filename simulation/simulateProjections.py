"""
Gijs' tilt series simulation - used for Ilja's contest paper
"""
# import numpy
# import sys
# from chimera import runCommand
# from pytom_volume import read
# from pytom_numpy import vol2npy
# import os
# import mrcfile
# import scipy
# from nufft.reconstruction import fourier_2d1d_iter_reconstruct
# from pytom.tools.script_helper import ScriptHelper, ScriptOption
# from pytom.tools.parse_script_options import parse_script_options
# from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
# from pytom.gui.mrcOperations import read_mrc
# from scipy.ndimage import rotate

# IO related modules
from pytom.gui.mrcOperations import *
from pytom.gui.guiFunctions import loadstar, datatype
import pytom.tompy.io
import configparser
import logging
import os
import mrcfile
import datetime

# Plotting
import matplotlib
from pylab import *
matplotlib.use('Qt5Agg')

# math
from scipy.ndimage import gaussian_filter
from pytom.reconstruction.reconstructionStructures import *
from pytom.basic.files import *
import numpy as xp
import random

phys_const_dict = {
    # Dictionary of physical constants rquired for calculation.
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

class ConfigLogger(object):
    '''Facilitates writing the conf file to a .log file in the outputFolder for reference of settings.'''
    def __init__(self, log):
        self.__log = log

    def __call__(self, config):
        self.__log.info("Config:")
        config.write(self)

    def write(self, data):
        # stripping the data makes the output nicer and avoids empty lines
        line = data.strip()
        self.__log.info(line)

def wavelength_eV2m(V):
    # OLD FUNCTION
    # h / sqrt( 2 * me * el * ev * 1)
    # return 6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * ev * 1) ** 0.5

    # NEW FUNCTION
    # function calculates the electron wavelength given a certain accelaration voltage V
    h = phys_const_dict["h"]
    e = phys_const_dict["el"]
    m = phys_const_dict["me"]
    c = phys_const_dict["c"]

    # matlab original: lambda = h/sqrt(e*V*m*(e/m*V/c^2 + 2 ));
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
        [x, y, z] = xp.meshgrid(xp.arange(-xp.floor(s1 / 2), -xp.floor(s1 / 2) + s1 - 1),
                                xp.arange(-xp.floor(s2 / 2), -xp.floor(s2 / 2) + s2 - 1),
                                xp.arange(-xp.floor(s3 / 2), -xp.floor(s3 / 2) + s3 - 1))

        r = xp.sqrt(x ** 2 + y ** 2 + z ** 2)
        lowp = (r <= hi)
        highp = (r >= low)
        image = xp.fft.fftshift(xp.fft.fftn(image))
        image = scf * xp.real(xp.fft.ifftn(xp.fft.fftshift(highp * lowp * image)))

    else:
        image = xp.fft.fftshift(xp.fft.fftn(image));
        if low > 0:
            image = xp.real(xp.fft.ifftn(xp.fft.fftshift(spheremask(image, hi, smooth) - spheremask(image, low, smooth))))
        else:
            image = xp.real(xp.fft.ifftn(xp.fft.fftshift(spheremask(image, hi, smooth))))

    return image


def tom_error(a, m=0, v=1):
    if len(a.shape) == 1:
        s1 = a.shape[0]
        c = a + [xp.sqrt(v) * xp.random.randn(s1) + m]
    elif len(a.shape) == 2:
        s1, s2 = a.shape
        c = a + [xp.sqrt(v) * xp.random.randn(s1, s2) + m]
    elif len(a.shape) == 3:
        s1, s2, s3 = a.shape
        c = a + [xp.sqrt(v) * xp.random.randn(s1, s2, s3) + m]
    return c


def spheremask(vol, radius, smooth=0):
    shape = vol.shape

    mask = xp.ones_like(vol).astype(xp.float32)
    if len(shape) == 2:
        x, y = xp.meshgrid(xp.arange(-shape[0] / 2, shape[0] / 2, 1), xp.arange(-shape[1] / 2, shape[1] / 2))
        r = xp.sqrt(x ** 2 + y ** 2)
        # imshow(mask[:,:])
        # show()
    if len(shape) == 3:
        x, y, z = xp.meshgrid(xp.arange(-shape[1] / 2, shape[1] / 2, 1), xp.arange(-shape[0] / 2, shape[0] / 2),
                           xp.arange(-shape[2] / 2, shape[2] / 2, 1))
        r = xp.sqrt(x ** 2 + y ** 2 + z ** 2)
        # imshow(mask[mask.shape[0]//2,:,:])
        # show()

    mask[r > radius] = 0

    if smooth > 0.01:
        mask = gaussian_filter(mask, smooth)
    return vol * mask


def tom_dev(image):
    return image.mean(), image.max(), image.min(), image.std(), image.std() ** 2

def create_ctf(Dz, vol, pix_size, voltage, Cs, sigma):
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

    # TODO remove: deprecated function!

    Lambda = wavelength_eV2m(voltage)

    Ny = 1 / (2 * pix_size)

    # print('in create_ctf the pixelsize = ', pix_size)
    # print('in create_ctf the Ny freq = ', Ny)

    if len(vol.shape) > 2:
        R, Y, Z = xp.meshgrid(xp.arange(-Ny, Ny, 2 * Ny / vol.shape[0]), xp.arange(-Ny, Ny, 2 * Ny / vol.shape[1]),
                           xp.arange(-Ny, Ny, 2 * Ny / vol.shape[2]))
        r = xp.sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[1])))
        r = xp.sqrt(R ** 2 + Y ** 2)

    # print(r)

    phase = xp.sin(pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2) ))
    amplitude = xp.cos(pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2) ))

    if sigma:
        phase = phase * xp.exp(-(r / (sigma * Ny)) ** 2)
        amplitude = amplitude * xp.exp(-(r / (sigma * Ny)) ** 2)

    return phase, amplitude

def calcCTF(Dz, vol, pix_size, voltage=200E3, Cs=2.7E-3, sigma_decay_ctf=0.4, amplitude_contrast=0.07):

    # calcCTF does not call create_ctf any longer!

    Lambda = wavelength_eV2m(voltage)

    Ny = 1 / (2 * pix_size)

    # print('in create_ctf the pixelsize = ', pix_size)
    # print('in create_ctf the Ny freq = ', Ny)

    if len(vol.shape) > 2:
        R, Y, Z = xp.meshgrid(xp.arange(-Ny, Ny, 2 * Ny / vol.shape[0]), xp.arange(-Ny, Ny, 2 * Ny / vol.shape[1]),
                              xp.arange(-Ny, Ny, 2 * Ny / vol.shape[2]))
        r = xp.sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[1])))
        r = xp.sqrt(R ** 2 + Y ** 2)

    # print(r)

    phase = xp.sin(pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2)))
    amplitude = xp.cos(pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2)))

    ctf = amplitude * amplitude_contrast - 1j * phase * (1-amplitude_contrast)

    if sigma_decay_ctf:
        ctf = ctf * xp.exp(-(r / (sigma_decay_ctf * Ny)) ** 2)

    # ctf1 = -ctf1
    # amplitude = -amplitude
    # ctf = (1-amplitude_contrast) * ctf1 + amplitude_contrast * amplitude
    return ctf

def simutomo(noisefree_projections, defocus, pixelsize, SNR, outputFolder='./', modelID=0,
             voltage=200e3, amplitude_contrast=0.07, Cs=2.7E-3, sigma_decay_ctf=0.4):
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

    # Assuming images are same size in x and y
    size = noisefree_projections.shape[0]
    n_images = noisefree_projections.shape[2]

    print('number of projections is ', n_images)

    # calculate CTF
    ctf = calcCTF(defocus, xp.zeros((size,size)), pixelsize, voltage=voltage, Cs=Cs, sigma_decay_ctf=sigma_decay_ctf,
                  amplitude_contrast=amplitude_contrast)

    sigma = 0
    # pre-calculate average stdv for a whole tilt series
    for n in range(n_images):
        proj = noisefree_projections[:,:,n]
        [mv, mn, mx, stv, tmp] = tom_dev(proj)
        sigma = sigma + tmp

    sigma = sigma / n_images

    print('average sigma over all projections is ', sigma)

    # take care of signal reduction due to ctf
    # ratio of white noise with and without mult by ctf
    # tom_dev() returns single values, no arrays
    [mv, mn, mx, stv, corrfac] = tom_dev(xp.real(xp.fft.ifftn(xp.fft.fftshift(
        ctf * spheremask(xp.fft.fftshift(xp.fft.fftn(tom_error(xp.zeros((size, size)), 0, 1.)[0, :, :])), 10, 0)))))
    [mv, mn, mx, stv, whitenoise] = tom_dev(xp.real(xp.fft.ifftn(
        xp.fft.fftshift(spheremask(xp.fft.fftshift(xp.fft.fftn(tom_error(xp.zeros((size, size)), 0, 1.)[0, :, :])), 10, 0)))))
    corrfac = whitenoise / corrfac
    sigma = sigma * corrfac

    print('sigma after reduction by ctf is ', sigma)

    # generate weighting function for reconstruction # WHY GENERATE WEIGHTING???
    # s1, s2, s3 = rotvol.shape
    [x, y] = xp.meshgrid(xp.arange(-size / 2, size / 2), xp.arange(-size / 2, size / 2))
    r = xp.sqrt(x ** 2 + y ** 2)
    mask = xp.ones((size, size))
    mask[r >= (size // 2) - 1] = 0

    if not os.path.exists(f'{outputFolder}/model_{modelID}/noisyProjections/'):
        os.mkdir(f'{outputFolder}/model_{modelID}/noisyProjections/')

    projections = xp.zeros((size, size, n_images), dtype=xp.float32)

    for n in range(n_images):
        proj = noisefree_projections[:,:,n]

        # ctf dependent contribution of noise
        # flux = 80 # electrons per square A
        # flux_per_tilt = flux / 60 # electrons per square A per tilt
        # flux_per_pixel = flux_per_tilt * 100
        # projFlux = proj * flux_per_pixel
        # noisy = xp.random.poisson(lam=(projFlux/SNR))
        proj_fft = xp.fft.fftn(xp.fft.ifftshift(proj))
        noisy = tom_error(proj_fft, 0, 0.5 / SNR * sigma)[0, :, :] # 0.5 / SNR * sigma

        print('noise added is samples from a standard normal distribution multiplied by sqrt( 0.5/SNR * sigma) = ',
              xp.sqrt(0.5/SNR*sigma))

        # noisy = noisy[nx/2-sx/2:nx/2+sx/2,nx/2-sx/2:nx/2+sx/2]

        # ctf independent contribution of noise
        bb = tom_error(xp.zeros(proj.shape), 0, 0.5 / SNR * sigma)[0, :, :] # 0.5 / SNR * sigma
        bg = tom_bandpass(bb, 0, 1, 0.2 * size)

        plot = False
        if plot:
            fig, ax = subplots(1, 3, figsize=(12, 4))
            ax[0].imshow(proj)
            ax[1].imshow(noisy)
            ax[2].imshow(bg)
            show()

        # add both contributions (ratio 1:1) in Fourier space
        tmp = noisy + xp.fft.fftn(xp.fft.ifftshift(bg)) # = fftshift(fftn(noisy)) + fftshift(fftn(bg))
        tmp = tmp * xp.fft.ifftshift(mask)
        noisy = xp.real(xp.fft.fftshift(xp.fft.ifftn(tmp)))

        out = mrcfile.new(
            f'{outputFolder}/model_{modelID}/noisyProjections/simulated_proj_model{modelID}_{n+1}.mrc',
            noisy.astype(xp.float32), overwrite=True)
        out.close()

        plot = False
        if plot:
            fig, ax = subplots(1, 2, figsize=(8, 4))
            ax[0].imshow(bg)
            ax[1].imshow(noisy)
            show()

        projections[:,:,n] = noisy

    # sys.exit()

    # simtomo =  fourier_2d1d_iter_reconstruct(projections, list(range(tiltrange[0],tiltrange[1]+1,tiltincr)), iter)
    # simtomo =  fourier_2d1d_iter_reconstruct(nfp, list(range(tiltrange[0],tiltrange[1]+1,tiltincr)), iter)
    # normalize to stdv
    # [mnv, maxv, minv, stv, dummy] = tom_dev(simtomo)
    # simtomo = (simtomo-mnv)/stv
    # return simtomo
    return projections


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
    return xp.array([w, x, y, z])


def quat_vec_mult(q, v):
    r"""
    Return the product of a quaternion and a vector

    Args:
       :q (array): Length-4 array [:math:`w`, :math:`x`, :math:`y`, :math:`z`] that represents the quaternion

       :v (array): Length-3 array [:math:`v_x`, :math:`v_y`, :math:`v_z`] that represents the vector
    """
    q2 = xp.array([0., v[0], v[1], v[2]])
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
    return xp.sum(list((a * b) for a, b in zip(v1, v2)))


def rand_quat():
    r"""
    Obtain a uniform random rotation in quaternion representation ([Shoemake1992]_ pages 129f)
    """
    x0, x1, x2 = xp.random.random(3)
    theta1 = 2. * xp.pi * x1
    theta2 = 2. * xp.pi * x2
    s1 = xp.sin(theta1)
    s2 = xp.sin(theta2)
    c1 = xp.cos(theta1)
    c2 = xp.cos(theta2)
    r1 = xp.sqrt(1 - x0)
    r2 = xp.sqrt(x0)
    q = xp.array([s1 * r1, c1 * r1, s2 * r2, c2 * r2])
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
    v3 = xp.array([0., 0., 0.])
    v3[i3] = 1.
    v3r = rotate_quat(v3, q)
    if ((i1 == 0) and (i2 == 2) and (i3 == 0)) or \
            ((i1 == 1) and (i2 == 0) and (i3 == 1)) or \
            ((i1 == 2) and (i2 == 1) and (i3 == 2)):
        e0 = xp.arctan2(v3r[(i1 + 2) % 3], v3r[(i1 + 1) % 3])
        e1 = xp.arccos(v3r[i1])
    elif ((i1 == 0) and (i2 == 2) and (i3 == 1)) or \
            ((i1 == 1) and (i2 == 0) and (i3 == 2)) or \
            ((i1 == 2) and (i2 == 1) and (i3 == 0)):
        e0 = xp.arctan2(v3r[(i1 + 2) % 3], v3r[(i1 + 1) % 3])
        e1 = -xp.arcsin(v3r[i1])
    elif ((i1 == 0) and (i2 == 1) and (i3 == 0)) or \
            ((i1 == 1) and (i2 == 2) and (i3 == 1)) or \
            ((i1 == 2) and (i2 == 0) and (i3 == 2)):
        e0 = xp.arctan2(v3r[(i1 + 1) % 3], -v3r[(i1 + 2) % 3])
        e1 = xp.arccos(v3r[i1])
    else:
        e0 = xp.arctan2(-v3r[(i1 + 1) % 3], v3r[(i1 + 2) % 3])
        # The reference states this:
        # e1 = -xp.arcsin(v3r[i1])
        # The tests only pass with the inverse sign, so I guess this is a typo.
        e1 = xp.arcsin(v3r[i1])
    q1 = xp.array([xp.cos(e0 / 2.), 0., 0., 0.])
    q1[1 + i1] = xp.sin(e0 / 2.)
    q2 = xp.array([xp.cos(e1 / 2.), 0., 0., 0.])
    q2[1 + i2] = xp.sin(e1 / 2.)
    q12 = quat_mult(q1, q2)
    v3n = xp.array([0., 0., 0.])
    v3n[(i3 + 1) % 3] = 1.
    v3n12 = quat_vec_mult(q12, v3n)
    v3nG = quat_vec_mult(q, v3n)
    e2_mag = xp.arccos(dotproduct(v3n12, v3nG))
    vc = crossproduct(v3n12, v3nG)
    m = dotproduct(vc, v3r)
    e2 = xp.sign(m) * e2_mag
    return xp.array([e0, e1, e2])


def crossproduct(a, b):
    c = xp.array([a[1] * b[2] - a[2] * b[1],
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

def addStructuralNoise(model, water=.94, th=1.3):
    # What is the threshold value??

    dx, dy, dz = model.shape
    data2 = model.copy()

    # set all values larger than threshold to the threshold value
    data2[data2 > th] = th

    # values in data2 equaling th will be zero, others will be larger than zero
    # zero values will become 1
    data4 = (1 - data2 / th)

    data5 = xp.ones((dx, dy, dz))
    data5[:, :] = data4
    a2 = xp.zeros((dx, dy, dz))
    a2[:, :] = model
    # noise = addNoise(dx*9//10, dz)
    # noise = noise[:dx,:dy,:dz]
    noise = xp.normal(water, 0.05, (dx, dy, dz))

    final = a2 + data5 * noise
    # final = xp.ones((128,128,128))

    # final = data4

    # imshow(final[64,:,:],cmap='binary')
    # imshow(final.sum(axis=0))
    # show()
    print(xp.median(noise[final > th]))
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
    ft = xp.fft.fftshift(xp.fft.fftn(data))
    x, y, z = xp.array(ft.shape) // 2
    ff = factor
    ft = ft[int(x - x // ff):int(x + x // ff), int(y - y // ff):int(y + y // ff), int(z - z // ff):int(z + z // ff)]
    particle = abs(xp.fft.ifftn(xp.fft.fftshift(ft)))
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


def generate_model(particleFolder, outputFolder, modelID, listpdbs, size=1024, thickness=200, waterdensity=.94,
                   proteindensity=1.3, numberOfParticles=1000):

    from voltools import transform

    dims = []

    #================================== binning

    # factor = 10
    # skipped = 0
    # models = []
    # for n, pdbid in enumerate(listpdbs):
    #     if not os.path.exists('{}/{}_model.mrc'.format(outputFolder, pdbid)):
    #         fname = generate_map(pdbid, n + skipped)
    #     else:
    #         fname = '{}/{}_model.mrc'.format(outputFolder, pdbid)
    #         skipped += 1
    #     if 1 or not os.path.exists('{}/{}_model_binned.mrc'.format(outputFolder, pdbid)):
    #         m = read_mrc(fname)
    #         size = max(m.shape)
    #         m2 = xp.zeros((size, size, size))
    #         dx, dy, dz = m.shape
    #         sx, ex = (size - dx) // 2, size - int(ceil((size - dx) / 2.))
    #         sy, ey = (size - dy) // 2, size - int(ceil((size - dy) / 2.))
    #         sz, ez = (size - dz) // 2, size - int(ceil((size - dz) / 2.))
    #         m2[sx:ex, sy:ey, sz:ez] = m
    #         m_r = crop(m2, factor)
    #         # print m_r.sum(), m2.sum(), m_r.max(), m2.max()
    #         # m_str = addStructuralNoise(m_r,m_r.shape[0])
    #         # print median(m_str[m_str > 1.]),median(m_r[m_r > 1.])
    #         convert_xp_array3d_mrc(m_r, '{}/{}_model_binned.mrc'.format(outputFolder, pdbid))
    #
    #
    #     else:
    #         m_r = read_mrc('{}/{}_model_binned.mrc'.format(outputFolder, pdbid))
    #
    #     models.append(m_r)
    #     dims.append(m_r.shape[0])

    #==================================

    X, Y, Z = size, size, thickness
    cell = xp.zeros((X, Y, Z))

    mask = xp.zeros_like(cell)

    volumes = []
    for i in range(len(listpdbs)):
        try:
            vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{listpdbs[i]}_model_binned.mrc')
            dx,dy,dz = vol.shape
            vol2 =xp.zeros((dx*2,dy*2, dz*2))
            # print(dx, dy, dz, vol2.shape)
            vol2[dx//2:-dx//2,dy//2:-dy//2,dz//2:-dz//2] = vol
            volumes.append(vol2)
            dims.append(vol2.shape[0])
        except Exception as e:
            print(e)
            raise Exception('Could not open pdb ', listpdbs[i])


    total = [0, ] * len(volumes)
    particleNr = 0
    outline = ''
    for i in range(numberOfParticles):
        if i % 200 == 0:
            print(i)
        a = xp.random.randint(0, len(volumes))
        q = rand_quat()
        euler = euler_from_quat(q)
        euler *= 180. / pi
        mdl = transform(volumes[a], rotation=(float(euler[0]), float(euler[1]), float(euler[2])), rotation_order='szyx', interpolation='filt_bspline', device='cpu')

        # mdl_vol = rotate(volumes[a], float(euler[0]), x=float(euler[1]), z2=float(euler[2]))
        # mdl = vol2npy(mdl_vol)
        # mdl = rotate(models[a],randint(0,360),reshape=False)

        xx, yy, zz = mdl.shape

        for i in range(900):
            loc_x = xp.random.randint(xx // 2 + 1, X - xx // 2 - 1)
            loc_y = xp.random.randint(yy // 2 + 1, Y - yy // 2 - 1)
            loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            if mask[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
               loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
               loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2].sum() == 0:
                break
        if i > 898:
            continue
        # print ('place particle {} at: {},{},{}'.format(particleNr,loc_x,loc_y,loc_z))

        mask[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
        loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
        loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2] = (mdl > 0.01).astype(xp.float32)

        cell[loc_x - dims[a] // 2:loc_x + dims[a] // 2 + dims[a] % 2,
        loc_y - dims[a] // 2:loc_y + dims[a] // 2 + dims[a] % 2,
        loc_z - dims[a] // 2:loc_z + dims[a] // 2 + dims[a] % 2] = mdl

        if xp.abs(size - loc_x) < size//2 and xp.abs(size - loc_y) < size//2:
            outline += "{} {:4d} {:4d} {:4d} {:4.0f} {:4.0f} {:4.0f}\n".format(listpdbs[a], loc_x -size//2, loc_y - size//2,
                                                                               loc_z, euler[0], euler[1],
                                                                               euler[2])

        particleNr += 1

        total[a] += 1

    print(outline)

    # Noise free model
    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/model_{modelID}_noisefree.mrc', cell)

    # Add water desnity with structural noise
    cell = addStructuralNoise(cell, water=waterdensity, th=proteindensity)

    # Grand model
    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/grandmodel_{modelID}.mrc', cell)

    # Cropped version
    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/model_{modelID}_cropped.mrc',
                         cell[size // 4:-size // 4, size // 4:-size // 4, :])

    # Save particle locations
    outfile = open(f'{outputFolder}/model_{modelID}/particle_locations_model_{modelID}.txt', 'w')
    try:
        outfile.write(outline)
    except Exception as e:
        print(e)
        raise Exception('Could not write to outfile.')
    finally:
        outfile.close()

    return


def generate_projections(angles, outputFolder='./', modelID=0, pixelSize=1e-9, voltage = 200e3,
                         sphericalAberration=2.7E-3, multislice=None, msdz=5e-9, amplitudeContrast=0.07, defocus=2,
                         sigmaDecayCTF=0.4):

    grandcell = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{0}.mrc')

    SIZE = grandcell.shape[0]
    heightBox = grandcell.shape[2]

    noisefree_projections = xp.zeros((SIZE//2, SIZE//2, len(angles)), dtype=float32)

    #TODO allow for different defocus per tilt image. Now single defocus for all images
    ctf = calcCTF(defocus, xp.zeros((SIZE,SIZE)), pixelSize, voltage=voltage, Cs=sphericalAberration,
                  sigma_decay_ctf=sigmaDecayCTF, amplitude_contrast=amplitudeContrast)

    plot = False
    if plot:
        fig, ax = subplots(1, 2, figsize=(8, 4))
        ax[0].imshow(ctf.real)
        ax[1].imshow(ctf.imag)
        show()

    for n, angle in enumerate(angles):


        filename = f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{int(angle)}.mrc'
        rotated_volume = pytom.tompy.io.read_mrc(filename)  #* 1E12

        print(rotated_volume.shape)

        #transform(gpu_volume, rotation=(0, angle, 0), rotation_order='szyx', interpolation='linear', device='cpu')
        # rotated_volume = d_vol.transform(rotation=(0, angle, 0), rotation_order='szyx', rotation_units='deg') # if center = None: center = np.divide(volume.shape, 2, dtype=np.float32)

        if not multislice:
            print('simulating projection (without ms) from tilt angle ', angle)
            projected_tilt_image = rotated_volume.sum(axis=2)#.get() # remove .get() if fully cupy
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*ctf))**2
            #TODO is abs()**2 the way to go for exit wave field?

        else:
            print('simulating projection (with ms) from tilt angle ', angle)

            zheight = heightBox * pixelSize # thickness of volume in nm???

            # Determine the number of slices
            if msdz > zheight:
                n_slices = 1
                msdz = zheight
                print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
            elif msdz < pixelSize:
                n_slices = heightBox
                msdz = pixelSize
                print('The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
            else:
                n_slices = int(xp.ceil(xp.around(zheight/msdz,3)))

            print('number of slices: ', n_slices)
            # Allocate space for multislice projection
            projected_potent_ms = xp.zeros((SIZE, SIZE, n_slices), dtype=complex)

            px_per_slice = int(xp.ceil(xp.around(heightBox / n_slices,3))) # CORRECT
            # PROJECTED POTENTIAL whithin slices (phase grating)
            for ii in range(n_slices):
                projected_potent_ms[:, :, ii] = rotated_volume[:,:, ii*px_per_slice : (ii+1)*px_per_slice].mean(axis=2) #.get() # remove .get() if fully cupy

            # zheight of slices in nm for Fresnel propagator
            dzprop = px_per_slice * pixelSize # CORRECT
            # wavelength
            Lambda = wavelength_eV2m(voltage)
            # relative mass
            relmass = phys_const_dict["me"] + phys_const_dict["el"] * voltage / (phys_const_dict["c"]**2)
            # sigma_transfer
            sig_transfer = 2 * pi * relmass * phys_const_dict["el"] * Lambda / (phys_const_dict["h"]**2) # CORRECT

            # TRANSMISSON FUNCTION: the slice thickness is constant
            psi_t = xp.exp(1j * sig_transfer * projected_potent_ms * dzprop)


            # Get value of q in Fourier space because the Fresnel propagator needs them
            xwm = pixelSize * (SIZE)  # pixelsize for multislice * size sample
            q_true_pix_m = 1 / xwm
            q_m = rr(SIZE, SIZE) * q_true_pix_m  # frequencies in Fourier domain

            # FRESNEL PROPAGATOR
            P = xp.exp(-1j * pi * Lambda * (q_m**2) * dzprop) # CORRECT

            # MULTISLICE
            psi_multislice = xp.zeros((SIZE, SIZE), dtype=complex) + 1 # should be complex datatype

            num_px_last_slice = heightBox % px_per_slice # CORRECT

            for ii in range(n_slices-min(1,num_px_last_slice)):
                #TODO ADD APPLICATION OF PYTOM.TOMPY

                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )

                psi_multislice = xp.fft.fftshift( xp.fft.ifftn((waveField * xp.fft.ifftshift(P) )) )

                plot = False
                if (not (ii % 10) and ii) and plot:
                    fig, ax = subplots(1, 4, figsize=(16, 4))
                    ax[0].imshow(xp.abs(waveField) ** 2)
                    ax[1].imshow(xp.abs(P))
                    ax[2].imshow(abs(psi_multislice))
                    ax[3].imshow(abs(psi_t[:,:, ii]))
                    show()

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                dzprop_end = num_px_last_slice * pixelSize
                psi_t[:, :, -1] = xp.exp(1j * sig_transfer * projected_potent_ms[:, :, -1] * dzprop_end)
                P_end = xp.exp(-1j * pi * Lambda * (q_m ** 2) * dzprop_end)
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn( waveField * xp.fft.ifftshift(P_end) ) )

            # GET INTENSITIES IN IMAGE PLANE
            waveCTF = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
            projected_tilt_image = xp.abs(xp.fft.fftshift(xp.fft.ifftn(waveCTF))) ** 2

            plot = False
            if plot:
                fig, ax = subplots(1,3, figsize=(12, 4))
                ax[0].imshow(xp.abs(psi_multislice))
                ax[1].imshow(xp.abs(waveCTF))
                ax[2].imshow(xp.abs(projected_tilt_image))
                show()

        noisefree_projections[:,:,n] = projected_tilt_image[SIZE//4:-SIZE//4, SIZE//4:-SIZE//4]

        # outproj = f'{outputFolder}/model_{modelID}/projections_grandmodel_{modelID}_angle_{angle}.em'
        # pytom.tompy.io.write(outproj, noisefree_projections[:,:,n])

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections_noisefree_{modelID}.mrc',
                         noisefree_projections) # noisefree_projections.astype(xp.float32)?

    return

def add_effects_microscope(outputFolder, modelID, defocus, pixelsize, SNR, voltage, amplitude_contrast, Cs, sigma_decay_ctf):

    # TODO add try except for opening these files as I am not specifying paths
    noisefree_projections = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/projections_noisefree_{modelID}.mrc')

    projections = simutomo(noisefree_projections, defocus, pixelsize, SNR, outputFolder=outputFolder, modelID=modelID,
                           voltage=voltage, amplitude_contrast=amplitude_contrast, Cs=Cs, sigma_decay_ctf=sigma_decay_ctf)

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections_{modelID}.mrc', projections)

    return

def reconstruct_tomogram(prefix, suffix, start_idx, end_idx, vol_size, angles, outputFolder, modelID, weighting=-1):
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()

    for i in range(start_idx, end_idx+1):

        p = Projection(prefix+str(i)+suffix,tiltAngle=angles[i-1])
        projections.append(p)

    outputname = os.path.join(outputFolder, f'model_{modelID}/tomogram_model_{modelID}.em')

    vol = projections.reconstructVolume( dims=vol_size, reconstructionPosition=[0,0,0], binning=1, applyWeighting=weighting)
    vol.write(outputname)
    os.system('em2mrc.py -f {} -t {}'.format(outputname, os.path.dirname(outputname)) ) 
    os.system(f'rm {outputname}')

    return

if __name__ == '__main__':

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
        SEED = int(config['General']['Seed'])

        # meta file
        metadata = loadstar(config['General']['MetaFile'], dtype=datatype)
        angles = metadata['TiltAngle'] # specified in degrees
        defocus = metadata['DefocusU'][0] * 1E-6 # defocus in um
        voltage = metadata['Voltage'][0] * 1E3 # voltage in keV
        sphericalAberration = metadata['SphericalAberration'][0] * 1E-3 # spherical aberration in mm
        amplitudeContrast = metadata['AmplitudeContrast'][0] # fraction of amplitude contrast
        # pixelSize = metadata['PixelSpacing'][0] * 1E-9 # pixel size in nm
        pixelSize = 1E-9
        sigmaDecayCTF = float(config['General']['SigmaDecayCTF'])
    except Exception as e:
        print(e)
        raise Exception('No general parameters specified in the config file.')

    if 'GenerateModel' in config.sections():
        try:
            particleFolder = config['GenerateModel']['ParticleFolder']
            listpdbs = eval(config['GenerateModel']['Models'])
            size = int(config['GenerateModel']['Size'])
            thickness = int(config['GenerateModel']['Thickness'])
            waterdensity = float(config['GenerateModel']['WaterDensity'])
            proteindensity = float(config['GenerateModel']['ProteinDensity'])
            numberOfParticles = int(config['GenerateModel']['NumberOfParticles'])
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters.')

    if 'Rotation' in config.sections():
        try:
            heightBox = int(config['Rotation']['HeightBox'])
            nodes = int(config['Rotation']['Nodes'])
        except Exception as e:
            print(e)
            raise Exception('Missing rotation parameters.')

    if 'GenerateProjections' in config.sections():
        try:
            multislice = config['GenerateProjections'].getboolean('MultiSlice')
            if multislice:
                msdz = float(config['GenerateProjections']['MultiSliceSize']) * 1E-9 # multislice step size in nm
        except Exception as e:
            print(e)
            raise Exception('Missing generate projection parameters.')

    if 'AddEffectsMicroscope' in config.sections():
        try:
            SNR = float(config['AddEffectsMicroscope']['SNR'])
        except Exception as e:
            print(e)
            raise Exception('Missing add effects microscope parameters.')

    if 'ReconstructTomogram' in config.sections():
        try:
            # start = int(config['ReconstructTomogram']['StartIdx'])
            # end = int(config['ReconstructTomogram']['EndIdx'])
            weighting = int(config['ReconstructTomogram']['Weighting'])
            sizeRecon = int(config['ReconstructTomogram']['SizeRecon'])
        except Exception as e:
            print(e)
            raise Exception('Missing reconstruct tomogram parameters.')

    if not os.path.exists(os.path.join(outputFolder, f'model_{modelID}')):
        os.mkdir(os.path.join(outputFolder, f'model_{modelID}'))

    if os.path.exists(f'{outputFolder}/model_{modelID}/simulator.log'):
        os.remove(f'{outputFolder}/model_{modelID}/simulator.log')

    logging.basicConfig(filename='{}/model_{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(outputFolder, modelID,
                                                                                             date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    # Generate or read a grand model
    if 'GenerateModel' in config.sections():
        # SET SEED for random number generation
        xp.random.seed(SEED)
        random.seed(SEED)

        print('Generating model')
        generate_model(particleFolder, outputFolder, modelID, listpdbs, size=size, thickness=thickness,
                                   waterdensity=waterdensity, proteindensity=proteindensity, numberOfParticles=numberOfParticles)

    if 'Rotation' in config.sections():
        print('Rotating model with ', nodes, ' nodes')
        dir = f'{outputFolder}/model_{modelID}/rotations'
        if not (os.path.exists(dir)):
            os.mkdir(dir)
        filename = f'{dir}/rotated_volume_0.mrc'
        if os.path.exists(filename):
            print('rotated volume 0 exists, so remove it...')
            os.remove(filename)
        os.system(f'mpiexec -n {nodes} pytom rotateMPI.py')

    # Generate or read noise-free projections
    if 'GenerateProjections' in config.sections():
        print('Simulating projections')
        generate_projections(angles, outputFolder=outputFolder, modelID=modelID,pixelSize=pixelSize,
                             voltage=voltage, sphericalAberration=sphericalAberration,multislice=multislice, msdz=msdz,
                             amplitudeContrast=amplitudeContrast, defocus=defocus, sigmaDecayCTF=sigmaDecayCTF)

    # Add effects of the microscope to the noise-free projections
    if 'AddEffectsMicroscope' in config.sections():
        # SET SEED for random number generation
        xp.random.seed(SEED)
        random.seed(SEED)

        print('Adding effects of microscope')
        add_effects_microscope(outputFolder, modelID, defocus, pixelSize, SNR, voltage, amplitudeContrast,
                               sphericalAberration, sigmaDecayCTF)

    # Reconstruct tomogram
    if 'ReconstructTomogram' in config.sections():
        print('Reconstructing tomogram')
        prefix  = os.path.join(outputFolder, f'model_{modelID}/noisyProjections/simulated_proj_model{modelID}_')
        suffix  = '.mrc'
        vol_size = [sizeRecon, sizeRecon, sizeRecon]
        reconstruct_tomogram(prefix, suffix, 13, 47, vol_size, angles, outputFolder, modelID, weighting=weighting)
