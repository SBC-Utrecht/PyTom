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
import datetime
import sys
from tqdm import tqdm
import time

# Plotting
# import matplotlib
# from pylab import *
# matplotlib.use('Qt5Agg')

# math
from pytom.reconstruction.reconstructionStructures import *
from pytom.basic.files import *
import scipy.ndimage
import numpy as xp
import random

phys_const_dict = {
    # Dictionary of physical constants required for calculation.
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
        [x, y, z] = xp.meshgrid(xp.arange(-xp.floor(s1 / 2), -xp.floor(s1 / 2) + s1 -1 ),
                                xp.arange(-xp.floor(s2 / 2), -xp.floor(s2 / 2) + s2 -1 ),
                                xp.arange(-xp.floor(s3 / 2), -xp.floor(s3 / 2) + s3 -1 ))

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


def spheremask(vol, radius, smooth=0, ellipsoid=0):
    """
    Create a spherical mask in vol of size radius. Smooth option can gaussian blur the mask at its edges. The ellipsoid
    option will instead create a ellipsoid mask. Radius should in this case correspond to radius of the ellipsoid in the
    smallest dimension of vol. Radii in other dimensions will be fit accordingly.

    @param vol: Volume
    @type vol: 2d or 3d array
    @param radius: Radius of the sphere
    @type radius: float
    @param smooth:
    @type smooth:
    @param ellipsoid: Flag for ellipsoid mask, 0 by default
    @type ellipsoid: bool

    @return:
    @rtype:

    @author: Marten Chaillet
    """
    shape = vol.shape
    mask = xp.ones_like(vol).astype(xp.float32)

    if ellipsoid==0:
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
    elif ellipsoid:
        shortest_side = min(shape)
        if len(shape) == 2:
            x, y = xp.meshgrid(xp.arange(-shortest_side / 2, shortest_side / 2, shortest_side/shape[0]),
                               xp.arange(-shortest_side / 2, shortest_side / 2, shortest_side/shape[1]))
            r = xp.sqrt(x ** 2 + y ** 2)
            # imshow(mask[:,:])
            # show()
        if len(shape) == 3:
            x, y, z = xp.meshgrid(xp.arange(-shortest_side / 2, shortest_side / 2, shortest_side/shape[1]),
                                  xp.arange(-shortest_side / 2, shortest_side / 2, shortest_side/shape[0]),
                                  xp.arange(-shortest_side / 2, shortest_side / 2, shortest_side/shape[2]))
            r = xp.sqrt(x ** 2 + y ** 2 + z ** 2)
            # imshow(mask[mask.shape[0]//2,:,:])
            # show()
    else:
        print('Non-valid option for ellipoid')

    mask[r > radius] = 0

    if smooth > 0.01:
        # gaussian_filter(input,sigma)
        gaussMask = scipy.ndimage.gaussian_filter(mask, smooth)
        return vol * gaussMask
    else:
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

    phase = xp.sin(xp.pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2) ))
    amplitude = xp.cos(xp.pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2) ))

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
    ctf = xp.exp( -1j * xp.pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2)) )

    ctf.real = ctf.real * amplitude_contrast
    ctf.imag = ctf.imag * (1-amplitude_contrast)

    if sigma_decay_ctf:
        decay = xp.exp(-(r / (sigma_decay_ctf * Ny)) ** 2)
        ctf = ctf * decay

    return - ctf

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
    # maybe at this point there should be a '-' for the CTF like in the old code
    ctf = - calcCTF(defocus, xp.zeros((size,size)), pixelsize, voltage=voltage, Cs=Cs, sigma_decay_ctf=sigma_decay_ctf,
                    amplitude_contrast=amplitude_contrast)

    sigma = 0
    # pre-calculate average stdv for a whole tilt series
    for n in range(n_images):
        proj = noisefree_projections[:,:,n]
        [mv, mn, mx, stv, tmp] = tom_dev(proj)
        sigma = sigma + tmp

    sigma = sigma / n_images

    # take care of signal reduction due to ctf
    # ratio of white noise with and without mult by ctf
    # tom_dev() returns single values, no arrays
    masked_noise = spheremask(xp.fft.fftshift(xp.fft.fftn(tom_error(xp.zeros((size, size)), 0, 1.)[0, :, :])), 10, 0)
    [mv, mn, mx, stv, corrfac] = tom_dev(xp.real(xp.fft.ifftn(xp.fft.fftshift(ctf * masked_noise))))
    [mv, mn, mx, stv, whitenoise] = tom_dev(xp.real(xp.fft.ifftn(xp.fft.fftshift(masked_noise))))
    corrfac = whitenoise / corrfac
    sigma = sigma * corrfac

    # generate weighting function for reconstruction # WHY GENERATE WEIGHTING???
    # s1, s2, s3 = rotvol.shape
    [x, y] = xp.meshgrid(xp.arange(-size / 2, size / 2), xp.arange(-size / 2, size / 2))
    r = xp.sqrt(x ** 2 + y ** 2)
    mask = xp.ones((size, size))
    mask[r >= (size // 2) - 1] = 0

    if not os.path.exists(f'{outputFolder}/model_{modelID}/noisyProjections/'):
        os.mkdir(f'{outputFolder}/model_{modelID}/noisyProjections/')

    projections = xp.zeros((size, size, n_images), dtype=xp.float32)

    print('noise will be scaled to sqrt( 0.5/SNR * sigma) = ',
          xp.sqrt(0.5 / SNR * sigma))

    for n in range(n_images):
        projection = noisefree_projections[:,:,n]

        # POISSONIAN NOISE!
        # flux = 80 # electrons per square A
        # flux_per_tilt = flux / 60 # electrons per square A per tilt
        # flux_per_pixel = flux_per_tilt * 100
        # projFlux = proj * flux_per_pixel
        # noisy = xp.random.poisson(lam=(projFlux/SNR))
        # ADD AT LATER STAGE

        # CTF dependent noise
        preNoise = tom_error(xp.zeros(projection.shape), 0, 0.5 / SNR * sigma)[0, :, :] # 0.5 / SNR * sigma
        ctfNoise = xp.fft.fftn(preNoise) * xp.fft.ifftshift(ctf) # for preNoise the ifftshift will have no effect

        # ctf independent contribution of noise
        # tom_bandpass(im, lo hi, smooth) ---> hi = 7 gives good results
        bg = tom_bandpass(preNoise, 0, 1, smooth=0.2 * size)

        # add both contributions (ratio 1:1) in Fourier space
        tmp = xp.fft.fftn(xp.fft.ifftshift(projection)) + xp.fft.fftn(xp.fft.ifftshift(bg)) # + ctfNoise
        tmp = tmp * xp.fft.ifftshift(mask)
        noisy = xp.real(xp.fft.fftshift(xp.fft.ifftn(tmp)))
        noisy = noisy.astype(xp.float32)

        pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/noisyProjections/simulated_proj_{n+1}.mrc', noisy)

        projections[:, :, n] = noisy

    return projections

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
    noise = xp.random.normal(water, 0.01, (dx, dy, dz))

    final = a2 + data5 * noise
    # final = xp.ones((128,128,128))

    # final = data4

    # imshow(final[64,:,:],cmap='binary')
    # imshow(final.sum(axis=0))
    # show()
    # print(xp.median(noise[final > th]))
    return final


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
                   proteindensity=1.3, numberOfParticles=1000, placementSize=512, retries=5000):

    from voltools import transform

    # outputs
    X, Y, Z = size, size, thickness
    cell = xp.zeros((X, Y, Z))

    occupancy_bbox_mask = xp.zeros_like(cell)
    occupancy_accurate_mask = xp.zeros_like(cell)
    class_bbox_mask = xp.zeros_like(cell)
    class_accurate_mask = xp.zeros_like(cell)
    ground_truth_txt_file = ''

    # load pdb volumes and pad them
    volumes = []
    for i in range(len(listpdbs)):
        try:
            # First find the voxel size of the map...
            # rotate
            # then bin
            vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{listpdbs[i]}_model_binned.mrc')
            dx, dy, dz = vol.shape
            vol2 = xp.zeros((dx*2, dy*2, dz*2), dtype=xp.float32)
            vol2[dx//2:-dx//2, dy//2:-dy//2, dz//2:-dz//2] = vol
            volumes.append(vol2)
        except Exception as ee:
            print(ee)
            raise Exception('Could not open pdb ', listpdbs[i])

    # attributes
    number_of_classes = len(listpdbs)
    save_path = f'{outputFolder}/model_{modelID}'
    dims = [v.shape for v in volumes]
    particles_by_class = [0, ] * number_of_classes
    particle_nr = 1
    default_tries_left = retries
    skipped_particles = 0
    loc_x_start = loc_y_start = (size - placementSize) - (placementSize // 2)
    loc_x_end = loc_y_end = (size - placementSize) + (placementSize // 2)

    for _ in tqdm(range(numberOfParticles), desc='Placing particles'):

        # select random class
        cls_id = xp.random.randint(0, number_of_classes)

        # generate random rotation
        # to be properly random, use uniform sphere sampling
        # https://math.stackexchange.com/a/442423/72032
        # https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
        # http://corysimon.github.io/articles/uniformdistn-on-sphere/
        u = xp.random.uniform(0.0, 1.0, (2,))
        theta = xp.arccos(2 * u[0] - 1)
        phi = 2 * xp.pi * u[1]
        psi = xp.random.uniform(0.0, 2 * xp.pi)
        p_angles = xp.rad2deg([theta, phi, psi])

        # rotate particle
        try:
            rotated_particle = transform(volumes[cls_id], rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
        except Exception as e:
            print(e)
            print('Something went wrong while rotating?')
            continue

        # remove particle rotation artifacts
        # threshold = min(volumes[cls_id][volumes[cls_id] > 0]) / 10
        # the approach above doesn't work well with PDB particles, there are often voxels with values of ^-11
        threshold = 0.001
        rotated_particle[rotated_particle < threshold] = 0

        # thresholded rotated particle
        accurate_particle_occupancy = rotated_particle > 0

        # find random location for the particle
        xx, yy, zz = rotated_particle.shape
        tries_left = default_tries_left
        while tries_left > 0:
            loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
            loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
            loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            tries_left -= 1

            # calculate coordinates of bbox for the newly rotated particle
            bbox_x = [loc_x - dims[cls_id][0] // 2, loc_x + dims[cls_id][0] // 2 + dims[cls_id][0] % 2]
            bbox_y = [loc_y - dims[cls_id][1] // 2, loc_y + dims[cls_id][1] // 2 + dims[cls_id][1] % 2]
            bbox_z = [loc_z - dims[cls_id][2] // 2, loc_z + dims[cls_id][2] // 2 + dims[cls_id][2] % 2]

            # create masked occupancy mask
            masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]]
            masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy

            # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
            if masked_occupancy_mask.sum() == 0:
                break

        # however if still can't fit, ignore this particle (also adds variance in how many particles are actually put)
        if tries_left < 1:
            skipped_particles += 1
            continue

        # populate occupancy volumes
        occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
        occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * particle_nr

        # populate class masks
        class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
        class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * (cls_id + 1)

        # populate density volume
        cell[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += rotated_particle

        # update stats
        particle_nr += 1
        particles_by_class[cls_id] += 1

        # update text
        ground_truth_txt_file += f'{listpdbs[cls_id]} {loc_x:.4f} {loc_y - 256:.4f} {loc_z - 256:.4f} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # add water density and structural noise
    print('Adding water and structural noise')
    noisy_cell = addStructuralNoise(cell, water=waterdensity, th=proteindensity)

    # save grandmodels
    print('Saving grandmodels')
    pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree.mrc', cell)
    pytom.tompy.io.write(f'{save_path}/grandmodel.mrc', noisy_cell)

    # save class masks
    print('Saving class volumes')
    pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', class_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/class_mask.mrc', class_accurate_mask)

    # save occupancy masks
    print('Saving occupancy volumes')
    pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', occupancy_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/occupancy_mask.mrc', occupancy_accurate_mask)

    # save particle text file
    with open(f'{save_path}/particle_locations.txt', 'w') as f:
        f.write(ground_truth_txt_file)

    # reporting
    print(f'Total number of particles in the tomogram: {particle_nr - 1}\n'
          f'Skipped {skipped_particles} particles ({default_tries_left} random location retries)\n'
          f'Particles by class: {particles_by_class}')

    # # debug: inspect all generated volumes
    # import napari
    # with napari.gui_qt():
    #     viewer = napari.Viewer()
    #     viewer.add_image(occupancy_bbox_mask, name='occupancy bbox mask', interpolation='bicubic')
    #     viewer.add_image(occupancy_accurate_mask, name='occupancy accurate mask', interpolation='bicubic')
    #     viewer.add_image(class_bbox_mask, name='class bbox mask', interpolation='bicubic')
    #     viewer.add_image(class_accurate_mask, name='class accurate mask', interpolation='bicubic')
    #     viewer.add_image(cell, name='cell', interpolation='bicubic')

    # trying to reduce intermediate memory usage
    del cell, noisy_cell, class_accurate_mask, class_bbox_mask, occupancy_accurate_mask, occupancy_bbox_mask

    return


def generate_projections(angles, outputFolder='./', modelID=0, pixelSize=1e-9, voltage = 200e3,
                         sphericalAberration=2.7E-3, multislice=None, msdz=5e-9, amplitudeContrast=0.07, defocus=2,
                         sigmaDecayCTF=0.4, noise_sigma = 0.04):

    grandcell = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{0}.mrc')

    SIZE = grandcell.shape[0]
    imageSize = SIZE//2
    heightBox = grandcell.shape[2]

    noisefree_projections = xp.zeros((imageSize, imageSize, len(angles)), dtype=float32)

    #TODO allow for different defocus per tilt image. Now single defocus for all images
    ctf = calcCTF(defocus, xp.zeros((imageSize,imageSize)), pixelSize, voltage=voltage, Cs=sphericalAberration,
                  sigma_decay_ctf=sigmaDecayCTF, amplitude_contrast=amplitudeContrast)


    for n, angle in enumerate(angles):

        filename = f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{int(angle)}.mrc'
        rotated_volume = pytom.tompy.io.read_mrc(filename)  #* 1E12

        # add structural noise with sigma dependent on the structural SNR
        # default value structSNR = 1.4 from Baxter et al. (2009)

        print(f'Adding structural noise to rotated volume with sigma = {noise_sigma:5.3f}')

        rotated_volume += xp.random.normal(0, noise_sigma, rotated_volume.shape) * (rotated_volume > 0)

        print(f'Volumes\' shape after structural noise addition is {rotated_volume.shape}')

        if not multislice:
            print('Simulating projection (without ms) from tilt angle ', angle)
            projected_tilt_image = rotated_volume[imageSize//2:-imageSize//2,imageSize//2:-imageSize//2,:].sum(axis=2)#.get() # remove .get() if fully cupy
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*ctf))**2
            #TODO is abs()**2 the way to go for exit wave field?

        else:
            print('Simulating projection (with ms) from tilt angle ', angle)

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

            print('Number of slices: ', n_slices)
            # Allocate space for multislice projection
            projected_potent_ms = xp.zeros((imageSize,imageSize, n_slices), dtype=complex)

            px_per_slice = int(xp.ceil(xp.around(heightBox / n_slices,3))) # CORRECT
            # PROJECTED POTENTIAL whithin slices (phase grating)
            for ii in range(n_slices):
                projected_potent_ms[:, :, ii] = rotated_volume[imageSize//2:-imageSize//2,imageSize//2:-imageSize//2,
                                                ii*px_per_slice : (ii+1)*px_per_slice].mean(axis=2) #.get() # remove .get() if fully cupy

            # zheight of slices in nm for Fresnel propagator
            dzprop = px_per_slice * pixelSize # CORRECT
            # wavelength
            Lambda = wavelength_eV2m(voltage)
            # relative mass
            relmass = phys_const_dict["me"] + phys_const_dict["el"] * voltage / (phys_const_dict["c"]**2)
            # sigma_transfer
            sig_transfer = 2 * xp.pi * relmass * phys_const_dict["el"] * Lambda / (phys_const_dict["h"]**2) # CORRECT

            # TRANSMISSON FUNCTION: the slice thickness is constant
            psi_t = xp.exp(1j * sig_transfer * projected_potent_ms * dzprop)


            # Get value of q in Fourier space because the Fresnel propagator needs them
            xwm = pixelSize * (imageSize)  # pixelsize for multislice * size sample
            q_true_pix_m = 1 / xwm
            q_m = rr(imageSize,imageSize) * q_true_pix_m  # frequencies in Fourier domain

            # FRESNEL PROPAGATOR
            P = xp.exp(-1j * xp.pi * Lambda * (q_m**2) * dzprop) # CORRECT

            # MULTISLICE
            psi_multislice = xp.zeros((imageSize,imageSize), dtype=complex) + 1 # should be complex datatype

            num_px_last_slice = heightBox % px_per_slice # CORRECT

            for ii in range(n_slices-min(1,num_px_last_slice)):
                #TODO ADD APPLICATION OF PYTOM.TOMPY

                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )

                psi_multislice = xp.fft.fftshift( xp.fft.ifftn((waveField * xp.fft.ifftshift(P) )) )

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                dzprop_end = num_px_last_slice * pixelSize
                psi_t[:, :, -1] = xp.exp(1j * sig_transfer * projected_potent_ms[:, :, -1] * dzprop_end)
                P_end = xp.exp(-1j * xp.pi * Lambda * (q_m ** 2) * dzprop_end)
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn( waveField * xp.fft.ifftshift(P_end) ) )

            # GET INTENSITIES IN IMAGE PLANE
            waveCTF = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
            projected_tilt_image = xp.abs(xp.fft.fftshift(xp.fft.ifftn(waveCTF))) ** 2

        noisefree_projections[:,:,n] = projected_tilt_image

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections_noisefree.mrc',
                         noisefree_projections) # noisefree_projections.astype(xp.float32)?

    del grandcell, noisefree_projections, ctf

    return

def add_effects_microscope(outputFolder, modelID, defocus, pixelsize, SNR, voltage, amplitude_contrast, Cs, sigma_decay_ctf):

    # TODO add try except for opening these files as I am not specifying paths
    noisefree_projections = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/projections_noisefree.mrc')

    projections = simutomo(noisefree_projections, defocus, pixelsize, SNR, outputFolder=outputFolder, modelID=modelID,
                           voltage=voltage, amplitude_contrast=amplitude_contrast, Cs=Cs, sigma_decay_ctf=sigma_decay_ctf)

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections.mrc', projections)

    return

def reconstruct_tomogram(prefix, suffix, start_idx, end_idx, vol_size, angles, outputFolder, modelID, weighting=-1):
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()

    for i in range(start_idx, end_idx+1):

        p = Projection(prefix+str(i)+suffix, tiltAngle=angles[i-1])
        projections.append(p)

    outputname = os.path.join(outputFolder, f'model_{modelID}/reconstruction.em')

    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0,0,0], binning=1, applyWeighting=weighting)
    vol.write(outputname)
    os.system(f'em2mrc.py -f {outputname} -t {os.path.dirname(outputname)}')
    os.system(f'rm {outputname}')

    return

def parallel_rotate_model(volume, outname, angle):
    print(f'Starting rotation process for angle {angle}')
    sys.stdout.flush()
    from voltools import transform
    # volume = pytom.tompy.io.read_mrc(filename)
    rotated_volume = transform(volume, rotation=(0, angle, 0), rotation_order='sxyz', interpolation='filt_bspline', device='cpu')
    pytom.tompy.io.write(outname, rotated_volume)
    print(f'Process for angle {angle} is finished ({outname})')
    sys.stdout.flush()
    return True

def create_rotation_model(outputFolder, modelID):

    # Load grandmodel
    save_path = f'{outputFolder}/model_{modelID}'
    if os.path.exists(f'{save_path}/grandmodel.mrc'):
        grandcell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel.mrc')
    else:
        raise Exception(f'create_rotation_model expects grandmodel be created before ({save_path}/grandmodel.mrc)')

    print(f'Successfully loaded grandcell (shape: {grandcell.shape})')

    size = grandcell.shape[0]
    height = grandcell.shape[2]
    if grandcell.shape[2] >= heightBox:
        raise Exception('Your model is larger than than the box that we will rotate (heightBox parameter)')

    volume = xp.zeros((size, size, heightBox), dtype=xp.float32)
    offset = (heightBox - height) // 2

    if (heightBox - height) % 2:
        volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
    else:
        volume[:, :, offset:-offset] = grandcell[:, :, :]

    dir = f'{save_path}/rotations'
    filename = f'{dir}/rotated_volume_0.mrc'
    pytom.tompy.io.write(filename, volume)

    print(f'Saved initial rotation volume at {filename}')
    return volume

if __name__ == '__main__':

    # Read config
    config = configparser.ConfigParser()
    try:
        if len(sys.argv) > 1:
            config_given = sys.argv[1]
            if config_given and os.path.exists(config_given):
                print(f'\nLoading a given configuration file: {config_given}')
                config.read_file(open(config_given))
        else:
            print(f'\nLoading default configuration file: pytom/simulation/simulation.conf')
            config.read_file(open('simulation.conf'))
    except Exception as e:
        print(e)
        raise Exception('Could not open config file.')

    print('Configuration sections:', config.sections())

    # Set simulation parameters
    try:
        outputFolder = config['General']['OutputFolder']
        modelID = int(config['General']['ModelID'])
        seed = int(config['General']['Seed'])

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

        print(f'Generating model {modelID} in folder {outputFolder}')
    except Exception as e:
        print(e)
        raise Exception('No general parameters specified in the config file.')

    if 'GenerateModel' in config.sections():
        try:
            particleFolder = config['GenerateModel']['ParticleFolder']
            listpdbs = eval(config['GenerateModel']['Models'])
            placementSize = int(config['GenerateModel']['PlacementSize'])
            size = int(config['GenerateModel']['Size'])
            thickness = int(config['GenerateModel']['Thickness'])
            waterdensity = float(config['GenerateModel']['WaterDensity'])
            proteindensity = float(config['GenerateModel']['ProteinDensity'])

            # parse range of number of particles
            p_range = config['GenerateModel']['NumberOfParticles'].split('-')
            if len(p_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                numberOfParticles = np.random.randint(int(p_range[0]), int(p_range[1]))
            else:
                numberOfParticles = int(p_range[0])
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
            start = int(config['ReconstructTomogram']['StartIdx'])
            end = int(config['ReconstructTomogram']['EndIdx'])
            weighting = int(config['ReconstructTomogram']['Weighting'])
            sizeRecon = int(config['ReconstructTomogram']['SizeRecon'])
        except Exception as e:
            print(e)
            raise Exception('Missing reconstruct tomogram parameters.')

    # Create directories and logger
    if not os.path.exists(os.path.join(outputFolder, f'model_{modelID}')):
        os.mkdir(os.path.join(outputFolder, f'model_{modelID}'))

    if os.path.exists(f'{outputFolder}/model_{modelID}/simulator.log'):
        os.remove(f'{outputFolder}/model_{modelID}/simulator.log')

    logging.basicConfig(filename='{}/model_{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(outputFolder, modelID,
                                                                                             date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    # Generate a grand model
    if 'GenerateModel' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating grand model')
        generate_model(particleFolder, outputFolder, modelID, listpdbs, size=size, thickness=thickness,
                       placementSize=placementSize, waterdensity=waterdensity, proteindensity=proteindensity,
                       numberOfParticles=numberOfParticles)

    # Generated rotated grand model versions
    if 'Rotation' in config.sections():
        print(f'\n- Rotating model')

        # If needed, creating rotations directory, removing leftover initial rotation volume
        dir = f'{outputFolder}/model_{modelID}/rotations'
        if not (os.path.exists(dir)):
            os.mkdir(dir)
        filename = f'{dir}/rotated_volume_0.mrc'
        if os.path.exists(filename):
            print(f'Found leftover rotated volume 0, removing it (path: {filename}')
            os.remove(filename)

        #  Create initial rotation volume (0 degree rotation)
        volume = create_rotation_model(outputFolder, modelID)

        # parallel rotate
        sys.stdout.flush()
        time.sleep(5)
        # sleep helps with memory usage
        # my theory (ilja): it gives time for garbage collector to clean up

        from joblib import Parallel, delayed
        verbosity = 55 # set to 55 for debugging, 11 to see progress, 0 to turn off output

        # joblib automatically memory maps a numpy array to child processes
        print(f'Rotating the model with {nodes} processes')
        results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads")\
            (delayed(parallel_rotate_model)(volume, f'{dir}/rotated_volume_{int(ang)}.mrc', ang)
             for ang in angles if ang != 0.)

        sys.stdout.flush()
        if all(results):
            print('All rotation processes finished successfully')
        else:
            print(f'{results.count(None)} rotation processes did not finish successfully')

    # Generate noise-free projections
    if 'GenerateProjections' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating projections')
        generate_projections(angles, outputFolder=outputFolder, modelID=modelID,pixelSize=pixelSize,
                             voltage=voltage, sphericalAberration=sphericalAberration,multislice=multislice, msdz=msdz,
                             amplitudeContrast=amplitudeContrast, defocus=defocus, sigmaDecayCTF=sigmaDecayCTF)

    # Add effects of the microscope to the noise-free projections
    if 'AddEffectsMicroscope' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Adding effects of microscope')
        add_effects_microscope(outputFolder, modelID, defocus, pixelSize, SNR, voltage, amplitudeContrast,
                               sphericalAberration, sigmaDecayCTF)

    # Reconstruct tomogram
    if 'ReconstructTomogram' in config.sections():
        print('\n- Reconstructing tomogram')
        prefix = os.path.join(outputFolder, f'model_{modelID}/noisyProjections/simulated_proj_')
        suffix = '.mrc'
        vol_size = [sizeRecon, sizeRecon, sizeRecon]
        reconstruct_tomogram(prefix, suffix, start, end, vol_size, angles, outputFolder, modelID, weighting=weighting)

