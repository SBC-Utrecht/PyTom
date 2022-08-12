"""
Gijs' tilt series simulation - used for Ilja's contest paper: SHREC19
"""
import numpy
import sys
# from chimera import runCommand
from pytom.gui.mrcOperations import *
import matplotlib

try:
    matplotlib.use('Qt5Agg')
except:
    pass
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


# from scipy.ndimage import rotate


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


def wavelength_eV2nm(ev):
    return 6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * ev * 1) ** 0.5


def create_ctf(Dz, vol, pix_size, voltage=200E3, Cs=2.7, sigma=0):
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

    Lambda = wavelength_eV2nm(voltage)

    Ny = 1 / (2 * pix_size)
    nyqvist = 2 * pix_size * 10 ** 9

    if len(vol.shape) > 2:
        R, Y, Z = meshgrid(arange(-Ny, Ny, 2 * Ny / vol.shape[0]), arange(-Ny, Ny, 2 * Ny / vol.shape[1]),
                           arange(-Ny, Ny, 2 * Ny / vol.shape[2]))
        r = sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = meshgrid(arange(-Ny, Ny, 2. * Ny / (vol.shape[0])), arange(-Ny, Ny, 2. * Ny / (vol.shape[1])))
        r = sqrt(R ** 2 + Y ** 2)

    vol = sin(pi / 2 * (Cs * Lambda ** 3 * r ** 4 - 2 * Dz * Lambda * r ** 2))
    amplitude = cos(pi / 2 * (Cs * Lambda ** 3 * r ** 4 - 2 * Dz * Lambda * r ** 2))

    if sigma:
        vol = vol * exp(-(r / (sigma * Ny)) ** 2)
        amplitude = amplitude * exp(-(r / (sigma * Ny)) ** 2)

    print(vol.shape, amplitude.shape)
    return vol, amplitude


def tom_dev(image):
    return image.mean(), image.max(), image.min(), image.std(), image.std() ** 2


def simutomo(dens, defocus, pixelsize, tiltrange, tiltincr, SNR, the, psi, phi, globalsigma=0, gpu=False, pytom=False,
             outputFolder='./', modelID=0):
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

    # rotate volume according to specified euler angles
    # rotvol = rotate(dens,phi)

    rotvol = dens

    # calculate CTF
    ctf1, amplitude = create_ctf(defocus, (dens.sum(axis=0)), pixelsize / 10., 200, 2.7, 0.4)

    ctf1 = -ctf1
    amplitude = -amplitude
    ctf = 0.93 * ctf1 + 0.07 * amplitude

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
        tmp = fftshift(fftn(noisy)) * ctf + fftshift(fftn(bg))
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


def generate_model(outputFolder, modelID, listpdbs, waterdensity=1, numberOfParticles=1000):
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

    # if addStructNoise:
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


def generate_projections(outputFolder, modelID, SIZE, angles, grandcell):
    from voltools.volume import Volume

    volume = zeros((1200, SIZE, SIZE), dtype=float32)
    volume[225:-225, :, :] = grandcell[:, :, :]
    d_vol = Volume(volume, cuda_warmup=False)

    dx, dy, dz = grandcell.shape
    angles = numpy.arange(angles[0], int(angles[1]) + angleIncrement / 2, angleIncrement)

    noisefree_projections = zeros((len(angles), SIZE // 2, SIZE // 2), dtype=float32)

    for n, angle in enumerate(angles):
        gpu_proj = d_vol.transform(rotation=(0, angle, 0), rotation_order='szyx', rotation_units='deg',
                                   around_center=True).project(cpu=True)
        # noisefree_projections[n] = gpu_proj[300:-300,300:-300]
        noisefree_projections[n] = gpu_proj[256:-256, 256:-256]

        # =  create_projections(dens, range(-60, 61, 3))
        # noisefree_projections[n] /= noisefree_projections[n][:,140:160].mean()

    if not os.path.exists('{}/model_{}'.format(outputFolder, modelID)): os.mkdir(
        '{}/model_{}'.format(outputFolder, modelID))

    a = mrcfile.new('{}/model_{}/projections_grandmodel_{}.mrc'.format(outputFolder, modelID, modelID),
                    data=noisefree_projections, overwrite=True)
    a.close()
    del d_vol
    return noisefree_projections


def add_effects_microscope(outputFolder, modelID, grandcell, noisefree_projections):
    iter = 20

    for n in range(len(noisefree_projections)):
        noisefree_projections[n] = rot90(noisefree_projections[n], 3)

    diff = (grandcell.shape[-1] - noisefree_projections.shape[-1])
    simtomo = simutomo(grandcell[:, diff // 2:-diff // 2, diff // 2:-diff // 2], defocus, pixelSize,
                       [angles[0], angles[1]], angleIncrement, SNR,
                       0, 0, 0, outputFolder=outputFolder, modelID=modelID)
    if simtomo:
        convert_numpy_array3d_mrc(simtomo, f"{outputFolder}/simutomo_{modelID}.mrc")


def reconstruct_tomogram(prefix, suffix, start_idx, end_idx, volsize, angles, outputFolder, modelID, weighting=-1):
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()

    for i in range(start_idx, end_idx + 1):
        p = Projection(prefix + str(i) + suffix, tiltAngle=angles[i - 1])
        projections.append(p)

    outputname = os.path.join(outputFolder, f'tomogram_model_{modelID}.em')

    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0, 0, 0], binning=1,
                                        weighting=weighting)
    vol.write(outputname)
    os.system('em2mrc.py -f {} -t {}'.format(outputname, os.path.dirname(outputname)))
    os.system(f'rm {outputname}')


if __name__ == '__main__':

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['--generateModel'], 'generate 3D object', False, True),
               ScriptOption(['--generateProjections'], 'generate 2D projections from 3D model.', False, True),
               ScriptOption(['--addEffectsMicroscope'], 'Simulate effects of microscope on top of projections', False,
                            True),
               ScriptOption(['--reconstructTomogram'], 'Reconstruct a tomogram from projections', False, True),

               ScriptOption(['-s', '--size'], 'Size of the model. Default 1024', True, True),
               ScriptOption(['--sizeRecon'], 'Size of the reconstruction. Default 512', True, True),
               ScriptOption(['-m', '--models'], 'PDBIDs, seperated by a comma', True, True),
               ScriptOption(['-w', '--waterdensity'], 'Water density. 1 by default.', True, True),
               ScriptOption(['-m', '--modelID'], 'Model ID.', True, True),
               ScriptOption(['-g', '--grandmodelFile'], 'Path to grand model', True, True),
               ScriptOption(['--projectionsFile'], 'Stack of projections', True, True),
               ScriptOption(['-n', '--numberOfParticles'], 'Number of particles. 10.000 by default', True, True),
               ScriptOption(['-o', '--outputFolder'], 'Output folder, Current by default.', True, True),
               ScriptOption(['-a', '--angles'], 'Output folder, "-60,60" by default.', True, True),
               ScriptOption(['-i', '--angleIncrement'], 'Output folder, "3" by default.', True, True),
               ScriptOption(['--SNR'], 'Signal to Noise Ratio, 0.02 by default.', True, True),
               ScriptOption(['-d', '--defocus'], 'Defocus (um), "-2" by default.', True, True),
               ScriptOption(['-p', '--pixelSize'], 'Size of a detector pixel (A), "10" by default.', True, True),
               ScriptOption(['--start'], 'Start index of projections used for reconstruction, "2" by default.', True,
                            True),
               ScriptOption(['--end'], 'End index of projections used for reconstruction, "40" by default.', True,
                            True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Simulate Cell, Projections and Noisy projections.',
                          authors='GS',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        generateModel, generateProjections, addEffectsMicroscope, reconstructTomogram, \
        SIZE, sizeRecon, models, waterdensity, modelID, grandmodelFile, \
        projectionsFile, numberOfParticles, outputFolder, angles, angleIncrement, SNR, defocus, pixelSize, start, end, help = \
            parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if models is None:
        listpdbs = ['3cf3', '1s3x', '1u6g', '4b4t', '1qvr', '3h84', '2cg9', '3qm1', '3gl1', '3d2f', '4d8q', '1bxn']
    else:
        listpdbs = models.split(',')

    if SIZE is None:
        SIZE = 1024
    else:
        SIZE = int(SIZE)

    if sizeRecon is None:
        sizeRecon = 512
    else:
        sizeRecon = int(sizeRecon)

    if waterdensity is None:
        waterdensity = 1.
    else:
        waterdensity = float(waterdensity)

    if numberOfParticles is None:
        numberOfParticles = 1000
    else:
        numberOfParticles = int(numberOfParticles)

    if angles is None:
        angles = [-60, 60]
    else:
        angles = map(float, angles.split(','))

    if angleIncrement is None:
        angleIncrement = 3
    else:
        angleIncrement = float(angleIncrement)

    if SNR is None:
        SNR = 0.02
    else:
        SNR = float(SNR)

    if defocus is None:
        defocus = -2
    else:
        defocus = float(defocus)

    if pixelSize is None:
        pixelSize = 10
    else:
        pixelSize = float(pixelSize)

    if modelID == None:
        modelID = 0
    else:
        modelID = int(modelID)

    if outputFolder is None:
        outputFolder = './'

    if start is None:
        start = 2
    else:
        start = int(start)

    if end is None:
        end = 40
    else:
        end = int(start)

    if not os.path.exists(outputFolder):
        raise Exception('Please make sure the output directory exists.')

    # Generate or read a grand model
    if generateModel:
        grandcell = generate_model(outputFolder, modelID, listpdbs, waterdensity, numberOfParticles)
    elif grandmodelFile:
        from pytom.gui.mrcOperations import read_mrc

        grandcell = read_mrc(grandmodelFile)
    else:
        grandcell = None

    # Generate or read noise-free projections
    if generateProjections and not (grandcell is None):
        noisefree_projections = generate_projections(outputFolder, modelID, SIZE, angles, grandcell)
    elif projectionsFile:
        noisefree_projections = read_mrc(projectionsFile)
    else:
        noisefree_projections = None

    # Add effects of the microscope to the noise-free projections
    if addEffectsMicroscope and not (grandcell is None) and not (noisefree_projections is None):
        print('Adding Effects of Microscope')
        add_effects_microscope(outputFolder, modelID, grandcell, noisefree_projections)

    if reconstructTomogram:
        prefix = os.path.join(outputFolder, f'model_{modelID}/noisyProjections/simulated_proj_model{modelID}_')
        suffix = '.mrc'
        vol_size = [sizeRecon, sizeRecon, sizeRecon]
        angles = range(angles[0], angles[1] + 1, angleIncrement)
        reconstruct_tomogram(prefix, suffix, start, end, vol_size, angles, outputFolder, modelID, weighting=True)