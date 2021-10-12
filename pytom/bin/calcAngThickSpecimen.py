#!/usr/bin/env pytom

"""
Created on January, 2021

@author: GvdS

Estimate the thickness and angle of the specimen.
"""


def calculate_angle_thickness_specimen(fname, radius=5, verbose=True, binning=1):
    import matplotlib
    matplotlib.use('Qt5Agg')
    from pylab import imshow, show, subplots, savefig
    from pytom.agnostic.tools import create_sphere
    from pytom.agnostic.io import read, write
    from pytom.voltools import transform
    from pytom.gpu.initialize import xp, device
    from pytom.agnostic.correlation import meanVolUnderMask, stdVolUnderMask
    from pytom.agnostic.transform import resize
    import sys
    import os

    if verbose:
        print(fname)
    folder = os.path.dirname(os.popen(f'ls -alrt {fname}').read().split()[-1])
    print(folder)
    data = read(fname)
    if binning != 1:
        data = resize(data, 1/binning)

    mask2 = create_sphere(data.shape, radius=data.shape[0]*2//5)
    mask = create_sphere(data.shape, radius=radius)

    volume = data #* mask2

    meanV = meanVolUnderMask(volume, mask)
    stdV = stdVolUnderMask(volume, mask, meanV)

    return stdV



def calculate_angle_thickness_specimen(fname, radius=5, verbose=True, binning=6, radius2=2):
    import matplotlib
    matplotlib.use('Qt5Agg')
    from pylab import imshow, show, subplots, savefig
    from pytom.agnostic.tools import create_sphere
    from pytom.agnostic.io import read, write
    from pytom.voltools import transform
    from pytom.gpu.initialize import xp, device
    from pytom.agnostic.correlation import meanVolUnderMask, stdVolUnderMask
    from pytom.agnostic.transform import resize
    import sys
    import os


    if verbose:
        print(fname)
    folder = os.path.dirname(os.popen(f'ls -alrt {fname}').read().split()[-1])
    print(folder)
    data = read(fname)

    data[data > data.mean() + 5 * data.std()] = xp.mean(data)
    if binning != 1:
        data = resize(data, 1/binning)

    mask2 = create_sphere(data.shape, radius=data.shape[0]*2.2//5)
    mask = create_sphere(data.shape, radius=radius)

    volume = data #* mask2

    meanV = meanVolUnderMask(volume, mask)
    stdV = stdVolUnderMask(volume, mask, meanV)

    if 0 and verbose:
        if 1:
            outname = os.path.join(folder, 'stdV.mrc')
            if 0:
                r=6
                import scipy
                import numpy as np
                stdV[stdV < stdV.mean() + stdV.std()*2] = 0
                stdV /= stdV.max() / 100
                stdV = stdV.get()
                stdV[:, :r, :] = 0
                stdV[:r,:, :] = 0
                stdV[-r:, :, :] = 0
                stdV[:, -r:, :] = 0

                a = np.zeros((stdV.shape[0], stdV.shape[1])).astype(np.uint8)
                xm, ym = np.ogrid[0:stdV.shape[0]:10, 0:stdV.shape[1]:10]
                markers = np.zeros_like(a).astype(np.int16)
                markers[xm, ym] = np.arange(xm.size * ym.size).reshape((xm.size, ym.size))
                print(markers.shape, stdV.shape)
                res2 = scipy.ndimage.watershed_ift(stdV.astype(np.uint8).squeeze(), markers)
                res2[xm, ym] = res2[xm - 1, ym - 1]  # remove the isolate seeds
                res2[res2 > res2.max()-1] = -1

                for n, i in enumerate(np.unique(res2)):
                    res2[res2 == i] = n

                imshow(res2, cmap='viridis')
                show()

            if 0:
                import numpy as np
                import matplotlib.pyplot as plt
                from scipy import ndimage as ndi

                from skimage.segmentation import watershed
                from skimage.feature import peak_local_max

                image = stdV.get().squeeze()

                # Now we want to separate the two objects in image
                # Generate the markers as local maxima of the distance to the background
                distance = ndi.distance_transform_edt(image)
                coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=image)
                mask = np.zeros(distance.shape, dtype=bool)
                mask[tuple(coords.T)] = True
                markers, _ = ndi.label(mask)
                labels = watershed(-distance, markers, mask=image)

                fig, axes = plt.subplots(ncols=3, figsize=(9, 3), sharex=True, sharey=True)
                ax = axes.ravel()

                ax[0].imshow(image, cmap=plt.cm.gray)
                ax[0].set_title('Overlapping objects')
                ax[1].imshow(-distance, cmap=plt.cm.gray)
                ax[1].set_title('Distances')
                ax[2].imshow(labels, cmap=plt.cm.nipy_spectral)
                ax[2].set_title('Separated objects')

                for a in ax:
                    a.set_axis_off()

                fig.tight_layout()
                plt.show()

            if 1:
                r=8
                from skimage import restoration
                from skimage.morphology import watershed, label
                import scipy, skimage
                from skimage.feature import peak_local_max
                import numpy as np

                stdV[:, :r, :] = 0
                stdV[:r, :, :] = 0
                stdV[-r:, :, :] = 0
                stdV[:, -r:, :] = 0
                stdV[stdV < stdV.mean() + stdV.std() * 4] = 0
                stdV /= stdV.max() / 100


                footprint = np.ones((int(radius), int(radius)))

                image = restoration.denoise_tv_chambolle(stdV.get().squeeze(), weight=0.1)


                distance = scipy.ndimage.distance_transform_edt(image)
                local_maxi = peak_local_max(image, indices=False, footprint=footprint)

                markers = skimage.morphology.label(local_maxi)

                labels_ws = watershed(-distance, markers, mask=image)
                #labels_ws = label(image > 0)
                x, y = np.indices((stdV.shape[0], stdV.shape[1]))
                out = np.zeros((stdV.shape[0], stdV.shape[1]))
                particles = np.unique(labels_ws)
                for particle in particles[particles > 0]:
                    x1,y1 = scipy.ndimage.center_of_mass(labels_ws*(labels_ws == particle))
                    mask_circle1 = (x - x1) ** 2 + (y - y1) ** 2 < radius2 ** 2
                    out += mask_circle1.astype(np.float32)
                print(out)

                fig,ax = subplots(1,3,figsize=(15,5))
                ax[0].imshow(label(image.T > 0))
                ax[1].imshow(out.T)
                ax[2].imshow(data.squeeze().T.get())
                show()

            write(outname, stdV)
        # except Exception as e:
        #     print(e)
        #     pass

    stdV *= mask2
    stdV[mask2<0.5] = stdV.min()*0

    proj = xp.expand_dims((stdV).sum(axis=1), 2)
    rot = xp.zeros_like(proj, dtype=xp.float32)

    axis = 0

    bestA = 0
    top = 0
    curve = []
    curve_norm = []

    reference_line = (stdV[stdV>0].min()*mask2).sum(axis=2).sum(axis=1)

    if verbose:
        fig, ax = subplots(1, 3, figsize=(15, 5))

    for a in range(-12,12):
        transform(proj, output=rot, rotation=[0,0,a], rotation_order='rxyz', interpolation='filt_bspline',device=device)
        rot2 = rot.squeeze().sum(axis=axis)
        if top < rot2.max():
            top = rot2.max()
            bestA = a
        if verbose:
            try:
                ax[1].plot(rot2)
            except Exception as e:
                ax[1].plot(rot2.get())

    rA = bestA


    for a in xp.arange(bestA-1, bestA+1, 0.1):
        transform(proj, output=rot, rotation=[0,0,float(a)], rotation_order='rxyz', interpolation='filt_bspline',device=device)
        rot2 = rot.squeeze().sum(axis=axis)
        if top < rot2.max():
            top = rot2.max()
            bestA = a
            curve_norm = rot2-reference_line
            curve = rot2
            vol = rot.copy()
        if verbose:
            try:
                ax[2].plot(rot2-reference_line)
            except Exception as e:
                ax[2].plot((rot2-reference_line).get())


    maxid = -2
    ref = -2

    for offset in range(5,15):
        aa = curve.get()
        bb = reference_line.get()

        cc = bb * 0 + 1
        cc[150:-150] = 0

        from numpy.fft import ifftn, fftn, fftshift

        top = ifftn(fftn(aa*cc)*fftn(cc*bb*offset/10.).conj())
        print(offset/10, top.max())
        if top.max() > maxid:
            maxid = top.max()
            ref = offset/10.



    left2, right2 = 0, len(curve_norm)
    left = int(curve_norm[:curve.argmax()].argmin())

    right = int(curve_norm[curve.argmax():].argmin() + curve.argmax())

    for n, v in enumerate(curve_norm[left:curve.argmax()]):
        print(n, v)
        if v > 0:
            left2 = int(left+n)
            break

    for n, v in enumerate(curve_norm[curve.argmax():right][::-1]):
        if v >0:
            right2 = int(right-n)
            break






    try:
        sra_file = f'{folder}/z_limits.txt'

        out = open(sra_file, 'w')
        out.write(f'{left} {right}')
        out.close()

        if os.path.exists(os.path.join(os.path.dirname(sra_file), 'WBP_Reconstruction.sh')):
            d = open(os.path.join(os.path.dirname(sra_file), 'WBP_Reconstruction.sh'), 'r').read()[:-1]
            try:
                used_angle = float(d.split('specimenAngle ')[-1].split()[0])
            except Exception as e:
                used_angle = 0
        else:
            used_angle = 0

        out = open(f'{folder}/specimen_rotation_angle.txt', 'w')
        out.write(f'{float(bestA) + used_angle :.1f}\n')
        out.close()


        if verbose:
            print(f'written files in {folder}')
            print(f'best angle: {float(bestA+used_angle):.1f}\n\tconservative estimate:\t {left:4d} {right:4d}\n\tdetermined thickness:\t {left2:4d} {right2:4d}')
    except:
        pass

    if verbose:
        ax[1].set_title(f'{rA}')
        ax[1].plot(reference_line.get())
        ax[2].set_title(f'best angle: {float(bestA):.1f}\n thickness: {left} {right}\n thickness: {left2} {right2}')
        ax[2].vlines(left, curve_norm.min(), curve_norm.max(), lw=1, colors='#54d777', label='Cut off left')
        ax[2].vlines(left2, curve_norm.min(), curve_norm.max(), lw=1, colors='#54d777', label='Cut off left2')

        ax[2].vlines(right, curve_norm.min(), curve_norm.max(), lw=1, colors='#54d777', label='Cut off right')
        ax[2].vlines(right2, curve_norm.min(), curve_norm.max(), lw=1, colors='#54d777', label='Cut off right2')

        ax[2].legend()
        print(vol.shape)
        ax[0].imshow(vol.squeeze().get(),cmap='gray')
        ax[0].vlines(left, 0, stdV.shape[2], lw=1, colors='#54d777', label='Cut off left')
        ax[0].vlines(right,0,stdV.shape[2], lw=1, colors='#54d777', label='Cut off left')

        try:
            savefig(f'{folder}/anglethickness.png')
        except:
            pass
        show()


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='This script calculated the angle and the thickness of the sample',
        authors='GvdS',
        options=[ScriptOption2(
                     ['-f', '--file'], 'Filename will convert this single file to the target format', 'file', 'required'),
                 ScriptOption2(
                     ['-r', '--radius'], 'Raius of the mask used to calculate the stdv of volume', 'float', 'optional', 5),
                 ScriptOption2(
                     ['-v', '--verbose'], 'save and display the resuls in a graph', 'no arguments', 'optional'),
                 ScriptOption2(
                     ['-b', '--binning'], 'binning factor of image', 'int', 'optional', 1),
                 ScriptOption2(
                     ['-g', '--gpuID'], 'Run Job on single gpu', 'int', 'optional')])

    #TODO write --filter to filter input files maybe on (glob) pattern or else on extension or similar

    filename, radius, verbose, binning, gpuID = parse_script_options2(sys.argv[1:], helper)

    calculate_angle_thickness_specimen(filename, radius, verbose=verbose, binning=binning)