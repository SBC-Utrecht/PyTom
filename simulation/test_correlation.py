import numpy as xp
import sys

from pytom.tompy.io import read_mrc, read, write


V_WATER = 4.5301 # potential value of low density amorphous ice (from Vulovic et al., 2013)

def open_as_vol(npy_volume):
    # read can be used for reading mrc files as pytom_volumes
    # Essential the function first stores numpy array as mrc and then reopens it as a pytom volume. Function vol2npy
    # does not function correctly.
    import os
    from pytom.basic.files import read
    from pytom.tompy.io import write
    tmp_file = '/data2/mchaillet/structures/correlation/tmp.mrc'
    write(tmp_file, npy_volume)
    new_vol = read(tmp_file)
    os.remove(tmp_file)
    return new_vol

def max_correlation(volume, template, mask=None):
    """
    @param volume: pytom_volume
    @param template: pytom_volume
    @param mask: pytom_volume
    @return: score, [xshift, yshift, zshift]
    @rtype: float, [float,float,float]
    """
    from pytom.basic.correlation import subPixelPeak
    from pytom.score.score import FLCFScore
    from pytom_volume import vol, peak

    assert len(set([template.sizeX(), template.sizeY(), template.sizeZ()])) == 1, 'dimensions of template not equal'

    # generate the mask
    if mask.__class__ != vol:
        from pytom_volume import initSphere
        # cubic volume??
        mask = vol(template.sizeX(), template.sizeY(), template.sizeZ())
        mask.setAll(0)
        radius = template.sizeX() / 2 - 5
        sigma = 4
        centerX, centerY, centerZ = template.sizeX() / 2, template.sizeY() / 2, template.sizeZ() / 2
        # syntax = initSphere(volume, radius, sigma, maxradius, centerX, centerY, centerZ)
        initSphere(mask, radius, sigma, 0, centerX, centerY, centerZ)
    else:
        # print(template.sizeX(), mask.sizeX(), template.sizeY(), mask.sizeY(), template.sizeZ(), mask.sizeZ())
        if template.sizeX() != mask.sizeX() and template.sizeY() != mask.sizeY() and template.sizeZ() != mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')

    # Center position of volume
    centerX, centerY, centerZ = (volume.sizeX() / 2), (volume.sizeY() / 2), (volume.sizeZ() / 2)

    # Fast local correlation function does cross correlation in Fourier space
    scoreObject = FLCFScore()
    scoringResult = scoreObject.score(volume, template, mask)
    pk = peak(scoringResult)
    # subPixelPeak with spline interpolation for higher accuracy
    [peakValue,peakPosition] = subPixelPeak(scoreVolume=scoringResult, coordinates=pk,
                                           interpolation='Spline', verbose=False)

    # determine shift relative to center
    shiftX = (peakPosition[0] - centerX)
    shiftY = (peakPosition[1] - centerY)
    shiftZ = (peakPosition[2] - centerZ)
    # Return max correlation (and shifts relative to center??)
    return peakValue, (shiftX, shiftY, shiftZ)


def max_correlation_nxcf(volume, template, mask=None):
    """
    @param volume: pytom_volume
    @param template: pytom_volume
    @param mask: pytom_volume
    @return: score, [xshift, yshift, zshift]
    @rtype: float, [float,float,float]
    """
    from pytom.basic.correlation import subPixelPeak
    from pytom.score.score import nxcfScore
    from pytom_volume import vol, peak

    assert len(set([template.sizeX(), template.sizeY(), template.sizeZ()])) == 1, 'dimensions of template not equal'

    # generate the mask
    if mask.__class__ != vol:
        from pytom_volume import initSphere
        # cubic volume??
        mask = vol(template.sizeX(), template.sizeY(), template.sizeZ())
        mask.setAll(0)
        radius = template.sizeX() / 2 - 5
        sigma = 4
        centerX, centerY, centerZ = template.sizeX() / 2, template.sizeY() / 2, template.sizeZ() / 2
        # syntax = initSphere(volume, radius, sigma, maxradius, centerX, centerY, centerZ)
        initSphere(mask, radius, sigma, 0, centerX, centerY, centerZ)
    else:
        # print(template.sizeX(), mask.sizeX(), template.sizeY(), mask.sizeY(), template.sizeZ(), mask.sizeZ())
        if template.sizeX() != mask.sizeX() and template.sizeY() != mask.sizeY() and template.sizeZ() != mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')

    # Center position of volume
    centerX, centerY, centerZ = (volume.sizeX() / 2), (volume.sizeY() / 2), (volume.sizeZ() / 2)

    # Fast local correlation function does cross correlation in Fourier space
    scoreObject = nxcfScore()
    scoringResult = scoreObject.score(volume, template, mask)
    pk = peak(scoringResult)
    # subPixelPeak with spline interpolation for higher accuracy
    [peakValue,peakPosition] = subPixelPeak(scoreVolume=scoringResult, coordinates=pk,
                                           interpolation='Spline', verbose=False)

    # determine shift relative to center
    shiftX = (peakPosition[0] - centerX)
    shiftY = (peakPosition[1] - centerY)
    shiftZ = (peakPosition[2] - centerZ)
    # Return max correlation (and shifts relative to center??)
    return peakValue, (shiftX, shiftY, shiftZ)


def exclude_solvent(volume, solvent_potential=4.5301):
    excluded_volume = volume - solvent_potential
    excluded_volume[excluded_volume<0] = 0
    return excluded_volume


def correlation_at(volume, template, voxel_size=None, at=None, mask=None, filter_map=True, sol=4.5301):
    if at is None:
        template = exclude_solvent(template, solvent_potential=sol)
        if mask is None:
            score, shift = max_correlation(open_as_vol(volume), open_as_vol(template), mask)
        else:
            score, shift = max_correlation(open_as_vol(volume), open_as_vol(template), open_as_vol(mask))
        print(f'correlation score at {voxel_size}')
        print(f'{score:.3f}')
    elif (type(at) is list) and (type(voxel_size) is float):
        from potential import reduce_resolution
        for resolution in at:
            if filter_map:
                new_volume = reduce_resolution(volume, voxel_size, resolution)
            else:
                new_volume = volume
            new_template = exclude_solvent(template, solvent_potential=sol)
            new_template = reduce_resolution(new_template, voxel_size, resolution)
            if mask is None:
                score, shift = max_correlation(open_as_vol(new_volume), open_as_vol(new_template), mask)
            else:
                score, shift = max_correlation(open_as_vol(new_volume), open_as_vol(new_template), open_as_vol(mask))
            print(f'correlation score at {resolution}')
            print(f'{score:.3f}')
    else:
        print('invalid argument type for keyword at')
    return


def center_potential(emdb, potential, shift):
    from potential import extend_volume
    import scipy.ndimage
    em_size = emdb.shape[0]
    pt_size = potential.shape[0]
    # These steps perfectly center the em map and the constructed potential maps!
    potential = extend_volume(potential, [em_size - pt_size] * 3, symmetrically=True, true_center=False)
    return scipy.ndimage.interpolation.shift(potential, shift, order=3)

if __name__ == '__main__':

    em_map = read_mrc('/data2/mchaillet/structures/em_maps/emd_4877.map')
    # volume = read('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-3.0A.mrc')
    volume = read('/data2/mchaillet/structures/potential_map/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
    mask = read_mrc('/data2/mchaillet/structures/correlation/structural_mask.mrc')

    if 1:
        # iasa = read('/data2/mchaillet/structures/potential_iasa/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
        from potential import reduce_resolution
        bench = reduce_resolution(volume, 0.81, 2.6)
        bench = exclude_solvent(bench)

        # BENCHMARK

        score, shift = max_correlation(open_as_vol(em_map), open_as_vol(bench), open_as_vol(mask))
        print(f'{score:.3f}')
        print(shift)

        # bench_shift = center_potential(em_map, bench, shift)
        mask_shift = center_potential(em_map, mask, shift)

        # write('/data2/mchaillet/structures/correlation/6RGQ_full_shifted.mrc', bench_shift)
        # write('/data2/mchaillet/structures/correlation/6RGQ_mask_shifted.mrc', mask_shift)

        volume_shift = center_potential(em_map, volume, shift)
        # iasa_shift = exclude_solvent(iasa_shift)
        # iasa_shift = reduce_resolution(iasa_shift, 0.81, 2.6)

        correlation_at(em_map, volume_shift, mask=mask_shift, sol=V_WATER)

        # res = [3, 5, 7, 10, 20, 30]
        # correlation_at(em_map, volume_shift, voxel_size=0.81, at=res, mask=mask_shift, filter_map=True, sol=V_WATER)

        # Fourier shell correlation
        from pytom.basic.correlation import FSC
        # from pytom.basic.plot import plotFSC
        number_bands = 100
        fsc_result = FSC(open_as_vol(em_map), open_as_vol(volume_shift), number_bands, mask=open_as_vol(mask_shift), verbose=True,
                         filename='/data2/mchaillet/structures/correlation/6RGQ_full_FSC_results_reverse_sol.txt')
        # plot_filename = '/data2/mchaillet/structures/correlation/fsc_test'
        # plotFSC(fsc_result, plot_filename)

    if 0:

        # Read updated volume
        # v_atom = read_mrc('/data2/mchaillet/structures/potential_iasa/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
        mask = read_mrc('/data2/mchaillet/structures/correlation/structural_mask.mrc')

        # res = [2, 2.2, 2.4, 2.6, 2.8, 3]
        # correlation_at(em_map, v_atom, voxel_size=0.81, at=res, mask=mask, filter_map=False)
        # res = [3, 5, 7, 10, 20, 30]
        # correlation_at(em_map, v_atom, voxel_size=0.81, at=res, mask=mask, filter_map=True)
        #
        # v_int = read('/data2/mchaillet/structures/potential_map/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
        #
        # res = [2, 2.2, 2.4, 2.6, 2.8, 3]
        # correlation_at(em_map, v_int, voxel_size=0.81, at=res, mask=mask, filter_map=False)
        # res = [3, 5, 7, 10, 20, 30]
        # correlation_at(em_map, v_int, voxel_size=0.81, at=res, mask=mask, filter_map=True)

        v_pdb = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-3.0A.mrc')
        res = [3, 5, 7, 10, 20, 30]
        correlation_at(em_map, v_pdb, voxel_size=0.81, at=res, mask=mask, filter_map=True)
#
# diff = [a - template_size for a in v_pdb.shape] # a - tempalate_size because the chimera density is 400x400x400
# offset1 = [a // 2 for a in diff]
# offset2 = [a-b for a,b in zip(diff,offset1)]
# v_pdb = v_pdb[offset1[0]:-offset2[0], offset1[1]:-offset2[1], offset1[2]:-offset2[2]]
#
# score, shift = max_correlation(open_as_vol(em_map), open_as_vol(v_pdb))
#
# print('Correlation score of pdb electron density 5.0A (chimera)')
# print(score)
# print(shift)
#
# v_pdb = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-2.60A.mrc')
#
# diff = [abs(a - template_size) for a in v_pdb.shape] # a - tempalate_size because the chimera density is 400x400x400
# v_pdb = extend_volume(v_pdb, diff, pad_value=0, symmetrically=False)
#
# score, shift = max_correlation(open_as_vol(em_map), open_as_vol(v_pdb))
#
# print('Correlation score of pdb electron density 2.6A (chimera)')
# print(score)
# print(shift)