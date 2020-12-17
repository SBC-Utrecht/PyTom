import numpy as xp
import sys

from pytom.tompy.io import read_mrc, read, write
from potential import extend_volume
from pytom_numpy import vol2npy


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
        if template.sizeX() != mask.sizeX() or template.sizeY() != mask.sizeY() or template.sizeZ() != mask.sizeZ():
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
    return peakValue, (shiftX, shiftY, shiftZ), mask


def max_correlation_nxcf(volume, template, mask=None):
    """
    @param volume: pytom_volume
    @param template: pytom_volume
    @param mask: pytom_volume
    @return: score, [xshift, yshift, zshift]
    @rtype: float, [float,float,float], pytom_volume
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
        if template.sizeX() != mask.sizeX() or template.sizeY() != mask.sizeY() or template.sizeZ() != mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')

    # Center position of volume
    centerX, centerY, centerZ = (volume.sizeX() / 2), (volume.sizeY() / 2), (volume.sizeZ() / 2)

    # normalised cross correlation function
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
    return peakValue, (shiftX, shiftY, shiftZ), mask


def exclude_solvent(volume, solvent_potential=4.5301):
    excluded_volume = volume - solvent_potential
    excluded_volume[excluded_volume<0] = 0
    return excluded_volume


def correlation_at(volume, template, voxel_size=None, at=None, mask=None, filter_map=True, sol=4.5301):
    if at is None:
        # template = exclude_solvent(template, solvent_potential=sol)
        if mask is None:
            score, shift, _ = max_correlation(open_as_vol(volume), open_as_vol(template), mask)
        else:
            score, shift, _ = max_correlation(open_as_vol(volume), open_as_vol(template), open_as_vol(mask))
        print(voxel_size,score)
        return score
    elif (type(at) is list) and (type(voxel_size) is float):
        from potential import reduce_resolution
        scores = []
        for resolution in at:
            if filter_map:
                new_volume = reduce_resolution(volume, voxel_size, resolution)
            else:
                new_volume = volume
            # new_template = exclude_solvent(template, solvent_potential=sol)
            new_template = reduce_resolution(template, voxel_size, resolution)
            if mask is None:
                score, shift, _ = max_correlation(open_as_vol(new_volume), open_as_vol(new_template), mask)
            else:
                score, shift, _ = max_correlation(open_as_vol(new_volume), open_as_vol(new_template), open_as_vol(mask))
            print(resolution, score)
            scores.append(score)
        return scores
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
    return scipy.ndimage.interpolation.shift(potential, shift, order=2)


if __name__ == '__main__':
    '''
    Should do circle_filter from pytom.tompy.tools instead of low pass gaussian for stronger suppresion of signal
    beyond the filter radius?
    '''
    from pytom.tompy.correlation import nxcc
    from potential import bin
    from potential import reduce_resolution #_2 as reduce_resolution # reduce resolution 2 contains stronger filtering
    # more appropriate for the experimental map as there the obtained resolution is often overestimated and not
    # homogeneous throughout the map.
    voxel_size = 1.08
    template_size = 180
    resolution = 3.2

    pdb = '6m54'

    em_map = read('/data2/mchaillet/structures/em_maps/apoFer_unmasked.mrc')
    # em_map = read('/data2/mchaillet/structures/em_maps/emd_4877.em')
    # em_map = covid
    # em_map = read_mrc('/data2/mchaillet/structures/em_maps/B_gal_unmasked.mrc')
    # em_map_filtered = reduce_resolution(em_map, voxel_size, 5)
    # write('/data2/mchaillet/structures/correlation/6m54/fitted/em_map_filtered_5A.mrc', em_map_filtered)

    # em_map = reduce_resolution(em_map, voxel_size/2, voxel_size)
    # em_map = bin(em_map, 2)

    em_map_filtered = reduce_resolution(em_map, voxel_size, resolution)
    write(f'/data2/mchaillet/structures/correlation/{pdb}/fitted/em_map_filtered_{resolution:.1f}A.mrc', em_map_filtered)
    # em_map_filtered=em_map

    # resolution = 4.0

    names = []
    scores = []
    res = [] # [5, 10, 20, 30]

    # mask = read_mrc('/data2/mchaillet/structures/potential/6m54/structural_mask.mrc')

    print('pdb chimera correlation through flcf scores')

    volume = read(f'/data2/mchaillet/structures/correlation/{pdb}/{pdb}_molmap_{resolution:.1f}A.mrc')
    # assert len(set(volume.shape)) == 1
    if any([s != template_size for s in volume.shape]):
        difference = [template_size - x for x in volume.shape]
        volume = extend_volume(volume, difference, symmetrically=True, true_center=False)
    score, shift, mask = max_correlation(open_as_vol(em_map_filtered), open_as_vol(volume))

    print(resolution, score)

    sub_scores = [score]
    for r in res:
        score, _,_ = max_correlation(open_as_vol(reduce_resolution(em_map, voxel_size, r)),
                                     open_as_vol(reduce_resolution(volume, voxel_size, r)), mask)
        sub_scores.append(score)
        print(r, score)

    names.append('chimera')
    scores.append(sub_scores)

    mask_shift = center_potential(em_map, vol2npy(mask), shift)

    volume_shift = center_potential(em_map, volume, shift)

    score = nxcc(em_map_filtered, volume_shift, mask_shift)
    print(score)

    # Store volume_shift and mask_shift
    write(f'/data2/mchaillet/structures/correlation/{pdb}/fitted/pdb_{resolution}A.mrc', volume_shift)
    write(f'/data2/mchaillet/structures/correlation/{pdb}/fitted/mask.mrc', mask_shift)

    print(f'FSC correlation with pdb chimera {resolution}A')
    from pytom.basic.correlation import FSC

    # from pytom.basic.plot import plotFSC
    number_bands = 50
    fsc_result = FSC(open_as_vol(em_map_filtered), open_as_vol(volume_shift), number_bands, mask=open_as_vol(mask_shift),
                     verbose=True,
                     filename=f'/data2/mchaillet/structures/correlation/{pdb}/FSC_pdb_{resolution}A.txt')

    # IASA map with solvent masking
    # pot = 'iasa_solvent_masked'
    # # volume = read('/data2/mchaillet/structures/potential/6m54/6m54_1.08A_solvent-7.000V_real.mrc')
    # volume = read('/data2/mchaillet/structures/potential/6m54/6m54_1.08A_solvent-4.530V_real.mrc')
    #
    # assert len(set(volume.shape)) == 1
    # if volume.shape[0] != template_size:
    #     difference = [template_size - x for x in volume.shape]
    #     volume = extend_volume(volume, difference, symmetrically=True, true_center=False)
    #
    # score, shift, _ = max_correlation(open_as_vol(em_map_filtered),
    #                                   open_as_vol(reduce_resolution(volume, voxel_size, resolution)),
    #                                   mask)
    # print('flcf correlation with solvent masked map')
    # print(resolution,score)
    #
    # sub_scores = [score]
    # for r in res:
    #     score, _, _ = max_correlation(open_as_vol(reduce_resolution(em_map, voxel_size, r)),
    #                                   open_as_vol(reduce_resolution(volume, voxel_size, r)), mask)
    #     sub_scores.append(score)
    #     print(r, score)
    #
    # names.append(pot)
    # scores.append(sub_scores)
    #
    # volume_shift = center_potential(em_map, volume, shift)
    #
    # # Store shifted volume for grayscale maps
    # write(f'/data2/mchaillet/structures/correlation/6m54/fitted/iasa_{resolution}A_smask.mrc', volume_shift)

    # print(f'FSC correlation with v_atom {resolution}A')
    # # from pytom.basic.plot import plotFSC
    # number_bands = 50
    # fsc_result = FSC(open_as_vol(em_map_filtered), open_as_vol(reduce_resolution(volume_shift, voxel_size, resolution)),
    #                  number_bands, mask=open_as_vol(mask_shift), verbose=True,
    #                  filename=f'/data2/mchaillet/structures/correlation/6m54/FSC_{pot}_{resolution}A.txt')

    # Now the IASA and IASA+PBE map
    for pot in [f'{pdb}_{voxel_size:.2f}A_real', f'{pdb}_{voxel_size:.2f}A_solvent-4.530V-gauss_real',
                f'{pdb}_{voxel_size:.2f}A_solvent-4.530V-mask_real']:

        volume = read(f'/data2/mchaillet/structures/potential/{pdb}/{pot}.mrc')
        # mask = read_mrc('/data2/mchaillet/structures/potential/6m54/structural_mask.mrc')

        assert len(set(volume.shape)) == 1
        if volume.shape[0] != template_size:
            difference = [template_size - x for x in volume.shape]
            volume = extend_volume(volume, difference, symmetrically=True, true_center=False)

        score, shift, _ = max_correlation(open_as_vol(em_map_filtered),
                                          open_as_vol(reduce_resolution(volume, voxel_size, resolution)),
                                          mask)

        print(f'correlation with flcf {pot} and map filtered')
        print(resolution,score)

        volume_shift = center_potential(em_map, reduce_resolution(volume, voxel_size, resolution), shift)

        score = nxcc(em_map_filtered, volume_shift, mask_shift)
        print(score)

        # Store shifted volume for grayscale maps
        write(f'/data2/mchaillet/structures/correlation/{pdb}/fitted/{pot}_{resolution}A.mrc', volume_shift)

        sub_scores = [score]
        for r in res:
            score, _, _ = max_correlation(open_as_vol(reduce_resolution(em_map, voxel_size, r)),
                                          open_as_vol(reduce_resolution(volume, voxel_size, r)), mask)
            sub_scores.append(score)
            print(r, score)

        names.append(pot)
        scores.append(sub_scores)

        print(f'FSC correlation with v_atom {resolution}A')
        # from pytom.basic.plot import plotFSC
        number_bands = 50
        fsc_result = FSC(open_as_vol(em_map_filtered), open_as_vol(reduce_resolution(volume_shift, voxel_size, resolution)),
                         number_bands, mask=open_as_vol(mask_shift), verbose=True,
                         filename=f'/data2/mchaillet/structures/correlation/{pdb}/FSC_{pot}_{resolution}A.txt')

    res = [resolution] + res
    scores = xp.ndarray.tolist(xp.transpose(xp.array(scores)))

    file = f'/data2/mchaillet/structures/correlation/{pdb}/correlation_flcf_results.txt'
    with open(file, 'w') as f:
        f.write('resolution\n')
        for n in names:
            f.write(f'#_{n}\n')

        for r, s in zip(res, scores):
            f.write(f'{r}')
            for s_i in s:
                f.write(f', {s_i}')
            f.write('\n')
