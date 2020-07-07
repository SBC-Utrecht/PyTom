import numpy as xp
import sys

from pytom.tompy.io import read_mrc, read, write
from potential import extend_volume


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
    from potential import reduce_resolution
    voxel_size = 1.08
    template_size = 150
    resolution = 3.2

    em_map = read_mrc('/data2/mchaillet/structures/em_maps/apoFer_masked.mrc')
    em_map_filtered = reduce_resolution(em_map, voxel_size, 5)
    write('/data2/mchaillet/structures/correlation/6m54/fitted/em_map_filtered_5A.mrc', em_map_filtered)
    em_map_filtered = reduce_resolution(em_map, voxel_size, resolution)
    write('/data2/mchaillet/structures/correlation/6m54/fitted/em_map_filtered_3.2A.mrc', em_map_filtered)

    mask = read_mrc('/data2/mchaillet/structures/potential/6m54/structural_mask.mrc')

    map_type = []
    scores_flcf = []

    # from pytom.tompy.score import FLCF

    # print('correlation with FLCF from tompy without peak interpolation')
    # template = read(f'/data2/mchaillet/structures/potential/6m54/vatom_ph7.5_1.08A.mrc')
    # score = FLCF(em_map, reduce_resolution(template, voxel_size, 3), mask=mask)
    # print(score.max())
    # template = read('/data2/mchaillet/structures/potential/6m54/6m54_pdb_chimera_3A.mrc')
    # difference = [template_size - x for x in template.shape]
    # template = extend_volume(template, difference, symmetrically=True, true_center=False)
    # score = FLCF(em_map, reduce_resolution(template, voxel_size, 3), mask=mask)
    # print(score.max())

    volume = read(f'/data2/mchaillet/structures/potential/6m54/6m54_pdb_chimera_{resolution}A.mrc')
    assert len(set(volume.shape)) == 1
    if volume.shape[0] != template_size:
        difference = [template_size - x for x in volume.shape]
        volume = extend_volume(volume, difference, symmetrically=True, true_center=False)
    score, shift, _ = max_correlation(open_as_vol(em_map_filtered), open_as_vol(volume),
                                      open_as_vol(mask))
    map_type.append('chimera')
    scores_flcf.append(score)

    mask_shift = center_potential(em_map, mask, shift)

    volume_shift = center_potential(em_map, volume, shift)

    # Store volume_shift and mask_shift
    write(f'/data2/mchaillet/structures/correlation/6m54/fitted/pdb_{resolution}A.mrc', volume_shift)
    write('/data2/mchaillet/structures/correlation/6m54/fitted/mask.mrc', mask_shift)

    # print(f'FSC correlation with pdb chimera {resolution}A')
    # from pytom.basic.correlation import FSC
    #
    # # from pytom.basic.plot import plotFSC
    # number_bands = 50
    # fsc_result = FSC(open_as_vol(em_map_filtered), open_as_vol(volume_shift), number_bands, mask=open_as_vol(mask_shift),
    #                  verbose=True,
    #                  filename=f'/data2/mchaillet/structures/correlation/6m54/FSC_pdb_{resolution}A.txt')

    from pytom.tompy.correlation import nxcc

    print('normalised cross correlation with pdb chimera and filtered map')
    res = [resolution, 5, 7, 10, 20, 30]
    scores = []
    for r in res:
        volume = read(f'/data2/mchaillet/structures/potential/6m54/6m54_pdb_chimera_{r}A.mrc')
        assert len(set(volume.shape)) == 1
        # if volume.shape[0] < template_size:
        #     difference = [template_size - x for x in volume.shape]
        #     volume = extend_volume(volume, difference, symmetrically=True, true_center=False)
        if volume.shape[0] > template_size:
            difference = [x - template_size for x in volume.shape]
            volume = volume[difference[0] // 2:-(difference[0] // 2 + difference[0] % 2),
                     difference[1] // 2:-(difference[1] // 2 + difference[1] % 2),
                     difference[2] // 2:-(difference[2] // 2 + difference[2] % 2)]
        volume_shift = center_potential(em_map, volume, shift)
        if r == 5:
            write(f'/data2/mchaillet/structures/correlation/6m54/fitted/pdb_{r}A.mrc', volume_shift)
        score = nxcc(reduce_resolution(em_map, voxel_size, r), volume_shift,
                     mask=mask_shift)
        print(r,score)
        scores.append(score)

    file = f'/data2/mchaillet/structures/correlation/6m54/chimera_correlation_map-filter.txt'
    with open(file, 'w') as f:
        for r, s in zip(res, scores):
            f.write(f'{r}, {s}\n')

    # Now the IASA and IASA+PBE map

    for pot in ['full_integration', 'full_integration_solvent_subtracted',
                'iasa_integration', 'iasa_integration_solvent_excluded']:

        volume = read(f'/data2/mchaillet/structures/potential/6m54/{pot}_1.08A.mrc')
        # mask = read_mrc('/data2/mchaillet/structures/potential/6m54/structural_mask.mrc')

        assert len(set(volume.shape)) == 1
        if volume.shape[0] != template_size:
            difference = [template_size - x for x in volume.shape]
            volume = extend_volume(volume, difference, symmetrically=True, true_center=False)

        score, shift, _ = max_correlation(open_as_vol(em_map_filtered),
                                          open_as_vol(reduce_resolution(volume, voxel_size, resolution)),
                                          open_as_vol(mask))
        map_type.append(pot)
        scores_flcf.append(score)

        volume_shift = center_potential(em_map, volume, shift)

        # Store shifted volume for grayscale maps
        write(f'/data2/mchaillet/structures/correlation/6m54/fitted/{pot}_{resolution}A.mrc', volume_shift)

        # print(f'FSC correlation with v_atom {resolution}A')
        # from pytom.basic.correlation import FSC
        #
        # # from pytom.basic.plot import plotFSC
        # number_bands = 50
        # fsc_result = FSC(open_as_vol(em_map_filtered), open_as_vol(reduce_resolution(volume_shift, voxel_size, resolution)),
        #                  number_bands, mask=open_as_vol(mask_shift), verbose=True,
        #                  filename=f'/data2/mchaillet/structures/correlation/6m54/FSC_{pot}_{resolution}A.txt')

        print(f'correlation with nxcc {pot} and map filtered')
        res = [resolution, 5, 7, 10, 20, 30]
        scores = []
        for r in res:
            if r == 5:
                write(f'/data2/mchaillet/structures/correlation/6m54/fitted/{pot}_{r}A.mrc', reduce_resolution(volume_shift, voxel_size, r))
            score = nxcc(reduce_resolution(em_map, voxel_size, r), reduce_resolution(volume_shift, voxel_size, r),
                         mask=mask_shift)
            print(r, score)
            scores.append(score)

        file = f'/data2/mchaillet/structures/correlation/6m54/{pot}_correlation_map-filter.txt'
        with open(file, 'w') as f:
            for r, s in zip(res, scores):
                f.write(f'{r}, {s}\n')

    # overall scores
    print('FLCF score for each input map')
    print(map_type)
    print(scores_flcf)