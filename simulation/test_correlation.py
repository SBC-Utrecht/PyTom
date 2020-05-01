import numpy as xp

from potential import reduce_resolution, extend_volume
from pytom.tompy.io import read_mrc, write

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
    return peakValue, [shiftX, shiftY, shiftZ]


def correlation_at(volume, template, voxel_size=None, at=None, mask=None):
    if at is None:
        template -= V_WATER
        template[template<0] = 0
        score, shift = max_correlation(open_as_vol(volume), open_as_vol(template), mask)
        print(f'correlation score at {voxel_size}')
        print(score)
    elif (type(at) is list) and (type(voxel_size) is float):
        for resolution in at:
            new_volume = volume
            # new_volume = reduce_resolution(volume, voxel_size, resolution)
            new_template = reduce_resolution(volume, voxel_size, resolution)
            new_template -= V_WATER
            new_template[new_template<0] = 0
            score, shift = max_correlation(open_as_vol(new_volume), open_as_vol(new_template), mask)
            print(f'correlation score at {resolution}')
            print(score)
    else:
        print('invalid argument type for keyword at')
    return



em_map = read_mrc('/data2/mchaillet/structures/em_maps/emd_4877.map')

# BENCHMARK
v_atom = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_atom_ph7_grid-0.81A_res-2.60A.mrc')
v_atom = extend_volume(v_atom, [10,10,10])
# store template size for use later
template_size = v_atom.shape[0]

correlation_at(em_map, v_atom)

# Read updated volume
v_atom = read_mrc('/data2/mchaillet/structures/potential_iasa/6RGQ_rem-solvent_ph7.0_0.81A.mrc')


diff = [template_size - a for a in v_atom.shape]
v_atom = extend_volume(v_atom, diff, pad_value=0, symmetrically=False)

# correlate at given resolution
correlation_at(em_map, v_atom)

# reduce resolution needs cubic volume
res = [1, 1.5, 2, 2.5, 3, 3.5, 4]
# print(type(res), type(res) is list)
correlation_at(em_map, v_atom, voxel_size=0.81, at=res)



# v_pdb = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-5.0A.mrc')
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