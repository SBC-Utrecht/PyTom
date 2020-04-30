import numpy as xp
import scipy.ndimage

from pytom.tompy.io import read_mrc, write
from pytom.tompy.correlation import subPixelPeakParabolic
from pytom.tompy.score import FLCF, find_peak_position
from pytom_numpy import npy2vol, vol2npy
from pytom.tompy.tools import create_sphere

from potential import iasa_potential, resample_apbs, combine_potential, reduce_resolution, extend_volume

# plot
# import matplotlib
# matplotlib.use('Qt5Agg')
# from pylab import *

V_WATER = 4.5301

def correlate(volume, template):
    print('Start correlation for differently sized volumes.')
    iter_xyz = [a-b for a,b in zip(volume.shape, template.shape)]
    correlation = xp.zeros(tuple(iter_xyz))
    print(f'number of iterations = {xp.prod(iter_xyz)}')
    count = 1
    for x in range(iter_xyz[0]):
        for y in range(iter_xyz[1]):
            for z in range(iter_xyz[2]):
                if xp.mod(count,1000) == 0:
                    print(f'iteration #{count}')
                coef = nxcc(npy2vol(volume[x:x+template.shape[0], y:y+template.shape[1], z:z+template.shape[2]],3),
                           npy2vol(template,3))
                correlation[x,y,z] = coef
                count += 1
    return correlation

# Reduce resolution
# v_atom = read_mrc('/data2/mchaillet/structures/potential_iasa/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
# v_atom = reduce_resolution(v_atom, 0.81, 2.6)
# write('/data2/mchaillet/structures/correlation/6RGQ_atom_ph7_grid-0.81A_res-2.60A.mrc', v_atom)

# v_bond = read_mrc('/data2/mchaillet/structures/potential_bond/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
# v_bond = reduce_resolution(v_bond, 0.81, 2.6)
# write('/data2/mchaillet/structures/correlation/6RGQ_bond_ph7_grid-0.81A_res-2.60A.mrc', v_bond)

# v_full = read_mrc('/data2/mchaillet/structures/potential_map/6RGQ_rem-solvent_ph7.0_0.81A.mrc')
# v_full = reduce_resolution(v_full, 0.81, 2.6)
# write('/data2/mchaillet/structures/correlation/6RGQ_full_ph7_grid-0.81A_res-2.60A.mrc', v_full)

# chim = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-0.81A.mrc')
# diff = [max(chim.shape) - x for x in chim.shape]
# chim = reduce_resolution(extend_volume(chim, diff, pad_value=0, symmetrically=True), 0.81, 2.6)
# write('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-2.60A_pytom.mrc', chim)



resolution_em_map = 2.6
voxel_size_em_map = 0.81
# mask_radius = 130 #?

em_map = read_mrc('/data2/mchaillet/structures/em_maps/emd_4877.map')
write('/data2/mchaillet/structures/correlation/em_map_regrid.mrc', em_map)

chimera = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-2.60A.mrc')
write('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-2.60A_regrid.mrc', chimera)

'''

v_atom = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_atom_ph7_grid-0.81A_res-2.60A.mrc')
# v_atom = reduce_resolution(v_atom, 0.81, 4)
v_atom = extend_volume(v_atom, [10,10,10], pad_value=0, symmetrically=True)
mask = create_sphere(v_atom.shape, radius=v_atom.shape[0]/2-5, sigma=4, gpu=False)


correlation = FLCF(em_map, v_atom, mask)
ind  = find_peak_position(correlation)
print('Fast local correlation with v_atom:')
print(correlation[ind], ind)

score, coordinates = subPixelPeakParabolic(correlation, ind)
print('subpixel peak:')
print(score, coordinates)

# ========== Make same size volume as em map using paste_in_center as also used by FLCF
# from pytom.tompy.tools import paste_in_center
# new = xp.zeros((size,size,size))
# v_atom = paste_in_center(v_atom, new)
# v_atom = scipy.ndimage.interpolation.shift(v_atom, center_shift, order=3)
# write('/data2/mchaillet/structures/correlation/atom_interpolated.mrc', v_atom)
'''
v_atom = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_atom_ph7_grid-0.81A_res-2.60A.mrc')
v_atom = reduce_resolution(v_atom, 0.81, 3.5)
v_atom -= V_WATER
v_atom[v_atom<0] = 0
write('/data2/mchaillet/structures/correlation/6RGQ_atom_ph7_grid-0.81A_res-3.5A_excluded_solvent.mrc', v_atom)
# v_atom = extend_volume(v_atom, [10,10,10], pad_value=0, symmetrically=True)
'''
correlation = FLCF(em_map, v_atom, mask)
ind  = find_peak_position(correlation)
print('Fast local correlation with v_atom - V_WATER:')
print(correlation[ind], ind)

score, coordinates = subPixelPeakParabolic(correlation, ind)
print('subpixel peak:')
print(score, coordinates)


# v_atom2 = iasa_potential('/data2/mchaillet/structures/apbs/6RGQ/6RGQ_rem-solvent.pqr', voxel_size=0.81, filetype='pqr')
# v_atom2 = reduce_resolution(v_atom2, 0.81, 4)
# v_atom2 -= V_WATER
# v_atom2[v_atom2<0] = 0
#
# diff = [a-b for a,b in zip(v_atom.shape,v_atom2.shape)]
# v_atom2 = extend_volume(v_atom2, diff, pad_value=0, symmetrically=True)
#
# correlation = FLCF(em_map, v_atom2, mask)
# ind  = find_peak_position(correlation)
# print('Fast local correlation with v_atom2 - V_WATER:')
# print(correlation[ind], ind)
#
# score, coordinates = subPixelPeakParabolic(correlation, ind)
# print('subpixel peak:')
# print(score, coordinates)


v_full = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_full_ph7_grid-0.81A_res-2.60A.mrc')
diff = [a-b for a,b in zip(v_atom.shape,v_full.shape)]
v_full = extend_volume(v_full, diff, pad_value=0, symmetrically=True)

correlation = FLCF(em_map, v_full, mask)
ind  = find_peak_position(correlation)
print('Fast local correlation with v_full:')
print(correlation[ind], ind)

score, coordinates = subPixelPeakParabolic(correlation, ind)
print('subpixel peak:')
print(score, coordinates)
'''
v_full = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_full_ph7_grid-0.81A_res-2.60A.mrc')
v_full = reduce_resolution(v_full, 0.81, 3.5)
v_full -= V_WATER
v_full[v_full<0] = 0
write('/data2/mchaillet/structures/correlation/6RGQ_full_ph7_grid-0.81A_res-3.5A_excluded_solvent.mrc', v_full)
'''
correlation = FLCF(em_map, v_full, mask)
ind  = find_peak_position(correlation)
print('Fast local correlation with v_full - V_WATER:')
print(correlation[ind], ind)

score, coordinates = subPixelPeakParabolic(correlation, ind)
print('subpixel peak:')
print(score, coordinates)

v_chim = read_mrc('/data2/mchaillet/structures/correlation/6RGQ_pdb_chimera_grid-0.81A_res-2.60A.mrc')
diff = [a-b for a,b in zip(v_atom.shape,v_chim.shape)]
v_chim = extend_volume(v_chim, diff, pad_value=0, symmetrically=True)

correlation = FLCF(em_map, v_chim, mask)
ind  = find_peak_position(correlation)
print('Fast local correlation with v_chim at 2.6A:')
print(correlation[ind], ind)

score, coordinates = subPixelPeakParabolic(correlation, ind)
print('subpixel peak:')
print(score, coordinates)

if 0:
    three_angstrom = reduce_resolution(v_chim, resolution_em_map, 3)

    correlation = FLCF(em_map, three_angstrom, mask)
    ind  = find_peak_position(correlation)
    print('Fast local correlation with v_chim at 3A:')
    print(correlation[ind], ind)

    score, coordinates = subPixelPeakParabolic(correlation, ind)
    print('subpixel peak:')
    print(score, coordinates)

    four_angstrom = reduce_resolution(v_chim, resolution_em_map, 4)

    correlation = FLCF(em_map, four_angstrom, mask)
    ind  = find_peak_position(correlation)
    print('Fast local correlation with v_chim at 4A:')
    print(correlation[ind], ind)

    score, coordinates = subPixelPeakParabolic(correlation, ind)
    print('subpixel peak:')
    print(score, coordinates)

    five_angstrom = reduce_resolution(v_chim, resolution_em_map, 5)

    correlation = FLCF(em_map, five_angstrom, mask)
    ind  = find_peak_position(correlation)
    print('Fast local correlation with v_chim at 5A:')
    print(correlation[ind], ind)

    score, coordinates = subPixelPeakParabolic(correlation, ind)
    print('subpixel peak:')
    print(score, coordinates)

    six_angstrom = reduce_resolution(v_chim, resolution_em_map, 6)

    correlation = FLCF(em_map, six_angstrom, mask)
    ind  = find_peak_position(correlation)
    print('Fast local correlation with v_chim at 6A:')
    print(correlation[ind], ind)

    score, coordinates = subPixelPeakParabolic(correlation, ind)
    print('subpixel peak:')
    print(score, coordinates) '''
