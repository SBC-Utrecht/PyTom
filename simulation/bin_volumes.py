from potential import reduce_resolution, bin, extend_volume
from pytom.tompy.io import read_mrc, write

if __name__ == '__main__':
    voxel_size = 1
    bin_factor = 10
    folder = '/data2/mchaillet/structures/potential'
    pdb_ids = ['3cf3', '1s3x', '1u6g', '4cr2', '1qvr', '3h84', '2cg9', '3qm1', '3gl1', '3d2f', '4d8q', '1bxn']
    pdb_ids = [id.upper() for id in pdb_ids]
    output = '/data2/mchaillet/structures/particles'

    for id in pdb_ids:
        print(f'binning {id}')
        volume = read_mrc(f'{folder}/{id}/vatom_ph7.0_1.00A.mrc')
        # extend volume by 10 angstrom to give a full ring of empty space after binning
        vol_ext = extend_volume(volume, [10,10,10], pad_value=0, symmetrically=True, true_center=False)
        write(f'{output}/{id}_ph7_10A.mrc',
              bin(reduce_resolution(vol_ext, voxel_size, voxel_size*bin_factor), bin_factor))
