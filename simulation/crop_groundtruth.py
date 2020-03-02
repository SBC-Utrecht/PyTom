import pytom.tompy.io as io
import shutil
import os

gt_size = (512, 512)

for folder in range(10):
    print(f'\nProcessing model {folder}\n')

    for file in ['grandmodel', 'grandmodel_noisefree',
                 'class_bbox', 'class_mask',
                 'occupancy_bbox', 'occupancy_mask']:

        path = f'model_{folder}/{file}.mrc'
        print(f'Loading {path}')

        vol = io.read_mrc(path)

        # backup original
        shutil.copy(path, f'model_{folder}/{file}_original.mrc')

        middle_coordinates = (vol.shape[0] // 2, vol.shape[1] // 2)
        print(f'\told shape: {vol.shape}')
        print(f'\tmiddle coordinates: {middle_coordinates}')

        vol = vol[middle_coordinates[0] - (gt_size[0] // 2) : middle_coordinates[0] + (gt_size[0] // 2),
                  middle_coordinates[1] - (gt_size[1] // 2) : middle_coordinates[1] + (gt_size[1] // 2),
                  :]

        print(f'\tnew shape: {vol.shape}')

        io.write(f'model_{folder}/{file}.mrc', vol)

    # backup
    p_txt = f'model_{folder}/particle_locations.txt'
    p_txt_backup = f'model_{folder}/particle_locations_original.txt'
    shutil.copy(p_txt, p_txt_backup)
    os.remove(p_txt)

    with open(p_txt_backup, 'r') as f:
        with open(p_txt, 'a') as f_cropped:
            for line in f:
                pdb_id, X, Y, Z, rot_Z1, rot_X, rot_Z2 = line.rstrip('\n').split()
                # there was an error first with ground truth text file writing, it was assumed X, Y - 256, Z - 256
                # but it should have been X - 256, Y - 256, Z
                # to accomodate, we use X-256, Y, Z+256
                f_cropped.write(f'{pdb_id} {int(float(X) - (gt_size[0] // 2))} {int(float(Y) - (gt_size[1] // 2))} {int(float(Z))} {rot_Z1} {rot_X} {rot_Z2}\n')

