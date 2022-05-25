#!/bin/bash

# emd_2938.map was downloaded from pdb and resampled on 2.62A grid
# emd_5592.map gives very similar result
# NOTE: create_template.py can also take care of the resampling with the '--map_spacing' argument
create_template.py -f emd_2938_2.62A.mrc -d . -o emd_2938_21A.mrc -s 2.62 -b 8 -c -z 3 -a 0.1 --cut_first_zero -l 50 -g 0 -x 20

# create the mask
create_mask.py -o mask_20_8_1.mrc -b 20 -r 8 -s 1