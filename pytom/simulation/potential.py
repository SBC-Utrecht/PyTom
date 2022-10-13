#!/usr/bin/env pytom

import multiprocessing as mp
import numpy as np
import pytom.simulation.physics as physics
import os, sys
# from numba import jit


# // N should be len(atoms) // 3
# // <<<N/TPB, TPB>>>
#
# // #include "math.h"  => can probably be skipped as cupy math.h probably includes erf()

iasa_integrate_text = '''
#define M_PI 3.141592654

extern "C" __global__ 
void iasa_integrate(float3 *atoms, unsigned char *elements, float *b_factors, float *occupancies, 
                    float *potential, float *solvent, float *scattering_factors, unsigned int *potential_dims,
                    float *displaced_volume, float voxel_size, unsigned int n_atoms, unsigned char exclude_solvent) {
    // exclude_solvent is used a Boolean value
    
    // get the atom index                                                                                                                                      
    unsigned int i = (blockIdx.x*(blockDim.x) + threadIdx.x); // correct for each atom having 3 coordinates
        
    if (i < n_atoms) {
    
        // make atoms a float3
        // potential_dims can also be a float3
        
        unsigned int j, l, m, n, potent_idx;
        int3 ind_min, ind_max;
        float atom_voxel_pot, sqrt_b, pi2_sqrt_b, pi2, sqrt_pi, factor3;
        float integral_x, integral_y, integral_z, integral_voxel;
        float3 voxel_bound_min, voxel_bound_max;
    
        // get all the atom information
        //atom.x = atoms[(i * 3) + 0];
        //atom.y = atoms[(i * 3) + 1];
        //atom.z = atoms[(i * 3) + 2];
        int elem = elements[i];  // its much easier if elements are ints, because dict wont work in C
        float b_factor = b_factors[i];
        float occupancy = occupancies[i];
        
        // scattering factors for this atom
        float *a = &scattering_factors[elem * 10];
        float *b = &scattering_factors[elem * 10 + 5];
        
        // get the radius of the volume displacement
        float r0 = cbrt(displaced_volume[elem] / sqrt(M_PI * M_PI * M_PI));
        
        // find the max radius over all gaussians with a certain cutoff
        float r2 = 0;
        for (j = 0; j < 5; j++) {
            r2 = max(r2, 15 / (4 * M_PI * M_PI / b[j]));
        };
        float r = sqrt(r2 / 3);
        
        // set indices in potential box
        ind_min.x = (int)floor((atoms[i].x - r) / voxel_size); // use floor and int casting for explicitness
        ind_max.x = (int)floor((atoms[i].x + r) / voxel_size);
        ind_min.y = (int)floor((atoms[i].y - r) / voxel_size);
        ind_max.y = (int)floor((atoms[i].y + r) / voxel_size);
        ind_min.z = (int)floor((atoms[i].z - r) / voxel_size);
        ind_max.z = (int)floor((atoms[i].z + r) / voxel_size);
        
        // enforce ind_min can never be smaller than 0 and ind_max can never be larger than potential_dims
        // ind min should also always be smaller than the dims...
        ind_min.x = max(ind_min.x, 0);
        ind_min.y = max(ind_min.y, 0);
        ind_min.z = max(ind_min.z, 0);
        ind_max.x = min(ind_max.x, potential_dims[0]);
        ind_max.y = min(ind_max.y, potential_dims[1]);
        ind_max.z = min(ind_max.z, potential_dims[2]);
        
        // precalc sqrt of pi
        sqrt_pi = sqrt(M_PI);
        // two times pi
        pi2 = 2 * M_PI;
        // loop over coordinates where this atom is present
        for (l = ind_min.x; l < (ind_max.x + 2); l++) {
        
            voxel_bound_min.x = l * voxel_size - atoms[i].x;
            voxel_bound_max.x = (l + 1) * voxel_size - atoms[i].x;
            
            for (m = ind_min.y; m < (ind_max.y + 2); m++) {
            
                voxel_bound_min.y = m * voxel_size - atoms[i].y;
                voxel_bound_max.y = (m + 1) * voxel_size - atoms[i].y;
                
                for (n = ind_min.z; n < (ind_max.z + 2); n++) {
                
                    voxel_bound_min.z = n * voxel_size - atoms[i].z;
                    voxel_bound_max.z = (n + 1) * voxel_size - atoms[i].z;
                    
                    // initialize to zero for this voxel
                    atom_voxel_pot = 0;
                    
                    for (j = 0; j < 5; j++) {
                        sqrt_b = sqrt(b[j]);
                        pi2_sqrt_b = pi2 / sqrt_b;
                        factor3 = pow(sqrt_b / (4 * sqrt_pi), 3);
                        
                        integral_x = (erf(voxel_bound_max.x * pi2_sqrt_b) - erf(voxel_bound_min.x * pi2_sqrt_b));
                        integral_y = (erf(voxel_bound_max.y * pi2_sqrt_b) - erf(voxel_bound_min.y * pi2_sqrt_b));
                        integral_z = (erf(voxel_bound_max.z * pi2_sqrt_b) - erf(voxel_bound_min.z * pi2_sqrt_b));
                        integral_voxel = integral_x * integral_y * integral_z * factor3;
                        
                        atom_voxel_pot += (a[j] / pow(b[j], (float)3 / 2)) * integral_voxel;
                    };
                    
                    potent_idx = l * potential_dims[1] * potential_dims[2] + m * potential_dims[2] + n;
                    atomicAdd( potential + potent_idx, atom_voxel_pot );
                    
                    if (exclude_solvent == 1) {
                        factor3 = pow(sqrt_pi * r0 / 2, 3);
                        
                        integral_x = erf(voxel_bound_max.x / r0) - erf(voxel_bound_min.x / r0);
                        integral_y = erf(voxel_bound_max.y / r0) - erf(voxel_bound_min.y / r0);
                        integral_z = erf(voxel_bound_max.z / r0) - erf(voxel_bound_min.z / r0);
                        
                        atomicAdd( solvent + potent_idx, factor3 * integral_x * integral_y * integral_z );
                    
                    };
                };
            };
        };
    };
};
'''


def extend_volume(vol, increment, pad_value=0, symmetrically=False, true_center=False, interpolation='filt_bspline'):
    """
    Increase volume by parameter increment ([x,y,z]). Options for changing the padding value, extending symmetrically on
    both sides of the input volume, shifting the original volume to the true new center if the increment/2 is
    non-integer.

    @param vol: input volume to be extended, 3d array
    @type  vol: L{np.ndarray}
    @param increment: list with increment value for each dimension
    @type  increment: L{list} - [L{int},] * 3
    @param pad_value: value to use as padding
    @type  pad_value: L{float}
    @param symmetrically: If False (default) the volume is just padded at the end of each dimension.
    @type  symmetrically: L{bool}
    @param true_center: If True interpolate the volume to true center in case increment is uneven in one direction.
    @type  true_center: L{bool}
    @param interpolation: voltools options ('filt_bspline', 'linear'), needed when shift is done to true_center
    @type  interpolation: L{string}

    @return: Extended volume, 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    if symmetrically:
        # condition = any([x%2 for x in increment])
        if true_center:
            from pytom.voltools import transform
            new = np.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)
            return transform(new, translation=tuple([x/2 for x in increment]), interpolation=interpolation)
        else:
            from pytom.agnostic.tools import paste_in_center
            new = np.zeros([a + b for a, b in zip(vol.shape, increment)])
            if pad_value:
                new += pad_value
            return paste_in_center(vol, new)
    else:
        return np.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)


def call_chimera(filepath, output_folder):
    """
    Run chimera for pdb file in order to add hydrogens and add symmetry units. The resulting structure is stored as a
    new pdb file with the extension {id}_symmetry.pdb. The function calls chimera on command line to execute.

    Reference Chimera: UCSF Chimera--a visualization system for exploratory research and analysis. Pettersen EF,
    Goddard TD, Huang CC, Couch GS, Greenblatt DM, Meng EC, Ferrin TE. J Comput Chem. 2004 Oct;25(13):1605-12.

    @param pdb_folder: path to a folder where pdb files are stored
    @type  pdb_folder: L{string}
    @param pdb_id: ID of pdb file
    @type  pdb_id: L{string}

    @return: name of file where the new pdb stucture is stored in, this can differ for pdb depending on
    symmetry thus we need to return it
    @rtype: L{string}

    @author: Marten Chaillet
    """
    print(f' - Calling chimera for {filepath}')

    input_folder, file = os.path.split(filepath)
    pdb_id, extension = os.path.splitext(file)
    # pdb_id = file.split('.')[0]
    # extension = file.split('.')[-1]

    # Skip if output files from this function already exist
    if os.path.exists(os.path.join(output_folder, f'{pdb_id}_rem-solvent_sym_addh.pdb')):
        output_filepath = os.path.join(output_folder, f'{pdb_id}_rem-solvent_sym_addh.pdb')
        print(f'File already exists: {output_filepath}')
        # return f'{pdb_id}_rem-solvent_sym_addh'
        return output_filepath
    elif os.path.exists(os.path.join(output_folder, f'{pdb_id}_rem-solvent_addh.pdb')):
        output_filepath = os.path.join(output_folder, f'{pdb_id}_rem-solvent_addh.pdb')
        print(f'File already exists: {output_filepath}')
        # return  f'{pdb_id}_rem-solvent_addh'
        return output_filepath

    if extension == '.pdb':
        # Command 'sym' in chimera crashes when there is no BIOMT symmetry in the pdb file. We need to make sure sym is only
        # executed when the BIOMT information specifies symmetrical units.
        symmetry = []
        try:
            with open(filepath,'r') as pdb:
                line = pdb.readline().split()
                while line[0] != 'REMARK':
                    line = pdb.readline().split()
                while line[0] == 'REMARK':
                    if line[1] == '350' and len(line) > 3:
                        if 'BIOMT' in line[2]:
                            symmetry.append(int(line[3]))
                    line = pdb.readline().split()
        except Exception as e:
            print(e)
            raise Exception('Could not read pdb file.')

        print(f'{pdb_id} has {len(set(symmetry))} symmetrical {"unit" if len(set(symmetry)) == 1 else "units"}.')

        if len(set(symmetry)) > 1:
            scriptpath = os.path.join(output_folder, f'_rem-solvent_sym_addh_{pdb_id}.py')
            output_filepath = os.path.join(output_folder, f'{pdb_id}_rem-solvent_sym_addh.pdb')
            try:
                with open(scriptpath, 'w') as chimera_script:
                    chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                         f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                         f'from chimera import runCommand as rc\n'
                                         f'rc("open {filepath}")\n'
                                         f'rc("delete solvent")\n'
                                         f'rc("delete ions")\n'
                                         f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                         f'rc("sym group biomt")\n'             # group biomt is also the default
                                         f'rc("combine all modelId 10")\n'
                                         f'rc("write format pdb #10 {output_filepath}")\n'
                                         f'rc("stop")\n')
            except Exception as e:
                print(e)
                raise Exception('Could not create chimera script.')
        else:
            scriptpath = os.path.join(output_folder, f'_rem-solvent_addh_{pdb_id}.py')
            output_filepath = os.path.join(output_folder, f'{pdb_id}_rem-solvent_addh.pdb')
            try:
                with open(scriptpath, 'w') as chimera_script:
                    chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                         f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                         f'from chimera import runCommand as rc\n'
                                         f'rc("open {filepath}")\n'
                                         f'rc("delete solvent")\n'
                                         f'rc("delete ions")\n'
                                         f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                         f'rc("write format pdb #0 {output_filepath}")\n'
                                         f'rc("stop")\n')
            except Exception as e:
                print(e)
                raise Exception('Could not create chimera script.')
        # module chimera should be loaded here...
        try:
            os.system(f'chimera --nogui --script {scriptpath}')
        except Exception as e:
            print(e)
            raise Exception('Chimera is likely not on your current path.')

        if len(set(symmetry)) > 1:
            return output_filepath
            # return f'{pdb_id}_rem-solvent_sym_addh' # returns new pdb name
        else:
            return output_filepath
            # return f'{pdb_id}_rem-solvent_addh'

    elif extension == '.cif':
        # If cif, assume the file does not contain any structural symmetry information as this is usually not the case
        # for large complexes.
        # Do transformations with chimera, and store as pdb. Chimera cannot write mmCIF files. ChimeraX can.
        # Also not that opening mmCIF files in chimera takes significantly longer.
        scriptpath = os.path.join(output_folder,f'_rem-solvent_addh_{pdb_id}.py')
        output_filepath = os.path.join(output_folder, f'{pdb_id}_rem-solvent_addh.pdb')
        try:
            with open(scriptpath, 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {filepath}")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                     f'rc("write format pdb #0 {output_filepath}")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
        # module chimera should be loaded here...
        try:
            os.system(f'chimera --nogui --script {scriptpath}')
        except Exception as e:
            print(e)
            raise Exception('Chimera is likely not on your current path.')
        # return f'{pdb_id}_rem-solvent_addh'
        return output_filepath
    else:
        print('non-valid structure file extension')
        return 0


def modify_structure_file(filepath, pattern, replacement, line_start=''):
    """
    Function required to make pqr files with large negative coordinates readible for APBS. APBS can only parse
    columns in the file when they are properly separated by white spaces. Input file will be overwritten.

    @param filepath: file path of pqr type file (with extension)
    @type  filepath: L{string}
    @param pattern: pattern to be replaced
    @type  pattern: L{string}
    @param replacement: replacement string for pattern
    @type  replacement: L{string}
    @param line_start: keyword argument, only modify line starting with line_start string
    @type  line_start: L{string}

    @return: - (input file overwritten)
    @rtype:  None

    @author: Marten Chaillet
    """
    from tempfile import mkstemp
    from shutil import move, copymode

    try:
        # Create temp file
        fh, abs_path = mkstemp()
        with os.fdopen(fh, 'w') as new_file:
            with open(filepath) as old_file:
                for line in old_file:
                    if line_start=='':
                        new_file.write(line.replace(pattern,replacement))
                    elif line.split()[0] == line_start:
                        new_file.write(line.replace(pattern,replacement))
                    else:
                        new_file.write(line)
        # Copy the file permissions from the old file to the new file
        copymode(filepath, abs_path)
        # Remove original file
        os.remove(filepath)
        # Move new file
        move(abs_path, filepath)
    except Exception as e:
        print(e)
        raise Exception('Unsuccessful in modifying pqr file with white space delimiters.')
    return


def call_apbs(pdb_filepath, force_field='amber', ph=7.):
    """
    Calls external programs pdb2pqr and apbs to execute on pdb structure. Both programs need to be on the path for this
    function to run. This function puts commands to the terminal to execute both programs.

    References
      - PDB2PQR: Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA. PDB2PQR: an automated pipeline for the setup,
        execution, and analysis of Poisson-Boltzmann electrostatics calculations. Nucleic Acids Research 32 W665-W667
        (2004).
      - APBS: Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of nanosystems: application to
        microtubules and the ribosome. Proc. Natl. Acad. Sci. USA 98, 10037-10041 2001.

    @param pdb_filepath: path to pdb file, this will be the folder where .pqr and .in file are stored
    @type  pdb_filepath: L{string}
    @param force_field: force field for parameterizing atoms (option: amber, ...)
    @type  force_field: L{string}
    @param ph: pH value of solvent surrounding the protein
    @type  ph: L{float}

    @return: - (output of programs called is stored folder of pdb_filepath)
    @rtype:  None

    @author: Marten Chaillet
    """
    folder, file = os.path.split(pdb_filepath)
    pdb_id, _ = os.path.splitext(file)

    print(f' - Running pdb2pqr and APBS on {pdb_filepath}')
    cwd = os.getcwd()

    pqr_filepath = os.path.join(folder, f'{pdb_id}.pqr')
    apbs_config  = os.path.join(folder, f'{pdb_id}.in')
    try:
        # Also PDB2PKA ph calculation method. Requires PARSE force field, can take very long for large proteins.
        os.system(f'pdb2pqr.py --ff={force_field} --ph-calc-method=propka --with-ph={ph} --apbs-input {pdb_filepath} {pqr_filepath}')
        print(' - Add white space delimiters to pqr file.')
        modify_structure_file(pqr_filepath, '-', ' -', line_start='ATOM')
        # APBS needs to execute from the folder where the structure is present
        os.chdir(folder)
        os.system(f'apbs {apbs_config}')
        # os.system(f'module load apbs_mc/1.5; apbs {apbs_in}')
        # ? subprocess.run(['apbs', f'{structure}.in'], cwd=output_folder) # executes in cwd
    except Exception as e:
        print(e)
        raise Exception('pdb2pqr or APBS potentially not on path.')

    # Change back to original directory
    os.chdir(cwd)
    return


def read_structure(filepath):
    """
    Read pdb, cif, or pqr file and return atom data in lists.

    todo move to basic.files or agnostic.io ??


    @param filepath: full path to the file, either .pdb, .cif, or .pqr
    @type  filepath: L{str}

    @return: a tuple of 6 lists (x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies)
    @rtype: L{tuple} -> (L{list},) * 6 with types (float, float, float, str, float, float)

    @author: Marten Chaillet
    """
    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = [], [], [], [], [], []

    _, extension = os.path.splitext(filepath)

    if extension == '.pdb':
        try:
            with open(filepath, 'r') as pdb:
                lines = pdb.readlines()
                atoms = [line for line in lines if line[:4] == 'ATOM']
                for line in atoms:
                    '''
        PDB example
        ATOM   4366  OXT SER I 456      10.602  32.380  -1.590  1.00 53.05           O
                    '''
                    x_coordinates.append(float(line[30:38]))
                    y_coordinates.append(float(line[38:46]))
                    z_coordinates.append(float(line[46:54]))
                    elements.append(line[76:78].strip())
                    b_factors.append(float(line[60:66]))
                    occupancies.append(float(line[54:60]))
                hetatms = [line for line in lines if line[:6] == 'HETATM']
                for line in hetatms:
                    '''
        PDB example
        HETATM13897 MG    MG A 501     120.846  94.563  17.347  1.00 79.97          MG
                    '''
                    x_coordinates.append(float(line[30:38]))
                    y_coordinates.append(float(line[38:46]))
                    z_coordinates.append(float(line[46:54]))
                    elements.append(line[76:78].strip())
                    b_factors.append(float(line[60:66]))
                    occupancies.append(float(line[54:60]))

        except Exception as e:
            print(e)
            raise Exception('Could not read pdb file.')
    elif extension == '.cif':
        try:
            with open(filepath, 'r') as pdb:
                lines = pdb.readlines()
                for line in lines:
                    if line.strip():
                        split_line = line.split()
                        if split_line[0] == 'ATOM':
                            '''
            PDBx/mmCIF example
            ATOM   171293 O  OP1   . G   WB 75 255  ? 252.783 279.861 251.593 1.00 50.94  ? 255  G   aa OP1   1
                            '''
                            x_coordinates.append(float(split_line[10]))
                            y_coordinates.append(float(split_line[11]))
                            z_coordinates.append(float(split_line[12]))
                            elements.append(split_line[2].strip())
                            b_factors.append(float(split_line[14]))
                            occupancies.append(float(split_line[13]))
                        elif split_line[0] == 'HETATM':
                            '''
            PDBx/mmCIF example
            HETATM 201164 MG MG    . MG  FD 79 .    ? 290.730 254.190 214.591 1.00 30.13  ? 3332 MG  A  MG    1
                            '''
                            x_coordinates.append(float(split_line[10]))
                            y_coordinates.append(float(split_line[11]))
                            z_coordinates.append(float(split_line[12]))
                            elements.append(split_line[2].strip())
                            b_factors.append(float(split_line[14]))
                            occupancies.append(float(split_line[13]))
        except Exception as e:
            print(e)
            raise Exception('Could not read cif file.')
    elif extension == '.pqr':
        try:
            with open(filepath, 'r') as pqr:
                lines = pqr.readlines()
                for line in lines:
                    if not line.strip():
                        split_line = line.split()
                        # TODO Whay about HETATM lines?
                        if split_line[0] == 'ATOM':
                            '''
                PQR example
                ATOM   5860  HA  ILE   379      26.536  13.128  -3.443  0.0869 1.3870
                            '''
                            x_coordinates.append(float(split_line[5]))
                            y_coordinates.append(float(split_line[6]))
                            z_coordinates.append(float(split_line[7]))
                            elements.append(split_line[2][0])  # first letter of long atom id is the element
                            b_factors.append(0.0) # not avalaible in PQR format
                            occupancies.append(1.0) # not avalaible in PQR format
                        # HETATM not working here because extracting element type from double letter elements, like MG, does
                        # not work properly. Should be tested though.
        except Exception as e:
            print(e)
            raise Exception('Could not read pqr file.')
    else:
        print(f'invalid filetype in read_structure() for {filepath}, return 0')
        return 0

    return x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies


def create_gold_marker(voxel_size, solvent_potential, oversampling=1, solvent_factor=1.0, imaginary=False, voltage=300E3):
    """
    From Rahman 2018 (International Journal of Biosensors and Bioelectronics).
    Volume of unit cell gold is 0.0679 nm^3 with 4 atoms per unit cell.
    Volume of gold bead is 4/3 pi r^3.

    @param voxel_size: voxel size of the box where gold marker is generated, in A
    @type  voxel_size: L{float}
    @param solvent_potential: solvent background potential
    @type  solvent_potential: L{float}
    @param oversampling: number of times to oversample the voxel size for more accurate generation
    @type  oversampling: L{int}
    @param solvent_factor: factor for denser solvent
    @type  solvent_factor: L{float}
    @param imaginary: flag for generating imaginary part of the potential
    @type  imaginary: L{bool}
    @param voltage: voltage of electron beam in eV, default 300E3
    @type  voltage: L{float}

    @return: if imaginary is True, return tuple (real, imaginary), if false return only real. boxes real and imag are
    3d arrays.
    @rtype: L{tuple} -> (L{np.ndarray},) * 2 or L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.agnostic.tools import create_sphere
    from pytom.simulation.support import reduce_resolution_real, create_ellipse, add_correlated_noise
    from pytom.agnostic.transform import resize

    assert (type(oversampling) is int) and (oversampling >= 1), print('Stop gold marker creation oversampling factor is not a positive'
                                                            ' integer.')

    # select a random size for the gold marker in nm
    diameter = np.random.uniform(low=4.0, high=10.0)

    # constants
    unit_cell_volume = 0.0679  # nm^3
    atoms_per_unit_cell = 4
    C = 2 * np.pi * physics.constants['h_bar']**2 / (physics.constants['el'] * physics.constants['me']) * 1E20  # nm^2
    voxel_size_nm = (voxel_size/10) / oversampling
    voxel_volume = voxel_size_nm**3

    # dimension of gold box, always add 5 nm to the sides
    dimension = int(np.ceil(diameter / voxel_size_nm)) * 3
    # sigma half of radius?
    r = 0.8 * ((diameter * 0.5) / voxel_size_nm)  # fraction of radius due to extension with exponential smoothing
    ellipse = True
    if ellipse:
        r2 = r * np.random.uniform(0.8, 1.2)
        r3 = r * np.random.uniform(0.8, 1.2)
        bead = create_ellipse(dimension, r, r2, r3, smooth=2)
    else:
        bead = create_sphere((dimension,)*3, radius=r)

    bead *= add_correlated_noise(int(r*0.75), dimension) * add_correlated_noise(int(r*0.25), dimension)
    # SIGMA DEPENDENT ON VOXEL SIZE
    # rounded_sphere = gaussian3d(sphere, sigma=(1 * 0.25 / voxel_size_nm))
    bead[bead < 0.9] = 0 # remove too small values
    # add random noise to gold particle to prevent perfect CTF ringing around the particle.
    # random noise also dependent on voxel size maybe?
    # rounded_sphere = (rounded_sphere > 0) * (rounded_sphere * np.random.normal(1, 0.3, rounded_sphere.shape))
    # rounded_sphere[rounded_sphere < 0] = 0

    if imaginary:
        solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY,
                                                        physics.WATER_MW, voltage) * solvent_factor
        gold_amplitude = physics.potential_amplitude(physics.GOLD_DENSITY, physics.GOLD_MW, voltage)
        gold_imaginary = bead * (gold_amplitude - solvent_amplitude)
        # filter and bin
        gold_imaginary = resize(reduce_resolution_real(gold_imaginary, 1, 2*oversampling), 1/oversampling,
                                interpolation='Spline')

    # values transformed to occupied volume per voxel from 1 nm**3 per voxel to actual voxel size
    solvent_correction = bead * (solvent_potential * solvent_factor)
    unit_cells_per_voxel = (bead * voxel_volume / unit_cell_volume)
    gold_atoms = unit_cells_per_voxel * atoms_per_unit_cell

    # interaction potential
    gold_scattering_factors = np.array(physics.scattering_factors['AU']['g'])
    # gold_scattering_factors[0:5].sum() == 10.57
    # C and scattering factor are in A units thus divided by 1000 A^3 = 1 nm^3 to convert
    gold_potential = gold_atoms * gold_scattering_factors[0:5].sum() * C / voxel_volume / 1000
    gold_real = gold_potential - solvent_correction
    # filter and bin
    gold_real = resize(reduce_resolution_real(gold_real, 1, 2*oversampling), 1/oversampling, interpolation='Spline')

    if imaginary:
        return gold_real, gold_imaginary
    else:
        return gold_real


def split_data(data_length, cores):
    if data_length == cores:
        indices = []
        for i in range(cores):
            indices.append((i,i+1))
        return indices
    elif data_length > cores:
        indices = [None] * cores
        n, N = 0, data_length
        for i in range(cores):
            l = N // cores + (N % cores > i)
            indices[i] = (n, n + l)
            n += l
        return indices
    else:
        indices = []
        for i in range(data_length):
            indices.append((i, i+1))
        return indices


def init(shared_data_, potential_shared_, solvent_shared_):
    global shared_data
    shared_data = shared_data_  # must be inherited, not passed as an argument to workers
    global potential_shared
    global solvent_shared
    potential_shared = potential_shared_
    solvent_shared = solvent_shared_


def tonumpyarray(mp_arr, shape, dt):
    return np.frombuffer(mp_arr.get_obj(), dtype=dt).reshape(shape)


def parallel_integrate(index, size, solvent_exclusion, voxel_size, dtype):
    from scipy.special import erf

    # return np.ndarray view of mp.Array
    potential = tonumpyarray(potential_shared, size, dtype)
    solvent = tonumpyarray(solvent_shared, size, dtype)

    print(f' --- process {mp.current_process().name} calculating {index[1]-index[0]} atoms')

    for i in range(index[0], index[1]):
        # order = x y z e b o
        x, y, z, element, b_factor, occupancy = shared_data[0][i], shared_data[1][i], shared_data[2][i], \
                                                shared_data[3][i], shared_data[4][i], shared_data[5][i]

        # atom type
        atom = element.upper()
        # atom center
        rc = [x, y, z]

        sf = np.array(physics.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b += (b_factor) # units in A

        if atom in list(physics.volume_displaced):
            r_0 = np.cbrt(physics.volume_displaced[atom] / (np.pi ** (3 / 2)))
        else:  # If not H,C,O,N we assume the same volume displacement as for carbon
            r_0 = np.cbrt(physics.volume_displaced['C'] / (np.pi ** (3 / 2)))

        r2 = 0  # it was 15 / (1 / r_0 ** 2) before but this gives very high values for carbon
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            r2 = np.maximum(r2, 15 / (4 * np.pi ** 2 / b[j]))
        # Radius of gaussian sphere
        r = np.sqrt(r2 / 3)

        ind_min = [int((c - r) // voxel_size) for c in rc]  # Smallest index to contain relevant potential x,y,z
        ind_max = [int((c + r) // voxel_size) for c in rc]  # Largest index to contain relevant potential x,y,z
        # Explicit real space coordinates for the max and min boundary of each voxel
        x_min_bound = np.arange(ind_min[0], ind_max[0] + 1, 1) * voxel_size - rc[0]
        x_max_bound = np.arange(ind_min[0] + 1, ind_max[0] + 2, 1) * voxel_size - rc[0]
        y_min_bound = np.arange(ind_min[1], ind_max[1] + 1, 1) * voxel_size - rc[1]
        y_max_bound = np.arange(ind_min[1] + 1, ind_max[1] + 2, 1) * voxel_size - rc[1]
        z_min_bound = np.arange(ind_min[2], ind_max[2] + 1, 1) * voxel_size - rc[2]
        z_max_bound = np.arange(ind_min[2] + 1, ind_max[2] + 2, 1) * voxel_size - rc[2]

        # x_min_bound, y_min_bound, z_min_bound = np.meshgrid(x_min_bound, y_min_bound, z_min_bound, indexing='ij')
        # x_max_bound, y_max_bound, z_max_bound = np.meshgrid(x_max_bound, y_max_bound, z_max_bound, indexing='ij')

        atom_potential = 0

        for j in range(5):
            sqrt_b = np.sqrt(b[j])  # calculate only once
            # Difference of error function == integrate over Gaussian
            int_x = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(x_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(x_min_bound * 2 * np.pi / sqrt_b))
            x_matrix = np.tile(int_x[:, np.newaxis, np.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(y_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(y_min_bound * 2 * np.pi / sqrt_b))
            y_matrix = np.tile(int_y[np.newaxis, :, np.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(z_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(z_min_bound * 2 * np.pi / sqrt_b))
            z_matrix = np.tile(int_z[np.newaxis, np.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])

            atom_potential += a[j] / b[j] ** (3 / 2) * x_matrix * y_matrix * z_matrix

        potential[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += \
            atom_potential

        if solvent_exclusion == 'gaussian':
            # excluded solvent potential
            int_x = np.sqrt(np.pi) * r_0 / 2 * (erf(x_max_bound / r_0) - erf(x_min_bound / r_0))
            x_matrix = np.tile(int_x[:, np.newaxis, np.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = np.sqrt(np.pi) * r_0 / 2 * (erf(y_max_bound / r_0) - erf(y_min_bound / r_0))
            y_matrix = np.tile(int_y[np.newaxis, :, np.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = np.sqrt(np.pi) * r_0 / 2 * (erf(z_max_bound / r_0) - erf(z_min_bound / r_0))
            z_matrix = np.tile(int_z[np.newaxis, np.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])

            solvent[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += (
                    x_matrix * y_matrix * z_matrix)

    print(f' --- process {mp.current_process().name} finished')


def iasa_integration_parallel(filepath, voxel_size=1., oversampling=1, solvent_exclusion=None,
                     V_sol=physics.V_WATER, absorption_contrast=False, voltage=300E3, density=physics.PROTEIN_DENSITY,
                     molecular_weight=physics.PROTEIN_MW, structure_tuple=None, cores=1):
    """
    Calculates interaction potential map to 1 A volume as described initially by Rullgard et al. (2011) in TEM
    simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013). This function applies averaging of
    the potential over the voxels to obtain precise results without oversampling.

    @param filepath: full filepath to pdb file
    @type  filepath: L{string}
    @param voxel_size: size of voxel in output map, default 1 A
    @type  voxel_size: L{float}
    @param oversampling: number of times to oversample final voxel size
    @type  oversampling: L{int}
    @param solvent_exclusion: flag to execute solvent exclusion using gaussian spheres. Solvent exclusion can be set
    with a string, either 'gaussian' or 'masking'. Default is None.
    @type  solvent_exclusion: L{str}
    @param solvent_masking: flag to do solvent exclusion using smoothed occupation mask (considered more accurate)
    @type  solvent_masking: L{bool}
    @param V_sol: average solvent background potential (V/A^3)
    @type  V_sol: L{float}
    @param absorption_contrast: flag to generate absorption factor for imaginary part of potential
    @type  absorption_contrast: L{bool}
    @param voltage: electron beam voltage, absorption factor depends on voltage, default 300e3
    @type  voltage: L{float}
    @param density: average density of molecule that is generated, default 1.35 (protein)
    @type  density: L{float}
    @param molecular_weight: average molecular weight of the molecule that is generated, default protein MW
    @type  molecular_weight: L{float}
    @param structure_tuple: structure information as a tuple (x_coordinates, y_coordinates, z_coordinates, elements,
    b_factors, occupancies), if provided this overrides file reading
    @type  structure_tuple: L{tuple} - (L{list},) * 6 with types (float, float, float, str, float, float)

    @return: A volume with interaction potentials, either tuple of (real, imag) or single real, both real and imag
    are 3d arrays.
    @rtype: L{tuple} -> (L{np.ndarray},) * 2 or L{np.ndarray}

    @author: Marten Chaillet

    TODO: This function could be way more elegant using multiprocessing.shared_memory.SharedMemory. However this
    TODO: functionality is only available from python3.8 onwards.
    """
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier
    from functools import partial, reduce
    from contextlib import closing
    import operator
    import ctypes
    import time

    assert (type(oversampling) is int) and (oversampling >= 1), print('oversampling parameter is not an integer')

    if structure_tuple is None:
        print(f' - Calculating IASA potential from {filepath}')
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filepath)
    else:
        print(f' - Calculating IASA potential from structure tuple')
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = structure_tuple

    # fix voxel size by oversampling
    if oversampling > 1:
        voxel_size /= oversampling
    # extend volume by 30 A in all directions to
    extra_space = 30
    # dV volume of a single voxel, needed for integration
    dV = voxel_size ** 3
    # conversion of electrostatic potential to correct units
    C = 4 * np.sqrt(np.pi) * physics.constants['h'] ** 2 / (
            physics.constants['el'] * physics.constants['me']) * 1E20  # angstrom**2

    # make coordinates start from origin
    x_max = np.max(x_coordinates - np.min(x_coordinates))
    y_max = np.max(y_coordinates - np.min(y_coordinates))
    z_max = np.max(z_coordinates - np.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)

    # give some extra space for atoms at the edges
    x_coordinates = x_coordinates - np.min(x_coordinates) + extra_space  # + difference[0] / 2
    y_coordinates = y_coordinates - np.min(y_coordinates) + extra_space  # + difference[1] / 2
    z_coordinates = z_coordinates - np.min(z_coordinates) + extra_space  # + difference[2] / 2
    # Define the volume of the protein
    sz_final = (int((largest_dimension + 2 * extra_space) / voxel_size),) * 3
    sz_initial = tuple([int((d + 2 * extra_space) / voxel_size) for d in dimensions])
    difference = tuple([f - i for f, i in zip(sz_final, sz_initial)])

    # split the data into fractions over the nodes
    indices = split_data(len(x_coordinates), cores)  # adjust nodes in case n_atoms is smaller than n_cores

    x_shared, y_shared, z_shared = mp.Array('d', x_coordinates, lock=False), \
                                   mp.Array('d', y_coordinates, lock=False), \
                                   mp.Array('d', z_coordinates, lock=False)
    b_shared, o_shared = mp.Array('d', b_factors, lock=False), mp.Array('d', occupancies, lock=False)
    e_shared = mp.Array(ctypes.c_wchar_p, elements, lock=False)
    # element_to_id = {k: i for i, k in enumerate(physics.scattering_factors.keys())}
    # e_shared = mp.Array('i', [element_to_id[e] for e in elements], lock=False)

    print(f'Number of atoms to go over is {len(x_coordinates)} spread over {len(indices)} processes')

    start = time.time()

    # create shared arrays
    potential_shared = mp.Array(ctypes.c_float, reduce(operator.mul, sz_initial))
    solvent_shared = mp.Array(ctypes.c_float, reduce(operator.mul, sz_initial))
    # initialize them via numpy
    dtype = np.float32
    potential, solvent = tonumpyarray(potential_shared, sz_initial, dtype), \
                                      tonumpyarray(solvent_shared, sz_initial, dtype)
    potential[:], solvent[:] = np.zeros(sz_initial, dtype=dtype), \
                                        np.zeros(sz_initial, dtype=dtype)

    with closing(mp.Pool(len(indices), initializer=init, initargs=((x_shared, y_shared, z_shared,
                                                                    e_shared, b_shared, o_shared),
                                                                   potential_shared, solvent_shared))) as p:
        p.map_async(partial(parallel_integrate,
                                          size=sz_initial,
                                          solvent_exclusion=solvent_exclusion,
                                          voxel_size=voxel_size,
                                          dtype=dtype), indices)
    p.join()

    end = time.time()
    print(f'shared mp.Array elapsed {end - start}')

    # convert potential to correct units and correct for solvent exclusion
    if solvent_exclusion == 'gaussian':
        # Correct for solvent and convert both the solvent and potential array to the correct units.
        real = (potential / dV * C) - (solvent / dV * V_sol)
    elif solvent_exclusion == 'masking':  # only if voxel size is small enough for accurate determination of mask
        solvent_mask = (potential > 1E-5) * 1.0
        # construct solvent mask
        # gaussian decay of mask
        if oversampling == 1:
            solvent_mask = reduce_resolution_fourier(solvent_mask, voxel_size, voxel_size * 2)
            solvent_mask[solvent_mask < 0.001] = 0
        # subtract solvent from the protein electrostatic potential
        real = (potential / dV * C) - (solvent_mask * V_sol)
    else:
        real = potential / dV * C

    # determine absorption contrast if set
    if absorption_contrast:
        # voltage by default 300 keV
        molecule_absorption = physics.potential_amplitude(density, molecular_weight, voltage)
        solvent_absorption = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY,
                                                         physics.WATER_MW, voltage) * (V_sol / physics.V_WATER)
        print('Calculating absorption contrast')
        print(f'Molecule absorption = {molecule_absorption:.3f}')
        print(f'Solvent absorption = {solvent_absorption:.3f}')

        if solvent_exclusion == 'masking':
            imaginary = solvent_mask * (molecule_absorption - solvent_absorption)
        elif solvent_exclusion == 'gaussian':
            imaginary = solvent / dV * (molecule_absorption - solvent_absorption)
        else:
            print('ERROR: Absorption contrast cannot be generated if solvent exclusion is not set to either gaussian '
                  'or masking.')
            sys.exit(0)

        real = np.pad(real, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                      mode='constant', constant_values=0)
        imaginary = np.pad(imaginary, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                           mode='constant', constant_values=0)

        if oversampling > 1:
            print('Rescaling after oversampling')
            real = reduce_resolution_fourier(real, voxel_size, voxel_size * 2 * oversampling)
            # TODO Bug with resizing!!
            real = resize(real, 1 / oversampling, interpolation='Spline')
            imaginary = resize(reduce_resolution_fourier(imaginary, voxel_size, voxel_size * 2 * oversampling),
                               1 / oversampling, interpolation='Spline')
        return real + 1j * imaginary
    else:
        # extend volume to a box
        real = np.pad(real, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                      mode='constant', constant_values=0)
        if oversampling > 1:
            print('Rescaling after oversampling')
            real = resize(reduce_resolution_fourier(real, voxel_size, voxel_size * 2 * oversampling),
                          1 / oversampling, interpolation='Spline')
        return real


def iasa_integration_gpu(filepath, voxel_size=1., oversampling=1, solvent_exclusion=None,
                     V_sol=physics.V_WATER, absorption_contrast=False, voltage=300E3, density=physics.PROTEIN_DENSITY,
                     molecular_weight=physics.PROTEIN_MW, structure_tuple=None, gpu_id=0):
    """
    Calculates interaction potential map to 1 A volume as described initially by Rullgard et al. (2011) in TEM
    simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013). This function applies averaging of
    the potential over the voxels to obtain precise results without oversampling.

    @param filepath: full filepath to pdb file
    @type  filepath: L{string}
    @param voxel_size: size of voxel in output map, default 1 A
    @type  voxel_size: L{float}
    @param oversampling: number of times to oversample final voxel size
    @type  oversampling: L{int}
    @param solvent_exclusion: flag to execute solvent exclusion using gaussian spheres. Solvent exclusion can be set
    with a string, either 'gaussian' or 'masking'. Default is None.
    @type  solvent_exclusion: L{str}
    @param V_sol: average solvent background potential (V/A^3)
    @type  V_sol: L{float}
    @param absorption_contrast: flag to generate absorption factor for imaginary part of potential
    @type  absorption_contrast: L{bool}
    @param voltage: electron beam voltage, absorption factor depends on voltage, default 300e3
    @type  voltage: L{float}
    @param density: average density of molecule that is generated, default 1.35 (protein)
    @type  density: L{float}
    @param molecular_weight: average molecular weight of the molecule that is generated, default protein MW
    @type  molecular_weight: L{float}
    @param structure_tuple: structure information as a tuple (x_coordinates, y_coordinates, z_coordinates, elements,
    b_factors, occupancies), if provided this overrides file reading
    @type  structure_tuple: L{tuple} - (L{list},) * 6 with types (float, float, float, str, float, float)

    @return: A volume with interaction potentials, either tuple of (real, imag) or single real, both real and imag
    are 3d arrays.
    @rtype: L{tuple} -> (L{np.ndarray},) * 2 or L{np.ndarray}

    @author: Marten Chaillet
    """
    import cupy as cp
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier

    # if (absorption_contrast) or (solvent_exclusion in ['gaussian', 'masking']) or (oversampling != 1):
    #     print('abs contrast, gaussian/masking solvent_exclusion, and oversampling still need to be implement for gpu')
    #     sys.exit(0)

    # cp set device ...
    cp.cuda.Device(gpu_id).use()

    assert (type(oversampling) is int) and (oversampling >= 1), print('oversampling parameter is not an integer')

    if structure_tuple is None:
        print(f' - Calculating IASA potential from {filepath}')
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filepath)
    else:
        print(f' - Calculating IASA potential from structure tuple')
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = structure_tuple

    # fix voxel size by oversampling
    if oversampling > 1:
        voxel_size /= oversampling
    # extend volume by 30 A in all directions to
    extra_space = 30
    # dV volume of a single voxel, needed for integration
    dV = voxel_size ** 3
    # conversion of electrostatic potential to correct units
    C = 4 * cp.sqrt(cp.pi) * physics.constants['h'] ** 2 / (
            physics.constants['el'] * physics.constants['me']) * 1E20  # angstrom**2

    # setup easier gpu indexing for scattering factors
    scattering = cp.array([v["g"] for k, v in physics.scattering_factors.items()], dtype=cp.float32)
    map_element_to_id = {k: i for i, k in enumerate(physics.scattering_factors.keys())}

    # setup gpu indexing for solvent displacement
    displacement = cp.array([physics.volume_displaced[k] if k in physics.volume_displaced.keys()
                            else physics.volume_displaced['C'] for k in map_element_to_id.keys()], dtype=cp.float32)

    # find dimensions
    n_atoms = len(x_coordinates)
    min_x, max_x, min_y, max_y, min_z, max_z = min(x_coordinates), max(x_coordinates), min(y_coordinates), \
                                        max(y_coordinates), min(z_coordinates), max(z_coordinates)
    x_size = max_x - min_x
    y_size = max_y - min_y
    z_size = max_z - min_z
    dimensions = (x_size, y_size, z_size)

    # Define the volume of the protein
    sz_potential = tuple([int((d + 2 * extra_space) / voxel_size) for d in dimensions])
    sz_potential_gpu = cp.array(sz_potential, dtype=cp.uint32)
    sz_final = (max(sz_potential),) * 3
    difference = tuple([a - b for a, b in zip(sz_final, sz_potential)])

    # initiate the final volume
    potential = cp.zeros(sz_potential, dtype=cp.float32)
    if solvent_exclusion == 'gaussian':
        solvent = cp.zeros(sz_potential, dtype=cp.float32)
    else:
        solvent = cp.array([0], dtype=cp.float32)

    # find the number of atoms that can be loaded to gpu each time
    # get current space on gpu
    free_space = cp.cuda.Device().mem_info[0]  # this value is in bytes
    # space per atom
    bytes_per_atom = 3 * 4 + 1 + 4 + 4  # 3 coordinate floats, element type uint, b_factor float, and occupancy float
    atoms_per_iter = free_space // bytes_per_atom
    niter = -(n_atoms // -atoms_per_iter)

    # print to user the number of atoms in system
    print(f'Number of atoms to go over is {n_atoms}')

    # create kernel
    iasa_integrate = cp.RawKernel(iasa_integrate_text, 'iasa_integrate')
    mempool = cp.get_default_memory_pool()

    # distribute the atoms in multiple runs if there are too many
    for i in range(niter):
        #TODO check that large systems are correctly split in batches
        # set the subset index for this iteration
        index = [i * atoms_per_iter, (i + 1) * atoms_per_iter if (i + 1) < niter else n_atoms]

        # move the selected atoms to gpu!
        atoms = cp.ascontiguousarray(cp.array([x_coordinates[index[0]:index[1]],
                                               y_coordinates[index[0]:index[1]],
                                               z_coordinates[index[0]:index[1]]]).T,
                                     dtype=cp.float32)
        n_atoms_iter = atoms.shape[0]
        elements = cp.array([map_element_to_id[e.upper()] for e in elements[index[0]:index[1]]], dtype=cp.uint8)
        b_factors = cp.array(b_factors[index[0]:index[1]], dtype=cp.float32)
        occupancies = cp.array(occupancies[index[0]:index[1]], dtype=cp.float32)

        print(f'{n_atoms_iter} on iter {i} ---- using index {index}')

        # give some extra space for atoms at the edges
        atoms[:, 0] += (- min_x + extra_space)  # + difference[0] / 2
        atoms[:, 1] += (- min_y + extra_space)  # + difference[1] / 2
        atoms[:, 2] += (- min_z + extra_space)  # + difference[2] / 2

        # threads and blocks
        n_threads = 1024
        n_blocks = int(cp.ceil(n_atoms_iter / n_threads).get())

        # call gpu integration either with or without gaussian solvent exclusion
        if solvent_exclusion == 'gaussian':
            iasa_integrate((n_blocks, 1, 1,), (n_threads, 1, 1), (atoms, elements, b_factors, occupancies, potential,
                                                                  solvent, scattering, sz_potential_gpu, displacement,
                                                                  cp.float32(voxel_size), cp.uint32(n_atoms_iter),
                                                                  cp.uint8(1)))
            # last argument is a Boolean for using solvent exclusion
        else:
            iasa_integrate((n_blocks, 1, 1,), (n_threads, 1, 1), (atoms, elements, b_factors, occupancies, potential,
                                                                  solvent, scattering, sz_potential_gpu, displacement,
                                                                  cp.float32(voxel_size), cp.uint32(n_atoms_iter),
                                                                  cp.uint8(0)))

        del atoms, elements, b_factors, occupancies
        mempool.free_all_blocks()

    # convert potential to correct units and correct for solvent exclusion
    if solvent_exclusion == 'gaussian':
        # Correct for solvent and convert both the solvent and potential array to the correct units.
        real = (potential * (C / dV)) - (solvent / dV * V_sol)
    elif solvent_exclusion == 'masking':  # only if voxel size is small enough for accurate determination of mask
        solvent_mask = (potential > 1E-5) * 1.0
        # construct solvent mask
        # gaussian decay of mask
        if oversampling == 1:
            solvent_mask = reduce_resolution_fourier(solvent_mask, voxel_size, voxel_size * 2)
            solvent_mask[solvent_mask < 0.001] = 0
        # subtract solvent from the protein electrostatic potential
        real = (potential * (C / dV)) - (solvent_mask * V_sol)
    else:
        real = potential * (C / dV)

    # determine absorption contrast if set
    if absorption_contrast:
        # voltage by default 300 keV
        molecule_absorption = physics.potential_amplitude(density, molecular_weight, voltage)
        solvent_absorption = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY,
                                                         physics.WATER_MW, voltage) * (V_sol / physics.V_WATER)
        print('Calculating absorption contrast')
        print(f'Molecule absorption = {molecule_absorption:.3f}')
        print(f'Solvent absorption = {solvent_absorption:.3f}')

        if solvent_exclusion == 'masking':
            imaginary = solvent_mask * (molecule_absorption - solvent_absorption)
        elif solvent_exclusion == 'gaussian':
            imaginary = solvent / dV * (molecule_absorption - solvent_absorption)
        else:
            print('ERROR: Absorption contrast cannot be generated if solvent exclusion is not set to either gaussian '
                  'or masking.')
            sys.exit(0)

        real = cp.pad(real, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                      mode='constant', constant_values=0)
        imaginary = cp.pad(imaginary, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                           mode='constant', constant_values=0)

        if oversampling > 1:
            print('Rescaling after oversampling')
            real = resize(reduce_resolution_fourier(real, voxel_size, voxel_size * 2 * oversampling),
                          1 / oversampling, interpolation='Spline')
            imaginary = resize(reduce_resolution_fourier(imaginary, voxel_size, voxel_size * 2 * oversampling),
                               1 / oversampling, interpolation='Spline')
        return real + 1j * imaginary
    else:
        # extend volume to a box
        real = cp.pad(real, ((difference[0] // 2, difference[0] // 2 + difference[0] % 2),
                             (difference[1] // 2, difference[1] // 2 + difference[1] % 2),
                             (difference[2] // 2, difference[2] // 2 + difference[2] % 2)),
                      mode='constant', constant_values=0)
        if oversampling > 1:
            print('Rescaling after oversampling')
            real = resize(reduce_resolution_fourier(real, voxel_size, voxel_size * 2 * oversampling),
                          1 / oversampling, interpolation='Spline')
        return real


# @jit(nopython=True) should use this? too many other function calls probably
def iasa_integration(filepath, voxel_size=1., oversampling=1, solvent_exclusion=None,
                     V_sol=physics.V_WATER, absorption_contrast=False, voltage=300E3, density=physics.PROTEIN_DENSITY,
                     molecular_weight=physics.PROTEIN_MW, structure_tuple=None):
    """
    Calculates interaction potential map to 1 A volume as described initially by Rullgard et al. (2011) in TEM
    simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013). This function applies averaging of
    the potential over the voxels to obtain precise results without oversampling.

    @param filepath: full filepath to pdb file
    @type  filepath: L{string}
    @param voxel_size: size of voxel in output map, default 1 A
    @type  voxel_size: L{float}
    @param oversampling: number of times to oversample final voxel size
    @type  oversampling: L{int}
    @param solvent_exclusion: flag to execute solvent exclusion using gaussian spheres. Solvent exclusion can be set
    with a string, either 'gaussian' or 'masking'. Default is None.
    @type  solvent_exclusion: L{str}
    @param V_sol: average solvent background potential (V/A^3)
    @type  V_sol: L{float}
    @param absorption_contrast: flag to generate absorption factor for imaginary part of potential
    @type  absorption_contrast: L{bool}
    @param voltage: electron beam voltage, absorption factor depends on voltage, default 300e3
    @type  voltage: L{float}
    @param density: average density of molecule that is generated, default 1.35 (protein)
    @type  density: L{float}
    @param molecular_weight: average molecular weight of the molecule that is generated, default protein MW
    @type  molecular_weight: L{float}
    @param structure_tuple: structure information as a tuple (x_coordinates, y_coordinates, z_coordinates, elements,
    b_factors, occupancies), if provided this overrides file reading
    @type  structure_tuple: L{tuple} - (L{list},) * 6 with types (float, float, float, str, float, float)

    @return: A volume with interaction potentials, either tuple of (real, imag) or single real, both real and imag
    are 3d arrays.
    @rtype: L{tuple} -> (L{np.ndarray},) * 2 or L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier
    from scipy.special import erf

    assert (type(oversampling) is int) and (oversampling >= 1), print('oversampling parameter is not an integer')
    if oversampling > 1:
        voxel_size /= oversampling

    extra_space = 30  # extend volume by 30 A in all directions

    print(f' - Calculating IASA potential from {filepath}')

    if structure_tuple is None:
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filepath)
    else:
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = structure_tuple

    x_max = np.max(x_coordinates - np.min(x_coordinates))
    y_max = np.max(y_coordinates - np.min(y_coordinates))
    z_max = np.max(z_coordinates - np.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)
    difference = [largest_dimension - a for a in dimensions]

    x_coordinates = x_coordinates - np.min(x_coordinates) + extra_space + difference[0]/2
    y_coordinates = y_coordinates - np.min(y_coordinates) + extra_space + difference[1]/2
    z_coordinates = z_coordinates - np.min(z_coordinates) + extra_space + difference[2]/2
    # Define the volume of the protein
    sz = (int((largest_dimension + 2 * extra_space) / voxel_size), ) * 3

    potential = np.zeros(sz)
    if solvent_exclusion:
        solvent = np.zeros(sz)

    print(f'Number of atoms to go over is {len(x_coordinates)}')

    for i in range(len(elements)):
        if np.mod(i, 5000) == 0:
            print(f'Calculating atom {i}.')

        atom = elements[i].upper()
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = np.array(physics.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b += (b_factor) # units in A

        if atom in list(physics.volume_displaced):
            r_0 = np.cbrt(physics.volume_displaced[atom] / (np.pi ** (3 / 2)))
        else:  # If not H,C,O,N we assume the same volume displacement as for carbon
            r_0 = np.cbrt(physics.volume_displaced['C'] / (np.pi ** (3 / 2)))

        r2 = 0  # this was 15 / (1 / r_0 ** 2) before but this gives very high values for carbon
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            r2 = np.maximum(r2, 15 / (4 * np.pi ** 2 / b[j]))
        # Radius of gaussian sphere
        r = np.sqrt(r2 / 3)

        rc = [x_coordinates[i], y_coordinates[i], z_coordinates[i]]  # atom center
        ind_min = [int((c - r) // voxel_size) for c in rc]  # Smallest index to contain relevant potential x,y,z
        ind_max = [int((c + r) // voxel_size) for c in rc]  # Largest index to contain relevant potential x,y,z
        # Explicit real space coordinates for the max and min boundary of each voxel
        x_min_bound = np.arange(ind_min[0], ind_max[0] + 1, 1) * voxel_size - rc[0]
        x_max_bound = np.arange(ind_min[0] + 1, ind_max[0] + 2, 1) * voxel_size - rc[0]
        y_min_bound = np.arange(ind_min[1], ind_max[1] + 1, 1) * voxel_size - rc[1]
        y_max_bound = np.arange(ind_min[1] + 1, ind_max[1] + 2, 1) * voxel_size - rc[1]
        z_min_bound = np.arange(ind_min[2], ind_max[2] + 1, 1) * voxel_size - rc[2]
        z_max_bound = np.arange(ind_min[2] + 1, ind_max[2] + 2, 1) * voxel_size - rc[2]

        atom_potential = 0

        for j in range(5):
            sqrt_b = np.sqrt(b[j])  # calculate only once
            # Difference of error function == integrate over Gaussian
            int_x = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(x_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(x_min_bound * 2 * np.pi / sqrt_b))
            x_matrix = np.tile(int_x[:, np.newaxis, np.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(y_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(y_min_bound * 2 * np.pi / sqrt_b))
            y_matrix = np.tile(int_y[np.newaxis, :, np.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = sqrt_b / (4 * np.sqrt(np.pi)) * (erf(z_max_bound * 2 * np.pi / sqrt_b) -
                                                     erf(z_min_bound * 2 * np.pi / sqrt_b))
            z_matrix = np.tile(int_z[np.newaxis, np.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])

            atom_potential += a[j] / b[j] ** (3 / 2) * x_matrix * y_matrix * z_matrix

        # scatter_add instead of +=
        potential[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += atom_potential

        if solvent_exclusion == 'gaussian':
            # excluded solvent potential
            int_x = np.sqrt(np.pi) * r_0 / 2 * (erf(x_max_bound / r_0) - erf(x_min_bound / r_0))
            x_matrix = np.tile(int_x[:, np.newaxis, np.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = np.sqrt(np.pi) * r_0 / 2 * (erf(y_max_bound / r_0) - erf(y_min_bound / r_0))
            y_matrix = np.tile(int_y[np.newaxis, :, np.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = np.sqrt(np.pi) * r_0 / 2 * (erf(z_max_bound / r_0) - erf(z_min_bound / r_0))
            z_matrix = np.tile(int_z[np.newaxis, np.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])

            solvent[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += (x_matrix *
                                                                                                         y_matrix *
                                                                                                         z_matrix)

    # Voxel volume
    dV = voxel_size ** 3
    # Convert to correct units
    C = 4 * np.sqrt(np.pi) * physics.constants['h'] ** 2 / (physics.constants['el'] * physics.constants['me']) * 1E20  # angstrom**2

    if solvent_exclusion == 'gaussian':
        # Correct for solvent and convert both the solvent and potential array to the correct units.
        real = (potential / dV * C) - (solvent / dV * V_sol)
    elif solvent_exclusion == 'masking': # only if voxel size is small enough for accurate determination of mask
        solvent_mask = (potential > 1E-5) * 1.0
        # construct solvent mask
        # gaussian decay of mask
        if oversampling == 1:
            solvent_mask = reduce_resolution_fourier(solvent_mask, voxel_size, voxel_size * 2)
            solvent_mask[solvent_mask < 0.001] = 0
        # fig, (ax1, ax2) = plt.subplots(1, 2)
        # slice = int(solvent_mask.shape[2] // 2)
        # ax1.imshow(potential[:, :, slice])
        # ax2.imshow(solvent_mask[:, :, slice])
        # show()
        real = (potential / dV * C) - (solvent_mask * V_sol)
    else:
        real = potential / dV * C

    if absorption_contrast:
        # voltage by default 300 keV
        molecule_absorption = physics.potential_amplitude(density, molecular_weight, voltage)
        solvent_absorption = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY,
                                                         physics.WATER_MW, voltage) * (V_sol / physics.V_WATER)
        print(f'molecule absorption = {molecule_absorption:.3f}')
        print(f'solvent absorption = {solvent_absorption:.3f}')

        if solvent_exclusion == 'masking':
            imaginary = solvent_mask * (molecule_absorption - solvent_absorption)
        elif solvent_exclusion == 'gaussian':
            imaginary = solvent/dV * (molecule_absorption - solvent_absorption)
        else:
            print('ERROR: Absorption contrast cannot be generated if the solvent masking or solvent exclusion method '
                  'are not used.')
            sys.exit(0)

        real = reduce_resolution_fourier(real, voxel_size, voxel_size*2*oversampling)
        real = resize(real, 1/oversampling, interpolation='Spline')
        imaginary = reduce_resolution_fourier(imaginary, voxel_size, voxel_size*2*oversampling)
        imaginary = resize(imaginary, 1/oversampling, interpolation='Spline')
        return real + 1j * imaginary
    else:
        real = reduce_resolution_fourier(real, voxel_size, voxel_size*2*oversampling)
        return resize(real, 1/oversampling, interpolation='Spline')


def iasa_rough(filepath, voxel_size=10, oversampling=1, solvent_exclusion=False, V_sol=physics.V_WATER):
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).
    This function applies averaging of the potential over the voxels BUT in a coarse way. Use only for a desired voxel
    spacing of 10 A.

    @param filepath: full path to pdb file
    @type  filepath: L{string}
    @param voxel_size: size (A) of voxel in output map, default 10 A
    @type  voxel_size: L{float}
    @param oversampling: number of times to oversample voxel size
    @type  oversampling: L{int}
    @param solvent_exclusion: flag to correct potential on each atom for solvent exclusion
    @type  solvent_exclusion: L{bool}
    @param V_sol: average solvent background potential (V/A^3)
    @type  V_sol: L{float}

    @return: a volume with electrostatic potential, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier

    assert (type(oversampling) is int) and (oversampling >= 1), print('oversampling parameter is not an integer')
    if oversampling > 1:
        voxel_size /= oversampling

    extra_space = 30  # extend volume by 10 A in all directions

    print(f' - Calculating IASA potential from {filepath}')

    x_coordinates, y_coordinates, z_coordinates, elements, _, _ = read_structure(filepath)

    x_max = np.max(x_coordinates - np.min(x_coordinates))
    y_max = np.max(y_coordinates - np.min(y_coordinates))
    z_max = np.max(z_coordinates - np.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)
    difference = [largest_dimension - a for a in dimensions]

    x_coordinates = x_coordinates - np.min(x_coordinates) + extra_space + difference[0]/2
    y_coordinates = y_coordinates - np.min(y_coordinates) + extra_space + difference[1]/2
    z_coordinates = z_coordinates - np.min(z_coordinates) + extra_space + difference[2]/2
    # Define the volume of the protein
    szx = np.abs(np.max(x_coordinates) - np.min(x_coordinates)) + 2 * extra_space + difference[0]
    szy = np.abs(np.max(y_coordinates) - np.min(y_coordinates)) + 2 * extra_space + difference[1]
    szz = np.abs(np.max(z_coordinates) - np.min(z_coordinates)) + 2 * extra_space + difference[2]
    sz = np.round(np.array([szx, szy, szz]) / voxel_size).astype(int)

    potential = np.zeros(sz)

    C = 2 * np.pi * physics.constants['h_bar'] ** 2 / (physics.constants['el'] * physics.constants['me']) * 1E20  # angstrom**2
    dV = voxel_size ** 3

    print(f'Number of atoms to go over is {len(x_coordinates)}.')

    for i in range(len(elements)):
        if np.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i].upper()

        rc = [x_coordinates[i], y_coordinates[i], z_coordinates[i]]  # atom center
        ind = tuple([int((c) // (voxel_size/oversampling)) for c in rc])  # Indexes to contain the potential

        potential[ind] += (np.array(physics.scattering_factors[atom]['g'])[0:5].sum() * C / dV)

        if solvent_exclusion:
            if atom in list(physics.volume_displaced):
                potential[ind] -= physics.volume_displaced[atom] * V_sol / dV
            else:  # If not H,C,O,N we assume the same volume displacement as for carbon
                potential[ind] -= physics.volume_displaced['C'] * V_sol / dV
    # Apply Gaussian filter
    potential = reduce_resolution_fourier(potential, 1, 2*oversampling)
    # Bin the volume
    potential = resize(potential, 1/oversampling, interpolation='Spline')

    return potential


def iasa_potential(filepath, voxel_size=1., oversampling=1): # add params voxel_size, oversampling?
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).
    In this algorithm the electrostatic potential is sampled from the sum of Gaussians function described on
    each atom. Can lead to incorrect sampling for large voxel sizes, requires voxel spacings of 0.25 A and smaller
    for correct result.

    @param filepath: full path to pdb file
    @type  filepath: L{string}
    @param voxel_size: size (A) of voxel in output map, default 1 A
    @type  voxel_size: L{float}
    @param oversampling: oversample potential function (multiple of 1), default 1 i.e. no oversampling
    @type  oversampling: L{int}

    @return: A volume with the electrostatic potential, 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier

    assert (type(oversampling) is int) and (oversampling >= 1), print('oversampling parameter is not an integer')
    if oversampling > 1:
        voxel_size /= oversampling
    print(f'IASA: volume will be sampled at {spacing}A')

    extra_space = 30  # extend volume by 10 A

    print(f' - Calculating IASA potential from {filepath}')

    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filepath)

    x_max = np.max(x_coordinates - np.min(x_coordinates))
    y_max = np.max(y_coordinates - np.min(y_coordinates))
    z_max = np.max(z_coordinates - np.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)
    difference = [largest_dimension - a for a in dimensions]

    x_coordinates = x_coordinates - np.min(x_coordinates) + extra_space + difference[0] / 2
    y_coordinates = y_coordinates - np.min(y_coordinates) + extra_space + difference[1] / 2
    z_coordinates = z_coordinates - np.min(z_coordinates) + extra_space + difference[2] / 2
    # Define the volume of the protein
    szx = np.abs(np.max(x_coordinates) - np.min(x_coordinates)) + 2 * extra_space + difference[0]
    szy = np.abs(np.max(y_coordinates) - np.min(y_coordinates)) + 2 * extra_space + difference[1]
    szz = np.abs(np.max(z_coordinates) - np.min(z_coordinates)) + 2 * extra_space + difference[2]

    sz = np.round(np.array([szx, szy, szz]) / voxel_size).astype(int)

    potential = np.zeros(sz)

    # C = 2132.8 A^2 * V; 1E20 is a conversion factor for Angstrom^2,
    # h and not h_bar, which is why we multiply with 4 * sqrt(pi) instead of 16*pi^(5/2)
    C = 4 * np.sqrt(np.pi) * physics.constants['h']**2 / (physics.constants['el'] * physics.constants['me']) * 1E20

    print(f'#atoms to go over is {len(x_coordinates)}.')

    for i in range(len(elements)):
        if np.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i].upper()
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = np.array(physics.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b is used for calculating the radius of the potential. See Rullgard et al. (2011) for addition of 16 R^2.
        # Incorporates low pass filtering by extending the gaussian radius. NOTE: b_factor is zero anyway
        b += (b_factor) # + 16 * spacing**2)

        r2 = 0
        b1 = np.zeros(5)
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            b1[j] = 4 * np.pi ** 2 / b[j] * voxel_size ** 2
            r2 = np.maximum(r2, 10/b1[j])
        # Radius of gaussian sphere
        r = np.sqrt(r2 / 3)
        xc1, yc1, zc1 = x_coordinates[i] / voxel_size, y_coordinates[i] / voxel_size, z_coordinates[i] / voxel_size
        # Center of the gaussian sphere.
        rc = [xc1, yc1, zc1]
        # Calculate the absolute indexes for the potential matrix.
        kmin = [np.maximum(0,x).astype(int) for x in np.ceil(rc-r)]
        kmax = [np.minimum(np.floor(x)-1,np.floor(y+r)).astype(int) for x,y in zip(sz,rc)]
        kmm = max([x-y for x,y in zip(kmax,kmin)])
        # Determine the coordinates for sampling from the sum of Gaussians.
        x = xc1 - np.arange(kmin[0], kmin[0]+kmm+1, 1)
        y = yc1 - np.arange(kmin[1], kmin[1]+kmm+1, 1)
        z = zc1 - np.arange(kmin[2], kmin[2]+kmm+1, 1)

        atom_potential = 0
        for j in range(5):
            x_matrix = np.tile(np.exp(-b1[j] * x**2)[:,np.newaxis,np.newaxis], [1, kmm+1, kmm+1])
            y_matrix = np.tile(np.exp(-b1[j] * y**2)[np.newaxis,:,np.newaxis], [kmm+1, 1, kmm+1])
            z_matrix = np.tile(np.exp(-b1[j] * z**2)[np.newaxis,np.newaxis,:], [kmm+1, kmm+1, 1])
            tmp = a[j] / b[j]**(3/2) * C * x_matrix * y_matrix * z_matrix
            atom_potential += tmp
        atom_potential *= occupancy
        # add the potential of element i to the full potential map
        potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] = \
            potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] + atom_potential

    # subtract the average background of amorphous ice
    potential -= physics.V_WATER

    # using spline interpolation here does not decrease our accuracy as the volume is already oversampled!
    potential = reduce_resolution_fourier(potential, 1, 2*oversampling)
    potential = resize(potential, 1/oversampling, interpolation='Spline')

    return potential


def parse_apbs_output(filepath):
    """
    Parses output file from APBS to a 3d volume, and returns the volume along the voxel spacing in x, y, and z.

    @param filepath: full path to .pqr.dx file produced by APBS software
    @type  filepath: L{str}

    @return: electrostatic shifts calculated by APBS as a 3d array, and the voxel spacing dx, dy, dz in A
    @rtype:  L{tuple} - (L{np.ndarray}, L{float}, L{float}, L{float})

    @author: Marten Chaillet
    """
    # Read file header
    try:
        with open(filepath, 'r') as apbs_file:
            grid, origin, spacing = [],[],[]
            delta_count = 0
            for i in range(11):
                line = apbs_file.readline()
                if 'object 1 class gridpositions counts' in line:
                    grid.extend([int(x) for x in line.split()[5:8]])
                elif 'origin' in line:
                    origin.extend([float(x) for x in line.split()[1:4]])
                elif 'delta' in line:
                    spacing.append(float(line.split()[1+delta_count]))
                    delta_count+=1
            print(f'grid\t= {grid}')
            print(f'origin\t= {origin}')
            print(f'spacing\t= {spacing}')

            line = apbs_file.readline().split()
            data = []
            while line[0] != 'attribute':
                data.extend(line)
                line = apbs_file.readline().split()
    except Exception as e:
        print(e)
        raise Exception('Could not open APBS data file.')

    # Reshape to sequence
    data = np.array(data,dtype=float)
    # Form to 3D grid
    data = np.reshape(data,(grid[2],grid[1],grid[0]),order='F').transpose(2,1,0) # x, y, z needs to be done as such to correspond to iasa
    dx = spacing[0]
    dy = spacing[1]
    dz = spacing[2]
    return data, dx, dy, dz


def resample_apbs(filepath, voxel_size=1.0):
    """
    First calls parse_abps_output to read an apbs output file, then scales voxels to voxel_size and
    refactors the values to volt.

    @param filepath: full filepath
    @type  filepath: L{str}
    @param voxel_size: desired voxel size after resampling in A, default is 1 A
    @type  voxel_size: L{float}

    @return: Resampled and scaled V_bond potential, 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from skimage.transform import rescale
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier

    print(f' - Parsing and resampling APBS file {filepath}')

    # Parse APBS data file
    potential, dxnew, dynew, dznew = parse_apbs_output(filepath)
    # Check if the voxel size is allowed, and adjust if neccessary
    smallest_possible = max([dxnew, dynew, dznew])
    if not(voxel_size >= smallest_possible):
        print(f'Requested voxel size is smaller than voxel size of the apbs map. Adjust to smallest possible voxel size'
              f' of {smallest_possible}.')
        voxel_size = smallest_possible
    # Make voxels same size by scaling to dx
    factor = (dxnew/dxnew, dynew/dxnew, dznew/dxnew)
    # Use skimage rescale for different factors along x, y, and z
    # order 3 for bi-cubic (splines), preserve_range otherwise image is returned as float between -1 and 1
    potential = rescale(potential, factor, mode='constant', order=3, preserve_range=True, multichannel=False,
                        anti_aliasing=False)
    # filter volume before downsampling to voxel_size
    potential = reduce_resolution_fourier(potential, dxnew, 2*voxel_size)
    # Scale the volume to voxels with voxel_size
    potential = resize(potential, dxnew/voxel_size, interpolation='Spline')

    print(f'Data after reshaping to {voxel_size} A voxels: {potential.shape}')

    # Convert to from kT/e to volts
    temperature = 291 # [K] = 18 C (room temperature)
    convert_to_volts = physics.constants['kb'] * temperature / physics.constants['el']
    potential = convert_to_volts * potential
    return potential


def combine_potential(potential1, potential2):
    """
    Combine V_atom and V_bond potential. Potentials do not need to have the same shape, but do need to have the same
    voxel spacing. Also the potential should be centered in the volume identically.

    @param potential1: first potential, 3d array
    @type  potential1: L{np.ndarray}
    @param potential2: second potential, 3d array
    @type  potential2: L{np.ndarray}

    @return: Combined V_atom and V_bond, 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.simulation.support import reduce_resolution_fourier

    print(' - Combining iasa and bond potential')

    size1 = potential1.shape
    size2 = potential2.shape

    # determine the extension for both volumes
    ext1 = [(s2 - s1 if s2 > s1 else 0) for s1, s2 in zip(size1, size2)]
    ext2 = [(s1 - s2 if s1 > s2 else 0) for s1, s2 in zip(size1, size2)]

    # extend volumes and filter before interpolating
    potential1_ext = extend_volume(reduce_resolution_fourier(potential1, 1, 2),
                                   ext1, symmetrically=True, true_center=True)
    potential2_ext = extend_volume(reduce_resolution_fourier(potential2, 1, 2),
                                   ext2, symmetrically=True, true_center=True)

    potential_combined = potential1_ext + potential2_ext

    return potential1_ext, potential2_ext, potential_combined


def wrapper(filepath, output_folder, voxel_size, oversampling=1, binning=1, solvent_exclusion=None,
            solvent_potential=physics.V_WATER, absorption_contrast=False, voltage=300E3, solvent_factor=1.0, cores=1,
            gpu_id=None):
    """
    Execution of generating an electrostatic potential (and absorption potential) from a pdb/cif file. Process
    includes preprocessing with chimera to add hydrogens and symmetry, then passing to IASA_intergration method to
    correctly sample the electrostatic potential to a 3d array. Two options can be provided for solvent correction.

    @param filepath: full path to pdb or cif filed
    @type  filepath: L{str}
    @param output_folder: folder to write all output to
    @type  output_folder: L{str}
    @param voxel_size: voxel size in A to sample the interaction potential to.
    @type  voxel_size: L{float}
    @param oversampling: number of times to oversample the interaction potential for better accuracy, multiple of 1
    @type  oversampling: L{int}
    @param binning: number of times to bin the volume after sampling, this file will be saved separately
    @type  binning: L{int}
    @param exclude_solvent: flag to exclude solvent with a Gaussian sphere
    @type  exclude_solvent: L{bool}
    @param solvent_masking: flag to excluded solvent by masking (thresholding method)
    @type  solvent_masking: L{bool}
    @param solvent_potential: background solvent potential, default 4.5301
    @type  solvent_potential: L{float}
    @param absorption_contrast: flag for generating absorption potential
    @type  absorption_contrast: L{bool}
    @param voltage: electron beam voltage in eV, parameter for absorption contrast, default 300e3
    @type  voltage: L{float}
    @param solvent_factor: solvent factor to increase background potential
    @type  solvent_factor: L{float}

    @return: - (files are written to output_folder)
    @rtype:  Nonee

    @author: Marten Chaillet
    """
    from pytom.agnostic.io import write
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier

    # Id does not makes sense to apply absorption contrast if solvent exclusion is not turned on
    if absorption_contrast:
        assert solvent_exclusion is not None, print('absorption contrast can only be applied if solvent exclusion is '
                                                    'used.')

    _, file = os.path.split(filepath)
    pdb_id, _ = os.path.splitext(file)

    # output_folder = os.path.join(output_folder, pdb_id)
    # if not os.path.exists(output_folder):
    #     print(f'Making folder {output_folder}...')
    #     os.mkdir(output_folder)

    # Call external programs for structure preparation and PB-solver
    # call_apbs(folder, structure, ph=ph)
    try:
        filepath = call_chimera(filepath, os.path.split(filepath)[0])  # output structure name is dependent on
        # modification by chimera
    except Exception as e:
        print(e)

    assert filepath != 0, 'something went wrong with chimera'

    # Calculate atom and bond potential, and store them
    # 4 times oversampling of IASA yields accurate potentials
    # Could be {structure}.{extension}, but currently chimera is only able to produce .pdb files, so the extended
    # structure file created by call chimera has a .pdb extension.
    if gpu_id is not None:
        v_atom = iasa_integration_gpu(filepath, voxel_size=voxel_size, oversampling=oversampling,
                                  solvent_exclusion=solvent_exclusion,
                                  V_sol=solvent_potential * solvent_factor, absorption_contrast=absorption_contrast,
                                  voltage=voltage, gpu_id=gpu_id)
    else:
        v_atom = iasa_integration_parallel(filepath, voxel_size=voxel_size, oversampling=oversampling,
                                  solvent_exclusion=solvent_exclusion,
                                  V_sol=solvent_potential * solvent_factor, absorption_contrast= absorption_contrast,
                                  voltage=voltage, cores=cores)

    # Absorption contrast map generated here will look blocky when generated at 2.5A and above!
    if np.iscomplexobj(v_atom):
        output_name = f'{pdb_id}_{voxel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
        print(f'writing real and imaginary part with name {output_name}')
        write(os.path.join(output_folder, f'{output_name}_real.mrc'), v_atom.real)
        write(os.path.join(output_folder, f'{output_name}_imag_{voltage*1E-3:.0f}V.mrc'), v_atom.imag)
    else:
        if solvent_exclusion:
            output_name = f'{pdb_id}_{voxel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
        else:
            output_name = f'{pdb_id}_{voxel_size:.2f}A'
        print(f'writing real part with name {output_name}')
        write(os.path.join(output_folder, f'{output_name}_real.mrc'), v_atom)

    if binning > 1:
        # v_atom_binned = iasa_integration(f'{output_folder}/{structure}.pdb', voxel_size=voxel_size*binning,
        #                       solvent_exclusion=exclude_solvent, V_sol=solvent_potential)
        # first filter the volume!
        if np.iscomplexobj(v_atom):
            print(' - Binning volume')
            filtered = [reduce_resolution_fourier(v_atom.real, voxel_size, voxel_size * 2 * binning),
                        reduce_resolution_fourier(v_atom.imag, voxel_size, voxel_size * 2 * binning)]
            binned = [resize(filtered.real, 1/binning, interpolation='Spline'),
                      resize(filtered.imag, 1/binning, interpolation='Spline')]

            output_name = f'{pdb_id}_{voxel_size*binning:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            print(f'writing real and imaginary part with name {output_name}')
            write(os.path.join(output_folder, f'{output_name}_real.mrc'), binned[0])
            write(os.path.join(output_folder, f'{output_name}_imag_{voltage*1E-3:.0f}V.mrc'), binned[1])

        else:
            print(' - Binning volume')
            filtered = reduce_resolution_fourier(v_atom, voxel_size, 2 * voxel_size * binning)
            binned = resize(filtered, 1/binning, interpolation='Spline')
            if solvent_exclusion:
                output_name = f'{pdb_id}_{voxel_size*binning:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            else:
                output_name = f'{pdb_id}_{voxel_size*binning:.2f}A'
            print(f'writing real part with name {output_name}')
            write(os.path.join(output_folder, f'{output_name}_real.mrc'), binned)
    return


if __name__ == '__main__':
    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA, PDB2PQR (modified), APBS

    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2

    # syntax is ScriptOption([short, long], description, requires argument, is optional)
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Calculate electrostatic potential from protein structure file.\n Script has dependencies on '
                    'pytom, chimera. (pdb2pqr apbs)',
        authors='Marten Chaillet',
        options=[ScriptOption2(['-f', '--file'], 'File path with protein structure, either pdb or cif.', 'file',
                               'required'),
               ScriptOption2(['-d', '--destination'], 'Folder to store the files produced by potential.py.',
                             'directory', 'optional', './'),
               ScriptOption2(['-s', '--spacing'], 'The size of the voxels of the output volume. 1A by default.',
                             'float', 'required'),
               ScriptOption2(['-n', '--oversampling'], 'n times pixel size oversampling.', 'int', 'optional', 1),
               ScriptOption2(['-b', '--binning'], 'Number of times to bin. Additional storage of binned volume.',
                             'int', 'optional', 1),
               ScriptOption2(['-x', '--exclude_solvent'],
                             'Whether to exclude solvent around each atom as a correction of the potential, '
                             'either "gaussian" or "masking".', 'string',
                             'optional'),
               ScriptOption2(['-p', '--solvent_potential'],
                             f'Value for the solvent potential. By default amorphous ice, {physics.V_WATER} V.',
                             'float', 'optional', physics.V_WATER),
               ScriptOption2(['-c', '--percentile'], 'Multiplication for solvent potential and absorption contrast of '
                                                    'solvent to decrease/increase contrast. Value between 0 and 3 (could'
                                                    ' be higher but would not make sense).', 'float', 'optional', 1),
               ScriptOption2(['-a', '--absorption_contrast'],
                             'Whether to add imaginary part of absorption contrast, can only be done if solvent'
                             'is excluded.', 'no arguments', 'optional'),
               ScriptOption2(['-v', '--voltage'],
                             'Value for the electron acceleration voltage. Need for calculating the inelastic mean '
                             'free '
                             'path in case of absorption contrast calculation. By default 300 (keV).', 'float',
                             'optional', 300),
               ScriptOption2(['--cores'], 'Number of cpu cores to use for the calculation.', 'int', 'optional', 1),
               ScriptOption2(['-g', '--gpuID'], 'GPU index to run the program on.', 'int', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    filepath, output_folder, voxel_size, oversampling, binning, solvent_exclusion, solvent_potential, solvent_factor, \
        absorption_contrast, voltage, cores, gpuID = options

    voltage *= 1e3

    wrapper(filepath, output_folder, voxel_size, oversampling=oversampling, binning=binning,
            solvent_exclusion=solvent_exclusion, solvent_potential=solvent_potential,
            absorption_contrast=absorption_contrast, voltage=voltage, solvent_factor=solvent_factor, cores=cores,
            gpu_id=gpuID)


