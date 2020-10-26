import numpy as xp
import math
import constant_dictionaries as phys
import os
# import pytom.basic.functions

# image display
import matplotlib as plt
plt.use('Qt5Agg')
from pylab import *

V_WATER = 4.5301 # potential value of low density amorphous ice (from Vulovic et al., 2013)

def extend_volume(vol, increment, pad_value=0, symmetrically=False, true_center=False):
    """
    Increase volume by parameter increment ([x,y,z]). Options for changing the padding value, extending symmetrically on
    both sides of the input volume, shifting the original volume to the true new center if the increment/2 is
    non-integer.

    @param vol: 3D matrix
    @type vol: numpy ndarray
    @param increment: list with increment value for each dimension, only integers
    @type increment: [int, int, int]
    @param pad_value: value to use as padding
    @type pad_value: float
    @param symmetrically: If False (default) the volume is just padded with zeros.
    @type symmetrically: Bool
    @param true_center: If True interpolate to true center
    @type true_center: Bool

    @return: Extended volume
    @rtype: numpy ndarray

    @author: Marten Chaillet
    """
    if symmetrically:
        # condition = any([x%2 for x in increment])
        if true_center:
            import scipy.ndimage
            vol = xp.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)
            # TODO Use voltools for this as well!
            # Filter removes potential artifacts from the interpolation
            return scipy.ndimage.interpolation.shift(vol, tuple([x / 2 for x in increment]), order=3)
        else:
            from pytom.tompy.tools import paste_in_center
            new = xp.zeros([a + b for a, b in zip(vol.shape, increment)])
            if pad_value:
                new += pad_value
            return paste_in_center(vol, new)
    else:
        return xp.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)


# def create_gaussian_low_pass(size, radius, center=None):
#     """
#     Create a 3D mask with gaussian edges.
#
#     @param size: size of the resulting volume.
#     @param cutoff: radius of the sphere inside the volume.
#     @param center: sigma of the Gaussian.
#     @param gpu: center of the sphere.
#
#     @return: sphere inside a volume.
#     """
#     assert len(size) == 3
#
#     if center is None:
#         center = [size[0]/2, size[1]/2, size[2]/2]
#
#     [x,y,z] = xp.mgrid[0:size[0], 0:size[1], 0:size[2]]
#     r = xp.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
#
#     filter = xp.exp(-(r ** 2) / (2 * radius ** 2))
#
#     return filter


def create_gaussian_low_pass(size, hwhm, center=None):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER
    Create a 3D mask with gaussian edges.

    @param size: size of the resulting volume.
    @param hwhm: cutoff value for the filter.
    @param center: center of the Gaussian. If not provided will be calculated.

    @return: sphere inside a volume.
    """
    assert len(size) == 3

    # full width at half maximum is two times the half width at half maximum (or cutoff)
    sigma = 2*hwhm / 2.355 # obtain sigma by dividing by 2 * sqrt( 2 * ln(2) )
    sigma *= 0.75 # factor for obtaining functioning filter? 0.8 gives the best results

    if center is None:
        center = [(size[0]-1)/2, (size[1]-1)/2, (size[2]-1)/2]

    [x,y,z] = xp.mgrid[0:size[0], 0:size[1], 0:size[2]]
    r = xp.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)

    filter = xp.exp(-r ** 2 / (2 * sigma ** 2))

    return filter


def reduce_resolution(map, voxel_size, resolution):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER. Tompy does not contain gaussian low pass in fourier space.
    Accepts only square input maps to easy low-pass filtering

    @param map: Volume
    @type map: 3d numpy array
    @param voxel_size: Voxel size of the volume in A
    @type voxel_size: float
    @param to_resolution: Desired resolution after applying low-pass filter
    @type to_resolution: float
    @param sharpness: ratio
    @type sharpness: float

    @return: Filtered volume
    @rtype: 3d numpy array

    @author: Marten Chaillet
    """
    # from pytom.tompy.transform import fourier_filter

    assert resolution >= voxel_size, "the requested resolution is non-valid as it is smaller than the voxel size"
    assert len(set(map.shape)) == 1, "dimensions of input are not equal"

    # resolution reduction factor
    nr_pixels_cutoff = ( map.shape[0] * voxel_size ) / resolution
    # NOTE: if resolution == 2*voxel_size, filter radius will be 0.5 * map.shape[0]. I.E. no filtering. This is correct
    # because the highest possible resolution (nyquist) equals 2*voxel_size.
    # create full gaussian mask
    mask = create_gaussian_low_pass( map.shape, nr_pixels_cutoff )
    # apply mask in fourier space
    # result = fourier_filter(map, mask, human=True)
    ft_map = xp.fft.fftn( xp.fft.ifftshift( map ) ) # shift center to origin
    result = xp.fft.fftshift( xp.fft.ifftn( ft_map * xp.fft.ifftshift( mask) ) ) # shift center of mask to origin
    return result.real


def reduce_resolution_2(map, voxel_size, resolution):

    from scipy.ndimage import fourier_gaussian

    input = xp.fft.fftn(map)
    result = fourier_gaussian(input, sigma=(resolution/(2*voxel_size))/2.35)
    return xp.fft.ifftn(result).real


def low_pass_filter(map, voxel_size, resolution):
    """
    Accepts only square input maps to easy low-pass filtering

    @param map: Volume
    @type map: 3d numpy array
    @param voxel_size: Voxel size of the volume in A
    @type voxel_size: float
    @param to_resolution: Desired resolution after applying low-pass filter
    @type to_resolution: float
    @param sharpness: ratio
    @type sharpness: float

    @return: Filtered volume
    @rtype: 3d numpy array

    @author: Marten Chaillet
    """
    from pytom.tompy.transform import fourier_filter
    from pytom.tompy.tools import create_sphere

    assert resolution > voxel_size, "the requested resolution is non-valid as it is smaller than the voxel size"
    assert len(set(map.shape)) == 1, "dimensions of input are not equal"

    # calculate radius of mask from pixels needed in fourier space for the desired resolution
    nr_pixels_fourier = (map.shape[0] * voxel_size) / resolution
    mask = create_sphere(map.shape, radius=nr_pixels_fourier, sigma=15, num_sigma=10)
    # apply mask in fourier space
    result = fourier_filter(map, mask, human=True)
    return result


def scale(volume, factor):
    """
    Scale volumes with skimage.
    @param potential: ndimage, numpy array
    @param factor: float, tuple of floats. Single float for same scaling in each dimension. Tuple for different factors
    @return: scaled multidimensional array (numpy)
    """
    # Skimage scale could be better for downsampling than scipy zoom
    import skimage
    # order 1 would be bi-linear, maybe better?
    # order 3 for bi-cubic (splines), preserve_range otherwise image is returned as float between -1 and 1
    return skimage.transform.rescale(volume, factor, mode='constant', order=3, preserve_range=True, multichannel=False,
                                     anti_aliasing=False)


def bin(potential, factor):
    """

    @param potential:
    @param factor: integer value
    @return:
    """
    assert type(factor) is int and factor >= 1, print('non-valid binning factor, should be integer above 1')

    if factor == 1:
        return potential

    s = (potential.shape[0] % factor) // 2
    d = (potential.shape[0] % factor) % 2

    potential = potential[s:potential.shape[0] - s - d, s:potential.shape[0] - s - d, s:potential.shape[0] - s - d]

    ds = int(potential.shape[0] // factor)
    image_size = potential.shape[0]

    binned = potential.reshape(ds, image_size // ds,
                               ds, image_size // ds, ds, image_size // ds).mean(-1).mean(1).mean(-2)

    return binned


def call_chimera(file, input_folder, output_folder):
    """
    Run chimera for pdb file in order to add hydrogens and add symmetry units. The resulting structure is stored as a
    new pdb file with the extension {id}_symmetry.pdb

    @param pdb_folder: Path to a folder where pdb files are stored
    @type pdb_folder: string
    @param pdb_id: ID of pdb file
    @type pdb_id: string

    @return: Returns the name of file where the new pdb stucture is stored in. This can differ for pdb depending on
    symmetry thus we need to return it.
    @rtype: string

    @author: Marten Chaillet
    """
    print(f' - Calling chimera for {file}')

    pdb_id = file.split('.')[0]
    extension = file.split('.')[-1]

    # Skip if output files from this function already exist
    if os.path.exists(f'{output_folder}/{pdb_id}_rem-solvent_sym_addh.pdb'):
        print(f'File already exists: {pdb_id}_rem-solvent_sym_addh.pdb')
        return f'{pdb_id}_rem-solvent_sym_addh'
    elif os.path.exists(f'{output_folder}/{pdb_id}_rem-solvent_addh.pdb'):
        print(f'File already exists: {pdb_id}_rem-solvent_addh.pdb')
        return  f'{pdb_id}_rem-solvent_addh'

    if extension == 'pdb':
        # Command 'sym' in chimera crashes when there is no BIOMT symmetry in the pdb file. We need to make sure sym is only
        # executed when the BIOMT information specifies symmetrical units.
        symmetry = []
        try:
            with open(f'{input_folder}/{file}','r') as pdb:
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
            scriptname = f'_rem-solvent_sym_addh_{pdb_id}'
            try:
                with open(f'{output_folder}/{scriptname}.py', 'w') as chimera_script:
                    chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                         f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                         f'from chimera import runCommand as rc\n'
                                         f'rc("open {input_folder}/{file}")\n'
                                         f'rc("delete solvent")\n'
                                         f'rc("delete ions")\n'
                                         f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                         f'rc("sym group biomt")\n'             # group biomt is also the default
                                         f'rc("combine all modelId 10")\n'
                                         f'rc("write format pdb #10 {output_folder}/{pdb_id}_rem-solvent_sym_addh.pdb")\n'
                                         f'rc("stop")\n')
            except Exception as e:
                print(e)
                raise Exception('Could not create chimera script.')
        else:
            scriptname = f'_rem-solvent_addh_{pdb_id}'
            try:
                with open(f'{output_folder}/{scriptname}.py', 'w') as chimera_script:
                    chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                         f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                         f'from chimera import runCommand as rc\n'
                                         f'rc("open {input_folder}/{file}")\n'
                                         f'rc("delete solvent")\n'
                                         f'rc("delete ions")\n'
                                         f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                         f'rc("write format pdb #0 {output_folder}/{pdb_id}_rem-solvent_addh.pdb")\n'
                                         f'rc("stop")\n')
            except Exception as e:
                print(e)
                raise Exception('Could not create chimera script.')
        # module chimera should be loaded here...
        try:
            os.system(f'chimera --nogui --script {output_folder}/{scriptname}.py')
        except Exception as e:
            print(e)
            raise Exception('Chimera is likely not on your current path.')

        if len(set(symmetry)) > 1:
            return f'{pdb_id}_rem-solvent_sym_addh' # returns new pdb name
        else:
            return f'{pdb_id}_rem-solvent_addh'

    elif extension == 'cif':
        # If cif, assume the file does not contain any structural symmetry information as this is usually not the case
        # for large complexes.
        # Do transformations with chimera, and store as pdb. Chimera cannot write mmCIF files. ChimeraX can.
        # Also not that opening mmCIF files in chimera takes significantly longer.
        scriptname = f'_rem-solvent_addh_{pdb_id}'
        try:
            with open(f'{output_folder}/{scriptname}.py', 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {input_folder}/{file}")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     f'rc("addh")\n' # If using pdb2pqr do not add hydrogens here.
                                     f'rc("write format pdb #0 {output_folder}/{pdb_id}_rem-solvent_addh.pdb")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
        # module chimera should be loaded here...
        try:
            os.system(f'chimera --nogui --script {output_folder}/{scriptname}.py')
        except Exception as e:
            print(e)
            raise Exception('Chimera is likely not on your current path.')
        return f'{pdb_id}_rem-solvent_addh'
    else:
        print('non-valid structure file extension')
        return 0


def modify_structure_file(filename, pattern, replacement, line_start=''):
    """
    Function required to make pqr files with large negative coordinates readible for APBS. APBS can only parse
    columns in the file when they are properly separated by white spaces. Input file will be overwritten.

    @param filename: File path of pqr type file (with extension)
    @type filename: string
    @param pattern: Pattern to be replaced
    @type pattern: string
    @param replacement: Replacement string for pattern
    @type replacement: string
    @param line_start: Keyword argument, only modify line starting with line_start
    #type line_start: string

    @return: empty
    @rtype:

    @author: Marten Chaillet
    """
    from tempfile import mkstemp
    from shutil import move, copymode

    try:
        # Create temp file
        fh, abs_path = mkstemp()
        with os.fdopen(fh, 'w') as new_file:
            with open(filename) as old_file:
                for line in old_file:
                    if line_start=='':
                        new_file.write(line.replace(pattern,replacement))
                    elif line.split()[0] == line_start:
                        new_file.write(line.replace(pattern,replacement))
                    else:
                        new_file.write(line)
        # Copy the file permissions from the old file to the new file
        copymode(filename, abs_path)
        # Remove original file
        os.remove(filename)
        # Move new file
        move(abs_path, filename)
    except Exception as e:
        print(e)
        raise Exception('Unsuccessful in modifying pqr file with white space delimiters.')
    return


def call_apbs(folder, structure, force_field='amber', ph=7.):
    """
    Calls external programs pdb2pqr and apbs to execute on pdb structure. References:

    @param pdb_folder: Folder where pdb structures are stored
    @type pdb_folder: string
    @param structure: Name of pdb file with coordinates of the protein structure
    @type structure: string
    @param apbs_folder: Folder to store output of pdb2pqr and apbs
    @type apbs_folder: string
    @param force_field: Force field for parameterizing atoms (option: amber, ...)
    @type force_field: string
    @param ph: pH value of solvent surrounding the protein
    @type ph: float

    @return: empty, output of programs called is stored in apbs_folder
    @rtype:

    @author: Marten Chaillet
    """
    # pdb2pqr and apbs should be on path for this function to run
    print(f' - Running pdb2pqr and APBS on {folder}/{structure}.pdb')
    cwd = os.getcwd()
    input_file = f'{folder}/{structure}.pdb'

    output_file = f'{folder}/{structure}.pqr'
    apbs_in = f'{folder}/{structure}.in'
    try:
        # Also PDB2PKA ph calculation method. Requires PARSE force field, can take very long for large proteins.
        os.system(f'pdb2pqr.py --ff={force_field} --ph-calc-method=propka --with-ph={ph} --apbs-input {input_file} {output_file}')
        print(' - Add white space delimiters to pqr file.')
        modify_structure_file(output_file, '-', ' -', line_start='ATOM')
        # APBS needs to execute from the folder where the structure is present
        os.chdir(folder)
        os.system(f'apbs {apbs_in}')
        # os.system(f'module load apbs_mc/1.5; apbs {apbs_in}')
        # ? subprocess.run(['apbs', f'{structure}.in'], cwd=output_folder) # executes in cwd
    except Exception as e:
        print(e)
        raise Exception('pdb2pqr or APBS potentially not on path.')

    # Change back to original directory
    os.chdir(cwd)
    return


def read_structure(filename):
    """
    Read pdb, cif, or pqr file and return atom data in lists.
    @param filename:
    @return:
    """
    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = [], [], [], [], [], []

    filetype = filename.split('.')[-1]

    if filetype == 'pdb':
        try:
            with open(filename, 'r') as pdb:
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
    elif filetype == 'cif':
        try:
            with open(filename, 'r') as pdb:
                lines = pdb.readlines()
                for line in lines:
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
            raise Exception('Could not read pdb file.')
    elif filetype == 'pqr':
        try:
            with open(filename, 'r') as pqr:
                lines = pqr.readlines()
                for line in lines:
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
            raise Exception('Could not read pdb file.')
    else:
        print('invalid filetype in iasa_potential, return 0')
        return 0

    return x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies


def iasa_integration(filename, voxel_size=1., solvent_exclusion=False, solvent_masking=False, V_sol=4.5301,
                     absorption_contrast=False, voltage=300E3, density=1.35, molecular_weight=7.2, solvent_factor=1.0,
                     structure_tuple=None):
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).
    This function applies averaging of the potential over the voxels to obtain precise results without oversampling.

    @param filename: ID of pdb file as present in pdb folder
    @type filename: string
    @param voxel_size: Size (A) of voxel in output map, default 1 A
    @type voxel_size: float
    @param solvent_exclusion: model the potential taking into account solvent exclusion
    @type solvent_exclusion: Bool
    @param V_sol: average solvent background potential (V/A^3)
    @type V_sol: float

    @return: A volume with interaction potentials
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    from scipy.special import erf

    extra_space = 20  # extend volume by 30 A in all directions

    print(f' - Calculating IASA potential from {filename}')

    if structure_tuple is None:
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filename)
    else:
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = structure_tuple

    x_max = xp.max(x_coordinates - xp.min(x_coordinates))
    y_max = xp.max(y_coordinates - xp.min(y_coordinates))
    z_max = xp.max(z_coordinates - xp.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)
    difference = [largest_dimension - a for a in dimensions]

    x_coordinates = x_coordinates - xp.min(x_coordinates) + extra_space + difference[0]/2
    y_coordinates = y_coordinates - xp.min(y_coordinates) + extra_space + difference[1]/2
    z_coordinates = z_coordinates - xp.min(z_coordinates) + extra_space + difference[2]/2
    # Define the volume of the protein
    szx = xp.abs(xp.max(x_coordinates) - xp.min(x_coordinates)) + 2 * extra_space + difference[0]
    szy = xp.abs(xp.max(y_coordinates) - xp.min(y_coordinates)) + 2 * extra_space + difference[1]
    szz = xp.abs(xp.max(z_coordinates) - xp.min(z_coordinates)) + 2 * extra_space + difference[2]
    sz = xp.round(xp.array([szx, szy, szz]) / voxel_size).astype(int)

    potential = xp.zeros(sz)
    if solvent_exclusion:
        solvent = xp.zeros(sz)

    for i in range(len(elements)):
        if xp.mod(i, 5000) == 0:
            print(f'Calculating atom {i}.')

        atom = elements[i].upper()
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = xp.array(phys.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b += (b_factor) # units in A

        if atom in list(phys.volume_displaced):
            r_0 = xp.cbrt(phys.volume_displaced[atom] / (xp.pi ** (3 / 2)))
        else:  # If not H,C,O,N we assume the same volume displacement as for carbon
            r_0 = xp.cbrt(phys.volume_displaced['C'] / (xp.pi ** (3 / 2)))

        r2 = 15 / (1 / r_0 ** 2)
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            r2 = xp.maximum(r2, 15 / (4 * xp.pi ** 2 / b[j]))
        # Radius of gaussian sphere
        r = xp.sqrt(r2 / 3)

        rc = [x_coordinates[i], y_coordinates[i], z_coordinates[i]]  # atom center
        ind_min = [int((c - r) // voxel_size) for c in rc]  # Smallest index to contain relevant potential x,y,z
        ind_max = [int((c + r) // voxel_size) for c in rc]  # Largest index to contain relevant potential x,y,z
        # Explicit real space coordinates for the max and min boundary of each voxel
        x_min_bound = xp.arange(ind_min[0], ind_max[0] + 1, 1) * voxel_size - rc[0]
        x_max_bound = xp.arange(ind_min[0] + 1, ind_max[0] + 2, 1) * voxel_size - rc[0]
        y_min_bound = xp.arange(ind_min[1], ind_max[1] + 1, 1) * voxel_size - rc[1]
        y_max_bound = xp.arange(ind_min[1] + 1, ind_max[1] + 2, 1) * voxel_size - rc[1]
        z_min_bound = xp.arange(ind_min[2], ind_max[2] + 1, 1) * voxel_size - rc[2]
        z_max_bound = xp.arange(ind_min[2] + 1, ind_max[2] + 2, 1) * voxel_size - rc[2]

        atom_potential = 0


        for j in range(5):
            sqrt_b = xp.sqrt(b[j])  # calculate only once
            # Difference of error function == integrate over Gaussian
            int_x = sqrt_b / (4 * xp.sqrt(xp.pi)) * (erf(x_max_bound * 2 * xp.pi / sqrt_b) -
                                                     erf(x_min_bound * 2 * xp.pi / sqrt_b))
            x_matrix = xp.tile(int_x[:, xp.newaxis, xp.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = sqrt_b / (4 * xp.sqrt(xp.pi)) * (erf(y_max_bound * 2 * xp.pi / sqrt_b) -
                                                     erf(y_min_bound * 2 * xp.pi / sqrt_b))
            y_matrix = xp.tile(int_y[xp.newaxis, :, xp.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = sqrt_b / (4 * xp.sqrt(xp.pi)) * (erf(z_max_bound * 2 * xp.pi / sqrt_b) -
                                                     erf(z_min_bound * 2 * xp.pi / sqrt_b))
            z_matrix = xp.tile(int_z[xp.newaxis, xp.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])
            temp = a[j] / b[j] ** (3 / 2) * x_matrix * y_matrix * z_matrix
            atom_potential += temp

        potential[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += atom_potential

        if solvent_exclusion:
            # excluded solvent potential
            int_x = xp.sqrt(xp.pi) * r_0 / 2 * (erf(x_max_bound / r_0) - erf(x_min_bound / r_0))
            x_matrix = xp.tile(int_x[:, xp.newaxis, xp.newaxis], [1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_y = xp.sqrt(xp.pi) * r_0 / 2 * (erf(y_max_bound / r_0) - erf(y_min_bound / r_0))
            y_matrix = xp.tile(int_y[xp.newaxis, :, xp.newaxis], [ind_max[0] - ind_min[0] + 1,
                                                                  1,
                                                                  ind_max[2] - ind_min[2] + 1])
            int_z = xp.sqrt(xp.pi) * r_0 / 2 * (erf(z_max_bound / r_0) - erf(z_min_bound / r_0))
            z_matrix = xp.tile(int_z[xp.newaxis, xp.newaxis, :], [ind_max[0] - ind_min[0] + 1,
                                                                  ind_max[1] - ind_min[1] + 1,
                                                                  1])

            solvent[ind_min[0]:ind_max[0] + 1, ind_min[1]:ind_max[1] + 1, ind_min[2]:ind_max[2] + 1] += ( x_matrix *
                                                                                                          y_matrix *
                                                                                                          z_matrix)

    # Voxel volume
    dV = voxel_size ** 3
    # Convert to correct units
    C = 4 * xp.sqrt(xp.pi) * phys.constants['h'] ** 2 / (phys.constants['el'] * phys.constants['me']) * 1E20  # angstrom**2

    if solvent_exclusion:
        # Correct for solvent and convert both the solvent and potential array to the correct units.
        real = (potential / dV * C) - (solvent / dV * (V_sol * solvent_factor))
    elif solvent_masking:
        # construct solvent mask
        # solvent_mask = xp.zeros_like(potential)
        cutoff = potential[potential > 0].mean()
        solvent_mask = (potential > 0.7 * cutoff) * 1.6
        # print(solvent_mask.shape)
        # gaussian decay of mask
        smoothed_mask = reduce_resolution(solvent_mask, voxel_size, voxel_size * 4)
        smoothed_mask[smoothed_mask < 0.001] = 0
        # fig, (ax1, ax2) = plt.subplots(1, 2)
        # slice = int(solvent_mask.shape[2] // 2)
        # ax1.imshow(solvent_mask[:, :, slice ])
        # ax2.imshow(smoothed_mask[:, :, slice ])
        # show()
        real = (potential / dV * C) - (smoothed_mask * (V_sol * solvent_factor))
    else:
        real = potential / dV * C

    if absorption_contrast:
        from simulateProjections import potential_amplitude
        # voltage by default 300 keV
        molecule_absorption = potential_amplitude(density, molecular_weight, voltage)
        solvent_absorption = potential_amplitude(0.93, 18, voltage) * solvent_factor
        # gold_amplitude = potential_amplitude(19.3, 197, voltage)
        print(f'molecule absorption = {molecule_absorption:.3f}')
        print(f'solvent absorption = {solvent_absorption:.3f}')

        if solvent_masking:
            imaginary = smoothed_mask * (molecule_absorption - solvent_absorption)
        else:
            imaginary = (real > 0) * (molecule_absorption - solvent_absorption)

        return [real, imaginary]
    else:
        return [real]


def iasa_rough(filename, voxel_size=10, solvent_exclusion=False, V_sol=4.5301):
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).
    This function applies averaging of the potential over the voxels BUT in a coarse way. Use only for desired foxel
    spacing of 10 A.

    @param filename: ID of pdb file as present in pdb folder
    @type filename: string
    @param voxel_size: Size (A) of voxel in output map, default 1 A
    @type voxel_size: float
    @param solvent_exclusion: model the potential taking into account solvent exclusion
    @type solvent_exclusion: Bool
    @param V_sol: average solvent background potential (V/A^3)
    @type V_sol: float

    @return: A volume with interaction potentials
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    oversampling = 2 # Adjust if neccesary. Voxel size of this function should never be below 5A.
    extra_space = 10  # extend volume by 10 A in all directions

    print(f' - Calculating IASA potential from {filename}')

    x_coordinates, y_coordinates, z_coordinates, elements, _, _ = read_structure(filename)

    x_max = xp.max(x_coordinates - xp.min(x_coordinates))
    y_max = xp.max(y_coordinates - xp.min(y_coordinates))
    z_max = xp.max(z_coordinates - xp.min(z_coordinates))
    dimensions = [x_max, y_max, z_max]
    largest_dimension = max(dimensions)
    difference = [largest_dimension - a for a in dimensions]

    x_coordinates = x_coordinates - xp.min(x_coordinates) + extra_space + difference[0]/2
    y_coordinates = y_coordinates - xp.min(y_coordinates) + extra_space + difference[1]/2
    z_coordinates = z_coordinates - xp.min(z_coordinates) + extra_space + difference[2]/2
    # Define the volume of the protein
    szx = xp.abs(xp.max(x_coordinates) - xp.min(x_coordinates)) + 2 * extra_space + difference[0]
    szy = xp.abs(xp.max(y_coordinates) - xp.min(y_coordinates)) + 2 * extra_space + difference[1]
    szz = xp.abs(xp.max(z_coordinates) - xp.min(z_coordinates)) + 2 * extra_space + difference[2]
    sz = xp.round(xp.array([szx, szy, szz]) / (voxel_size/oversampling)).astype(int)

    potential = xp.zeros(sz)

    C = 2 * xp.pi * phys.constants['h_bar'] ** 2 / (phys.constants['el'] * phys.constants['me']) * 1E20  # angstrom**2
    dV = (voxel_size/oversampling) ** 3

    for i in range(len(elements)):
        if xp.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i].upper()

        rc = [x_coordinates[i], y_coordinates[i], z_coordinates[i]]  # atom center
        ind = tuple([int((c) // (voxel_size/oversampling)) for c in rc])  # Indexes to contain the potential

        potential[ind] += (xp.array(phys.scattering_factors[atom]['g'])[0:5].sum() * C / dV)

        if solvent_exclusion:
            if atom in list(phys.volume_displaced):
                potential[ind] -= phys.volume_displaced[atom] * V_sol / dV
            else:  # If not H,C,O,N we assume the same volume displacement as for carbon
                potential[ind] -= phys.volume_displaced['C'] * V_sol / dV
    # Apply Gaussian filter
    potential = reduce_resolution(potential, (voxel_size/oversampling), voxel_size)
    # Bin the volume
    potential = bin(potential, oversampling)

    return potential


def iasa_potential(filename, voxel_size=1., oversampling=0, low_pass_filter=True): # add params voxel_size, oversampling?
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).
    This function using sampling to determine the potential.

    @param filename: ID of pdb file as present in pdb folder
    @type filename: string
    @param voxel_size: Size (A) of voxel in output map, default 1 A
    @type voxel_size: float
    @param oversampling: Increased sampling of potential (multiple of 1), default 1 i.e. no oversampling
    @type oversampling: int

    @return: A volume with interaction potentials
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    extra_pixels = 40  # extend volume by 10 A

    print(f' - Calculating IASA potential from {filename}')

    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filename)

    x_coordinates = x_coordinates - xp.min(x_coordinates) + extra_pixels
    y_coordinates = y_coordinates - xp.min(y_coordinates) + extra_pixels
    z_coordinates = z_coordinates - xp.min(z_coordinates) + extra_pixels
    # Define the volume of the protein
    szx = xp.abs(xp.max(x_coordinates) - xp.min(x_coordinates)) + 2 * extra_pixels
    szy = xp.abs(xp.max(y_coordinates) - xp.min(y_coordinates)) + 2 * extra_pixels
    szz = xp.abs(xp.max(z_coordinates) - xp.min(z_coordinates)) + 2 * extra_pixels
    shape = [szx,szy,szz]
    print(f'pre adjustment shape = {shape}')
    # Make the size of the box a result of 2^n to make fourier transforms more efficient.
    exponent = math.log(max(shape),2)
    if (exponent % 1) != 0:
        new_size = 2**xp.round(exponent)
        difference = [new_size - s for s in shape]
        x_coordinates += (difference[0] / 2)
        y_coordinates += (difference[1] / 2)
        z_coordinates += (difference[2] / 2)
        shape = [s+d for s,d in zip(shape,difference)]
    print(f'after adjustment = {shape}')

    if (oversampling > 1) and ( type(oversampling) is int ):
        spacing = voxel_size / oversampling
    else:
        spacing = voxel_size
    print(f'IASA: volume will be sampled at {spacing}A')

    sz = [ s / spacing for s in shape] # increase for oversampling
    # C = 2132.8 A^2 * V; 1E20 is a conversion factor for Angstrom^2,
    # h and not h_bar, which is why we multiply with 4 * sqrt(pi) instead of 16*pi^(5/2)
    C = 4 * xp.sqrt(xp.pi) * phys.constants['h']**2 / (phys.constants['el'] * phys.constants['me']) * 1E20

    potential = xp.zeros(tuple(xp.round(sz).astype(int)))
    print(f'#atoms to go over is {len(x_coordinates)}.')

    for i in range(len(elements)):
        if xp.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i].upper()
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = xp.array(phys.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b is used for calculating the radius of the potential. See Rullgard et al. (2011) for addition of 16 R^2.
        # Incorporates low pass filtering by extending the gaussian radius. NOTE: b_factor is zero anyway
        b += (b_factor) # + 16 * spacing**2)

        r2 = 0
        b1 = xp.zeros(5)
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            b1[j] = 4 * xp.pi ** 2 / b[j] * spacing ** 2
            r2 = xp.maximum(r2, 10/b1[j])
        # Radius of gaussian sphere
        r = xp.sqrt(r2 / 3)
        xc1, yc1, zc1 = x_coordinates[i] / spacing, y_coordinates[i] / spacing, z_coordinates[i] / spacing
        # Center of the gaussian sphere.
        rc = [xc1, yc1, zc1]
        # Calculate the absolute indexes for the potential matrix.
        kmin = [xp.maximum(0,x).astype(int) for x in xp.ceil(rc-r)]
        kmax = [xp.minimum(xp.floor(x)-1,xp.floor(y+r)).astype(int) for x,y in zip(sz,rc)]
        kmm = max([x-y for x,y in zip(kmax,kmin)])
        # Determine the coordinates for sampling from the sum of Gaussians.
        x = xc1 - xp.arange(kmin[0], kmin[0]+kmm+1, 1)
        y = yc1 - xp.arange(kmin[1], kmin[1]+kmm+1, 1)
        z = zc1 - xp.arange(kmin[2], kmin[2]+kmm+1, 1)

        atom_potential = 0
        for j in range(5):
            x_matrix = xp.tile(xp.exp(-b1[j] * x**2)[:,xp.newaxis,xp.newaxis], [1, kmm+1, kmm+1])
            y_matrix = xp.tile(xp.exp(-b1[j] * y**2)[xp.newaxis,:,xp.newaxis], [kmm+1, 1, kmm+1])
            z_matrix = xp.tile(xp.exp(-b1[j] * z**2)[xp.newaxis,xp.newaxis,:], [kmm+1, kmm+1, 1])
            tmp = a[j] / b[j]**(3/2) * C * x_matrix * y_matrix * z_matrix
            atom_potential += tmp
        atom_potential *= occupancy
        # add the potential of element i to the full potential map
        potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] = \
            potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] + atom_potential

    # potential = subtract_solvent(potential)
    print(f'oversampled potential shape {potential.shape}')
    # Only need to downscale if the volume was oversampled.
    if spacing != voxel_size:
        # Volume needs to be cubic before applying a low pass filter and binning (??)
        # size = max(potential.shape)
        # increment = [size - a for a in potential.shape]
        # # using spline interpolation here does not decrease our accuracy as the volume is already oversampled!
        # potential = extend_volume(potential, increment, pad_value=0, symmetrically=True, true_center=True)
        print('low pass filter')
        if low_pass_filter:
            potential = reduce_resolution(potential, spacing, voxel_size)
        print('binning')
        potential = bin(potential, oversampling)
        # Other option is to use scale instead of binning, which is proper downsampling.
        # potential = scale(potential, 1/oversampling)

    return potential


def parse_apbs_output(filename):
    """
    parse_apbs_output: Parses output file from APBS.
    @param filename: Path and name of .pqr.dx file produced by APBS software
    @return: data, dxnew, dynew, dznew are in A, thickness in m, data is the reshaped apbs output
    @author: Marten Chaillet
    """
    # Read file header
    try:
        with open(filename, 'r') as apbs_file:
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
    data = xp.array(data,dtype=float)
    # Form to 3D grid
    data = xp.reshape(data,(grid[2],grid[1],grid[0]),order='F').transpose(2,1,0) # x, y, z needs to be done as such to correspond to iasa
    dx = spacing[0]
    dy = spacing[1]
    dz = spacing[2]
    return data, dx, dy, dz


def resample_apbs(filename, voxel_size=1.0, low_pass_filter=False):
    """
    resample_APBS: First calls parse_abps_output to read an apbs output file, then scales voxels to 1 A voxels and
    refactors the values to volts.

    @param filename: Full filenamepath
    @type filename: string

    @return: Resampled and scaled V_bond potential
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    # from pytom_volume import vol
    # from pytom_numpy import vol2npy
    from voltools import transform

    print(f' - Parsing and resampling APBS file {filename}')

    # Parse APBS data file
    potential, dxnew, dynew, dznew = parse_apbs_output(filename)
    # Check if the voxel size is allowed, and adjust if neccessary
    smallest_possible = max([dxnew, dynew, dznew])
    if not(voxel_size >= smallest_possible):
        print(f'Requested voxel size is smaller than voxel size of the apbs map. Adjust to smallest possible voxel size '
              f'of {smallest_possible}.')
        voxel_size = smallest_possible
    # Make voxels same size by scaling to dx
    factor = (dxnew/dxnew, dynew/dxnew, dznew/dxnew)
    potential = scale(potential, factor)
    if low_pass_filter:
        size = max(potential.shape)
        increment = [size-a for a in potential.shape]
        potential = extend_volume(potential, increment, pad_value=0, symmetrically=True, true_center=True)
        potential = reduce_resolution(potential, dxnew, voxel_size)
    # Scale the volume to voxels with voxel_size
    potential = scale(potential, dxnew/voxel_size)

    # TODO use voltools

    print(f'Data after reshaping to {voxel_size} A voxels: {potential.shape}')

    # Convert to from kT/e to volts
    temperature = 291 # [K] = 18 C (room temperature)
    convert_to_volts = phys.constants['kb'] * temperature / phys.constants['el']
    potential = convert_to_volts * potential
    return potential


def combine_potential(iasa_potential, bond_potential, voxel_size):
    """
    Combine isolated atom and bond potential.

    @param iasa_potential:
    @type iasa_potential: 3d numpy/cupy array[x,y,z], float
    @param bond_potential:
    @type bond_potential: 3d numpy/cupy array[x,y,z], float

    @return: Combined IASA and bond potential (APBS)
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    print(' - Combining iasa and bond potential')
    # Reduce interpolation! Only shift one of the volumes to true center. Preferably V_bond.
    # This code could be more elegant.
    assert max(iasa_potential.shape) >= max(bond_potential.shape), print('Break! Attempting to extend IASA volume.')
    try:
        # Because volumes are shifted to true center (by interpolation) the precision is reduced slightly.
        difference = [max(iasa_potential.shape) - x for x in bond_potential.shape]
        bond_potential = extend_volume(bond_potential, difference, symmetrically=True, true_center=True)
        # Apply filter to remove the artifacts from interpolating
        bond_potential = reduce_resolution(bond_potential, voxel_size, voxel_size)
        full_potential = iasa_potential + bond_potential
    except Exception as e:
        print(e)
        raise Exception('Could not fit atom and bond potential together.')
    return iasa_potential, bond_potential, full_potential


def wrapper(file, input_folder, output_folder, voxel_size, binning=None, exclude_solvent=False, solvent_masking=False,
            solvent_potential=4.5301, absorption_contrast=False, voltage=300E3, solvent_factor=1.0):
    import pytom.tompy.io
    # TODO function can be executed in parallel for multiple structures

    # Id does not makes sense to apply absorption contrast if solvent exclusion is not turned on
    if absorption_contrast:
        assert exclude_solvent, print('absorption contrast can only be applied if solvent exclusion is used.')

    pdb_id = file.split('.')[0]

    # Call external programs for structure preparation and PB-solver
    structure = call_chimera(file, input_folder, output_folder) # output structure name is dependent on modification by chimera
    # call_apbs(folder, structure, ph=ph)

    assert structure != 0, 'something went wrong with chimera'

    # Calculate atom and bond potential, and store them
    # 4 times oversampling of IASA yields accurate potentials
    # Could be {structure}.{extension}, but currently chimera is only able to produce .pdb files, so the extended
    # structure file created by call chimera has a .pdb extension.
    v_atom = iasa_integration(f'{output_folder}/{structure}.pdb', voxel_size=voxel_size,
                              solvent_exclusion=exclude_solvent, solvent_masking=solvent_masking, V_sol=solvent_potential,
                              absorption_contrast= absorption_contrast, voltage=voltage, solvent_factor=solvent_factor)
    # Absorption contrast map generated here will look blocky when generated at 2.5A and above!
    if absorption_contrast:
        print(' - Filtering volume')
        print(v_atom[0].shape, v_atom[1].shape)
        filtered = [reduce_resolution_2(v_atom[0], voxel_size, 2*voxel_size),
                    reduce_resolution_2(v_atom[1], voxel_size, 2*voxel_size)]
        output_name = f'{pdb_id}_{voxel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
        pytom.tompy.io.write(f'{output_folder}/{output_name}_real.mrc', filtered[0])
        pytom.tompy.io.write(f'{output_folder}/{output_name}_imag.mrc', filtered[1])
    else:
        print(' - Filtering volume')
        filtered = reduce_resolution_2(v_atom[0], voxel_size, 2*voxel_size)
        if exclude_solvent or solvent_masking:
            output_name = f'{pdb_id}_{voxel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
        else:
            output_name = f'{pdb_id}_{voxel_size:.2f}A'
        pytom.tompy.io.write(f'{output_folder}/{output_name}.mrc', filtered)

    if binning is not None:
        # v_atom_binned = iasa_integration(f'{output_folder}/{structure}.pdb', voxel_size=voxel_size*binning,
        #                       solvent_exclusion=exclude_solvent, V_sol=solvent_potential)
        # first filter the volume!
        if absorption_contrast:
            print(' - Filtering volume')
            filtered = [reduce_resolution_2(v_atom[0], voxel_size, voxel_size * binning),
                        reduce_resolution_2(v_atom[1], voxel_size, voxel_size * binning)]
            print(' - Binning volume')
            binned = [bin(filtered[0], binning), bin(filtered[1], binning)]

            output_name = f'{pdb_id}_{voxel_size*binning:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            pytom.tompy.io.write(f'{output_folder}/{output_name}_real.mrc', binned[0])
            pytom.tompy.io.write(f'{output_folder}/{output_name}_imag.mrc', binned[1])

        else:
            print(' - Filtering volume')
            filtered = reduce_resolution_2(v_atom, voxel_size, voxel_size*binning)
            print(' - Binning volume')
            binned = bin(filtered, binning)
            if exclude_solvent or solvent_masking:
                output_name = f'{pdb_id}_{voxel_size*binning:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            else:
                output_name = f'{pdb_id}_{voxel_size*binning:.2f}A'
            pytom.tompy.io.write(f'{output_folder}/{output_name}.mrc', binned)
    return


if __name__ == '__main__':
    # Make into script with command line options, following options should suffice
    # pdb_id, mandatory
    # voxel_size, optional, default is 1 A?
    # resolution, optional, default is 2*voxel_size
    # pH, optional, default is 7

    # parameters: folder, pdb_id, ph, voxel_size, resolution?

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA, PDB2PQR (modified), APBS

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    # syntax is ScriptOption([short, long], description, requires argument, is optional)
    options = [ScriptOption(['-f', '--file'], 'Protein structure file, either pdb or cif.', True, False),
               ScriptOption(['-d', '--folder'], 'Folder where the structure file is located.', True, False),
               ScriptOption(['-o', '--outputfolder'], 'Folder to store the files produced by potential.py.', True, False),
               ScriptOption(['-s', '--spacing'], 'The size of the voxels of the output volume. 1A by default.', True, True),
               ScriptOption(['-b', '--binning'], 'Number of times to bin. Additional storage of binned volume.', True, True),
               ScriptOption(['-x', '--exclude_solvent'],
                            'Whether to exclude solvent around each atom as a correction of the potential.', False, True),
               ScriptOption(['-m', '--mask_solvent'], 'Whether to exclude solvent by masking.', False, True),
               ScriptOption(['-p', '--solvent_potential'],
                            'Value for the solvent potential. By default amorphous ice, 4.5301 V.', True, True),
               ScriptOption(['-c', '--percentile'], 'Multiplication for solvent potential and absorption contrast of '
                                                    'solvent to decrease/increase contrast. Value between 0 and 3 (could'
                                                    ' be higher but would no make sense).', True, True),
               ScriptOption(['-a', '--absorption_contrast'],
                            'Whether to add imaginary part of absorption contrast, can only be done if solvent'
                            'is excluded.', False, True),
               ScriptOption(['-v', '--voltage'],
                            'Value for the electron acceleration voltage. Need for calculating the inelastic mean free '
                            'path in case of absorption contrast calculation. By default 300 (keV).', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], description='Calculate electrostatic potential from protein structure file.\n'
                                                                  'Script has dependencies on pytom, chimera. (pdb2pqr apbs)',
                          authors='Marten Chaillet', options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    elif ('-h' in str(sys.argv)) or ('--help' in str(sys.argv)):
        # This needs to be added to correctly parse the help options. ScriptHelper does not do it correctly.
        # Better make a change to ScriptHelper.
        print(helper)
        sys.exit()
    try:
        file, input_folder, output_folder, voxel_size, binning, exclude_solvent, solvent_masking, solvent_potential, \
                    solvent_factor, absorption_contrast, voltage, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if not voxel_size: voxel_size=1
    else: voxel_size = float(voxel_size)

    if not solvent_potential: solvent_potential = 4.5301
    else: solvent_potential = float(solvent_potential)

    if not solvent_factor: solvent_factor = 1
    else: solvent_factor = float(solvent_factor)

    if not voltage: voltage = 300E3
    else: voltage = float(voltage)*1E3 # value should be given in keV

    if binning: binning = int(binning)

    if not os.path.exists(f'{input_folder}/{file}'):
        print('Protein structure file does not exist!')
        sys.exit()
    if not os.path.exists(output_folder):
        print('Output folder does not exist!')
        sys.exit()

    wrapper(file, input_folder, output_folder, voxel_size, binning=binning,
            exclude_solvent=exclude_solvent, solvent_masking=solvent_masking, solvent_potential=solvent_potential,
            absorption_contrast=absorption_contrast, voltage=voltage, solvent_factor=solvent_factor)


