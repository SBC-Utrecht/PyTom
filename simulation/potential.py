import numpy as xp
import constant_dictionaries as phys
import scipy.ndimage
import os
import pytom.basic.functions

V_WATER = 4.5301 # potential value of amorphous ice (from Vulovic et al., 2013)

def call_chimera(pdb_folder, pdb_id):
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
    print(f' - Calling chimera for {pdb_id}')

    # Command 'sym' in chimera crashes when there is no BIOMT symmetry in the pdb file. We need to make sure sym is only
    # executed when the BIOMT information specifies symmetrical units.
    symmetry = []
    try:
        with open(f'{pdb_folder}/{pdb_id}.pdb','r') as pdb:
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

    extension = 'py'
    if len(set(symmetry)) > 1:
        scriptname = f'_sym_addh_{pdb_id}'
        try:
            with open(f'{pdb_folder}/{scriptname}.{extension}', 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent and ions (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {pdb_folder}/{pdb_id}.pdb")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     f'rc("addh")\n'
                                     f'rc("sym group biomt")\n'             # group biomt is also the default
                                     f'rc("combine all modelId 10")\n'
                                     f'rc("write format pdb #10 {pdb_folder}/{pdb_id}_sym_addh.pdb")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
    else:
        scriptname = f'_addh_{pdb_id}'
        try:
            with open(f'{pdb_folder}/{scriptname}.{extension}', 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent and ions (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {pdb_folder}/{pdb_id}.pdb")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     f'rc("addh")\n'
                                     f'rc("write format pdb #0 {pdb_folder}/{pdb_id}_addh.pdb")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
    # module chimera should be loaded here...
    try:
        os.system(f'chimera --nogui --script {pdb_folder}/{scriptname}.{extension}')
    except Exception as e:
        print(e)
        raise Exception('Chimera is likely not on your current path.')

    if len(set(symmetry)) > 1:
        return f'{pdb_id}_sym_addh' # returns new pdb name
    else:
        return f'{pdb_id}_addh'


def modify_structure_file(filename, pattern, replacement, line_start=''):
    """
    Function required to make pqr files with large negative coordinates readible for APBS. This program can only parse
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


def call_apbs(pdb_folder, structure, apbs_folder, force_field='amber', ph=7.):
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
    print(f' - Running pdb2pqr and APBS on {pdb_folder}/{structure}.pdb')
    cwd = os.getcwd()
    input_file = f'{pdb_folder}/{structure}.pdb'
    output_folder = f'{apbs_folder}/{structure.split("_")[0]}'
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    output_file = f'{output_folder}/{structure}.pqr'
    apbs_in = f'{output_folder}/{structure}.in'
    try:
        # Also PDB2PKA ph calculation method. Requires PARSE force field, can take very long for large proteins.
        os.system(f'pdb2pqr.py --ff={force_field} --ph-calc-method=propka --with-ph={ph} --apbs-input {input_file} {output_file}')
        print('Add white space delimiters to pqr file.')
        modify_structure_file(output_file, '-', ' -', line_start='ATOM')
        # APBS needs to execute from the folder where the structure is present
        os.chdir(output_folder)
        os.system(f'module load apbs_mc/1.5; apbs {apbs_in}')
        # ? subprocess.run(['apbs', f'{structure}.in'], cwd=output_folder) # executes in cwd
    except Exception as e:
        print(e)
        raise Exception('pdb2pqr or APBS potentially not on path.')

    # Change back to original directory
    os.chdir(cwd)
    return


def iasa_potential(filename, voxel_size=1., oversampling=1): # add params voxel_size, oversampling?
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).

    @param pdb_folder: Path to folder with pdb files
    @type pdb_folder:
    @param pdb_id: ID of pdb file as present in pdb folder
    @type pdb_id:
    @param voxel_size: Size (A) of voxel in output map, default 1 A
    @type voxel_size:
    @param oversampling: Increased sampling of potential (multiple of 1), default 1 i.e. no oversampling
    @type oversampling:

    @return: A volume with interaction potentials
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """

    # oversampling = int(oversampling)   # oversampling to 0.25 A ==> Not used at the moment
    # if oversampling:
    #   voxel_size = voxel_size / oversampling
    extra_pixels = 10  # extend volume by 10 A

    print(f' - Calculating IASA potential from {filename}')

    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = [],[],[],[],[],[]

    try:
        with open(filename, 'r') as pdb:
            lines = pdb.readlines()
            atoms = [line for line in lines if line[:4]=='ATOM']
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
    except Exception as e:
        print(e)
        raise Exception('Could not read pdb file.')

    x_coordinates = x_coordinates - xp.min(x_coordinates) + extra_pixels
    y_coordinates = y_coordinates - xp.min(y_coordinates) + extra_pixels
    z_coordinates = z_coordinates - xp.min(z_coordinates) + extra_pixels
    # Define the volume of the protein
    szx = xp.abs(xp.max(x_coordinates) - xp.min(x_coordinates)) + 2 * extra_pixels
    szy = xp.abs(xp.max(y_coordinates) - xp.min(y_coordinates)) + 2 * extra_pixels
    szz = xp.abs(xp.max(z_coordinates) - xp.min(z_coordinates)) + 2 * extra_pixels

    sz = [szx,szy,szz] # [ x / voxel_size for x in [szx,szy,szz]] could increase size for oversampling here
    # C = 2132.8 A^2 * V; 1E20 is a conversion factor for Angstrom^2
    C = 4 * xp.sqrt(xp.pi) * phys.constants['h']**2 / (phys.constants['el'] *
                                                                phys.constants['me']) * 1E20

    potential = xp.zeros(tuple(xp.round(sz).astype(int)))
    print(f'#atoms to go over is {len(x_coordinates)}.')
    for i in range(len(elements)):
        if xp.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i]
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = xp.array(phys.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b is used for calculating the radius of the potential. See Rullgard et al. (2011) for addition of 16 R^2
        b += b_factor + 16 * voxel_size**2
        r2 = 0
        b1 = xp.zeros(5)
        for j in range(5):
            # calculate maximal radius assuming symmetrical potential
            b1[j] = 4 * xp.pi**2 / b[j] * voxel_size**2
            r2 = xp.maximum(r2, 10/b1[j])
        r = xp.sqrt(r2 / 3)
        xc1 = x_coordinates[i] # / voxel_size
        yc1 = y_coordinates[i] # / voxel_size
        zc1 = z_coordinates[i] # / voxel_size
        rc = [xc1, yc1, zc1]
        kmin = [xp.maximum(0,x).astype(int) for x in xp.ceil(rc-r)]
        kmax = [xp.minimum(xp.floor(x)-1,xp.floor(y+r)).astype(int) for x,y in zip(sz,rc)]
        kmm = max([x-y for x,y in zip(kmax,kmin)])

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

    potential -= V_WATER
    potential[potential<0] = 0

    # if oversampling:
    #     downsample with voltools transform

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


def resample_apbs(filename):
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
    # from voltools import transform

    print(f' - Parsing and resampling APBS file {filename}')

    # Parse APBS data file
    potential, dxnew, dynew, dznew = parse_apbs_output(filename)
    # Make the voxels cubic 1 A
    voxel_size = 1
    scaling = [dxnew/voxel_size, dynew/voxel_size, dznew/voxel_size] # dxnew, etc. are in A
    # TODO We want to use voltools here instead of scipy
    potential = scipy.ndimage.zoom(potential, tuple(scaling), order=5)

    # # Create a edge taper mask
    # taper_width = 5 # half of extra pixels added?
    # [x,y,z] = potential.shape
    # # taper_edges only takes pytom volume
    # _, taper_mask = pytom.basic.functions.taper_edges(vol(x,y,z),taper_width)
    # taper_mask = vol2npy(taper_mask)
    # # Apply the mask in Fourier space
    # print('Check')
    # potential = xp.real(xp.fft.fftshift(xp.fft.ifftn(xp.fft.fftn(xp.fft.ifftshift(potential)) * xp.fft.ifftshift(taper_mask))))

    print(f'Data after reshaping to 1A voxels: {potential.shape}')

    # Convert to from kT/e to volts
    temperature = 291 # [K] = 18 C (room temperature)
    convert_to_volts = phys.constants['kb'] * temperature / phys.constants['el']
    potential = convert_to_volts * potential
    return potential


def combine_potential(iasa_potential, bond_potential):
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
    try:
        iasa_size = iasa_potential.shape
        bond_size = bond_potential.shape
        if all(x >= y for x,y in zip(iasa_size,bond_size)):
            difference = [x-y for x,y in zip(iasa_size,bond_size)]
            bond_potential = xp.pad(bond_potential, tuple([(0,x) for x in difference]), 'constant', constant_values=0)
            bond_potential = scipy.ndimage.interpolation.shift(bond_potential,tuple([x/2 for x in difference]), order=5)
            print(f'Extended bond volume from {bond_size} to {bond_potential.shape}, and shifted values to new center.')
        elif all(x < y for x,y in zip(iasa_size,bond_size)):
            difference = [y-x for x, y in zip(iasa_size, bond_size)]
            iasa_potential = xp.pad(iasa_potential, tuple([(0,x) for x in difference]), 'constant', constant_values=0)
            iasa_potential = scipy.ndimage.interpolation.shift(iasa_potential, tuple([x/2 for x in difference]), order=5)
            print(f'Extended IASA volume from {iasa_size} to {iasa_potential.shape}, and shifted values to new center.')
        full_potential = iasa_potential + bond_potential
    except Exception as e:
        print(e)
        raise Exception('Could not fit atom and bond potential together.')
    return full_potential


def crop(data, factor):
    s = int((data.shape[0] % factor) // 2)
    d = int((data.shape[0] % factor) % 2)
    data = data[s:data.shape[0] - s - d, s:data.shape[0] - s - d, s:data.shape[0] - s - d]

    ds = int(data.shape[0] // factor)
    image_size = data.shape[0]
    binned = data.reshape(ds, image_size // ds, ds, image_size // ds, ds, image_size // ds).mean(-1).mean(1).mean(-2)

    # Apply edge tapering
    # ft = fftshift(fftn(data))
    # x, y, z = numpy.array(ft.shape) // 2
    # ff = factor
    # ft = ft[int(x - x // ff):int(x + x // ff), int(y - y // ff):int(y + y // ff), int(z - z // ff):int(z + z // ff)]
    # particle = abs(ifftn(fftshift(ft)))

    return binned


def scale(potential, resolution_in, resolution_out, order=5, taper_edges=0):

    size = max(potential.shape)
    m2 = xp.zeros((size, size, size))
    dx, dy, dz = potential.shape
    sx, ex = (size - dx) // 2, size - int(xp.ceil((size - dx) / 2.))
    sy, ey = (size - dy) // 2, size - int(xp.ceil((size - dy) / 2.))
    sz, ez = (size - dz) // 2, size - int(xp.ceil((size - dz) / 2.))
    m2[sx:ex, sy:ey, sz:ez] = potential
    factor = resolution_out/resolution_in

    return crop(m2, factor)


# def scale(potential, resolution_in, resolution_out, order=5, taper_edges=0):
#     """
#     Scale potential from resolution_in to resolution_out. Possibly with edge_tapering.
#     @param potential: 3d array (numpy)
#     @param resolution_in:
#     @param resolution_out:
#     @return:
#     """
#     from voltools import transform
#     scaling = resolution_in/resolution_out
#     potential = scipy.ndimage.zoom(potential, scaling, order=order)
#     if taper_edges:
#         print("Edge tapering has to be implemented.")
#         # from simulateProjections import spheremask
#         # sphere_mask = spheremask(xp.ones_like(potential), maxradius-10, smooth=10, ellipsoid=1)
#     return potential

def wrapper(pdb_id, pdb_folder, apbs_folder, iasa_folder, bond_folder, map_folder):
    import pytom.tompy.io
    # TODO function can be executed in parallel for multiple structures
    # pdb_folder = '/data2/mchaillet/structures/pdb'
    # apbs_folder = '/data2/mchaillet/structures/apbs'
    # iasa_folder = '/data2/mchaillet/structures/potential_iasa'
    # bond_folder = '/data2/mchaillet/structures/potential_bond'
    # map_folder = '/data2/mchaillet/structures/potential_map'

    # Call external programs for structure preparation and PB-solver
    structure = call_chimera(pdb_folder, pdb_id) # output structure name is dependent on modification by chimera
    call_apbs(pdb_folder, structure, apbs_folder)

    outfile = f'{structure}_1.0A.mrc'
    # Calculate atom and bond potential, and store them
    v_atom = iasa_potential(f'{pdb_folder}/{structure}.pdb')
    pytom.tompy.io.write(f'{iasa_folder}/{outfile}', v_atom)
    v_bond = resample_apbs(f'{apbs_folder}/{structure.split("_")[0]}/{structure}.pqr.dx')
    pytom.tompy.io.write(f'{bond_folder}/{outfile}', v_bond)
    map = combine_potential(v_atom, v_bond)
    pytom.tompy.io.write(f'{map_folder}/{outfile}', map)

    map = scale(map, 1, 10)
    map[map<0] = 0
    outfile = f'{structure.split("_")[0]}_10.0A.mrc'
    pytom.tompy.io.write(f'{map_folder}/{outfile}', map)
    return


if __name__ == '__main__':
    # DEFAULT FOLDERS FOR WRITING INPUT/OUTPUT
    pdb_folder = '/data2/mchaillet/structures/pdb'
    apbs_folder = '/data2/mchaillet/structures/apbs'
    iasa_folder = '/data2/mchaillet/structures/potential_iasa'
    bond_folder = '/data2/mchaillet/structures/potential_bond'
    map_folder = '/data2/mchaillet/structures/potential_map'

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA, PDB2PQR (modified), APBS

    # LIST OF PDBS TO EXECUTE ON
    pdb_ids = ['3cf3', '1s3x', '1u6g', '4cr2', '1qvr', '3h84', '2cg9', '3qm1', '3gl1', '3d2f', '4d8q', '1bxn']

    for id in [id.upper() for id in pdb_ids]:
        if not os.path.exists(f'{pdb_folder}/{id}.pdb'):
            print(f'Skipping {id} because the pdb file does not exist in {pdb_folder}.')
            continue
        elif os.path.exists(f'{map_folder}/{id}_addh_1.0A.mrc') or os.path.exists(f'{map_folder}/{id}_sym_addh_1.0A.mrc'):
            print(f'{id} already has a map in folder {map_folder}.')
            continue
        else:
            try:
                wrapper(id, pdb_folder, apbs_folder, iasa_folder, bond_folder, map_folder)
            except Exception as e:
                print(e)
                print(f'Something when wrong while creating map for {id}. Continuing with next pdb file in list.')
                continue
