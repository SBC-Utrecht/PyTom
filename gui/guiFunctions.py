import matplotlib
import matplotlib.backends.backend_qt5agg
try: matplotlib.use('Qt5Agg')
except: pass
import os
import numpy
import copy
import time
from numpy.fft import fft2, ifft2, fftshift
from numpy import int32, arange, conj, zeros, ones, bool8, sqrt, newaxis
from pytom_volume import read
from pytom.tools.files import checkFileExists, checkDirExists
from pytom.basic.files import read as read_pytom
from pytom.basic.files import read_em_header
from pytom.gui.mrcOperations import read_mrc, read_angle
from pytom_numpy import vol2npy
from pytom.gui.guiSupportCommands import multiple_alignment
from pylab import imread
import mrcfile
import copy
from pytom.basic.files import write_em

def kill_proc(runner):
    runner.terminate()

def initSphere(x, y, z, radius=-1, smooth=0, cent_x=-1, cent_y=-1, cent_z=-1, filename='', cutoff_SD=3):
    from scipy.ndimage import gaussian_filter
    if cent_x < 0 or cent_x > x: cent_x = x // 2
    if cent_y < 0 or cent_y > y: cent_y = y // 2
    if cent_z < 0 or cent_z > z: cent_z = z // 2

    if x < 1 or y < 1 or z < 1 or radius < 1:
        return 0

    Y, X, Z = numpy.meshgrid(numpy.arange(y), numpy.arange(x), numpy.arange(z))
    X = X.astype('float32')
    Y = Y.astype('float32')
    Z = Z.astype('float32')

    X -= x - cent_x - 0.5
    Y -= y - cent_y - 0.5
    Z -= z - cent_z - 0.5

    R = numpy.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    R2 = R.copy()
    sphere = numpy.zeros((x, y, z), dtype='float32')
    sphere[R <= radius] = 1.

    if smooth:
        R2[R <= radius] = radius
        sphere = numpy.exp(-1 * ((R2 - radius) / smooth) ** 2)
        sphere[sphere <= numpy.exp(-cutoff_SD**2/2.)] = 0
    if filename:
        from numpy import array
        from pytom_numpy import npy2vol
        #maskVol = npy2vol(array(sphere, dtype='float32', order='F'), 3)
        #write_em(filename, maskVol)
        import mrcfile
        mrcfile.new(filename.replace('.em', '.mrc'), sphere,overwrite=True)
        return 1
    else:
        return sphere

def read_markerfile(filename,tiltangles):
    if filename[-4:] == '.mrc':
        mark_frames = mrc2markerfile(filename, tiltangles)
    elif filename[-3:] == '.em':
        mark_frames = em2markerfile(filename, tiltangles)
    elif filename[-5:] == '.wimp':
        mark_frames = wimp2markerfile(filename, tiltangles)
    elif filename[-4:] == '.npy':
        mark_frames = npy2markerfile(filename, tiltangles)
    elif filename[-4:] == '.txt':
        mark_frames = txt2markerfile(filename,tiltangles)
    else:
        return 0

    return mark_frames

def readMarkerfile(filename, num_tilt_images):
    if filename.endswith('.em'):
        from pytom.basic.files import read
        from pytom_numpy import vol2npy
        markerfile = read(filename)
        markerdata = vol2npy(markerfile).copy()
        return markerdata

    elif filename.endswith('.txt'):
        
        data = loadstar(filename)
        datalen = data.shape[0]
        x, y = datalen // num_tilt_images, num_tilt_images
        markerdata = data.reshape(x, y, 4)[:, :, 1:].transpose(2, 1, 0)
        return markerdata

def txt2markerfile(filename,tiltangles):
    data = loadstar(filename)
    datalen = data.shape[0]
    x, y = datalen // len(tiltangles), len(tiltangles)
    mark_frames = data.reshape(x, y, 4)[:, :, 2:].transpose(1, 0, 2)
    print(mark_frames.shape)
    return mark_frames

def npy2markerfile(filename,tiltangles):
    return numpy.load(filename)

def wimp2markerfile(filename, tiltangles, write=False):
    data = [line for line in open(filename).readlines()]

    for i in range(10):
        if '# of object' in data[i]:
            num_markers = int(data[i].split()[-1])

    markerset = numpy.ones((len(tiltangles),num_markers,2))*-1.

    for line in data:
        if 'Object #:' in line:
            objectID = int(line.split()[-1])-1
            if objectID >= num_markers:
                raise Exception('Wimp file contains more tiltimages than you have loaded! Conversion Failed.')
        try:
            x,y,itilt = map(float, line.split()[1:4])
            markerset[int(itilt)][objectID] = [x,y]
        except:
            pass

    return markerset

def mrc2markerfile(filename, tiltangles):
    if 0:
        mf = read_mrc(filename)
    else:
        m = mrcfile.open(filename, permissive=True)
        mf = copy.deepcopy(m.data)
        m.close()
    num,angles,d = mf.shape
    markers = -1 * numpy.ones((len(tiltangles), mf.shape[0], 2))
    for i in range(num):
        for j in range(angles):
            for n, angle in enumerate(tiltangles):
                if abs(mf[2, j, i] - angle) < 0.1:
                    j = n
                    break

            markers[j, i, 0] = mf[i, j, 1]
            markers[j, i, 1] = mf[i, j, 2]

    return markers

def em2markerfile(filename,tiltangles):
    vol = read(filename)
    mf = copy.deepcopy(vol2npy(vol))

    (d, angles, num) = mf.shape
    # print mf.shape
    locX, locY = 1, 2

    mark_frames = -1 * numpy.ones((len(tiltangles), num, 2))

    # print mf[0,:,1:3].shape,self.coordinates[:,0,:].shape
    markers = []
    for i in range(num):
        for j in range(angles):
            m = -1

            for n, angle in enumerate(tiltangles):
                if abs(mf[0, j, i] - angle) < 1:
                    m = n
                    break
            if m == -1: continue

            mark_frames[m, i, 0] = mf[locX, j, i]
            mark_frames[m, i, 1] = mf[locY, j, i]

    return mark_frames

def conv_mrc2em(directory, output_folder):
    '''This function converts all mrc files in the specified
    directory to the .em format and saves the files in the
    specified output_folder.'''

    fileList = sorted(os.listdir(directory))

    for fname in [f for f in fileList if f[-4:] == '.mrc']:

        mrc2em(os.path.join(directory, fname), output_folder)

def renumber_gui2pytom(output_folder, prefix):
    # store em files under different index (old_index +1), marked with .temp
    folder = output_folder
    fileList = [f for f in os.listdir(output_folder) if (f[:len(prefix)] == prefix and f[-3:] == '.em')]
    fileList = sorted(fileList)
    for nn, fname in enumerate(fileList):

        num = int(fname[-len('??.em'):-len('.em')]) + 1
        finalname = '{}/{}_{}.em'.format(folder, prefix, nn + 1)
        os.system('mv {} {}.temp'.format(os.path.join(folder, fname), finalname))

    # move .em.temp to .em
    fileList = sorted([f for f in os.listdir(output_folder) if (f[:len(prefix)] == prefix and f[-5:] == '.temp')])
    for fname in fileList:
        fname = os.path.join(output_folder, fname)
        os.system('mv {} {}'.format(fname, fname[:-len('.temp')]))

def mrc2em(filename,destination):


    if not os.path.exists(filename):
        raise RuntimeError('MRC file not found! ',filename)

    if not os.path.exists(destination):
        raise RuntimeError('Destination directory not found! ', destination)

    emfile = read(filename)
    splitName = filename.split(os.sep)
    filename = splitName[len(splitName)-1]
    newFilename = destination + os.sep + filename[0:len(filename)-4] + '.em'
    try:
        emfile.write(newFilename,'em')
    except:
        pass

def slurm_command(name='TemplateMatch',folder='./', cmd='', num_nodes=1,
                  modules = ['python3/3.7', 'lib64/append', 'openmpi/2.1.1', 'pytom/dev/gui_devel']):

    module_load = ''
    if modules:
        module_load = 'module load '
    for module in modules:
        module_load += module+' '

    slurm_generic_command = '''#!/usr/bin/sh
#SBATCH --time        12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    {}                                                                       
#SBATCH --output      {}/%x-%j.out 

{}

{}'''.format(name,folder,module_load,cmd)
    return slurm_generic_command

def gen_queue_header(name='TemplateMatch', folder='./', cmd='', num_nodes=1, emailaddress='',id='',
                     modules=['openmpi/2.1.1', 'python3/3.7', 'lib64/append', 'pytom/dev/gui_devel'], suffix= '',
                     qtype='slurm', num_jobs_per_node=20, time=12, partition='defq', singleton=False):
    module_load = ''
    queue_command = ''
    if modules:
        module_load = 'module load '
    for module in modules:
        module_load += module + ' '

    if singleton:
        singletoncommand = '#SBATCH --dependency=singleton'
    else:
        singletoncommand = ''

    if partition == 'fastq':
        oversubscribe = '#SBATCH --oversubscribe'
    else: oversubscribe = ''

    if qtype == 'slurm':
        queue_command = '''#!/usr/bin/bash
#SBATCH --time        {}:00:00
#SBATCH -N {}
#SBATCH --partition {}
#SBATCH --ntasks-per-node {}
#SBATCH --job-name    {}                                                                       
#SBATCH --output      {}/%j-%x{}.out 
{}
{}

{}

{}'''.format(time, num_nodes, partition, num_jobs_per_node, name, folder, suffix, oversubscribe, singletoncommand, module_load, cmd)

    elif qtype == 'qsub':
        print ('qsub has not been defined.')
        return ''

    elif qtype == 'torque':
        if emailaddress:
            emailaddress = '#PBS -M {}\n        #PBS -m abe'.format(emailaddress)
        queue_command = '''#!/usr/bin/bash
#PBS -k o 
#PBS -l nodes={}:ppn={},walltime={}:00:00 
{}
#PBS -N {} 
#PBS -j oe
#PBS -q {}
'''.format(num_nodes, num_jobs_per_node, time, emailaddress, name, partition)

    return queue_command

def createGenericDict(fname='template',cmd='', folder='', partition='defq', num_jobs_per_node=20, time=12, suffix='',
                      num_nodes=1, modules=['openmpi/2.1.1', 'python3/3.7', 'lib64/append', 'pytom/dev/gui_devel'],
                      id=''):

    genericSbatchDict = {'fname':fname,'cmd':cmd,'folder':folder, 'modules':modules, 'time':time, 'partition':partition,
                         'num_jobs_per_node': num_jobs_per_node, 'suffix': suffix, 'num_nodes': num_nodes, 'id': id}

    return genericSbatchDict

def sort( obj, nrcol ):
    obj.sort(key=lambda i: float(i[nrcol]))

def sort_str( obj, nrcol ):
    obj.sort(key=lambda i: str(i[nrcol]))

def avail_gpu(cutoff_busy=.25, cutoff_space = 0.5):
    lines = [line.split() for line in os.popen('nvidia-smi').readlines()]

    list_gpu, available_gpu = [],[]
    busy_list = []
    for n, line in enumerate(lines):
        if len(line)>2 and line[2] == 'GeForce':
            list_gpu.append(int(line[1]))
            busy = float(lines[n+1][12].strip('%'))/100.
            busy_list.append(busy)
            space = float(lines[n+1][8].strip('MiB'))/float(lines[n+1][10].strip('MiB'))
            if busy < cutoff_busy and space < cutoff_space:
                available_gpu.append(int(line[1]))
        
    comb = list(zip(available_gpu,busy_list))

    sort(comb,1)

    try:
        av, b = zip(*comb)
    except:
        av = []
    return av

def delete( path):
    '''the function moves a file or folder to the .trash folder in the folder below. if .trash does not exists it will be create'''
    import shutil
    if path.endswith('/'):
        path = path[:-1]
    if len(path) == 0 or len(path.strip('/'))== 0:
        return
        
    try:
        trash = os.path.join( os.path.dirname(path), '.trash/')
    except:
        trash= '.trash/'
    
    if not os.path.exists( trash ): os.mkdir( trash )
    
    if not os.path.exists("{}{}".format(trash,path)):
        if os.path.isdir(path): shutil.move(path, trash)
        else: shutil.move(path,trash)
    else:
        for n in range(100):
            if not os.path.exists("{}{}_{}".format(trash,path,n) ):
                if os.path.isdir(path): shutil.move(path, "{}{}_{}".format(trash,path,n))
                else: shutil.move(path,"{}{}_{}".format(trash,os.path.basename(path),n) )
                break    

def Radial_average(image, mask=None):
    """Calculates the radial average of an array of any shape,
    the center is assumed to be at the physical center."""
    import pylab
    [cx,cy] = image.shape
    if mask == None:
        mask = ones((cx,cy), dtype='bool8')
    else:
        mask = bool8(mask)
    axis_values = [arange(0,l) - l/2. + 0.5 for l in image.shape]
    radius = zeros((image.shape[-1]))
    for i in range(len(image.shape)):
        radius = radius + ((axis_values[-(1+i)][(slice(0, None), ) + (newaxis, )*i])**2)
    radius = int32(sqrt(radius))
   
    number_of_bins = radius[mask].max() + 1
    radial_sum = zeros(number_of_bins)
    weight = zeros(number_of_bins)
    for value, this_radius in zip(image[mask], radius[mask]):
        radial_sum[this_radius] += value
        weight[this_radius] += 1.
    return radial_sum / weight

def detect_shift(arr0,arr1,image=[]):
    x,y = image.shape
    cross = abs(fftshift( ifft2(fftshift(fftshift(fft2(arr0))*conj(fftshift(fft2(arr1)))))))**2
    locx,locy =  (abs((cross))**2).flatten().argmax()%y, (abs((cross))**2).flatten().argmax()/y
    return cross, locx-y/2, locy-x/2, cross[int(locy)][int(locx)]
    
def create_folder(foldername):
    '''Checks if foldername exists, if not it will create it.'''
    if not os.path.exists(foldername): os.mkdir( foldername ) 

def write_text2file(text,fname,mode='a'):
    create_folder( os.path.dirname(fname) ) 
    out = open(fname,mode)
    out.write(text)
    out.close()

def batch_tilt_alignment( fnames_tomograms='', projectfolder='.', num_procs=[], num_procs_per_proc=1, deploy=False,
                          queue=False, tiltseriesname='sorted/sorted', markerfile='sorted/markerfile.em',
                          targets='alignment', qcommand='sbatch', logfolder = './'):
    '''BATCHMODE: tilt alignment. Submits a number of sbatch jobs to slurm queueing system.
    Each job calculates the tilt aligment for each marker in a markerfile.
    It divides the number or jobs with respect to the num_procs.'''

    pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
    num_submitted_jobs = 0
    for n in range(len(num_procs)-1):
        cmd = multiple_alignment.format( d=(projectfolder, pytompath, num_procs[n], num_procs[n+1], num_procs_per_proc,
                                            tiltseriesname, markerfile, targets, fnames_tomograms) )

        if queue:
            jobname = 'Alignment_BatchMode_Job_{:03d}'.format(num_submitted_jobs)
            cmd = gen_queue_header(name=jobname, folder=logfolder, cmd=cmd )

        write_text2file(cmd, '{}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n), 'w')

        if deploy:
            if queue:
                os.system('{} {}/jobscripts/alignment_{:03d}.job'.format(qcommand, projectfolder, n))
                num_submitted_jobs += 1
            else:
                os.system('bash {}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n))
                num_submitted_jobs += 1

def create_folderstructure(folderstructure, enter, projectdir='.'):
    for mainfolder in sorted(folderstructure):

        if mainfolder in ('copy_files', 'run_scripts'): continue

        if not os.path.exists('%s/%s' % (enter, mainfolder)): os.system('mkdir %s/%s' % (enter, mainfolder))

        if len(folderstructure[mainfolder]) and not type(folderstructure[mainfolder]) == type([]):
            create_folderstructure(folderstructure[mainfolder], "%s/%s" % (enter, mainfolder))

def copy_files(folderstructure, enter, projectdir='.'):
    for mainfolder in sorted(folderstructure):
        if 'copy_files' == mainfolder:
            for fname in folderstructure['copy_files']:
                if len(fname) == 0:
                    continue
                elif os.path.exists(fname):
                    os.system('cp %s %s/%s/' % (fname, projectdir, enter))
                else:
                    continue
                    # raise Exception("\n\nFile does not exists! \n\n\t Please check if you have inserted the correct path: \n\t %s \n\n" % fname)
            continue
        if len(folderstructure[mainfolder]) and not type(folderstructure[mainfolder]) == type([]):
            copy_files(folderstructure[mainfolder], "%s/%s" % (enter, mainfolder), projectdir)

def create_project_filestructure(projectdir='.'):
    # with open('config.json') as json_data_file:
    #    folderstructure = json.load(json_data_file)

    folderstructure = {
        "01_Raw_Nanographs": {
            "copy_files": ["/Users/gijs/Documents/PostDocUtrecht/ExperimentalData/180221_CPXV12_Strep/em-fscollect"],
            "run_scripts": ["em-fscollect"]
        },
        "02_Preprocessed_Nanographs": {
            "Motion_corrected": "",
            "CTF_corrected": "",
            "copy_files": [""],
            "run_scripts": [""]
        },
        "03_Tomographic_Reconstruction": {
            ".tomoname": {
                "raw": "",
                "stacks": "",
                "imod": "",
                "ctf": "",
                "sorted": {
                    "excluded": ""
                },
                "sorted_binned": "",
                "alignment": "",
                "reconstruction": {
                    "INFR": {
                        "temp_files_unweighted": "",
                        "temp_files_binned": "",
                        "backup": ""
                    },
                    "WBP": {
                        "temp_files_weighted": "",
                        "temp_files_binned": "",
                        "backup": ""
                    },
                    "copy_files": ["reconstruction.sh"]
                },
                "copy_files": [""]
            },
            "copy_files": ["em_tomoprep", "em_prepbatch_utrecht.sh", "validate_generated_tomograms.py",
                           "reconstruction_batch.py"],
            "run_scripts": [""]
        },
        "04_Particle_Picking": {
            "Tomograms": "",
            "Picked_Particles": "",
            "Template_Matching": {
                "cross_correlation": "",
                "motlfiles": "",
                "classification": ""
            },
            "copy_files": [""],
            "run_scripts": [""]
        },
        "05_Subtomogram_Analysis": {
            "Subtomograms": "",
            "Reconstruction": "",
            "ParticlePolishing": "",
            "Alignment": {
                "FRM": "",
                "GLocal": ""
            },
            "Classification": {
                "CPCA": "",
                "AutoFocus":""
            },
            "copy_files": [""],
            "run_scripts": [""]
        },
        "06_Segmentation": {
            "copy_files": [""],
            "run_scripts": [""]
        },
        "LogFiles": ""
    }

    if not os.path.exists(projectdir):
        os.mkdir(projectdir)
    create_folderstructure(folderstructure, projectdir)

datatype0 = [('DefocusU', 'f4'),
            ('DefocusV', 'f4'),
            ('DefocusAngle', 'f4'),
            ('Voltage', 'i4'),
            ('SphericalAberration', 'f4'),
            ('AmplitudeContrast', 'f4'),
            ('PhaseShift', 'f4'),
            ('PixelSpacing', 'f4'),
            ('MarkerDiameter', 'i4'),
            ('TiltAngle', 'f4'),
            ('RotationTheta', 'f4'),
            ('InPlaneRotation', 'f4'),
            ('TranslationX', 'f4'),
            ('TranslationY', 'f4'),
            ('TranslationZ', 'f4'),
            ('Magnification', 'f4'),
            ('Intensity', 'f4'),
            ('FileName', 'U1000')]

headerText0 = 'data_\nloop_\n'
units0 = ['um', 'um', 'deg', 'kV', 'mm', '', 'deg', 'A', 'A', 'deg', 'deg', 'deg', 'px', 'px', 'px', '', '', '' ]
fmt0='%11.6f %11.6f %6.2f %4d %6.2f %4.2f %11.6f %11.6f %4d %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %5.3f %5.3f %s'


datatype = [('DefocusU', 'f4'),
            ('DefocusV', 'f4'),
            ('DefocusAngle', 'f4'),
            ('Voltage', 'i4'),
            ('SphericalAberration', 'f4'),
            ('AmplitudeContrast', 'f4'),
            ('PhaseShift', 'f4'),
            ('PixelSpacing', 'f4'),
            ('MarkerDiameter', 'i4'),
            ('TiltAngle', 'f4'),
            ('RotationTheta', 'f4'),
            ('InPlaneRotation', 'f4'),
            ('TranslationX', 'f4'),
            ('TranslationY', 'f4'),
            ('TranslationZ', 'f4'),
            ('Magnification', 'f4'),
            ('Intensity', 'f4'),
            ('ImageSize', 'i4'),
            ('AcquisitionOrder', 'i4'),
            ('FileName', 'U1000')]

headerText = 'data_\nloop_\n'
units = ['um', 'um', 'deg', 'kV', 'mm', '', 'deg', 'A', 'A', 'deg', 'deg', 'deg', 'px', 'px', 'px', '', '','px', '', '' ]
fmt='%11.6f %11.6f %6.2f %4d %6.2f %4.2f %11.6f %11.6f %4d %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %5.3f %5.3f %4d %3d %s'

for n, h in enumerate(datatype):
    headerText += '_{} {}\n'.format(h[0], '({})'.format(units[n])*(units[n]!=''))

for n, h in enumerate(datatype):
    headerText0 += '_{} {}\n'.format(h[0], '({})'.format(units[n])*(units[n]!=''))

datatypeMR = [('MarkerIndex', 'i4'),
              ('OffsetX',     'f4'),
              ('OffsetY',     'f4'),
              ('OffsetZ',     'f4'),
              ('PositionX',   'f4'),
              ('PositionY',   'f4'),
              ('PositionZ',   'f4')]

headerMarkerResults = 'data_\nloop_\n'
unitsMR = ['', 'px', 'px', 'px', 'px', 'px', 'px']
fmtMR='%3d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f'
for n, h in enumerate(datatypeMR):
    headerMarkerResults += '_{} {}\n'.format(h[0], '({})'.format(unitsMR[n])*(unitsMR[n]!=''))


datatypeAR = [('AlignmentTransX', 'f4'),
              ('AlignmentTransY', 'f4'),
              ('TiltAngle',       'f4'),
              ('InPlaneRotation', 'f4'),
              ('Magnification',   'f4'),
              ('FileName', 'U1000')]

headerAlignmentResults = 'data_\nloop_\n'
unitsAR = ['px', 'px', 'degrees', 'degrees', '', '']
fmtAR='%15.10f %15.10f %15.10f %15.10f %15.10f %s'
for n, h in enumerate(datatypeAR):
    headerAlignmentResults += '_{} {}\n'.format(h[0], '({})'.format(unitsAR[n])*(unitsAR[n]!=''))


LOCAL_ALIGNMENT_RESULTS = [('ParticleIndex', 'i4'),
                           ('AlignmentTransX', 'f4'),
                           ('AlignmentTransY', 'f4'),
                           ('TiltAngle', 'f4'),
                           ('InPlaneRotation', 'f4'),
                           ('Magnification', 'f4'),
                           ('FileName', 'U1000')]

headerLocalAlignmentResults = ''
unitsLAR = ['', 'px', 'px', 'degrees', 'degrees', '', '']
fmtLAR = '%7d %15.10f %15.10f %15.10f %15.10f %15.10f %s'
for n, h in enumerate(LOCAL_ALIGNMENT_RESULTS):
    headerLocalAlignmentResults += '{} {}\n'.format(h[0], '({})'.format(unitsLAR[n]) * (unitsLAR[n] != ''))


def headerline(line):

    if line.startswith('data_') or line.startswith('loop_') or line.startswith('_') or line.startswith('#'):
        return False
    else:

        return True

def loadstar(filename, dtype='float32'):
    print(filename)
    with open(filename, 'r') as f:
        lines = [line for line in f if headerline(line)]

        arr = numpy.genfromtxt(lines, dtype=dtype)
    return arr

def savestar(filename, arr, header='', fmt='', comments='#'):
    numpy.savetxt(filename, arr, comments=comments, header=header,fmt=fmt)

def update_metafile(filename, columnID, values ):
    metadata= loadstar(filename,dtype=datatype)
    metadata[columnID] = values
    savestar(filename,metadata,fmt=fmt,header=headerText)

def createMetaDataFiles(nanographfolder, mdocfiles=[], target='', mdoc_only=False):

    if not mdocfiles:
        mdocfiles = [mdocfile for mdocfile in os.listdir(nanographfolder) if mdocfile.endswith('.mdoc')]
    datafiles = [fname for fname in os.listdir(nanographfolder) if fname.split('.')[-1] in ('tif','mrc','em') and not (fname[0] in '0123456789')]

    annotated = {}

    tomo_from_filename = {}

    for df in datafiles:
        annotated[os.path.basename(df)] = 0

    tiltAxis = 180
    pixelSpacing = 1.5
    voltage = 300


    for nr, mdocfile in enumerate(sorted(mdocfiles)):
        metadata = []
        mdocfilepath = os.path.join(nanographfolder, mdocfile)
        header = False
        datadict = {'TiltAngle':tiltAxis, 'Magnification': 79000, 'Intensity':0.0, 'PixelSpacing':pixelSpacing,
                    'Defocus':3, 'RotationAngle':270, 'SubFramePath':'', 'Voltage':voltage, 'MarkerDiameter': 100,
                    'SphericalAberration':2.7, 'AmplitudeContrast': 0.08, 'PhaseContrast':0.,}
        for description, dtype in datatype:
            if not description in datadict.keys(): datadict[description] = 0.

        mdocdata = [line.split() for line in open(mdocfilepath).readlines() if len(line.split()) >=3]
        nrTiltImages = 0
        for nn, line in enumerate(mdocdata):
            if '[ZValue' == line[0] and not header:
                header = True
                continue

            if not header:
                if line[0]   == 'PixelSpacing': datadict[line[0]] = float(line[2])
                elif line[0] == 'Voltage':      datadict[line[0]] = int(line[2])
                elif line[0] == 'ImageSize':
                    try: datadict[line[0]] = min( int(line[2]), int(line[3]) )
                    except: pass
                elif 'axis' in line:            datadict['InPlaneRotation'] = float(line[6].strip(','))-250
                continue

            if line[0] in datadict.keys():
                if line[0] == 'RotationAngle': line[0] = 'RotationTheta'
                try: datadict[line[0]] = float(line[2])
                except:
                    fname = os.path.basename( line[2].replace('\\','/') )
                    datadict['FileName'] = fname
                    if fname in annotated.keys():
                        annotated[fname] += 1

            if '[ZValue' == line[0] or nn+1 == len(mdocdata):
                data = [0.,]*len(datatype)
                for n, (description, dtype) in enumerate(datatype):
                    if description in datadict.keys():
                        data[n] = datadict[description]
                if 'Defocus' in datadict.keys():
                    data[0] = datadict['Defocus']
                    data[1] = datadict['Defocus']

                if 'AcquisitionOrder' == datatype[-2][0]:
                    data[-2] = nrTiltImages
                    nrTiltImages +=1

                metadata.append(tuple(data))

        if len(metadata) == 0: continue
        a = numpy.rec.array(metadata, dtype=datatype)
        a = numpy.sort(a,order='TiltAngle')

        outname = mdocfilepath.replace('mdoc','meta')
        if target: outname = os.path.join(target, mdocfile.replace('mdoc','meta'))
        savestar(outname, a, fmt=fmt, header=headerText)

    if mdoc_only: return


    acquisition_order = {}
    size = {}

    for k, v in annotated.items():
        if v ==0:
            tomoname = k.split('_')[0]
            if not tomoname in tomo_from_filename.keys():
                tomo_from_filename[tomoname] = []

                if k.split('.')[-1] in ('mrc', 'em'):
                    vol = read(os.path.join(nanographfolder,k))
                    data = copy.deepcopy(vol2npy(vol))
                    size[tomoname] = numpy.min(data.shape)
                else:

                    #try:
                    #    data = imread( os.path.join(nanographfolder,k) )
                    #    size[tomoname] = numpy.min(data.shape)
                    #except:
                    aa = os.popen('header {} | grep "Number of columns, rows, sections" '.format( os.path.join(nanographfolder,k) )).read()[:-1]
                    dd = list(map(int, aa.split()[-3:-1]))
                    size[tomoname] = numpy.min(dd)

            #Find tiltangle

            if k.endswith('em'):
                fileHeader = read_em_header(os.path.join(nanographfolder, k))
                tiltangle = fileHeader.get_tiltangle()

            elif k.endswith('mrc'):
                tiltangle = read_angle(os.path.join(nanographfolder,k))
            else:
                tiltangle = 0.


            tomo_from_filename[tomoname].append([k.replace('[','_').replace(']','_').split('_'), tiltangle, k])

    for NR, (k, v) in enumerate(tomo_from_filename.items()):

        neg = numpy.array( [0,]*len(v[0][0]), dtype=int)
        tiltangles_header = []
        for list2, angle, fname in v:
            tiltangles_header.append(angle)
            for n, part in enumerate(list2):
                if '-' in part:
                    neg[n] += 1


        loc = numpy.argmax(neg[neg < len(v)])

        if not neg[loc] > 0:
            loc = -1

        tiltangles_header = numpy.array(tiltangles_header)
        metadata = [[0., ] * len(datatype), ] * len(v)

        for NR, (d,t) in enumerate(datatype):
            if d == 'TiltAngle':

                break

        if len(v) < 10:
            return

        if loc > -0.5:
            for i in range(len(v)):
                metadata[i][NR] = float(v[i][0][loc])
                metadata[i][-1] = v[i][2]
                if datatype[-3][0] == 'ImageSize':  metadata[i][-3] = size[k]
                metadata[i] = tuple(metadata[i])

        elif len(numpy.unique(numpy.round(tiltangles_header).astype(int))) == len(tiltangles_header):
            for i in range(len(v)):
                metadata[i][NR] = float(tiltangles_header[i])
                metadata[i][-1] = v[i][2]
                if datatype[-3][0] == 'ImageSize':  metadata[i][-3] = size[k]
                metadata[i] = tuple(metadata[i])

        else:
            continue

        a = numpy.rec.array(metadata, dtype=datatype)
        a = numpy.sort(a, order='TiltAngle')

        outname = '{}/{}.meta'.format(nanographfolder, v[0][0][0])
        savestar(outname, a, fmt=fmt, header=headerText)

def update_metadata_from_defocusfile(metafile, defocusfile):
    metadata = loadstar(metafile, dtype=datatype)

    data = open(defocusfile,'r')
    header = list(map(float, data.readlines()[0].split()))
    skiprows = 0

    for h in header[2:-1]:
        if abs(h) < 0.001:
            skiprows=1

    data.close()

    if skiprows: print('First Line of {} is ignorned'.format(os.path.basename(defocusfile)))

    for dd in (10,9,8,7,6,5):
        try:
            defocusResults = loadstar(defocusfile,skiprows=skiprows, usecols=range(0,dd))
        except:
            continue
        break

    if len(defocusResults) > len(metadata):
        raise Exception('Defocus file is larger than meta file thus cannot uniquely determine all angles.')

    tiltImages, columns = defocusResults.shape

    resultsNames = ['ID_start', 'ID_end', 'AngleStart', 'AngleEnd', 'DefocusU', 'DefocusV', 'DefocusAngle',
                    'PhaseShift', 'CutOn']
    if columns < 6 or defocusResults[0][5] < 400:
        resultsNames[5] = 'Empty'
        resultsNames[6] = 'Empty'
        columns+=2
    resultsNames = resultsNames[:columns]

    for n, line in enumerate(defocusResults):
        angleN = line[2]

        NN = 999
        for index, j in enumerate( metadata['TiltAngle']):
            if abs(j-angleN) < 0.1:
                NN = index
        if NN > 900:
            angleN = (line[2] + line[3] ) /2.
            for index, j in enumerate(metadata['TiltAngle']):
                if abs(j - angleN) < 0.1:
                    NN = index
        if NN > 900:
            raise Exception('Angles from Defocus File does not correspond to angles in metafile.')

        for query in ('DefocusU', 'DefocusV', 'DefocusAngle', 'PhaseShift'):
            for nn, ii in enumerate(resultsNames):
                if query == ii:
                    metadata[query][NN] = line[nn]
                    if query in ('DefocusU', 'DefocusV'):
                        metadata[query][NN] /= 1000.
                    break
        if 'Empty' in resultsNames:
            metadata['DefocusV'][NN] = metadata['DefocusU'][NN]
            metadata['DefocusAngle'][NN] = 0.

    savestar(metafile, metadata, fmt=fmt, header=headerText)

def parseImodShiftFile(imodShiftFile):
   
    data = loadstar(imodShiftFile)
    return data[:, -2:]

def addShiftToMarkFrames(mark_frames, shifts, metadata, excluded):
    tiltAngles = metadata['TiltAngle'][excluded]
    zeroIndex = abs(tiltAngles).argmin()
    zerox, zeroy = shifts[zeroIndex]*0

    for i in range(len(mark_frames)):
        mark_frames[i, :, 0] -= shifts[i][1] - zeroy
        mark_frames[i, :, 1] -= shifts[i][0] - zerox



    return mark_frames, shifts
