
from numpy.fft import fft2, ifft2, fftshift
from numpy import int32, arange, conj, zeros, ones, bool8, sqrt, newaxis
import os
import numpy
import copy
import time
from pytom_volume import read
from pytom.tools.files import checkFileExists, checkDirExists
from pytom.basic.files import read as read_pytom
from pytom.gui.mrcOperations import read_mrc
from pytom_numpy import vol2npy
from pytom.gui.guiSupportCommands import multiple_alignment

def read_markerfile(filename,tiltangles):
    if filename[-4:] == '.mrc':
        mark_frames = mrc2markerfile(filename, tiltangles)
    elif filename[-3:] == '.em':
        mark_frames = em2markerfile(filename, tiltangles)
    elif filename[-5:] == '.wimp':
        mark_frames = wimp2markerfile(filename, tiltangles)
    else:
        return 0

    return mark_frames
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
    mf = read_mrc(filename)
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



def popup_messagebox(mtype, title, message):
    '''This function creates an Error message, info statement, warning message, question, ask ok cancel reques\
t'''
    d = {"Error":showerror, 
         "Info":showinfo, 
         "Warning":showwarning, 
         "Question":askquestion, 
         "Cancel":askokcancel}
    d[mtype](title,message)


def mrc2em(filename,destination):


    if not os.path.exists(filename):
        raise RuntimeError('MRC file not found! ',filename)

    if not os.path.exists(destination):
        raise RuntimeError('Destination directory not found! ', destination)

    emfile = read(filename)
    splitName = filename.split(os.sep)
    filename = splitName[len(splitName)-1]
    newFilename = destination + os.sep + filename[0:len(filename)-4] + '.em'
    emfile.write(newFilename,'em')


def slurm_command(name='TemplateMatch',folder='./', cmd='', modules = ['openmpi/2.1.1', 'pytom/0.971']):

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


def gen_queue_header(name='TemplateMatch', folder='./', cmd='',
                     modules=['openmpi/2.1.1', 'python3/3.7', 'lib64/append', 'pytom/dev/python3'],
                     qtype='slurm'):
    module_load = ''
    if modules:
        module_load = 'module load '
    for module in modules:
        module_load += module + ' '

    if qtype == 'slurm':
        queue_command = '''#!/usr/bin/bash
#SBATCH --time        12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    {}                                                                       
#SBATCH --output      {}/%x-%j.out 

{}

{}'''.format(name, folder, module_load, cmd)

    if qtype == 'qsub':
        print ('qsub has not been defined.')
        return ''

    return queue_command

def createGenericDict(fname='template',cmd='', folder='',
                      modules=['openmpi/2.1.1', 'python3/3.7', 'lib64/append', 'pytom/dev/python3']):
    genericSbatchDict = {'fname':fname,'cmd':cmd,'folder':folder, 'modules':modules}
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

    av, b = zip(*comb)
    return av

def delete( path):
    '''the function moves a file or folder to the .trash folder in the folder below. if .trash does not exists it will be create'''
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

'''
def browse_for_entry( entry_box, dialog_type, remote=1, initialdir='.',dir=True,pattern='*'):
    """Browse for a file and put that filename in an entry box.

    @param entry_box: The entry box to put the filename in.
    @type entry_box: tkinter.Entry
    @param dialog_type: The type of dialog to display. Should be 'open' or
        'save as'
    @type dialog_type: str
    """
    #filetypes=[('dm4 files', '.dm4')]
    dialog_type = dialog_type.lower()
    if dialog_type == 'file': name = askopenfilename( initialdir=initialdir,remote=remote,pattern=pattern)
    elif dialog_type == 'dir':  name = askdirectory(initialdir=initialdir,remote=remote,pattern=pattern)
    else: raise ValueError('%s not a valid dialog type.' % dialog_type)
    
    entry_box.delete(0, tk.END)
    entry_box.insert(0, str(name) )
'''

def browse_for_entry( entry_box, dialog_type, remote=1, initialdir='.',dir=True,pattern='*',function=None,p='',mod=0):
    """Browse for a file and put that filename in an entry box.

    @param entry_box: The entry box to put the filename in.
    @type entry_box: tkinter.Entry
    @param dialog_type: The type of dialog to display. Should be 'dir' or 'file'.
    @type dialog_type: str
    @param remote: Search a remote folder for file/folder.
    @type remote: int
    @initialdir: starting folder for the search.
    @type initialdir: '.'.
    @param dir: search for folder.
    @type dir: Boolean.
    @param pattern: search string for files in folders.
    @type: str
    @param function: command executed at end of browse entry.
    @type: function
    @param p: parameters passed to function.
    @type p: str
    @param mod: parameter passed to function, describing the mode.
    @type mod: int
    """
    dialog_type = dialog_type.lower()

    if dialog_type == 'file': name = askopenfilename( initialdir=initialdir,remote=remote,pattern=pattern)
    elif dialog_type == 'dir':  name = askdirectory(initialdir=initialdir,remote=remote,pattern=pattern)
    else: raise ValueError('%s not a valid dialog type.' % dialog_type)
    
    if not function:
        entry_box.delete(0, tk.END)
        entry_box.insert(0, str(name))
        entry_box.xview("end")
    else:
        function(name, p, mod)
    
def create_folder(foldername):
    '''Checks if foldername exists, if not it will create it.'''
    if not os.path.exists(foldername): os.mkdir( foldername ) 

def write_text2file(text,fname,mode='a'):
    create_folder( os.path.dirname(fname) ) 
    out = open(fname,mode)
    out.write(text)
    out.close()


def batch_tilt_alignment( number_tomonames, fnames_tomograms='', projectfolder='.', num_procs=20, num_procs_per_proc=1, tiltseriesname='sorted/sorted',
                         markerfile='sorted/markerfile.em',targets='alignment', firstindex=1, lastindex=21, refindex=11, weightingtype=0, deploy=False, queue=False):
    '''BATCHMODE: tilt alignment. Submits a number of sbatch jobs to slurm queueing system. Each job calculates the tilt aligment for each marker in a markerfile.  It divides the number or jobs with respect to the num_procs.'''

    pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])

    for n in range(number_tomonames):

        if not n % num_procs == 0:
            continue

        cmd = multiple_alignment.format( d=(projectfolder, pytompath, n, min(number_tomonames,num_procs+n),
                                          num_procs_per_proc, tiltseriesname, markerfile, targets, projectfolder,
                                          firstindex, lastindex, refindex, weightingtype, fnames_tomograms) )

        if queue:
            cmd = gen_queue_header(name='Alignment', folder=projectfolder, cmd=cmd )

        write_text2file(cmd,'{}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n), 'w' )

        if deploy:
            if queue:
                os.system('sbatch {}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n))
            else:
                os.system('bash {}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n))



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
                "sorted": {
                    "excluded": ""
                },
                "sorted_binned": "",
                "alignment": {
                    "unweighted_unbinned": ""
                },
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
                "ccf_out_mirr": "",
                "motlfiles": "",
                "classification": ""
            },
            "copy_files": [""],
            "run_scripts": [""]
        },
        "05_Subtomogram_Analysis": {
            "Subtomograms": "",
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
        }
    }

    if not os.path.exists(projectdir):
        os.mkdir(projectdir)
    create_folderstructure(folderstructure, projectdir)
