import os
import glob
import numpy
import time
import shutil
import atexit
from multiprocessing import Manager, Process, cpu_count

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# pytom gui imports
from pytom.gui.guiStyleSheets import DEFAULT_STYLE_PROGRESSBAR
from pytom.gui.fiducialAssignment import FiducialAssignment
from pytom.gui.guiStructures import GuiTabWidget, SimpleTable, Worker
from pytom.gui.guiFunctions import sort, axis_angle_from_ar_file, loadstar
from pytom.gui.guiSupportCommands import templateAlignment, templateWBP, templateCTFPlotter, \
    templateCTFCorrectionImod, ParamsFileCTFPlotter, multiple_alignment, templateINFR
import pytom.gui.guiFunctions as guiFunctions
from pytom.gui.mrcOperations import remove_hot_pixels
from pytom.agnostic.io import read_size
from pytom.basic.datatypes import DATATYPE_METAFILE, DATATYPE_ALIGNMENT_RESULTS_RO as ALIGNRESULTS_ORDER, \
    DATATYPE_ALIGNMENT_RESULTS as ALIGNRESULTS_OLD


class TomographReconstruct(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(TomographReconstruct, self).__init__(parent)

        self.stage = 'v02_'
        self.addGeneralVariables()

        # Select Tomograms
        headers = ["Select Tomograms", "Create Markerfile", "Alignment", 'CTF Correction', 'Reconstruction']
        subheaders = [[], [], ['Individual Alignment', 'Batch Alignment'],
                      ['CTF Determination', 'CTF Correction', 'Batch Correction'],
                      ['INFR', 'WBP', 'Batch Reconstruction']]
        tabUIs = [[self.tab1UI],
                  [self.tab2UI],
                  [self.tab31UI,self.tab32UI],
                  [self.tab41UI,self.tab42UI,self.tab43UI],
                  [self.tab51UI, self.tab52UI,self.tab53UI]]
        static_tabs = [[False], [True], [True, False], [True, True, False], [True, True, False]]
        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, tabUIs=tabUIs, tabs=self.tabs_dict,
                     tab_actions=self.tab_actions, static_tabs=static_tabs)

    def tab1UI(self,  id=''):
        self.filepath_tomodata = {}
        self.filepath_tomodata['Motion Corrected'] = self.motioncor_folder
        self.filepath_tomodata['Raw Nanographs']  = self.rawnanographs_folder

        headers = ["meta file", "create", 'square', 'pututive name',"input files", 'Name tomogram', 'delete', 'Bin Factor IMOD']
        types = ['txt2', 'checkbox','checkbox', 'txt', 'combobox', 'txt', 'checkbox', 'lineedit', 'txt']
        sizes = [0, 80, 80, 500, 0, 200, 80, 80]

        tooltip = ['Name of meta files. This file indicates which images belong together, as well as the tiltangles.',
                   'Create new tomogram folders.',
                   'Crop the images such that they become square',
                   'Putative name of tomogram',
                   'Which files do you want to use',
                   'Names of existing tomogram folder.',
                   'Select items to delete tomogram folders.',
                   'OPTIONAL: If your meta file originates from an imod stack, this binning factor can be used to create\n'
                   ' pytoms markerfile from imod defined alignment points. Binning used in IMOD should be given here.\n'
                   'Binning factor is ignored for other meta files.']

        processed = sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder)))
        processed_dict = {}
        metadata_dict = {}

        for i in processed:
            processed_dict[os.path.basename(i)] = i

        processed_fn = [os.path.basename(line) for line in processed]
        unprocessed = sorted(glob.glob('{}/*.meta'.format(self.rawnanographs_folder)))
        unprocessed = numpy.array(unprocessed + sorted(glob.glob('{}/import*/*.meta'.format(self.rawnanographs_folder))))

        ff = []

        for u_item in unprocessed:
            bn = os.path.basename(u_item)
            if bn in processed_fn:
                if not bn in metadata_dict.keys():
                    metadata_dict[bn] = loadstar(processed_dict[os.path.basename(u_item)],dtype=DATATYPE_METAFILE)
                c = loadstar(u_item,dtype=DATATYPE_METAFILE)

                if len(metadata_dict[bn]['TiltAngle']) != len(c['TiltAngle']):
                    print('double accepted', u_item)
                    ff.append(u_item)
                    continue

                if not abs(metadata_dict[bn]['TiltAngle'] - c['TiltAngle']).sum() < 0.001:
                    print('double accepted', u_item)
                    ff.append(u_item)
            else:
                ff.append(u_item)
                print('singlet accepted', u_item)

        unprocessed = ff
        values = []

        for t in processed:
            values.append([t, False, False, '', ['Raw Nanographs', 'Motion Corrected'], t.split('/')[-3], True, self.binningFactorIMOD])
        for t in unprocessed:
            values.append([t, True, True, '', ['Raw Nanographs','Motion Corrected'], '', False, self.binningFactorIMOD])

        if values:
            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, connect=self.update_create_tomoname)

        else:
            self.popup_messagebox('Warning','No meta files found', 'No meta files found. Have you downloaded your mdoc files, or are file names sufficient to determine the tiltangles?')
            return

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.create_tomogram_folders(pid, v))

    def tab2UI(self,  id=''):

        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = self.table_layouts[id]

        last, reftilt = 10, 5

        self.insert_pushbutton(parent,cstep=1, text='Create Markerfile', action=self.startFidAssignment, params=[parent],
                               wname='startFidAss',rstep=1)
        self.widgets['startFidAss'].setSizePolicy(self.sizePolicyC)
        self.insert_label(parent,sizepolicy=self.sizePolicyA)
        pass

    def tab31UI(self, id=''):

        self.row, self.column = 0,0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        h = mode = self.stage + 'SingleAlignment_'

        last, reftilt = -60, 60
        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted Tilt Images', mode +'FolderSorted',
                                    'Select the folder where the sorted tiltimages are located.\n')
        self.insert_label_spinbox(parent, mode+'FirstAngle', text='First Angle',
                                  tooltip='Lowest tilt angle of tilt images used for alignment.',
                                  value=-60,minimum=-90,maximum=90, stepsize=1)
        self.insert_label_spinbox(parent, mode +'LastAngle', text='Last Angle',
                                  tooltip='Highest tilt angle of tilt images used in alignment.',
                                  value=60, minimum=-90, maximum=90, stepsize=1)
        self.insert_label_spinbox(parent, mode +'RefTiltIndex', text='Reference Tilt Image', value=reftilt,minimum=1,
                                  tooltip='Index of reference tilt image. Typically zeros degree image.')
        self.insert_label_spinbox(parent, mode +'RefMarkerIndex', text='Reference Marker', value=1,
                                  tooltip='Index of reference marker. See previous step.')
        self.insert_label_spinbox(parent, mode + 'RotationTiltAxis', text='Angle Tilt Axis (degrees)',
                                  value=0, minimum=0, maximum=359,
                                  tooltip='Angle of the tilt axis (degrees). 0 degrees is facing north, '+
                                          '90 degrees is facing east.')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor',value=1, cstep=0,
                                  tooltip='Binning factor used for reconstruction')
        #self.insert_label(parent,'Orientation tilt axis', rstep=1, cstep=1,
        #                  tooltip='Orientation of the tiltaxis with respect to the orientation of the camera.'+
        #                 'The point of the arrow indicates the rotation direction, following the left hand rule.')
        #self.widgets[mode + 'tomofolder'] = QLineEdit()
        #self.widgets[mode + 'tomogramNR'] = QLineEdit()
        for name in ('tomofolder', 'tomogramNR','FirstIndex', 'LastIndex', 'Reduced', 'outFolder', 'markerfile',
                     'tiltSeriesName'):
            self.widgets[mode + name] = QLineEdit()

        self.widgets[mode + 'FolderSorted'].textChanged.connect(lambda dummy, m=mode: self.updateTomoFolder(m))
        self.widgets[mode + 'FirstAngle'].valueChanged.connect(lambda dummy, m=mode: self.updateIndex(m))
        self.widgets[mode + 'LastAngle'].valueChanged.connect(lambda dummy, m=mode: self.updateIndex(m))

        self.updateTomoFolder(mode)

        execfilename = [mode + 'tomofolder', 'alignment/alignment.sh']

        queue_name = 'Alignment'
        self.queue_job_names.append(queue_name)
        paramsSbatch = guiFunctions.createGenericDict(fname=queue_name,folder=self.logfolder, id='SingleAlignment')
        paramsCmd    = [mode + 'tomofolder', self.parent().pytompath, mode + 'tiltSeriesName', mode + 'FirstIndex',
                        mode + 'LastIndex', mode + 'RefTiltIndex', mode + 'RefMarkerIndex', mode + 'markerfile',
                        mode + 'outFolder', mode + 'BinningFactor', '0', mode+'RotationTiltAxis',
                        templateAlignment]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=execfilename, paramsSbatch = paramsSbatch,
                                 paramsCmd=paramsCmd, action=self.convert_em, paramsAction=[mode,'alignment','sorted'],
                                 mandatory_fill = [mode + 'FolderSorted'])
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)
        self.updateTomoFolder(mode)
        self.updateIndex(mode)

    def tab32UI(self, id=''):
        headers = ["name tomogram", "align", 'First Angle',"Last Angle", 'Ref. Image', 'Ref. Marker', 'Exp. Rot. Angle', 'Input Folder', 'Fix Marker Pos', 'Ref. Marker Tomo', '']
        types = ['txt', 'checkbox', 'lineedit', 'lineedit', 'lineedit', 'combobox', 'lineedit', 'combobox', 'checkbox', 'combobox', 'txt']
        sizes = [0, 80, 0, 0, 0, 0, 0, 0, 0, 0 ]

        tooltip = ['Names of existing tomogram folders.',
                   'Do alignment.',
                   'First angle of tiltimages.',
                   'Last angle of tiltimages.',
                   'Reference image number.',
                   'Reference marker index',
                   'Expected Rotation Angle',
                   'Input Folder for tilt images used in alignment.',
                   'Fix marker position during optimization step. Useful for local marker alignment.',
                   'Marker ID of reference marker in tomogram. Needed for Local Marker Refinement.']

        markerfiles = sorted(glob.glob('{}/tomogram_*/sorted/markerfile.txt'.format(self.tomogram_folder)))
        markerfilesEM = sorted(glob.glob('{}/tomogram_*/sorted/markerfile.em'.format(self.tomogram_folder)))
        novel = [markerfile for markerfile in markerfilesEM if not markerfile[:-3]+'.txt' in markerfiles]
        markerfiles += novel
        markerfiles = sorted(markerfiles)
        values = []
        self.mfiles = []
        for markerfile in markerfiles:
            qmarkerfile = os.path.join(os.path.dirname(markerfile), '*.meta')
            qsortedfiles = os.path.join(os.path.dirname(markerfile), 'sorted_*.mrc')

            #qmarkerfile = markerfile.replace('markerfile.txt','*.meta')
            #qsortedfiles = markerfile.replace('markerfile.txt','sorted_*.mrc')
            if 1:
                metafile = glob.glob( qmarkerfile )[0]
                metadata = loadstar(metafile,dtype=guiFunctions.datatype)
                tangs = metadata['TiltAngle']
                sortedfiles = sorted(glob.glob(qsortedfiles))
                last_frame = len(sortedfiles)
                index_zero_angle = 0
                mm = 9999

                for n, sortedfile in enumerate(sortedfiles):
                    #print(sortedfile)
                    index_s = int(sortedfile.split('_')[-1].split('.')[0])
                    if abs(tangs[index_s]) < mm:
                        mm = abs(tangs[index_s])
                        index_zero_angle = index_s


                data = guiFunctions.readMarkerfile(markerfile, len(sortedfiles))

                if len(data.shape) < 3: continue
                options_reference = list(map(str, range( data.shape[2] ))) + ['all']
                expect = int(float(metadata['InPlaneRotation'][0]))
                input_folders = ['sorted', 'ctf/sorted_ctf']
                values.append( [markerfile.split('/')[-3], True, numpy.floor(tangs.min()), numpy.ceil(tangs.max()),
                                index_zero_angle, options_reference, expect, input_folders, True, options_reference[:-1], ''] )
                self.mfiles.append(markerfile)
            else:continue
        if not values:
            return
        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, wname=self.stage + 'BatchAlignment_')
        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.run_multi_align(pid, v))

    def tab41UI(self, id=''):
        grid = self.table_layouts[id]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'CTFDetermination_'

        items += list(self.create_expandable_group(self.ctfDetermination, self.sizePolicyB, 'CTF Determination',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab42UI(self, id=''):
        grid = self.table_layouts[id]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleCTFCorrection_'

        items += list(self.create_expandable_group(self.ctfCorrection, self.sizePolicyB, 'CTF Correction',
                                                   mode=t0))
        items[-1].setVisible(False)
        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab43UI(self, id=''):

        headers = ['Name Tomogram', 'Correct', 'Index GPU', 'Defocus Tol', 'Width Interpol.', 'Pixel Size',
                   'Spherical Aberr.','Amp. Contrast', 'Voltage', 'Rot. Angle', 'Input Stack', 'Angle File',
                   'Defocus File', 'Output Folder', '']
        types = ['txt', 'checkbox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit',
                 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'txt']
        sizes = [0, 0, 30, 30, 30, 430, 30, 30, 30, 0, 0, 0, 0, 0, 0, ]

        tooltip = ['Names of existing tomogram folders.',
                   'Run checkbox',
                   'GPU ID (only one gpu allowed)',
                   'Defocus tolerance (nm)',
                   'Interpolation width (pixels)',
                   'Pixel Size (nm)',
                   'Spherical Aberration (mm)',
                   'Amplitude Contrast',
                   'Voltage (kV)',
                   'Rotation axis angle (degrees)',
                   'Input Stack.',
                   'Angle File',
                   'Defocus file (from CTF Determination).',
                   'Output Folder']

        defocusfiles = sorted(glob.glob('{}/tomogram_*/ctf/*.defocus'.format(self.tomogram_folder)))

        values = []

        for defocusfile in defocusfiles:
            tomogramName = os.path.dirname(os.path.dirname(defocusfile))
            ctffolder = os.path.dirname(defocusfile)
            ctfplotter = os.path.join(ctffolder, 'ctfplotter.com')

            try:

                if os.path.exists(ctfplotter):
                    d = {}
                    keys = ('InputStack', 'AngleFile', 'DefocusFile', 'AxisAngle', 'PixelSize', 'Voltage',
                            'SphericalAberration', 'AmplitudeContrast', 'DefocusTol')
                    for line in open(ctfplotter).readlines():
                        kvpair = line.split()
                        if len(kvpair) == 2 and kvpair[0] in keys:
                            k, v = kvpair
                            d[k] = v
                    origdir = os.path.dirname(d['InputStack'])
                    prefix = os.path.basename(origdir) + '_ctf'
                    outfolder = os.path.join(ctffolder, prefix)
                    if not os.path.exists(outfolder):
                        os.mkdir(outfolder)
                    values.append([tomogramName, True, '', d['DefocusTol'], 15, d['PixelSize'], d['SphericalAberration'],
                                   d['AmplitudeContrast'], d['Voltage'], d['AxisAngle'], d['InputStack'],
                                   d['AngleFile'], defocusfile, outfolder, ''])

            except Exception as e:
                print(e)
                continue

        if values:
            try:
                self.num_nodes[id].setParent(None)
            except:
                pass
            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, nn=True,
                          wname=self.stage + 'BatchCTFCorrection')
            self.tab43_widgets = self.tables[id].widgets
            self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.run_multi_ctf_correction(pid, v))

    def tab43UIOld(self, id=''):

        headers = ["Name Tomogram", 'Correct', "Reference Marker", "Defocus File", "Rot. Angle", 'Grid Spacing', 'Field Size', 'Bin Factor', '']
        types = ['txt', 'checkbox', 'combobox', 'combobox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'txt']
        sizes = [0, 80, 80, 0, 0, 0, 0, 0]

        tooltip = ['Names of existing tomogram folders.',
                   'CTF Correction.',
                   'Reference Marker.',
                   'Defocus file (from correction).',
                   'In-plane Rotation Angle. Consider which angle you have used last for the alignment.',
                   'Grid spacing [2,4,6,...] is the size of the area which is corrected with a' +
                   'constant ctf. The ctf value is taken from the center of this area. This ' +
                   'area builds up the corrected projection.',
                   'Field size[2, 4, 6, ...] & (fs >= gs) is the size of the area which is extracted ' +
                   'from the projection and corrected with a constant ctf. Fieldsize is also ' +
                   'the size of the modelled ctf.',
                   'Binning factor applied to images.']

        markerfiles = sorted(glob.glob('{}/tomogram_*/sorted/markerfile.em'.format(self.tomogram_folder)))

        values = []

        for markerfile in markerfiles:
            tomogramName = os.path.dirname( os.path.dirname(markerfile) )
            query = '{}/alignment/*marker*'.format(tomogramName)
            folders = [folder for folder in glob.glob(query) if os.path.isdir(folder)]
            folders = sorted( [folder for folder in folders if len(os.listdir(folder)) > 1] )
            if not folders: continue
            referenceMarkerOptions = [os.path.dirname(markerfile)] + folders + ['all']
            defocusFile = list(reversed(sorted(glob.glob( os.path.join(tomogramName, 'ctf/*.defocus') ))))
            if not defocusFile: continue
            values.append( [tomogramName, True, referenceMarkerOptions, defocusFile, 0, 16, 256, 1, ''] )

        if values:
            try:
                self.num_nodes[id].setParent(None)
            except:
                pass
            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, nn=True)
            self.tab43_widgets = self.tables[id].widgets
            self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.run_multi_ctf_correction(pid, v))

    def tab51UI(self, id=''):

        return  # TODO this options needs to be fixed

        alg = 'ReconstructINFR'
        mode = 'v02_{}_'.format(alg)
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]

        last, reftilt = 10, 5
        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB, width=400)
        self.insert_label_line_push(parent, 'Folder Sorted Tilt Images', mode + 'FolderSorted',
                                    'Select the folder where the sorted tiltimages are located.\n')
        self.insert_label_spinbox(parent, mode + 'FirstAngle', text='First Angle',
                                  tooltip='Tilt angle of first image (deg).',
                                  value=-60, minimum=-90, maximum=90, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'LastAngle', text='Last Angle',
                                  tooltip='Tilt angle of last image (deg).',
                                  value=60, minimum=-90, maximum=90, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'RefTiltIndex', text='Reference Tilt Image', value=reftilt, minimum=1,
                                  tooltip='Index of reference tilt image. Typically zeros degree image.')
        self.insert_label_spinbox(parent, mode + 'RefMarkerIndex', text='Reference Marker', value=1,
                                  tooltip='Index of reference marker. See previous step.')
        self.insert_label_spinbox(parent, mode + 'RotationTiltAxis', text='Angle Tilt Axis (degrees)',
                                  value=0, minimum=0, maximum=359,
                                  tooltip='Angle of the tilt axis (degrees). 0 degrees is facing norther, ' +
                                          '90 degrees is facing east.')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor', value=8, cstep=0,
                                  tooltip='Binning factor used for reconstruction')

        # self.insert_label(parent,'Orientation tilt axis', rstep=1, cstep=1,
        #                  tooltip='Orientation of the tiltaxis with respect to the orientation of the camera.'+
        #                 'The point of the arrow indicates the rotation direction, following the left hand rule.')

        for name in ('tomofolder', 'tomogramNR', 'FirstIndex', 'LastIndex', 'Reduced', 'tiltSeriesName', 'markerfile',
                     'specimenAngleFlag', 'DimY'):
            if mode + name in self.widgets.keys():
                continue
            self.widgets[mode + name] = QLineEdit()

        self.widgets[mode + 'FolderSorted'].textChanged.connect(lambda dummy, m=mode: self.updateTomoFolder(m))
        self.updateTomoFolder(mode)
        self.widgets[mode + 'FirstAngle'].valueChanged.connect(lambda dummy, m=mode: self.updateIndex(m))
        self.widgets[mode + 'LastAngle'].valueChanged.connect(lambda dummy, m=mode: self.updateIndex(m))

        execfilename = [mode + 'tomofolder', 'reconstruction/INFR/INFR_reconstruction.sh']
        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'ReconstructionINFR'
        paramsSbatch['folder'] = self.logfolder  # os.path.dirname(execfilename)
        #paramsSbatch['partition'] = 'fastq'
        paramsSbatch['time'] = 1
        paramsSbatch['num_jobs_per_node'] = 1
        paramsSbatch['id'] = 'ReconstructINFR'
        self.insert_gen_text_exe(parent, mode, jobfield=False, action=self.convert_em,
                                 paramsAction=[mode, 'reconstruction/INFR', 'sorted'],
                                 exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=[mode + 'tomofolder', self.parent().pytompath, mode + 'FirstIndex',
                                            mode + 'LastIndex',
                                            mode + 'RefTiltIndex', mode + 'RefMarkerIndex', mode + 'BinningFactor',
                                            self.parent().pytompath, mode + 'tomogramNR', mode + 'RotationTiltAxis',
                                            mode + 'specimenAngleFlag',
                                            templateINFR], mandatory_fill=[mode+'FolderSorted'])

        # self.insert_label_action_label(parent,'Generate command',cstep=-2, sizepolicy=self.sizePolicyB)
        # self.insert_textfield(parent, h+'CommandText', columnspan=3, rstep=1, cstep=2)
        # self.insert_label_action_label(parent,'Execute command', rstep=1)
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)

    def tab52UI(self, id=''):
        alg = 'ReconstructWBP'
        mode = 'v02_{}_'.format(alg)
        self.row, self.column = 0,0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]

        available_tomograms = [os.path.basename(d) for d in sorted(glob.glob('{}/tomogram_*'.format(
            self.tomogram_folder)))]

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB, width=400)
        self.insert_label_combobox(parent, 'Tomogram', mode + 'tomogram', labels=available_tomograms,
                                   tooltip='Select the tomogram to be reconstructed.')
        self.insert_label_combobox(parent, 'Alignment', mode + 'alignment', labels=[],
                                   tooltip='Select the type of alingment to be used.')
        self.insert_label_spinbox(parent, mode + 'firstAngle', 'First Angle',
                                  'Select the first angle of the angular range used for reconstruction.',
                                  minimum=-90, maximum=90, stepsize=1, value=0)
        self.insert_label_spinbox(parent, mode + 'lastAngle', 'Last Angle',
                                  'Select the last angle of the angular range used for reconstruction.',
                                  minimum=-90, maximum=90, stepsize=1, value=0)
        self.insert_label_combobox(parent, 'Ctf corrected', mode + 'ctfCorrChoice', labels=[],
                                   tooltip='Select original (sorted) or ctf corrected images (sorted_ctf).')
        # self.insert_label_spinbox(parent, mode + 'SpecimenAngle', text='Specimen angle (degrees)',
        #                           value=0, minimum=-90, maximum=90,
        #                           tooltip='Angle of the specimen, will be added to the tilt axis rotation. Can be '
        #                                   'used to place the sample horizontally in the tomogram.')
        self.insert_label_spinbox(parent, mode + 'WeightingType', text='Weighting Type',
                                  value=1, minimum=-1, maximum=1, stepsize=1,
                                  tooltip='Select weighting type:\n\t 0: no weighting\n\t-1: analytical weighting'+
                                          '\n\t 1: "exact" weighting')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor', value=8, minimum=1,
                                  tooltip='Binning factor used for reconstruction')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', validator=QIntValidator(),
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs "
                                       "separate them using a comma, e.g. 0,1,2 ")
        self.insert_label_line(parent, "Cores", mode + 'cores', cstep=0, validator=QIntValidator(),
                               tooltip="How many procs to use if run on cpus")

        # initialize backend info for reconstruction
        for name in ('tomofolder', 'SpecimenAngleFlag', 'DimY', 'DimX', 'DimZ', 'Voldims', 'RotationTiltAxis',
                     'alignmentResultsFile', 'projectionDir', 'gpuString', 'cores_flag'):
            self.widgets[mode + name] = QLineEdit('')

        # update all alignment options whenever tomogram is changed
        self.widgets[mode + 'tomogram'].currentIndexChanged.connect(lambda dummy, m=mode:
                                                                    self.update_tomogram_reconstruction(m))
        self.update_tomogram_reconstruction(mode)

        # set the specimen angle flag if a value is selected
        # self.widgets[mode + 'SpecimenAngle'].valueChanged.connect(lambda dummy, m=mode: self.updateSpecimenAngle(m))
        # self.updateSpecimenAngle(mode)

        # volume dims can change for different bin, rotation axis, and projections
        self.widgets[mode + 'BinningFactor'].valueChanged.connect(lambda dummy, m=mode: self.update_vol_dims(m))
        self.widgets[mode + 'RotationTiltAxis'].textChanged.connect(lambda dummy, m=mode: self.update_vol_dims(m))
        self.widgets[mode + 'projectionDir'].textChanged.connect(lambda dummy, m=mode: self.update_vol_dims(m))
        self.update_vol_dims(mode)

        # update gpu string
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, ids=self.widgets[mode + 'gpuID'],
                                                         flag=self.widgets[mode + 'gpuString']:
                                                         self.update_gpu_flag(ids, flag, single_gpu=True))
        self.widgets[mode + 'cores'].textChanged.connect(lambda d, inp=self.widgets[mode + 'cores'],
                                                         flag=self.widgets[mode + 'cores_flag']:
                                                         self.update_cores_flag(inp, flag, suppress_message=False))

        execfilename = [mode + 'tomofolder', 'reconstruction/WBP/WBP_Reconstruction.sh']

        # TODO Remove these parameters as these should not be used. Instead we should just stick with the qparams.
        paramsSbatch = guiFunctions.createGenericDict(fname='ReconstructionWBP', folder=self.logfolder,
                                                      id='ReconstructWBP')

        paramsCmd = [mode + 'tomofolder', self.parent().pytompath, mode + 'BinningFactor', mode + 'tomogram',
                     mode + 'Voldims', mode + 'WeightingType', mode + 'RotationTiltAxis',
                     mode + 'SpecimenAngleFlag', mode + 'DimY', mode + 'DimZ', mode + 'alignmentResultsFile',
                     mode + 'gpuString', mode + 'projectionDir', mode + 'firstAngle', mode + 'lastAngle',
                     mode + 'cores_flag', templateWBP]

        self.insert_gen_text_exe(parent, mode, exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, mandatory_fill=[mode + 'tomofolder',
                                                                      mode + 'alignmentResultsFile'],
                                 paramsAction=[mode, 'reconstruction/WBP', 'sorted'], jobfield=False)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)

    def tab53UI(self, key='tab53'):
        # TODO INFR option needs to be put back in
        # TODO specimen angle back in?
        headers = ['name tomogram', 'WBP', 'alignment', 'first angle', 'last angle', 'ctf', 'weighting for WBP',
                   'binning', 'gpu id', '']
        types = ['txt', 'checkbox', 'combobox', 'lineedit', 'lineedit', 'combobox', 'lineedit', 'lineedit',
                 'lineedit', '']
        sizes = [0, 80, 0, 0, 0, 0, 0, 0, 0, 0]
        tooltip = ['Names of existing tomogram folders.',
                   # 'Do INFR Reconstruction.',
                   'Do WBP Reconstruction.',
                   'Select alignment type',
                   'First angle of tiltimages',
                   'Last angle of tiltimages',
                   'Whether to use ctf corrected images',
                   'Binning factor applied to images.',
                   'Weighting type: -1 ramp weighting, 1 analytical weighting, 0 no weighting',
                   'GPU ID, only applicable for WPB currently (e.g 3)']  # 'Angle of specimen around tilt-axis',

        # ctf combobox can be filled from beginning as the tomogram does not change
        # alignment options can as well be filled from beginning
        # there should be a link between alignment and the first/last angle. Everytime alignment is changed the first
        # last angle should be updated

        # get all tomograms present in the project
        tomogram_folders = sorted(glob.glob('{}/*'.format(self.tomogram_folder)))
        tomogram_folders = [t for t in tomogram_folders if not 'jobscripts' in os.path.basename(t)]
        tomograms = [os.path.basename(d) for d in tomogram_folders]

        # initialize list for data per tomo
        values = []

        # get info per tomogram about alignment options and ctf options
        for tomo_name, tomo_folder in zip(tomograms, tomogram_folders):

            # find available alignments
            alignment_dir = os.path.join(tomo_folder, 'alignment')
            alignment_choices = []
            for alignment in os.listdir(alignment_dir):
                d = os.path.join(alignment_dir, alignment)
                if os.path.exists(d) and os.path.isdir(d):
                    for pointer in os.listdir(os.path.join(d, 'GlobalAlignment')):
                        f = os.path.join(d, 'GlobalAlignment', pointer, 'alignmentResults.txt')
                        if os.path.exists(f):
                            alignment_choices.append(alignment)
                            break
                    # we might have found an alignment file

            # get ctf options
            sorted_dir = os.path.join(tomo_folder, 'sorted')
            ctf_sorted_dir = os.path.join(tomo_folder, 'ctf', 'sorted_ctf')
            tilt_choices = []
            if len([f for f in os.listdir(sorted_dir) if f.endswith('.mrc')]) > 0:
                tilt_choices.append('sorted')
            if len([f for f in os.listdir(ctf_sorted_dir) if f.endswith('.mrc')]) > 0:  # at least
                tilt_choices.append('sorted_ctf')

            # add to table fill, angles will be set in a bit
            values.append([tomo_name, True, alignment_choices, 0, 0, tilt_choices, -1, 8, '', ''])

        # set the base values and the connect alignment type with angle reading
        self.fill_tab(key, headers, types, values, sizes, tooltip=tooltip, wname=self.stage + 'BatchReconstruct')
        self.widgets[key + 'gpu_flag'] = QLineEdit('')
        for row in range(len(values)):
            w = self.tables[key].widgets['widget_{}_2'.format(row)]  # alignment choice widget
            w.currentIndexChanged.connect(lambda d, r=row, k=key: self.update_alignment_choice_batch(r, k))
            self.update_alignment_choice_batch(row, key)

            # gpu flag is not used but just to check if the input is valid
            w = self.tables[key].widgets['widget_{}_8'.format(row)]
            w.textChanged.connect(lambda d, ids=w, flag=self.widgets[key + 'gpu_flag']:
                                  self.update_gpu_flag(ids, flag, single_gpu=True, suppress_message=True))

        self.pbs[key].clicked.connect(lambda dummy, pid=key: self.run_multi_reconstruction(pid))

    def fill_tab(self, id, headers, types, values, sizes, tooltip=[], connect=0, nn=False, nc=False, wname=''):
        try:
            self.tables[id].setParent(None)
            self.pbs[id].setParent(None)
            self.ends[id].setParent(None)
            self.checkbox[id].setParent(None)
        except:
            pass

        self.tables[id] = SimpleTable(headers, types, values, sizes, tooltip=tooltip, connect=connect)
        self.widgets['v02_batch_aligntable{}'.format(id)] = self.tables[id]
        self.pbs[id] = QPushButton('Run')
        self.pbs[id].setSizePolicy(self.sizePolicyC)
        self.ends[id] = QWidget()
        self.ends[id].setSizePolicy(self.sizePolicyA)
        self.widgets[f'{wname}ExecutePushButton'] = self.pbs[id]
        if nn:
            num_nodes = QSpinBox()
            num_nodes.setValue(2)
            num_nodes.setRange(1, 9)
            num_nodes.setPrefix('Numbre of nodes: ')

        try:
            self.num_nodes[id] = num_nodes
        except:
            self.num_nodes = {}
            self.num_nodes[id] = 0   # should this not be num_nodes ??

        for n, a in enumerate((self.tables[id], self.num_nodes[id], self.checkbox[id], self.pbs[id], self.ends[id])):
            if n==1 and not nn: continue
            if n==2:
                self.widgets['{}_queue'.format(id)] = a
                a.setEnabled(self.qtype != 'none')
            self.table_layouts[id].addWidget(a)

    def startFidAssignment(self,parent=None):
        self.fidass = FiducialAssignment(self)
        self.fidass.show()
        pass

    def ctfDetermination(self, mode, title=''):

        tooltip = 'Run IMOD ctfplotter.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted & Aligned Tilt Images', mode + 'FolderSortedAligned',
                                    'Select the folder where the sorted tiltimages are located.\n',
                                    initdir=self.tomogram_folder, mode='folder')
        self.insert_label_line_push(parent, 'Angle File', mode + 'AngleFile',mode='file',filetype='tlt',
                                    tooltip='File with a list of tilt angles. (.tlt extension)',
                                    initdir=self.tomogram_folder)
        self.insert_label_line_push(parent, 'Defocus File', mode + 'DefocusFile',mode='file',filetype='defocus',
                                    tooltip='Filename to store results (.defocus extension).',
                                    initdir=self.tomogram_folder)
        self.insert_label_line_push(parent, 'Config File (OPTIONAL)', mode + 'ConfigFile', mode='file',filetype='cfg',
                                    tooltip='File with a list of noise files used to estimate the noise (.cfg extension).',
                                    initdir=self.tomogram_folder)

        self.insert_label_spinbox(parent, mode +'AxisAngle', text='Axis Angle (deg)',
                                  tooltip='Specifies how much the tilt axis deviates from vertical (Y axis).',
                                  value=270, minimum=0, maximum=359, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'PixelSpacing', text='Pixel Spacing (nm)',
                                  tooltip='Size of pixels in nm.', wtype=QDoubleSpinBox, decimals=3,
                                  value=0.175, minimum=.1, maximum=10, stepsize=.1)
        self.insert_label_spinbox(parent, mode + 'ExpectedDefocus', text='Expected Defocus (nm)',
                                  value=6000, minimum=-60000, maximum=60000, wtype=QDoubleSpinBox,
                                  stepsize=100, tooltip='Expected defocus at the tilt axis in nanometers.')
        self.insert_label_spinbox(parent, mode + 'AngleRangeMin', text='Angle Range(min)',
                                  value=-60, minimum=-180, maximum=-2, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'AngleRangeMax', text='Angle Range (max)',
                                  tooltip='Maximal angle of tilt images.',
                                  value=60, minimum=2, maximum=180, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'Voltage', text='Voltage (kV)',
                                  tooltip='Specifies how much the tilt axis deviates from vertical (Y axis).',
                                  value=200, minimum=2, maximum=1300, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'SphericalAberration', text='Spherical Aberration',
                                  value=2.7, stepsize=0.1, minimum=0, maximum=10., wtype=QDoubleSpinBox,
                                  tooltip='Expected defocus at the tilt axis in nanometers.')
        self.insert_label_spinbox(parent, mode + 'AmplitudeContrast', text='Amplitude Contrast',
                                  value=0.8, stepsize=0.1, minimum=0, maximum=1, wtype=QDoubleSpinBox)
        self.insert_label_checkbox(parent, mode+'FindAstigmatism', 'Find Astigmatism')
        self.insert_label_checkbox(parent, mode + 'FindPhaseShift', 'Find PhaseShift')
        self.insert_label_checkbox(parent, mode + 'FindCutOn', 'Find Cut On Frequency',rstep=1, cstep=0)

        exefilename = [mode+'tomofolder', 'ctf/CTFPlotter.sh']
        xmlfilename = [mode+'tomofolder', 'ctf/ctfplotter.com']
        outfolder = [mode+'tomofolder', 'ctf' ]

        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'CTF Plotting'
        paramsSbatch[ 'folder' ] = self.logfolder #outfolder

        paramsXML = [mode + 'StackFile', mode + 'AngleFile', mode + 'DefocusFile', mode + 'ConfigFile',
                     mode + 'AxisAngle', mode + 'PixelSpacing', mode + 'ExpectedDefocus', mode + 'AngleRangeMin',
                     mode + 'AngleRangeMax', mode + 'Voltage', mode + 'SphericalAberration', mode + 'AmplitudeContrast',
                     mode + 'FindAstigmatism', mode + 'FindPhaseShift', mode + 'FindCutOn', mode + 'ConfigHeader',
                     ParamsFileCTFPlotter]

        paramsCmd = [mode + 'tomofolder', 'ctf/ctfplotter.com', templateCTFPlotter]

        self.insert_gen_text_exe(parent, mode, jobfield=True, exefilename=exefilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, paramsXML=paramsXML, xmlfilename=xmlfilename,
                                 mandatory_fill=[mode + 'FolderSortedAligned', mode+'AngleFile', mode + 'DefocusFile'])

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        parent.addWidget(label)

        self.widgets[mode+'tomofolder'] = QLineEdit()
        self.widgets[mode+'StackFile'] = QLineEdit()
        self.widgets[mode+'ConfigHeader'] = QLineEdit()
        self.widgets[mode + 'FolderSortedAligned'].textChanged.connect(lambda dummy, m=mode: self.updateCTFPlotter(m))
        self.widgets[mode + 'ConfigFile'].textChanged.connect(lambda  dummy, m=mode: self.updateConfigHeader(m))
        self.updateCTFPlotter(mode,newstack=False)
        self.updateConfigHeader(mode)

        setattr(self, mode + 'gb_CD', groupbox)
        return groupbox

    def ctfCorrection(self, mode, title=''):
        tooltip = 'CTF Correction for aligned tilt images.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB, width=300)

        self.insert_label_line_push(parent, 'Defocus File', mode + 'DefocusFile', mode='file', filetype='defocus',
                                    tooltip='Filename with results from CTF Estimation (.defocus extension).',
                                    initdir=self.tomogram_folder)
        self.insert_label_line_push(parent, 'Stack with Sorted (& Aligned Tilt Images)', mode + 'InputStack',
                                    'Select the stack of sorted and possibly aligned tiltimages.\n',
                                    initdir=self.tomogram_folder, mode='folder')
        self.insert_label_line_push(parent, 'IMOD tlt File', mode + 'AngleFile', mode='file', filetype='tlt',
                                    tooltip='Filename with tilt angles (.tlt extension).',
                                    initdir=self.tomogram_folder)
        self.insert_label_line(parent, 'Output File', mode + 'OutputFile', value='ctfCorrected.st',
                               tooltip='Filename with results from CTF Correction.')
        self.insert_label_spinbox(parent, mode + 'DefocusTol', text='Defocus Tolerance (nm)',
                                  tooltip='',
                                  value=200, stepsize=10, minimum=1, maximum=1000)
        self.insert_label_spinbox(parent, mode + 'InterpolationWidth', text='Interpolation Width',
                                  tooltip='number of pixels of corrected strip which are interpolated',
                                  value=15, stepsize=1, minimum=1, maximum=1000)
        self.insert_label_spinbox(parent, mode + 'PixelSize', text='Pixel Size (nm)',
                                  tooltip='Physical size corresponding to one pixel in image (nm)', decimals=3,
                                  value=0.262, stepsize=.1, minimum=.01, maximum=1000, wtype=QDoubleSpinBox)
        self.insert_label_spinbox(parent, mode + 'SphericalAberration', text='Spherical Abberration',
                                  tooltip='Physical size corresponding to one pixel in image (nm)',
                                  value=2.7, stepsize=.1, minimum=.01, maximum=1000, wtype=QDoubleSpinBox)
        self.insert_label_spinbox(parent, mode + 'AmplitudeContrast', text='Amplitude Contrast ',
                                  tooltip='Amplittude Contrast (fraction between 0 and 1)',
                                  value=0.08, stepsize=.01, minimum=0, maximum=1, wtype=QDoubleSpinBox)
        self.insert_label_spinbox(parent, mode + 'Voltage', text='Voltage (kV)',
                                  tooltip='Voltage (kV)', value=200, stepsize=.1, minimum=10, maximum=1000)
        self.insert_label_spinbox(parent, mode + 'AxisAngle', text='In-plane Rotation Angle of tilt axis',
                                  tooltip='In-plane rotation angle of tilt axis after alignment.',
                                  value=0, stepsize=1, minimum=0, maximum=359)
        self.insert_label_line(parent, 'Index of the GPU', mode + 'gpuID', cstep=0, validator=QIntValidator(),
                               tooltip="This box determine the GPU you want to use. Empty field means no gpu.\nMultiple GPU's is not possible.", )

        self.widgets[mode + 'ctffolder'] = QLineEdit()
        self.widgets[mode + 'OutputFolder'] = QLineEdit()
        self.widgets[mode + 'prefix'] = QLineEdit()
        self.widgets[mode + 'origdir'] = QLineEdit()
        self.widgets[mode + 'gpuString'] = QLineEdit()

        self.widgets[mode + 'DefocusFile'].textChanged.connect(lambda d, m=mode: self.updateCTFCorrectionImod(m))
        exefilename = [mode + 'ctffolder', 'CTFCorrectionImod.sh']
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m, 1))

        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'CTF_Correction'
        paramsSbatch['folder'] = self.logfolder  # [mode + 'tomofolder', 'ctf']
        paramsSbatch['id'] = 'SingleCTFCorrection'
        paramsSbatch['gpus'] = self.widgets[mode + 'gpuID']

        paramsCmd = [mode + 'ctffolder', mode + 'InputStack', mode + 'OutputFile', mode + 'AngleFile',
                     mode + 'DefocusFile', mode + 'DefocusTol', mode + 'InterpolationWidth', mode + 'PixelSize',
                     mode + 'SphericalAberration', mode + 'AmplitudeContrast', mode + 'Voltage', mode + 'AxisAngle',
                     mode + 'OutputFolder', mode + 'prefix', mode + 'origdir', mode + 'gpuString',
                     templateCTFCorrectionImod]

        self.insert_gen_text_exe(parent, mode, exefilename=exefilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, paramsAction=[mode],
                                 mandatory_fill=[mode+'DefocusFile', mode+'AngleFile', mode+'OutputFile', mode+'InputStack'])

        # self.updateGpuString(mode)
        self.updateCTFCorrectionImod(mode)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        parent.addWidget(label)

        setattr(self, mode + 'gb_CD', groupbox)
        return groupbox

    def create_tomogram_folders(self, id, values):
        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables['tab1'].table
        widgets = self.tables['tab1'].widgets
        print('Creating dirs')
        jobs = []

        num_folders_to_be_deleted = 0

        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row, 6)
            #wname2= 'widget_{}_{}'.format(row, 7)
            #for wname in (wname1, wname2):
            if wname in widgets.keys() and widgets[wname].isChecked():
                reply = QMessageBox.question(self, 'Delete Folders', 'Are you sure you want to permanantly delete the selected folders?',
                        QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if reply == QMessageBox.No:
                    self.tab1UI(id)
                    return

        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row, 1)

            # Create
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = widgets['widget_{}_{}'.format(row, 3)].text()
                metafile = values[row][0]
                crop = widgets['widget_{}_{}'.format(row, 2)].isChecked()
                folder = self.filepath_tomodata[widgets['widget_{}_{}'.format(row, 4)].currentText()]
                binning = float(widgets[f'widget_{row}_{7}'].text())
                print(binning)
                jobs.append([self.create_tomodir_instance, (tomofoldername, metafile, folder, crop, binning)])

                #self.create_tomodir_instance(tomofoldername,metafile,folder)

            # Redo
            # wname = 'widget_{}_{}'.format(row, 6)
            # if wname in widgets.keys() and widgets[wname].isChecked():
            #     tomofoldername = widgets['widget_{}_{}'.format(row, 5)].text()
            #     metafile = os.path.basename(values[row][0])
            #     folder = self.filepath_tomodata[widgets['widget_{}_{}'.format(row, 4)].currentText()]
            #     crop = widgets['widget_{}_{}'.format(row, 2)].isChecked()
            #
            #     error = False
            #     metaf = glob.glob(os.path.join(self.rawnanographs_folder, '*', metafile))
            #     if metaf:
            #         metafile = metaf[0]
            #     elif os.path.join(self.rawnanographs_folder, metafile):
            #         metafile = os.path.join(self.rawnanographs_folder, metafile)
            #     else:
            #         error = True
            #
            #     if os.path.exists(metafile) and not error:
            #         sfolder = os.listdir(os.path.join(self.tomogram_folder, tomofoldername, 'sorted'))
            #         sfolder = [os.path.join(self.tomogram_folder, tomofoldername, 'sorted', f) for f in sfolder]
            #         try:
            #             [os.system('unlink {}'.format(f)) for f in sfolder if f.endswith('.mrc') and f.startswith('sorted_')]
            #         except:
            #             pass
            #         if os.path.exists(os.path.join(self.tomogram_folder, tomofoldername)):
            #             shutil.rmtree(os.path.join(self.tomogram_folder, tomofoldername))
            #
            #         jobs.append( [self.create_tomodir_instance, (tomofoldername, metafile, folder, crop)])
            #
            #         num_folders_to_be_deleted += 1
            #     else:
            #         print(f'{metafile} not found, skipping redo')

            # Delete
            wname = 'widget_{}_{}'.format(row, 6)
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = widgets['widget_{}_{}'.format(row, 5)].text()
                sfolder = os.listdir(os.path.join(self.tomogram_folder, tomofoldername, 'sorted'))
                sfolder = [os.path.join(self.tomogram_folder, tomofoldername, 'sorted', f) for f in sfolder]
                [os.system('unlink {}'.format(f)) for f in sfolder if f.endswith('.mrc') and f.startswith('sorted_')]
                if os.path.exists(os.path.join(self.tomogram_folder, tomofoldername)):
                    shutil.rmtree( os.path.join(self.tomogram_folder, tomofoldername) )

                num_folders_to_be_deleted += 1

                # self.create_tomodir_instance(tomofoldername,mdocfile,folder)

        events = []
        self.jobs_nr_create_tomofolders = len(jobs)
        print(self.jobs_nr_create_tomofolders)

        if len(jobs):
            self.divide_jobs(jobs, id)
        else:
            self.tab1UI(id)

    def divide_jobs(self, jobs, id):


        try:
            del self.counters
        except:
            pass

        counters = [0,]*len(jobs)

        manager = Manager()

        self.counters = manager.list(counters)

        self.generate_statusbar(len(jobs))

        procs = []
        self.completion = [0,]*self.num_parallel_procs
        for i in range(min(len(jobs), self.num_parallel_procs)):
            proc = Process(target=self.run_jobs, args=(jobs[i::self.num_parallel_procs], i, self.counters) )
            procs.append(proc)
            proc.start()
            atexit.register(guiFunctions.kill_proc, proc)

            time.sleep(1)


        self.rerunid = id
        self.aa = len(jobs)
        proc = Worker(fn=self.check_run, args=(len(jobs), id))
        proc.signals.result1.connect(self.update_progress_generate_tomogramdir)
        proc.signals.finished_mcor.connect(self.delete_progressbar_tomogramdir)
        proc.start()

    def check_run(self, a, id, signals):
        import time
        num = 0
        while sum(self.counters) < a and num < 1000:
            total = sum(self.counters)
            signals.result1.emit(total)
            time.sleep(1)
            num += 1

        signals.finished_mcor.emit()

    def delete_progressbar_tomogramdir(self):
        self.statusBar.removeWidget(self.progressBar)
        if self.aa: self.popup_messagebox("Info", "Completion", 'Successfully generated tomogram directories.')
        self.tab1UI(self.rerunid)

    def run_jobs(self, jobs, procid, counters):

        for n, (target, args) in enumerate(jobs):
            tomofoldername, metafile, folder, crop, binning = args
            target(tomofoldername, metafile, folder, crop, binning, procid, counters)

    def generate_statusbar(self, nrJobs):
        widget = QWidget(self)
        layout = QHBoxLayout()
        widget.setLayout(layout)

        self.progressBar = QProgressBar()

        self.progressBar.setSizePolicy(self.sizePolicyA)
        self.progressBar.setStyleSheet(DEFAULT_STYLE_PROGRESSBAR)
        self.progressBar.setFormat('Prepare Tomograms: %v/%m')
        self.progressBar.setSizePolicy(self.sizePolicyB)
        self.progressBar.setMaximum(nrJobs)
        layout.addWidget(self.progressBar)
        self.statusBar.addPermanentWidget(self.progressBar, stretch=1)

    def preprocdata(self, dst_mcor, crop):
        from pytom.agnostic.io import read as readNPY, write
        data = readNPY(dst_mcor)
        if crop:
            size = min(data.shape[:2])
            print(size, data.shape)
            data = data[:size, :size]
        write(dst_mcor, remove_hot_pixels(data))

    def create_tomodir_instance(self, tomofoldername, metafile, folder, crop, binning=4, procid=0, counters=[]):
        src = os.path.join(self.tomogram_folder, '.tomoname')
        num_subprocesses = 5
        dst = os.path.join(self.tomogram_folder, tomofoldername)
        if os.path.exists(dst):
            print('{} exists: QUIT'.format(dst))
            return
        shutil.copytree(src, dst)
        meta_dst = os.path.join(dst, 'sorted', os.path.basename(metafile))
        os.system('cp {} {}'.format(metafile, meta_dst) )

        metadata = loadstar(metafile, dtype=guiFunctions.datatype)

        tif_files  = metadata['FileName']
        tiltangles = metadata['TiltAngle']

        print(tif_files, folder)

        if len(list(tif_files)) == 0:
            return

        for n in range(len(tif_files)):
            tif_files[n] = [line.replace('.tif', '.mrc') for line in repr(tif_files[n]).split('\\') if '.' in line][-1]
            print( tif_files[n])

        if '/IMODSTACK_' in tif_files[0]:
            import mrcfile
            from pytom.bin.mrcs2mrc import extract_single_image

            imodstack = True
            stack = metafile.replace('.meta', '.st')
            fid = metafile.replace('.meta', '.fid')
            if os.path.exists(metafile.replace('.meta', '.prexg')):
                prexg = metafile.replace(".meta",".prexg")
                infile = metafile.replace(".meta",".fid")
                fid = infile[:-4]+'_raw.fid'
                os.system(f'xfmodel -back -scale {1/binning} -prealign {prexg} -input {infile} -o {fid}')


            # cmd = "newstack --InputFile {} --OutputFile {} --ModeToOutput 0 --FloatDensities 2 --BinByFactor 4 --ImagesAreBinned 1.0"
            # stackout = os.path.join(os.path.dirname(os.path.dirname(meta_dst)), 'imod', os.path.basename(stack))
            # cmd.format(stack, stackout)
            #
            # if os.path.exists(metafile.replace('.meta', '.prexg')):
            #     cmd += f' --TransformFile {metafile.replace(".meta",".prexg") }'
            # os.system(cmd)

            imoddata = mrcfile.open(stack, permissive=True).data.copy()

        else:
            imodstack = False

        name_angle = list(zip(tif_files, tiltangles))

        sort(name_angle, 1)
        metadata = numpy.sort(metadata , order='TiltAngle')

        procs = []
        num_copied = 0
        for n, (tif, angle) in enumerate(name_angle):
            while len(procs) >= num_subprocesses:
                time.sleep(0.5)
                procs = [proc for proc in procs if proc.is_alive()]

            # Ugly resetting of origin folder as we are now organising files in 01_Raw_Nanographs by import index to avoid cluttering one folder with 1000s of files,
            # should happen earlier when creating joblist

            if folder.endswith('01_Raw_Nanographs'):
                folder = os.path.dirname(metafile)

            src_mcor = os.path.join(folder, tif[1:-1])
            dst_mcor = os.path.join(os.path.dirname(meta_dst), 'sorted_{:02d}.mrc'.format(n))

            if os.path.exists(src_mcor):
                # print('test', src_mcor)
                num_copied += 1
                os.system('cp {} {}'.format(src_mcor, dst_mcor))
                #proc = Process(target=square_mrc, args=([dst_mcor]))
                #procs.append(proc)
                #proc.start()
                self.preprocdata(dst_mcor, crop)


            elif imodstack:
                sid = int(tif.split('_')[-1].split('.')[0])
                extract_single_image(imoddata[sid,:,:].T, n, '', '', os.path.dirname(meta_dst), 'sorted_')
                out = self.preprocdata(dst_mcor, crop)
                metadata['FileName'][n] = os.path.join(os.path.dirname(meta_dst), f'sorted_{n:02d}.mrc')
                num_copied += 1

        if imodstack:
            from pytom.gui.guiFunctions import fmt, headerText, savestar
            from pytom.reconstruction.tiltAlignmentFunctions import getIMODpreshifts
            from pytom.gui.guiFunctions import read_markerfile, write_markerfile
            # Save adjusted metadata (new file names pointing to extracted images)
            savestar(meta_dst, metadata, fmt=fmt, header=headerText)


            wimp = meta_dst.replace('.meta', '.wimp')
            prexg = metafile.replace('.meta', '.prexg')

            if os.path.exists(fid):
                os.system(f'convertmod {fid} {wimp}')

                markerFileName = os.path.join(os.path.dirname(meta_dst), 'markerfile.txt')
                markerfile = read_markerfile(wimp, tiltangles)
                if 0 and os.path.exists(prexg):
                    print('applying imods preshifts')
                    shiftx, shifty = getIMODpreshifts(prexg)
                    for i in range(markerfile.shape[1]):
                        markerfile[:,i,0] -= shiftx
                        markerfile[:,i,1] -= shifty
                print('write markerfile')

                out = markerfile.copy()
                cntr = 0

                for i in range(markerfile.shape[1]):
                    xx = markerfile[:, i, 0].sum()
                    yy = markerfile[:, i, 1].sum()
                    if yy > 0 and xx > 0:
                        out[:,cntr,:] = markerfile[:,i,:]
                        cntr += 1

                if cntr > 0:
                    markerfile = out[:,:cntr,:]

                write_markerfile(markerFileName, markerfile, tiltangles, bin_factor=binning)


        if num_copied < 1:
            shutil.rmtree(dst)

        counters[procid] += 1

    def updateSpecimenAngle(self, mode):
        angle = self.widgets[mode + 'SpecimenAngle'].text()
        if angle != '':
            res = '' if float(angle) == 0 else f'--specimenAngle {angle} '
        else:
            res = ''
        self.widgets[mode + 'SpecimenAngleFlag'].setText(res)

    def update_progress_generate_tomogramdir(self, total):
        self.progressBar.setValue(total)

    def updateTomoFolder(self, mode):

        folderSorted = self.widgets[mode+'FolderSorted'].text()
        if not folderSorted: return
        # basename() gets the second element of split()
        prefix = os.path.basename(folderSorted)
        self.widgets[mode + 'tiltSeriesName'].setText(f'{prefix}/{prefix}')
        t = folderSorted.replace('/sorted', '')  # remove sorted
        t = t.split('/alignment')[0]
        self.widgets[mode + 'tomofolder'].setText(t)
        self.widgets[mode + 'tomogramNR'].setText(os.path.basename(t))
        self.updateMarkerFile(folderSorted, mode)

        files = [line for line in os.listdir(folderSorted) if line.startswith('sorted') and line.endswith('.mrc')]
        lastIndex = len(files)
        #oself.widgets[mode+'LastIndex'].setText(lastIndex)
        metafiles = [line for line in os.listdir(folderSorted) if line.endswith('.meta')]

        if len(metafiles) == 1:
            metafile = metafiles[0]

            try:
                metadata = loadstar(os.path.join(folderSorted, metafile), dtype=guiFunctions.datatype)
            except:
                metadata = loadstar(os.path.join(folderSorted,metafile),dtype=guiFunctions.datatype0)

            angles = metadata['TiltAngle'].copy()
            self.widgets['RotationTiltAxis'] = metadata
            for i in range(len(files)):
                if not 'sorted_{:02d}.mrc'.format(i) in files:
                    angles[i]+=10000

            refIndex = abs(angles).argmin()
            self.widgets[mode+'RefTiltIndex'].setValue(refIndex+1*('INFR' in mode))

            fi,li = angles.argmin(), angles[angles < 1000].argmax()
            self.widgets[mode+'FirstIndex'].setText(str(fi))
            self.widgets[mode + 'LastIndex'].setText(str(li))

            self.updateIndex(mode)

        if 'WBP' in mode: self.updateVoldims(mode)

        self.updateIndex(mode)

    def updateMarkerFile(self, foldername, mode):
        markerfile = 'alignment/markerfile.txt' if os.path.exists(f'{foldername}/markerfile.txt') else ''
        if markerfile == '' and os.path.exists(f'{foldername}/markerfile.em'): markerfile = 'alignment/markerfile.em'

        if markerfile:
            os.system(f'cp {foldername}/{os.path.basename(markerfile)} {os.path.dirname(foldername)}/{markerfile}')

        self.widgets[mode + 'markerfile'].setText(markerfile)

    def updateVoldims(self, mode, rotation_axis=None):
        folderSorted = self.widgets[mode+'FolderSorted'].text()
        if not folderSorted: return
        files = [line for line in os.listdir(folderSorted) if line.startswith('sorted') and line.endswith('.mrc')]
        if not files: return
        files = files[0]
        imdims = read_size(os.path.join(folderSorted, files))

        dimx = str(int(float(imdims[0])/float(self.widgets[mode+'BinningFactor'].text())+.5))
        dimy = str(int(float(imdims[1])/float(self.widgets[mode+'BinningFactor'].text())+.5))

        if self.widgets[mode + 'RotationTiltAxis'].text() != '':
            expectedRotation = float(self.widgets[mode+'RotationTiltAxis'].text()) if \
                (rotation_axis is None) else rotation_axis
        else:
            expectedRotation = 0

        if abs(90 - (expectedRotation % 180)) < 45:
            dx, dy, dz = dimy, dimx, dimy
        else:
            dx, dy, dz = dimx, dimy, dimx

        self.widgets[mode+'Voldims'].setText(f'{dx}')
        self.widgets[mode + 'DimX'].setText(f'{dx}')
        self.widgets[mode + 'DimY'].setText(f'{dy}')
        self.widgets[mode + 'DimZ'].setText(f'{dz}')

    def updateIndex(self, mode):
        if 'INFR' in mode: INFR=1
        else: INFR=0
        sorted_folder = self.widgets[mode + 'FolderSorted'].text()
        if sorted_folder == '':
            return
        metafiles = [os.path.join(sorted_folder, meta) for meta in os.listdir(sorted_folder) if meta.endswith('meta')]
        #print(metafile)
        if len(metafiles) == 1:
            metafile = metafiles[0]
            metadata = loadstar(metafile,dtype=guiFunctions.datatype)

            try:
                firstAngle = self.widgets[mode + 'FirstAngle'].value()
                lastAngle = self.widgets[mode + 'LastAngle'].value()
            except:
                firstAngle = float(self.widgets[mode + 'FirstAngle'].text())
                lastAngle = float(self.widgets[mode + 'LastAngle'].text())
            files = sorted([f for f in os.listdir(sorted_folder) if f.startswith('sorted') and f.endswith('mrc')])
            first,last = 0, len(files)-1
            self.widgets[mode + 'FirstIndex'].setText( str(0))
            self.widgets[mode + 'LastIndex'].setText( str( len(files)-1) )
            for n, tiltAngle in enumerate( metadata['TiltAngle']):
                mrcfile = 'sorted_{:02d}.mrc'.format(n)
                if tiltAngle - firstAngle < -0.1 and mrcfile in files:
                    self.widgets[mode + 'FirstIndex'].setText(str(n+1+INFR))

                if tiltAngle - lastAngle < 0.1 and mrcfile in files:
                    self.widgets[mode + 'LastIndex'].setText(str(n+INFR))

            fi, li = int(self.widgets[mode+'FirstIndex'].text()), int(self.widgets[mode+'LastIndex'].text())

            if fi > 0 or li-INFR < len(files)-1:
                self.widgets[mode + 'Reduced'].setText('_reduced_{},{}'.format(firstAngle, lastAngle))
            else:
                self.widgets[mode + 'Reduced'].setText('')

            self.widgets[mode + 'Reduced'].setText('_{:.1f},{:.1f}'.format(firstAngle, lastAngle))

            try:
                refid = int(self.widgets[mode + 'RefMarkerIndex'].text())
                angs = self.widgets[mode + 'Reduced'].text()
                uudir = 'marker_{:04d}{}'.format(refid, angs)
                alignment_file = os.path.join('alignment', uudir, 'GlobalAlignment', 'sorted',
                                                                      'alignmentResults.txt')
                self.widgets[mode + 'AlignFile'].setText(alignment_file)
                self.widgets[mode + 'RotationTiltAxis'].setText(str(axis_angle_from_ar_file(
                    os.path.join(self.widgets[mode + 'tomofolder'].text(), alignment_file))))

                if 'WBP' in mode: self.updateVoldims(mode)

            except Exception as e:
                print(e)
                pass

            try:
                refid = int(self.widgets[mode + 'RefMarkerIndex'].text())
                angs = self.widgets[mode + 'Reduced'].text()
                uudir = 'marker_{:04d}{}'.format(refid, angs)
                self.widgets[mode + 'outFolder'].setText( os.path.join('', 'alignment', uudir, 'GlobalAlignment',
                                                                os.path.basename(sorted_folder), ) )
                #print(self.widgets[mode + 'outFolder'].text())

            except Exception as e:
                print(e)
                pass

    def updateCTFCorrectionImod(self, mode):
        defocusfile = self.widgets[mode + 'DefocusFile'].text()
        ctffolder = os.path.dirname(defocusfile)
        ctfplotter = os.path.join(ctffolder, 'ctfplotter.com')
        if not defocusfile:
            return

        try:
            self.widgets[mode + 'ctffolder'].setText(ctffolder)
            if os.path.exists(ctfplotter):
                d = {}
                keys = ('InputStack', 'AngleFile', 'DefocusFile', 'AxisAngle', 'PixelSize', 'Voltage',
                        'SphericalAberration', 'AmplitudeContrast', 'DefocusTol')
                for line in open(ctfplotter).readlines():
                    kvpair = line.split()
                    if len(kvpair) == 2 and kvpair[0] in keys:
                        k, v = kvpair
                        try:
                            d[k] = float(v)
                        except:
                            d[k] = v

                for k in d.keys():
                    if d[k].__class__ == str:
                        self.widgets[mode + k].setText(d[k])
                    else:
                        self.widgets[mode + k].setValue(d[k])
                origdir = os.path.dirname(d['InputStack'])
                prefix = os.path.basename(origdir) + '_ctf'
                outfolder = os.path.join(ctffolder, prefix)
                if not os.path.exists(outfolder):
                    os.mkdir(outfolder)
                self.widgets[mode + 'origdir'].setText(origdir)
                self.widgets[mode + 'prefix'].setText(prefix)
                self.widgets[mode + 'OutputFolder'].setText(outfolder)

        except Exception as e:
            print(e)

    def updateCTFCorrection(self, mode):
        folder = self.widgets[mode + 'FolderSortedAligned'].text()
        if not folder:
            print('CTF Correction: empty folder sorted aligned')
            return
        try:
            tomogramID = folder.split('tomogram_')[-1][:3]
            tomoname = '{}/tomogram_{}'.format(self.tomogram_folder, tomogramID)
            cOut = '{}/ctf/ctf_{}'.format(tomoname, os.path.basename(folder))
            if not os.path.exists(cOut): os.mkdir(cOut)
            cPrefix = os.path.join(cOut, 'sorted_aligned_ctf_')
            uPrefix = os.path.join(folder, 'sorted_aligned_')
            sortedFolder = os.path.join(tomoname, 'sorted')
            metafile = glob.glob(sortedFolder + '/*.meta')[0]
            # defocusfile = glob.glob( os.path.join( os.path.dirname(cOut), '*.defocus'))
        except:
            self.popup_messagebox('Error', 'Setting values for CTF Correction Failed',
                                  'The setting of the parameters for CTF correction has failed. \n' +
                                  'Running ctf correction might result in errors. ' +
                                  'Try to select Folder Sorted Aligned again.')
            return

        self.widgets[mode + 'tomofolder'].setText(tomoname)
        self.widgets[mode + 'uPrefix'].setText(uPrefix)
        self.widgets[mode + 'cPrefix'].setText(cPrefix)
        self.widgets[mode + 'MetaFile'].setText(metafile)

    def updateGridAndFieldSize(self, mode):
        gridSpacing = float(self.widgets[mode + 'GridSpacing'].value())
        fieldSize = float(self.widgets[mode + 'FieldSize'].value())
        if gridSpacing > fieldSize:
            self.widgets[mode + 'GridSpacing'].setValue(fieldSize)
        self.widgets[mode + 'GridSpacing'].setMaximum(fieldSize)

    def updateMetaFile(self, params):
        mode = params[0]
        metafile = self.widgets[mode + 'MetaFile'].text()
        defocusfile = self.widgets[mode + 'DefocusFile'].text()
        try:
            guiFunctions.update_metadata_from_defocusfile(metafile, defocusfile)
        except:
            self.popup_messagebox('Error', 'Update MetaData Failed',
                                  'Update metadata has failed. Your job is not using the paramaters from the selected defocus file.')

    def updateGpuString(self, mode):
        gpu_id = self.widgets[mode + 'gpuID'].text()
        if gpu_id == '':
            self.widgets[mode + 'gpuString'].setText('')
            return
        try:
            gpu_id = int(gpu_id)
        except ValueError:
            self.widgets[mode + 'gpuID'].setText('')
            self.widgets[mode + 'gpuString'].setText('')
            self.popup_messagebox('Warning', 'Invalid value in field',
                                  'Reconstruction can only use a single GPU')
            return

        self.widgets[mode + 'gpuString'].setText(f'--gpuID {gpu_id} ')

    def updateConfigHeader(self, mode):
        if self.widgets[mode + 'ConfigFile'].text():
            self.widgets[mode + 'ConfigHeader'].setText('ConfigFile')
        else:
            self.widgets[mode + 'ConfigHeader'].setText('          ')

    def updateCTFPlotter(self, mode, newstack=True):
        folder = self.widgets[mode + 'FolderSortedAligned'].text()
        if not folder:
            # print('No Folder Selected.')
            return
        try:
            tomogramID = folder.split('tomogram_')[-1][:3]
            sortedFolder = '{}/tomogram_{}/sorted/'.format(self.tomogram_folder, tomogramID)
            tomoname = '{}/tomogram_{}'.format(self.tomogram_folder, tomogramID)

            outstack = '{}/tomogram_{}_{}.st'.format(folder,tomogramID, os.path.basename(folder))
        except:
            print('update CTF Plotter failed.')
            return
        files = [line for line  in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted_')]
        if not files:
            print('No files wuth MRC-format found in: {}'.format(folder))
            return
        dd = []
        for file in files:
            id = int(file.split('_')[-1].split('.')[0])
            dd.append([file,id])
        guiFunctions.sort(dd,1)
        files, ids = zip(*dd)
        if newstack:
            infile, outfile = open(os.path.join(folder, 'filein.txt'), 'w'), open(os.path.join(folder, 'fileout.txt'),'w')
            cmd = 'cd {}; newstack -filei filein.txt -fileo fileout.txt '.format(folder)
            outfile.write('{}\n{}\n{}\n'.format(1, outstack, len(files)))
            infile.write('{}\n'.format(len(files)))
            for fname in files:
                infile.write('{}\n0\n'.format(fname))
            infile.close()
            outfile.close()
            os.system(cmd)

        metafile = [os.path.join(sortedFolder, line) for line in os.listdir(sortedFolder) if line.endswith('.meta')][0]
        metadata = loadstar(metafile,dtype=guiFunctions.datatype)
        outAngle   = outstack.replace('.st','.tlt')
        out = open(outAngle,'w')
        for id in ids: out.write('{}\n'.format(metadata['TiltAngle'][id]))
        out.close()
        outDefocus = os.path.join(tomoname, 'ctf', 'resultsCTFPlotter.defocus')

        self.widgets[mode + 'tomofolder'].setText(tomoname)
        self.widgets[mode + 'StackFile'].setText(outstack)
        self.widgets[mode + 'AngleFile'].setText(outAngle)
        self.widgets[mode + 'DefocusFile'].setText(outDefocus)

        self.widgets[mode + 'AxisAngle'].setValue(metadata['InPlaneRotation'][0])
        self.widgets[mode + 'PixelSpacing'].setValue(metadata['PixelSpacing'][0]/10.)
        self.widgets[mode + 'ExpectedDefocus'].setValue(1000.*metadata['DefocusU'][abs(metadata['TiltAngle']).argmin()])
        self.widgets[mode + 'AngleRangeMin'].setValue(numpy.floor(metadata['TiltAngle'][ids[0]]))
        self.widgets[mode + 'AngleRangeMax'].setValue(numpy.ceil(metadata['TiltAngle'][ids[-1]]))
        self.widgets[mode + 'Voltage'].setValue(metadata['Voltage'][0])
        self.widgets[mode + 'SphericalAberration'].setValue(metadata['SphericalAberration'][0])
        self.widgets[mode + 'AmplitudeContrast'].setValue(metadata['AmplitudeContrast'][0])

    def update_create_tomoname(self, id, columnId=0, rowId=0):
        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables['tab1'].table
        widgets = self.tables['tab1'].widgets
        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row,1)
            if wname in widgets.keys() and widgets[wname]:
                if widgets[wname].isChecked():
                    for index in range(1000):
                        tt = 'tomogram_{:03d}'.format(n)
                        if not os.path.exists(os.path.join(self.tomogram_folder, tt)):
                            break
                        n += 1
                    widgets['widget_{}_{}'.format(row, 3)].setText('tomogram_{:03d}'.format(n))
                    n+=1
                else:
                    widgets['widget_{}_{}'.format(row, 3)].setText('')

        table.resizeColumnsToContents()

    def update_tomogram_reconstruction(self, mode):
        # set tomofolder
        self.widgets[mode + 'tomofolder'].setText(os.path.join(self.tomogram_folder,
                                                         self.widgets[mode + 'tomogram'].currentText()))

        # remove previous alignment choices
        self.widgets[mode + 'alignment'].disconnect()
        self.widgets[mode + 'alignment'].clear()

        # find available alignments
        alignment_dir = os.path.join(self.tomogram_folder, self.widgets[mode + 'tomogram'].currentText(), 'alignment')
        alignment_choices = []
        for alignment in os.listdir(alignment_dir):
            d = os.path.join(alignment_dir, alignment)
            if os.path.exists(d) and os.path.isdir(d):
                for pointer in os.listdir(os.path.join(d, 'GlobalAlignment')):
                    f = os.path.join(d, 'GlobalAlignment', pointer, 'alignmentResults.txt')
                    if os.path.exists(f):
                        alignment_choices.append(alignment)
                        break
                # we might have found an alignment file

        # set the choices for alignment in the GUI
        self.widgets[mode + 'alignment'].currentIndexChanged.connect(
            lambda d, m=mode: self.update_alignment_choice(m))
        self.widgets[mode + 'alignment'].addItems(alignment_choices)

        # get ctf options
        self.widgets[mode + 'ctfCorrChoice'].disconnect()
        self.widgets[mode + 'ctfCorrChoice'].clear()
        sorted_dir = os.path.join(self.tomogram_folder, self.widgets[mode + 'tomogram'].currentText(), 'sorted')
        ctf_sorted_dir = os.path.join(self.tomogram_folder, self.widgets[mode + 'tomogram'].currentText(), 'ctf',
                                      'sorted_ctf')
        tilt_choices = []
        if len([f for f in os.listdir(sorted_dir) if f.endswith('.mrc')]) > 0:
            tilt_choices.append('sorted')
        if len([f for f in os.listdir(ctf_sorted_dir) if f.endswith('.mrc')]) > 0:  # at least
            tilt_choices.append('sorted_ctf')
        # set choices in the widget
        self.widgets[mode + 'ctfCorrChoice'].currentIndexChanged.connect(
            lambda d, m=mode: self.update_ctf_corr_choice(m))
        self.widgets[mode + 'ctfCorrChoice'].addItems(tilt_choices)

    def update_alignment_choice(self, mode):
        folder = os.path.join(self.tomogram_folder, self.widgets[mode + 'tomogram'].currentText(), 'alignment',
                              self.widgets[mode + 'alignment'].currentText(), 'GlobalAlignment')
        try:
            pointer = [d for d in os.listdir(folder) if 'sorted' in d][0]
            self.widgets[mode + 'alignmentResultsFile'].setText(os.path.join(folder, pointer, 'alignmentResults.txt'))

            # if an alignment file has been selected, we can start looking for angle range
            if os.path.exists(self.widgets[mode + 'alignmentResultsFile'].text()):
                try:  # TODO dit werkt niet als verwacht
                    alignment = loadstar(self.widgets[mode + 'alignmentResultsFile'].text(), dtype=ALIGNRESULTS_ORDER)
                except:
                    alignment = loadstar(self.widgets[mode + 'alignmentResultsFile'].text(), dtype=ALIGNRESULTS_OLD)
            tilt_angles = alignment['TiltAngle']
            min_angle = int(round(tilt_angles.min()))
            max_angle = int(round(tilt_angles.max()))
            # set bounds of first angle and set to min value
            self.widgets[mode + 'firstAngle'].setMinimum(min_angle)
            self.widgets[mode + 'firstAngle'].setMaximum(max_angle)
            self.widgets[mode + 'firstAngle'].setValue(min_angle)
            # set bounds of last angle and set to max value
            self.widgets[mode + 'lastAngle'].setMinimum(min_angle)
            self.widgets[mode + 'lastAngle'].setMaximum(max_angle)
            self.widgets[mode + 'lastAngle'].setValue(max_angle)

            # set rotation axis for tomogram dims
            self.widgets[mode + 'RotationTiltAxis'].setText(str(alignment['InPlaneRotation'][0]))
        except IndexError:
            print('No sorted or sorted_ctf folder in the alignment directory.')

    def update_alignment_choice_batch(self, row_id, table_id):
        w_tomogram = self.tables[table_id].widgets['widget_{}_0'.format(row_id)]
        w_alignment = self.tables[table_id].widgets['widget_{}_2'.format(row_id)]
        w_first_angle = self.tables[table_id].widgets['widget_{}_3'.format(row_id)]
        w_last_angle = self.tables[table_id].widgets['widget_{}_4'.format(row_id)]

        # set the alignment file but if sorted exists, if not get sorted_ctf
        ar_file = os.path.join(self.tomogram_folder, w_tomogram.text(), 'alignment',
                               w_alignment.currentText(), 'GlobalAlignment', 'sorted', 'alignmentResults.txt')
        if not os.path.exists(ar_file):
            ar_file = os.path.join(self.tomogram_folder, w_tomogram.text(), 'alignment',
                                   w_alignment.currentText(), 'GlobalAlignment', 'sorted_ctf', 'alignmentResults.txt')
        # read old and new type of alignment-results
        try:  # TODO same, dit werkt niet
            alignment = loadstar(ar_file, dtype=ALIGNRESULTS_ORDER)
        except:
            alignment = loadstar(ar_file, dtype=ALIGNRESULTS_OLD)
        # get min max tilt angles from alignment
        tilt_angles = alignment['TiltAngle']
        min_angle = int(round(tilt_angles.min()))
        max_angle = int(round(tilt_angles.max()))
        w_first_angle.setText(str(min_angle))
        w_last_angle.setText(str(max_angle))

    def update_ctf_corr_choice(self, mode):
        sorted_choice = self.widgets[mode + 'ctfCorrChoice'].currentText()
        if 'ctf' in sorted_choice:
            self.widgets[mode + 'projectionDir'].setText(os.path.join(self.tomogram_folder,
                                                                      self.widgets[mode + 'tomogram'].currentText(),
                                                                      'ctf', 'sorted_ctf'))
        else:
            self.widgets[mode + 'projectionDir'].setText(os.path.join(self.tomogram_folder,
                                                                      self.widgets[mode + 'tomogram'].currentText(),
                                                                      'sorted'))

    def update_vol_dims(self, mode):
        # check if files exist
        projection_dir = self.widgets[mode+'projectionDir'].text()
        if not projection_dir:
            return  # exit if not projection dir provided
        files = [line for line in os.listdir(projection_dir) if line.startswith('sorted') and line.endswith('.mrc')]
        if len(files) == 0:
            return  # exit if no files found

        # run calculation of volume dims
        files = files[0]
        imdims = read_size(os.path.join(projection_dir, files))

        dimx = str(int(float(imdims[0])/float(self.widgets[mode+'BinningFactor'].text())+.5))
        dimy = str(int(float(imdims[1])/float(self.widgets[mode+'BinningFactor'].text())+.5))

        if self.widgets[mode+'RotationTiltAxis'].text() != '':
            expected_rotation = float(self.widgets[mode+'RotationTiltAxis'].text())
        else:
            expected_rotation = 0

        if abs(90 - (expected_rotation % 180)) < 45:
            dx, dy, dz = dimy, dimx, dimy
        else:
            dx, dy, dz = dimx, dimy, dimx

        self.widgets[mode+'Voldims'].setText(f'{dx}')  # needed ??
        self.widgets[mode + 'DimX'].setText(f'{dx}')
        self.widgets[mode + 'DimY'].setText(f'{dy}')
        self.widgets[mode + 'DimZ'].setText(f'{dz}')

    def update_cores(self, mode):
        n_cores = self.widgets[mode + 'cores'].text()
        if n_cores == '':
            self.widgets[mode + 'cores_flag'].setText('')
            return
        try:
            n_cores = int(n_cores)
        except ValueError:
            self.widgets[mode + 'cores'].setText('')
            self.widgets[mode + 'cores_flag'].setText('')
            self.popup_messagebox('Warning', 'Invalid value in field',
                                  'Can only fill in integer number of cpu cores')
            return

        self.widgets[mode + 'cores_flag'].setText(f'--numProcesses {n_cores} ')

    def convert_em(self,params):
        mode = params[0]
        directory = self.widgets[mode+'tomofolder'].text()
        output_folder = '{}/{}'.format(directory, params[1])
        prefix = params[2]
        sorted_folder = self.widgets[mode+'FolderSorted'].text()
        alignmentType = os.path.basename(params[1])

        if not os.path.exists(f'{output_folder}/temp_files_unweighted') and alignmentType in ('INFR', 'WBP'):
            os.mkdir(f'{output_folder}/temp_files_unweighted')

        if alignmentType == 'INFR':
            if len([line for line in os.listdir(output_folder) if line.startswith(prefix.split('/')[-1])]):
                os.system('rm {}/sorted*.em'.format(output_folder))
            guiFunctions.conv_mrc2em(self.widgets[mode+'FolderSorted'].text(), output_folder)
            guiFunctions.renumber_gui2pytom(output_folder, prefix.split('/')[-1])
        #print (self.pytompath)
        #print ('{}/bin/pytom rename_renumber.py {} {} {}'.format(self.pytompath, sorted_folder, output_folder, prefix))


        # print '{}/reconstruction/{}/reconstruction.sh'.format(tomofolderpath,p)
        if 1:
            if os.path.exists('{}/{}_reconstruction.sh'.format(output_folder,os.path.basename(output_folder))):
                for i in range(1000):
                    if not os.path.exists('{}/backup'.format(output_folder)): os.mkdir(
                        '{}/backup'.format(output_folder))
                    fname = "{}/backup/reconstruction_{:03d}".format(output_folder, i)
                    if os.path.exists(fname): continue

                    os.mkdir(fname)
                    for f in ('temp_files_binned', 'temp_files_unweighted', 'temp_files_weighted', 'tomogram.em',
                              '{}_reconstruction.sh'.format(os.path.basename(output_folder)), 'markerfile.em'):
                        if os.path.exists('{}/{}'.format(output_folder, f)):
                            os.system('mv {}/{} {}/'.format(output_folder, f, fname))
                        if not '.' in f and not os.path.exists("{}/{}".format(output_folder, f)):
                            os.mkdir("{}/{}".format(output_folder, f))
                    break
        else:
            pass

        if os.path.exists(f'{directory}/sorted/markerfile.txt'):
            os.system(f'cp {directory}/sorted/markerfile.txt {output_folder}/markerfile.txt')
        elif os.path.exists(f'{directory}/sorted/markerfile.em'):
            guiFunctions.convert_markerfile(f'{directory}/sorted/markerfile.em', f'{directory}/sorted/markerfile.txt')
            os.system(f'cp {directory}/sorted/markerfile.txt {output_folder}/markerfile.txt')
        else:
            pass

        self.activate_stage(2)

    def run_multi_align(self, id, values):
        """
        Multi tilt alignment will create jobs with a specified number of cores (n). Each job will then give the
        alignment of one tomogram to a single core. If there are more tomograms then cores, the jobs will go over n
        tomogram each.

        @param id: tab id (tab32)
        @type id: str
        @param values: fillable values from the tab as a list
        @type values: list
        @return:
        @rtype: None
        """
        print('multi_align', id)
        table = self.tables[id].table
        widgets = self.tables[id].widgets

        for i in range(10000):
            file_tomoname = os.path.join(self.tomogram_folder, 'jobscripts/.multi_alignment_{:04d}.txt'.format(i))
            if not os.path.exists(file_tomoname):
                break

        tomofolder_info = []
        total_number_markers = 0
        number_tomonames = 0
        firstindices, lastindices, expectedangles = [], [], []
        firstangles, lastangles = [], []
        mode = 'v02_ba_'
        for name in ('FirstAngle', 'LastAngle', 'FirstIndex', 'LastIndex', 'Reduced', 'FolderSorted'):
            self.widgets[mode + name] = QLineEdit()
        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row, 1)

            # Align
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = values[row][0]
                firstindex = widgets['widget_{}_{}'.format(row, 2)].text()
                lastindex  = widgets['widget_{}_{}'.format(row, 3)].text()
                refindex   = widgets['widget_{}_{}'.format(row, 4)].text()
                markindex  = widgets['widget_{}_{}'.format(row, 5)].currentText()
                expected   = widgets['widget_{}_{}'.format(row, 6)].text()
                inputfolder= values[row][7][widgets['widget_{}_{}'.format(row, 7)].currentIndex()]
                fixmarkers = widgets['widget_{}_{}'.format(row, 8)].isChecked()
                refmarkIDTomo = widgets['widget_{}_{}'.format(row, 9)].currentText()
                tiltseriesname = os.path.join(inputfolder, os.path.basename(inputfolder))
                self.widgets[mode + 'FolderSorted'].setText(os.path.join(self.tomogram_folder,
                                                                         tomofoldername, 'sorted'))
                number_tomonames += 1
                folder = os.path.join(self.tomogram_folder, tomofoldername)
                os.system('cp {} {}/alignment/'.format(self.mfiles[row],folder))

                self.widgets[mode + 'FirstAngle'].setText(firstindex)
                self.widgets[mode + 'LastAngle'].setText(lastindex)
                self.updateIndex(mode)
                fi, li = self.widgets[mode + 'FirstIndex'].text(), self.widgets[mode + 'LastIndex'].text()
                fa, la = self.widgets[mode + 'FirstAngle'].text(), self.widgets[mode + 'LastAngle'].text()
                firstindices.append(fi)
                lastindices.append(li)
                firstangles.append(fa)
                lastangles.append(la)
                expectedangles.append(expected)

                markerfile = '{}/alignment/markerfile.{}'.format(folder,self.mfiles[row].split('.')[-1])
                l = len(glob.glob(os.path.join(os.path.dirname(self.mfiles[row]), 'sorted_*.mrc')))

                markerdata = guiFunctions.readMarkerfile(markerfile, l)

                if markindex == 'all':
                    numMark = markerdata.shape[2]
                else:
                    numMark = 1

                total_number_markers += numMark

                if not refmarkIDTomo: refmarkIDTomo = '*'

                ll = '{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'
                ll = ll.format(tomofoldername, refindex, numMark, markindex, fa, la, fi, li, 0,
                               expected, tiltseriesname, markerfile, fixmarkers, refmarkIDTomo)
                tomofolder_info.append([numMark, ll])

        if not tomofolder_info:
            return

        if self.checkbox[id].isChecked():
            _, _, cores, _, _, _ = self.qparams['BatchAlignment'].values()
        else:
            cores = cpu_count()

        new_list = sorted(tomofolder_info, key=lambda l: l[0], reverse=True)

        new_info = []
        lprocs = [0]
        taken = [0, ] * number_tomonames
        while number_tomonames - sum(taken):
            take = []
            partial_sum = 0
            for n in range(len(new_list)):
                if taken[n]: continue
                # this piece assumes 20 cores, could easily be updated with a number of cores
                if not partial_sum or (partial_sum + new_list[n][0]) < cores + 1:
                    partial_sum += new_list[n][0] % cores
                    if not partial_sum: partial_sum += cores
                    take.append(n)
            for t in take:
                taken[t] = 1
                new_info.append(new_list[t])
            lprocs.append(sum(taken))

        tomofolder_file = open(file_tomoname, 'w')
        for x,y in new_info:
            tomofolder_file.write(y)
        tomofolder_file.close()

        num_submitted_jobs = 0

        submissionIDs = []

        for n in range(len(lprocs) - 1):

            input_params = (self.tomogram_folder, self.pytompath, lprocs[n], lprocs[n + 1],
                            tiltseriesname, 'alignment', file_tomoname)

            cmd = multiple_alignment.format(d=input_params)
            exefilename = '{}/jobscripts/alignment_{:03d}.sh'.format(self.tomogram_folder, n)

            if self.checkbox[id].isChecked():
                qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchAlignment'].values()

                jobname = 'Alignment_BatchMode_Job_{:03d}'.format(num_submitted_jobs)

                cmd = guiFunctions.gen_queue_header(name=jobname, folder=self.logfolder, partition=qname, time=time,
                                                    num_nodes=n_nodes, cmd=qcmd, modules=modules,
                                                    num_jobs_per_node=cores) + cmd

            ID, num = self.submitBatchJob(exefilename, id, cmd)
            num_submitted_jobs += 1
            submissionIDs.append(ID)

        if num_submitted_jobs > 0:
            self.popup_messagebox('Info', 'Submission Status', f'Submitted {num_submitted_jobs} jobs to the queue.')
            self.addProgressBarToStatusBar(submissionIDs, key='QJobs', job_description='Alignment Batch')

    def run_multi_ctf_correction(self, id, values):
        print('multi_ctf_corrections', id)
        num_nodes = self.tables[id].table.rowCount()
        num_submitted_jobs = 0
        submissionIDs = []
        try:
            num_nodes = self.num_nodes[id].value()
        except:
            pass
        for row in range(self.tables[id].table.rowCount()):
            widget = 'widget_{}_1'.format(row)
            if self.tab43_widgets[widget].isChecked():
                ctffolder     = os.path.join(values[row][0], 'ctf')
                tomofolder    = self.tab43_widgets['widget_{}_{}'.format(row, 0)].text()
                gpu           = self.tab43_widgets['widget_{}_{}'.format(row, 2)].text()
                defocustol    = self.tab43_widgets['widget_{}_{}'.format(row, 3)].text()
                interpolation = self.tab43_widgets['widget_{}_{}'.format(row, 4)].text()
                pixelsize     = self.tab43_widgets['widget_{}_{}'.format(row, 5)].text()
                cs            = self.tab43_widgets['widget_{}_{}'.format(row, 6)].text()
                ac            = self.tab43_widgets['widget_{}_{}'.format(row, 7)].text()
                voltage       = self.tab43_widgets['widget_{}_{}'.format(row, 8)].text()
                axisangle     = self.tab43_widgets['widget_{}_{}'.format(row, 9)].text()
                inputstack    = self.tab43_widgets['widget_{}_{}'.format(row,10)].text()
                anglefile     = self.tab43_widgets['widget_{}_{}'.format(row,11)].text()
                defocusfile   = self.tab43_widgets['widget_{}_{}'.format(row,12)].text()
                outputfolder  = self.tab43_widgets['widget_{}_{}'.format(row,13)].text()

                try:
                    origdir = os.path.dirname(inputstack)
                    prefix = os.path.basename(origdir) + '_ctf'
                    outfolder = os.path.join(ctffolder, prefix)
                    if not os.path.exists(outfolder):
                        os.mkdir(outfolder)

                    outputfile = 'ctfCorrected.st'

                    metafile = glob.glob(f'{os.path.dirname(inputstack)}/*.meta')
                    if len(metafile) > 0: metaFlag = f' --metaFile {metafile[0]}'
                    else: metaFlag = ''

                    if gpu:
                        gpuFlag = f' -gpu {int(gpu)} '
                    else:
                        gpuFlag = ''

                    jobParams = [ctffolder, inputstack, outputfile, anglefile, defocusfile, defocustol, interpolation,
                                 pixelsize, cs, ac, voltage, axisangle, outputfolder, prefix, origdir, gpuFlag, metaFlag]

                    jobscript = templateCTFCorrectionImod.format(d=jobParams)

                    fname = 'CTF_Batch_ID_{}'.format(num_submitted_jobs % num_nodes)
                    suffix = f"_{tomofolder}_"

                    qname,n_nodes,cores,time, modules, qcmd = self.qparams['BatchCTFCorrection'].values()


                    job = guiFunctions.gen_queue_header(folder=self.logfolder, cmd = qcmd, name=fname,
                                                        suffix=suffix, time=time,
                                                        partition=qname, num_nodes=n_nodes, singleton=True,
                                                        num_jobs_per_node=cores, modules=modules, gpus=gpu)*self.checkbox[id].isChecked() + jobscript

                    exefilename = os.path.join(outputfolder, 'ctfCorrectionBatch.sh')

                    ID, num = self.submitBatchJob(exefilename, id, job)
                    num_submitted_jobs += 1
                    submissionIDs.append(ID)

                except Exception as e:
                    print(e)
                    print('Failed correction for ' + defocusfile)
                    continue

        if num_submitted_jobs > 0:
            self.popup_messagebox('Info', 'Submission Status', f'Submitted {num_submitted_jobs} jobs to the queue.')
            self.addProgressBarToStatusBar(submissionIDs, key='QJobs', job_description='CTF Correction Batch')

    def run_multi_reconstruction(self, id):
        print('multi_reconstructions', id)

        table = self.tables[id].table
        widgets = self.tables[id].widgets
        mode = 'batch_recon_'
        dd = {1: ('INFR', 'reconstruction/INFR'), 2: ('WBP', 'reconstruction/WBP')}

        # initialize widgets for calculating tomogam dims
        for name in ('Voldims', 'DimX', 'DimY', 'DimZ', 'RotationTiltAxis', 'FolderSorted',
                     'BinningFactor', 'gpu_flag'):
            if mode + name in self.widgets.keys():
                continue
            self.widgets[mode + name] = QLineEdit()

        # initialize for job submission
        nsj = 0
        jobCode = {}
        execfilenames = {}

        try:
            for row in range(table.rowCount()):
                tomogram = widgets['widget_{}_0'.format(row)].text()
                tomofolder = os.path.join(self.tomogram_folder, tomogram)

                # TODO increment by 1 when INFR is added back in
                alignment       = widgets['widget_{}_{}'.format(row, 2)].currentText()
                first_angle     = int(widgets['widget_{}_{}'.format(row, 3)].text())
                last_angle      = int(widgets['widget_{}_{}'.format(row, 4)].text())
                projection_dir  = widgets['widget_{}_{}'.format(row, 5)].currentText()
                weighting       = int(widgets['widget_{}_{}'.format(row, 6)].text())
                binning         = int(widgets['widget_{}_{}'.format(row, 7)].text())
                specimen_angle  = ''  # if activate => float(widgets['widget_{}_{}'.format(row, 8)].text())
                gpu_id          = ''  # set as empty as INFR cannot use gpus
                self.update_gpu_flag(widgets['widget_{}_{}'.format(row, 8)],
                                     self.widgets[mode + 'gpu_flag'], single_gpu=True, suppress_message=True)

                sorted_folder = os.path.join(tomofolder, 'ctf', 'sorted_ctf') if 'ctf' in projection_dir else \
                    os.path.join(tomofolder, 'sorted')

                for dd_key, (method, method_folder) in dd.items():  # TODO add INFR back in by looping over 1, 2
                    # ======= skip filling for INFR
                    if method == 'INFR':
                        widget = 'widget_{}_{}'.format(row, dd_key)
                        continue
                    else:
                        widget = 'widget_{}_{}'.format(row, 1)  # should be dd_key if INFR is added back
                    # ========

                    if widgets[widget].isChecked():

                        alignment_folder = os.path.join(tomofolder, 'alignment', alignment, 'GlobalAlignment')
                        pointer = [d for d in os.listdir(alignment_folder) if 'sorted' in d][0]
                        ar_file = os.path.join(alignment_folder, pointer, 'alignmentResults.txt')
                        rotation_tilt_axis = axis_angle_from_ar_file(ar_file)
                        self.widgets[mode + 'RotationTiltAxis'].setText(str(rotation_tilt_axis))
                        self.widgets[mode + 'FolderSorted'].setText(sorted_folder)
                        self.widgets[mode + 'BinningFactor'].setText(str(binning))
                        # determine dims based on rotation axis (same as imod)
                        self.updateVoldims(mode)
                        dx, dy, dz = (int(self.widgets[mode + 'DimX'].text()),
                                      int(self.widgets[mode + 'DimY'].text()),
                                      int(self.widgets[mode + 'DimZ'].text()))

                        qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchReconstruct'].values()
                        execfilename = os.path.join(tomofolder, '{}/{}_Reconstruction.sh'.format(method_folder,
                                                                                                 method))

                        if method == 'INFR':
                            paramsCmd = [tomofolder, self.pytompath, first_angle, last_angle,
                                         binning, sorted_folder, ar_file, specimen_angle, os.path.basename(tomofolder)]
                            txt = templateINFR.format(d=paramsCmd)

                        elif method == 'WBP':
                            specimen_angle_flag = '' if specimen_angle == 0 or specimen_angle == '' else \
                                f'--specimenAngle {specimen_angle} '
                            gpu_flag = self.widgets[mode + 'gpu_flag'].text()
                            gpu_id = widgets['widget_{}_{}'.format(row, 8)].text()
                            cores_flag = '' if gpu_flag != '' else f'--numProcesses {cores} '

                            paramsCmd = [tomofolder, self.pytompath, binning, tomogram, dx, weighting,
                                         rotation_tilt_axis, specimen_angle_flag, dy, dz, ar_file, gpu_flag,
                                         sorted_folder, first_angle, last_angle, cores_flag]

                            txt = templateWBP.format(d=paramsCmd)

                        else:
                            print('No Batch Submission')
                            continue

                        # what about gpus in slurm???
                        job_name = 'TomoRecon_{}'.format(nsj % n_nodes)  # ??
                        cmd = txt

                        if gpu_id != '' and method == 'WBP':  # force that only WBP can use GPUs (might be changed in
                            # future if INFR is ported to GPU)
                            job_code_key = method + '_' + str(gpu_id)

                            if not job_code_key in jobCode.keys():
                                if self.checkbox[id].isChecked():
                                    job = guiFunctions.gen_queue_header(folder=self.logfolder, name=job_name,
                                                                        time=time, num_nodes=n_nodes, partition=qname,
                                                                        modules=modules,
                                                                        num_jobs_per_node=cores, gpus=gpu_id,
                                                                        cmd=qcmd) + cmd
                                else:
                                    job = cmd
                                execfilenames[job_code_key] = execfilename
                                jobCode[job_code_key] = job
                                nsj += 1

                            elif job_code_key in jobCode.keys():
                                jobCode[job_code_key] += f'\nwait\n\n{cmd}\n'

                        else:  # set up CPU jobs
                            job_code_key = method + '_' + f'noGPU_{nsj % n_nodes}'
                            if not job_code_key in jobCode.keys():
                                if self.checkbox[id].isChecked():
                                    job = guiFunctions.gen_queue_header(folder=self.logfolder, name=job_name,
                                                                        time=time, num_nodes=n_nodes, partition=qname,
                                                                        modules=modules, num_jobs_per_node=cores,
                                                                        cmd=qcmd) + cmd
                                else:
                                    job = cmd

                                execfilenames[job_code_key] = execfilename
                                jobCode[job_code_key] = job
                                nsj += 1

                            else:
                                jobCode[job_code_key] += f'\nwait\n\n{cmd}\n'

            wid = []
            todoList = {}

            for key in jobCode.keys():
                if not key in todoList:
                    todoList[key] = []
                todoList[key].append([execfilenames[key], id, jobCode[key]])

            from time import sleep
            for key in todoList.keys():
                print(f'starting {len(todoList[key])} jobs on device {key}')
                self.localJobs[self.workerID] = []
                wid.append(self.workerID)
                proc = Worker(fn=self.multiSeq, args=((self.submitBatchJob, todoList[key], self.workerID)))
                self.threadPool.start(proc)
                self.workerID += 1
                sleep(.01)

            if nsj > 0:
                self.popup_messagebox('Info', 'Submission Status', f'Submitted {nsj} jobs to the queue.')
                self.addProgressBarToStatusBar(wid, key='QJobs', job_description='Tom. Reconstr. Batch')

        except ValueError:
            self.popup_messagebox('Warning', 'Invalid value in field',
                                  'One of filled in values could not be parsed as an integer/float.')

    def submit_multi_recon_job(self, params):

        try:
            exefilename = params[0]
            exefile = open(exefilename, 'w')
            exefile.write(params[1])
            exefile.close()

            if len(params[1].split('SBATCH')) > 2:

                dd = os.popen('{} {}'.format(self.qcommand, exefilename))
                text = dd.read()[:-1]
                id = text.split()[-1]
                logcopy = os.path.join(self.projectname, f'LogFiles/{id}_{os.path.basename(exefilename)}')
                os.system(f'cp {exefilename} {logcopy}')


            else:
                os.system('sh {}'.format(params[0]))
        except:
            print ('Please check your input parameters. They might be incomplete.')

    def prep_value(self, params):

        mode = params[0]
        outtext = '''# command file to run ctfplotter\n#\n'''

        folder = self.widgets[mode +'FolderSortedAligned'].text()
        metafile = [line for line in os.listdir('{}/../../sorted/'.format(folder)) if line.endswith('.meta')][0]
        metadata = loadstar(metafile,dtype=guiFunctions.datatype)

        files = [os.path.join(folder,line) for line in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted_aligned')]

        files = sorted(files)

