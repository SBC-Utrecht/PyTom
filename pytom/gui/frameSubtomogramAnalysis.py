import os
import glob

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# pytom gui functions import
from pytom.gui.guiSupportCommands import extractParticles, extractParticlesClosestMarker, polishParticles, \
    templateFRMJob, templateAverageParticleList, templateCCC, templateAC, templateFSC, templateGLocal, templateCPCA, \
    templateFRMSlurm
from pytom.gui.guiStructures import GuiTabWidget, SelectFiles, Worker, CreateFSCMaskFile, CreateMaskFile
import pytom.gui.guiFunctions as guiFunctions
from pytom.bin.extractTomoNameFromXML import extractParticleListsByTomoNameFromXML
from pytom.basic.structures import ParticleList
from pytom.agnostic.io import read_size
from pytom.basic.files import read, loadstar
from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO as ALIGNRESULTS_ORDER, DATATYPE_ALIGNMENT_RESULTS as \
    ALIGNRESULTS_OLD


class SubtomoAnalysis(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(SubtomoAnalysis, self).__init__(parent)
        self.stage          = 'v04_'
        self.addGeneralVariables()

        # CHANGE THE NEXT FOUR VARIABLES WHEN ADDING OR REMOVING TABS
        # ONE NEEDS 1) UI FUNCTION TO SET UP THE GENERAL TAB (NUMBER OF DROPDOWN FIELDS VS DYNAMIC TABLE)
        #           2) FUNCTION TO FILL DROP DOWN MENU / FILL THE TABLE
        #           3) HELPER FUNCTIONS TO AUTOFILL OR UPDATE FIELDS BASED ON USER INPUT

        headers = ["Reconstruct Subtomograms", "Particle Polishing", "Average", "Align Subtomograms", "Classify Subtomograms", "Validation"]
        subheaders = [['Single Reconstruction','Batch Reconstruction'],['Single', 'Batch'], [], ['FRM Alignment','GLocal'],['CPCA','Auto Focus'], [] ]
        tabUIs = [[self.SubtomoReconstrSingleUI, self.SubtomoReconstrBatchUI],
                  [self.PolishSingleUI, self.PolishBatchUI],
                  [self.AveragePL],
                  [self.FRMUI,self.GLocalUI],
                  [self.CPCAUI,self.AC3DUI],
                  [self.FSCUI]]

        self.names = [[],[],[],[],[],[],[]]

        static_tabs = [[True, False], [True, False], [True],[True, True], [True, True], [True, True]]




        self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders,tabUIs=tabUIs,tabs=self.tabs_dict,
                     tab_actions=self.tab_actions, static_tabs=static_tabs)



    # General UI functions

    def dd(self, key):
        pass

    def SubtomoReconstrSingleUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleSubtomoReconstruct_'

        items += list(self.create_expandable_group(self.addSubtomogramReconstructionFields, self.sizePolicyB,
                                                   'Single Reconstruction', mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def SubtomoReconstrBatchUI(self, key=''):
        try: self.extractLists.text()
        except: self.extractLists = QLineEdit()

        self.mass_extract = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.xml'],
                                        outputline=self.extractLists, id=key,
                                        run_upon_complete=self.populateSubtomoReconBatchTable,
                                        title='Select particlLists')

    def PolishSingleUI(self, key=''):

        return

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleParticlePolish_'

        items += list(self.create_expandable_group(self.addParticlePolishFields, self.sizePolicyB, 'Single Polishing',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def PolishBatchUI(self, key=''):

        return

        try: self.polishLists.text()
        except: self.polishLists = QLineEdit()

        self.selectPLWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.xml'],
                                          outputline=self.polishLists, id=key,
                                          run_upon_complete=self.selectTemplates,
                                          title='Select particlLists')

    def AveragePL(self,key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'AverageParticleList_'

        items += list(self.create_expandable_group(self.addAverageFields, self.sizePolicyB, 'Calculate Average from Particle List',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def FRMUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 =  self.stage + 'FRMAlignment_'

        items += list(self.create_expandable_group(self.addFRMAlignmentFields, self.sizePolicyB, 'FRM Alignment',
                                                   mode=t0))
        items[-1].setVisible(False)


        for n, item in enumerate(items):
            grid.addWidget(item, n, 0,1,4)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def GLocalUI(self, key=''):


        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        items += list(self.create_expandable_group(self.addGLocalAlignmentFields, self.sizePolicyB, 'GLocal Alignment',
                                                   mode=self.stage + 'GLocalAlignment_'))
        items[-1].setVisible(False)
        # items += list( self.create_expandable_group(self.createSubtomograms, self.sizePolicyB, 'Extract Subtomograms',
        #                                            mode=self.stage+'extract_') )
        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def CPCAUI(self, key=''):

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'PairwiseCrossCorrelation_', self.stage + 'CPCA_'

        items += list(self.create_expandable_group(self.addCCCFields, self.sizePolicyB, 'Pairwise Cross Correlation',
                                                   mode=t0))
        items[-1].setVisible(False)

        items += list(self.create_expandable_group(self.addCPCAFields, self.sizePolicyB, 'CPCA',
                                                   mode=t1))
        items[-1].setVisible(False)
        # items += list(self.create_expandable_group(self.sampleInformation, self.sizePolicyB, 'Sample Information',
        #                                           mode=t2))
        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def AC3DUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'AutoFocusClassification_'

        items += list(self.create_expandable_group(self.addAC3DFields, self.sizePolicyB, 'Auto Focus Classification',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def FSCUI(self,key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'FourierShellCorrelation_'

        items += list(self.create_expandable_group(self.addFSCFields, self.sizePolicyB, 'Fourier Shell Correlation',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)


    # Methods to fill tabs with input fields or an input table

    def addSubtomogramReconstructionFields(self, mode='', title=''):

        #title = "Single Reconstruction"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.tomogram = None
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particlelist',initdir=self.pickpartdir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_combobox(parent, 'Alignment', mode + 'alignment', labels=[],
                                   tooltip='Select the alignment for subtomo reconstruction.')
        self.insert_label_combobox(parent, 'Ctf corrected', mode + 'ctfCorrChoice', labels=[],
                                   tooltip='Select original (sorted) or ctf corrected images (sorted_ctf).')
        self.insert_label_spinbox(parent, mode + 'firstAngle', 'First Angle',
                                  'Select the first angle of the angular range used for reconstruction.',
                                  minimum=-90, maximum=90, stepsize=1, value=0)
        self.insert_label_spinbox(parent, mode + 'lastAngle', 'Last Angle',
                                  'Select the last angle of the angular range used for reconstruction.',
                                  minimum=-90, maximum=90, stepsize=1, value=0)
        self.insert_label_spinbox(parent, mode + 'BinFactorReconstruction',
                                  'Binning factor used in the reconstruction.',
                                  'Defines the binning factor used in the reconstruction of the tomogram from which' +
                                  'the particles are selected.',
                                  minimum=1, stepsize=1, value=8)
        self.insert_label_spinbox(parent, mode + 'WeightingFactor', 'Apply Weighting (-1,0,1)',
                                  'Sets the weighting scheme applied to the tilt images, -1 is recommended. \n(-1 = '
                                  'ramp, 0 = no weighting, 1 = exact weighting)',
                                  minimum=-1, maximum=1, stepsize=1, value=-1)
        self.insert_label_spinbox(parent, mode + 'SizeSubtomos', 'Size subtomograms.',
                                  'Sets the size of the subtomograms.',
                                  minimum=10, maximum=1000, stepsize=1, value=128)
        self.insert_label_spinbox(parent, mode + 'BinFactorSubtomos', 'Binning Factor Subtomograms.',
                                  'Sets the binning factor of the subtomograms.', rstep=1,
                                  value=1, stepsize=1, minimum=1)
        self.insert_label_spinbox(parent, mode + 'OffsetX', 'Offset in x-dimension',
                                  'Has the tomogram been cropped in the x-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the x-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetY', 'Offset in y-dimension',
                                  'Has the tomogram been cropped in the y-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the y-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetZ', 'Offset in z-dimension',
                                  'Has the tomogram been cropped in the z-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the z-dimension.',
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)
        # self.insert_label_line_push(parent, 'Result file particle polishing', mode + 'polishFile', mode='file',
        #                             filetype='txt', initdir=self.polishfolder, # cstep=-1,
        #                             tooltip='Select a resultfile from particle polishing (Optional).')
        self.insert_label_line(parent, 'cores', mode + 'cores', validator=QIntValidator(),
                               tooltip="Integer number of cpu cores to use if running on cpu" )
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0, validator=QIntValidator(),
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs "
                                       "separate them using a comma, e.g. 0,1,2 ")

        # TODO have taken out polishing as this is currently not enabled in GUI
        # self.widgets[mode + 'polishFlag'] = QLineEdit('')
        self.widgets[mode + 'projectionDir'] = QLineEdit('')
        self.widgets[mode + 'alignmentResultsFile'] = QLineEdit('')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('')
        self.widgets[mode + 'cores_flag'] = QLineEdit('')

        # self.widgets[mode + 'particlelist'].textChanged.connect(lambda d, m=mode: self.updateMeta(m))
        # self.widgets[mode + 'polishFile'].textChanged.connect(lambda d, m=mode: self.updatePolishFlag(m))
        self.widgets[mode + 'particlelist'].textChanged.connect(lambda d, m=mode: self.update_alignment_options(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, ids=self.widgets[mode + 'gpuID'],
                                                                flag=self.widgets[mode+'gpuString']:
                                                         self.update_gpu_flag(ids, flag))
        self.widgets[mode + 'cores'].textChanged.connect(lambda d, cores=self.widgets[mode + 'cores'],
                                                         flag=self.widgets[mode + 'cores_flag']:
                                                         self.update_cores(cores, flag))

        execfilename = os.path.join(self.subtomodir, 'Reconstruction/reconstructSubtomograms.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='subtomoReconstr', folder=self.logfolder,
                                                      id='SingleSubtomoReconstruct')
        paramsCmd = [mode + 'particlelist', mode + 'alignmentResultsFile', mode + 'BinFactorReconstruction',
                     mode + 'SizeSubtomos', mode + 'BinFactorSubtomos', mode + 'OffsetX', mode + 'OffsetY',
                     mode + 'OffsetZ', self.subtomodir, mode + 'WeightingFactor',
                     mode + 'cores_flag', mode + 'projectionDir', mode + 'firstAngle', mode + 'lastAngle',
                     mode + 'gpuString', extractParticles]
        # mode + projectionDir  # mode + 'polishFlag'
        mandatory_fill = [mode + 'particlelist', mode + 'alignmentResultsFile']

        self.update_gpu_flag(self.widgets[mode + 'gpuID'], self.widgets[mode + 'gpuString'])

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 mandatory_fill=mandatory_fill)
        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def populateSubtomoReconBatchTable(self, key='tab12'):
        self.mass_extract.close()
        particle_files_start = sorted(self.extractLists.text().split('\n'))
        particle_files = []
        tomograms = []

        for particle_file in particle_files_start:
            pl = ParticleList()
            pl.fromXMLFile(particle_file)
            pls = pl.splitByProjectDir()
            # set a flag if there was a project split
            proj_split = True if len(pls) > 1 else False

            for pp in pls:
                if len(pp) == 0: continue

                ppts = pp.splitByTomoName()

                # if there is only one project and only one tomo,
                # then we dont need to save any split lists
                if len(ppts) == 1 and not proj_split:
                    particle_files.append(particle_file)
                    continue

                projdir = os.path.basename(os.path.normpath(pp[0].getInfoGUI().getProjectDir()))
                projdir = os.path.basename(os.path.normpath(self.projectname)) if projdir in ['', '.'] else projdir
                dirname = os.path.splitext(os.path.basename(particle_file))[0]

                # create a dir to store the split output lists
                output_folder = os.path.join(self.pickpartdir, dirname)
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)

                # save a list per tomogram (reconstruction algo works per on a per tomogram basis)
                # + get the tomogram from the particle list file
                for ppt in ppts:
                    picking_origin = os.path.basename(ppt[0].getPickPosition().getOriginFilename())
                    if 'WBP' in picking_origin or 'INFR' in picking_origin:
                        tomogram = picking_origin.split('_WBP')[0] if 'WBP' in picking_origin else picking_origin.split(
                            '_INFR')[0]
                    else:
                        tomogram = os.path.splitext(picking_origin)[0]
                    tomograms.append(tomogram)
                    ppt.toXMLFile(f'{output_folder}/particleList_{projdir}_{tomogram}.xml')
                    particle_files.append(f'{output_folder}/particleList_{projdir}_{tomogram}.xml')

        # organize based on tomogram name
        particle_files = [x for _, x in sorted(zip(tomograms, particle_files), key=lambda pair: pair[0])]
        tomograms = sorted(tomograms)

        headers = ['particle_list', 'tilt-series', 'run', 'alignment', 'first angle', 'last angle', 'ctf',
                   'bin factor tomo', 'weighting', 'subtomo size', 'subtomo bin', 'offset x', 'offset y', 'offset z',
                   'scale factor', 'gpu ids', 'polish result', '']
        types = ['txt', 'txt', 'checkbox', 'combobox', 'lineedit', 'lineedit', 'combobox', 'lineedit', 'lineedit',
                 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'comboboxF', '']
        tooltip = ['Particle list file', 'tilt-series for which particles will be reconstructed.',
                   'Check this box to run subtomogram reconstruction.',
                   'Alignment type selection.',
                   'Negative limit tilt angles.', 'Positive limit tilt angles.',
                   'Origin folder of projections, either sorted or ctf corrected.',
                   'Binning factor in the tomographic reconstruction where the particles were located.',
                   'Weighting Type.\n0: No Weighting,\n1: Analytical Weighting.\n-1: Ramp Weighting',
                   'Size of the subtomograms after binning.', 'Binning factor for subtomograms.',
                   'Offset in X-dimension', 'Offset in Y-dimension', 'Offset in Z-dimension',
                   'Scaling factor to compensate for magnification differences between datasets',
                   'One or more GPU IDs (separated by a comma, e.g. 0,4,5 )',
                   'Results files from particle polishing.']
        sizes = [0, ] * (len(headers) - 1)

        values = []
        project_tomograms = sorted([d for d in os.listdir(self.tomogram_folder) if (not 'jobscripts' in d and
                                                                                    not d.startswith('.'))])
        skip_pop_up_shown = False

        for n, (tomogram, particle_file) in enumerate(zip(tomograms, particle_files)):

            tomo_folder = os.path.join(self.tomogram_folder, tomogram)

            # check if the tomogram is available in this pytom project, otherwise skip and put pop up message
            if not tomogram in project_tomograms or not os.path.exists(tomo_folder):
                if not skip_pop_up_shown:
                    self.popup_messagebox('Warning', 'Invalid particles found',
                                          'Some of the provided particles were picked in tilt-series that are not in '
                                          'the current project.')
                    skip_pop_up_shown = True
                continue

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

            # get tomogram binning from particle list
            particles = ParticleList()
            particles.fromXMLFile(particle_file)
            tomogram_binning = int(particles[0].xpath('PickPosition')[0].get('Binning'))

            polishfiles = glob.glob(f'{self.polishfolder}/resultsPolish*{tomogram}*.txt')
            polishfiles += glob.glob(f'{self.polishfolder}/*/resultsPolish*{tomogram}*.txt')
            polishfiles += ['']

            values.append([particle_file, tomogram, True, alignment_choices, 0, 0, tilt_choices, tomogram_binning,
                           -1, 128, 1, 0, 0, 0, 1, '', polishfiles, ''])

        try:
            self.num_nodes[key].setParent(None)
        except:
            pass

        if values:
            self.fill_tab(key, headers, types, values, sizes, tooltip=tooltip, nn=True,
                          wname=self.stage + 'BatchSubtomoReconstruct_')
            self.tab12_widgets = self.tables[key].widgets
            self.tab12_widgets['gpu_flag'] = QLineEdit('')

            for row in range(len(values)):
                w = self.tab12_widgets['widget_{}_3'.format(row)]  # alignment choice widget
                w.currentIndexChanged.connect(lambda d, r=row, k=key: self.update_alignment_choice_batch(r, k))
                self.update_alignment_choice_batch(row, key)

                # gpu flag is not used but just to check if the input is valid
                w = self.tab12_widgets['widget_{}_15'.format(row)]
                w.textChanged.connect(lambda d, ids=w, flag=self.tab12_widgets['gpu_flag']:
                                        self.update_gpu_flag(ids, flag, suppress_message=True))

            self.pbs[key].clicked.connect(lambda dummy, pid=key, v=values: self.massExtractParticles(pid, v))
        else:
            return

    def massExtractParticles(self, pid, values):
        num_nodes = int(self.num_nodes[pid].value())
        jobCode = {}
        execfilenames = {}
        nsj = 0

        try:
            for row in range(self.tables[pid].table.rowCount()):
                if self.tab12_widgets['widget_{}_2'.format(row)].isChecked():
                    particleXML = values[row][0]
                    tomogram = values[row][1]
                    alignment = self.tab12_widgets['widget_{}_{}'.format(row, 3)].currentText()
                    # folder_align_type = self.tab12_widgets['widget_{}_{}'.format(row, 3)].currentText()
                    min_angle = int(self.tab12_widgets['widget_{}_{}'.format(row, 4)].text())
                    max_angle = int(self.tab12_widgets['widget_{}_{}'.format(row, 5)].text())
                    folder_origin = self.tab12_widgets['widget_{}_{}'.format(row, 6)].currentText()

                    # get the alignment results file
                    al_dir = os.path.join(self.tomogram_folder, tomogram, 'alignment', alignment, 'GlobalAlignment')
                    point = [o for o in os.listdir(al_dir) if 'sorted' in o][0]
                    ar_file = os.path.join(al_dir, point, 'alignmentResults.txt')
                    projection_dir = os.path.join(self.tomogram_folder, tomogram,
                                                  folder_origin if not 'ctf' in folder_origin else 'ctf/sorted_ctf')

                    bin_read = int(self.tab12_widgets['widget_{}_{}'.format(row, 7)].text())
                    weight = int(self.tab12_widgets['widget_{}_{}'.format(row, 8)].text())
                    size = int(self.tab12_widgets['widget_{}_{}'.format(row, 9)].text())
                    bin_subtomo = int(self.tab12_widgets['widget_{}_{}'.format(row, 10)].text())
                    offx = int(self.tab12_widgets['widget_{}_{}'.format(row, 11)].text())
                    offy = int(self.tab12_widgets['widget_{}_{}'.format(row, 12)].text())
                    offz = int(self.tab12_widgets['widget_{}_{}'.format(row, 13)].text())

                    # set gpu flag
                    self.update_gpu_flag(self.tab12_widgets['widget_{}_{}'.format(row, 15)],
                                         self.tab12_widgets['gpu_flag'])
                    gpu_ids_flag = self.tab12_widgets['gpu_flag'].text()
                    gpu_ids = self.tab12_widgets['widget_{}_{}'.format(row, 15)].text()
                    # set scale flag
                    scale = self.tab12_widgets['widget_{}_{}'.format(row, 14)].text()
                    scale = None if (scale == '1' or scale == '') else float(scale)
                    scale_flag = '' if scale is None else f'--scaleFactorParticle {scale} '
                    # add polish file
                    pp_file = self.tab12_widgets['widget_{}_{}'.format(row, 16)].currentText()
                    pp_flag = f'--particlePolishResultFile {ppfile} ' if pp_file else ''

                    # no closest marker reconstruction
                    if not 'closest' in alignment.lower():

                        outname = f'Reconstruction/reconstruct_subtomograms_{tomogram}.sh'
                        execfilename = os.path.join(self.subtomodir, outname)
                        qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchSubtomoReconstruct'].values()

                        paramsCmd = [particleXML, ar_file, bin_read, size, bin_subtomo, offx, offy, offz,
                                     self.subtomodir, weight, f'--numProcesses {cores} ', projection_dir,
                                     min_angle, max_angle, scale_flag, pp_flag, gpu_ids_flag]

                        txt = extractParticles.format(d=paramsCmd)

                    else:
                        # TODO remove this, closest marker reconstruction will no longer be supported
                        outname = 'Reconstruction/reconstruct_subtomograms_{:03d}_refmarker_{}.sh'
                        # outname = outname.format(int(tomoindex), refid)
                        # execfilename = os.path.join(self.subtomodir, outname)
                        #
                        # folder_origin = os.path.join(markerPath, folder_align_type, folder_origin)
                        #
                        #
                        # tomodir = os.path.dirname(os.path.dirname(metafile))
                        #
                        # reconAlg = None
                        # for alg in ('WBP', 'INFR'):
                        #     if alg in particleXML:
                        #         reconAlg = alg
                        #
                        # if not reconAlg:
                        #     print('FAIL: subtomogram reconstruction for {} failed. No INFR or WBP in xml path to particle.')
                        #     continue
                        #
                        # # Find the most recent markerLocation file
                        # end = 'reconstruction/{}/markerLocations*_irefmark_*.txt'
                        # end = end.format(reconAlg, tomoindex, reconAlg)
                        # logfilequery = os.path.join(tomodir, end)
                        # logfiles = glob.glob(logfilequery)
                        # logfiles.sort(key=os.path.getmtime)
                        # logfile = logfiles[-1]
                        #
                        # qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchSubtomoReconstruct'].values()
                        #
                        # paramsCmd = [particleXML, folder_origin, bin_read, size, bin_subtomo, offx, offy, offz,
                        #              self.subtomodir, weight, metafile, logfile, str(cores*n_nodes), 'sorted_aligned']
                        #
                        # txt = extractParticlesClosestMarker.format(d=paramsCmd)

                    job_name = 'SubtomoRecon_{}'.format(nsj % num_nodes)
                    timed = time
                    cmd = txt

                    if gpu_ids and not gpu_ids in jobCode.keys():
                        if self.checkbox[pid].isChecked():
                            job = guiFunctions.gen_queue_header(folder=self.logfolder, name=job_name,
                                                                time=timed, num_nodes=n_nodes, partition=qname,
                                                                modules=modules, num_jobs_per_node=cores, gpus=gpu_ids,
                                                                cmd=qcmd) + cmd
                        else:
                            job = cmd

                        execfilenames[gpu_ids] = execfilename
                        jobCode[gpu_ids] = job
                        nsj += 1

                    elif gpu_ids and gpu_ids in jobCode.keys():
                        jobCode[gpu_ids] += f'\nwait\n\n{cmd}\n'

                    elif not f'noGPU_{nsj % num_nodes}' in jobCode.keys():
                        if self.checkbox[pid].isChecked():
                            job = guiFunctions.gen_queue_header(folder=self.logfolder, name=job_name,
                                                                time=timed, num_nodes=n_nodes, partition=qname,
                                                                modules=modules, num_jobs_per_node=cores,
                                                                cmd=qcmd) + cmd
                        else:
                            job = cmd

                        execfilenames[f'noGPU_{nsj % num_nodes}'] = execfilename
                        jobCode[f'noGPU_{nsj % num_nodes}'] = job
                        nsj += 1

                    else:
                        jobCode[f'noGPU_{nsj % num_nodes}'] += f'\nwait\n\n{cmd}\n'

            wid = []
            todoList = {}

            for key in jobCode.keys():
                if not key in todoList:
                    todoList[key] = []
                todoList[key].append([execfilenames[key], pid, jobCode[key]])

            from time import sleep
            for key in todoList.keys():
                print(f'starting {len(todoList[key])} jobs on device {key}')
                self.localJobs[self.workerID] = []
                wid.append(self.workerID)
                proc = Worker(fn=self.multiSeq, args=((self.submitBatchJob, todoList[key], self.workerID)))
                self.threadPool.start(proc)
                self.workerID += 1
                sleep(.01)

            if nsj:
                self.popup_messagebox('Info', 'Submission Status', f'Submitted {nsj} jobs to the queue.')
                self.addProgressBarToStatusBar(wid, key='QJobs', job_description='Subtom Recon Batch')

        except ValueError:
            self.popup_messagebox('Warning', 'Invalid value in field',
                                  'One of filled in values could not be parsed as an integer/float.')

    def addParticlePolishFields(self, mode='', title=''):
        # title = "Particle Polish (Single)"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Folder with aligned tilt images', mode + 'AlignedTiltDir',
                                    'Select the folder with the aligned tilt images.')
        self.insert_label_line_push(parent, 'Template File', mode + 'Template',mode='file',filetype=['em', 'mrc'],
                                    tooltip='Select a template file.')
        self.insert_label_line_push(parent, 'FSC File', mode + 'FSCPath', mode='file', filetype=['dat'],
                                    tooltip='Select a FSC file (Optional).')
        self.insert_label_line_push(parent, 'Meta File', mode + 'MetaFile', mode='file', filetype=['meta'],
                                    tooltip='Select the corresponding meta file (in tomogram_???/sorted)')
        self.insert_label_line_push(parent, 'Alignment Results File', mode + 'AlignmentResultsFile', mode='file', filetype=['txt'],
                                    tooltip='Select the alignment results file (in alignment/???/???/)')
        self.insert_label_line_push(parent, 'Output Directory', mode + 'destination', mode='folder',
                                    tooltip='Select the destination directory.')

        self.insert_label_spinbox(parent,mode+'dimZ', 'Size tomogram in Z-dimension.',
                                  'Sets the size of z-dimension of the origin tomogram.',
                                  minimum=10,maximum=4000,stepsize=1,value=200)
        self.insert_label_spinbox(parent, mode + 'BinFactorReconstruction', 'Binning factor used in the reconstruction.',
                                  'Defines the binning factor used in the reconstruction of the tomogram from which'+
                                  'the particles are selected.',
                                  minimum=1,stepsize=1,value=8)
        self.insert_label_spinbox(parent, mode+'maxParticleShift', 'Max particle shift.',
                                  'Sets the max of the particle shift.',rstep=1,
                                  value=25, stepsize=1, minimum=1)
        self.insert_label_spinbox(parent,  mode + 'OffsetX', 'Offset in x-dimension',
                                  'Has the tomogram been cropped in the x-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the x-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1,minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetY', 'Offset in y-dimension',
                                  'Has the tomogram been cropped in the y-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the y-dimension.', cstep=-1,rstep=1,
                                  value=0, stepsize=1,minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetZ', 'Offset in z-dimension',
                                  'Has the tomogram been cropped in the z-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the z-dimension.', cstep=0, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'FSCFlag'] = QLineEdit('')
        self.widgets[mode + 'MetaFileFlag'] = QLineEdit('')
        self.widgets[mode + 'AlignmentResultsFileFlag'] = QLineEdit('')


        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updateMeta(m))
        self.widgets[mode + 'FSCPath'].textChanged.connect(lambda d, m=mode: self.updateFSCFlag(m))
        self.widgets[mode + 'MetaFile'].textChanged.connect(lambda d, m=mode: self.updateMetaFileFlag(m))
        self.widgets[mode + 'AlignmentResultsFile'].textChanged.connect(lambda d, m=mode: self.updateAlignmentResultsFlag(m))

        execfilename = os.path.join(self.subtomodir, 'ParticlePolishing/reconstructSubtomograms.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='particlePolishSingle', folder=self.logfolder,
                                                      id='SingleParticlePolish')
        paramsCmd = [self.subtomodir, mode + 'numberMpiCores', self.pytompath, mode + 'particleList',
                     mode + 'AlignedTiltDir', mode + 'Template', mode + 'destination', mode + 'BinFactorReconstruction',
                     mode + 'maxParticleShift', mode + 'OffsetX', mode + 'OffsetY', mode + 'OffsetZ', mode + 'FSCFlag',
                     mode + 'MetaFileFlag', mode + 'AlignmentResultsFileFlag',polishParticles]
        mandatory_fill = [mode + 'particleList', mode + 'AlignedTiltDir', mode+'Template', mode + 'MetaFile',
                          mode + 'AlignmentResultsFile', mode + 'destination']



        [func(mode) for func in (self.updateMeta, self.updateFSCFlag, self.updateMetaFileFlag, self.updateAlignmentResultsFlag)]

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 mandatory_fill=mandatory_fill)
        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def populatePartPolishingBatchTable(self, key=''):
        self.selectFSCFilesWindow.close()
        particleFilesStart = sorted(self.polishLists.text().split('\n'))
        templates = sorted(self.templateList.text().split('\n'))
        fscfiles = [''] + sorted(self.fscList.text().split('\n'))
        particleFiles = []

        for particleFile in particleFilesStart:
            if '_tomogram_' in particleFile:
                particleFiles.append(particleFile)
            else:
                pLs = extractParticleListsByTomoNameFromXML(particleFile, directory=os.path.dirname(particleFile))
                for pl in pLs:
                    particleFiles.append(pl)

        particleFiles = sorted(particleFiles)

        headers = ["Filename particleList", "Run", "Origin", "Tilt Images",  'Template', 'dimZ', 'Bin factor recon',
                   "Max Shift", "Offset X", "Offset Y", "Offset Z", 'Meta File', 'Align Results', 'FSC Filter', '']
        types = ['txt', 'checkbox', 'combobox', 'combobox', 'combobox', 'lineedit', 'lineedit', 'lineedit', 'lineedit',
                 'lineedit', 'lineedit', 'combobox', 'comboboxF', 'combobox', 'txt']
        a = 40
        sizes = [0, 0, 80, 80, a, a, a, a, a, a, a, a, a, a]

        tooltip = ['Names of the particleList files',
                   'Check this box to run particle polishing.',
                   'From which step do you want to use tilt-images?',
                   'Folder with aligned results.',
                   'Template File',
                   'Z-dimension of tomogram from which particles are picked.',
                   'Binning factor used for the reconstruction.',
                   'Max particle shift (pixel). Set the area used to locate maximum.'
                   'Offset in X-dimension', 'Offset in Y-dimension', 'Offset in Z-dimension',
                   'Meta file from which angles are extracted',
                   'FSC File used for filtering the cross-correlation']

        values = []
        refmarkindices = []

        for n, particleFile in enumerate(particleFiles):
            if not particleFile: continue

            path = os.path.basename(particleFile)
            for i in ('particleList_', 'coords_','_flipped'):
                path = path.replace(i,'')

            base, ext = os.path.splitext(path)

            if '_tomogram_' in base: base = 'tomogram_' + base.split('_tomogram_')[1]

            for t in ('WBP', 'INFR'):
                if t in base: base = base.split(t)[0] + t

            if base + '.mrc' in os.listdir(self.tomogramfolder) or base + '.em' in os.listdir(self.tomogramfolder):
                if os.path.exists(os.path.join(self.tomogramfolder, base + '.mrc')):
                    folder = os.popen('ls -alrt {}.mrc'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                elif os.path.exists(os.path.join(self.tomogramfolder, base + '.em')):
                    folder = os.popen('ls -alrt {}.em'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                else:
                    folder = ''
                if not folder: continue


                dimZ = read_size(folder.split()[-1], 'z')
                folder = os.path.dirname(folder.split()[-1])

                # markerfile  = os.path.join(folder, 'markerfile.txt')
                # markerdata = readMarkerfile(markerfile,61)
                sor = os.path.join(os.path.dirname(os.path.dirname(folder)), 'sorted')
                al = os.path.join(os.path.dirname(os.path.dirname(folder)), 'alignment')
                ctf = os.path.join(os.path.dirname(os.path.dirname(folder)), 'ctf')
                choices = [al + '/' + f for f in os.listdir(al) if 'marker_' in f and os.path.isdir(al + '/' + f)]
                # choices = list(map(str,range(markerdata.sizeZ()))) # + ['closest']
                # a = sorted(glob.glob('{}/Reconstruction*-*.out'.format(folder)))[-1]

                try:
                    from lxml import etree
                    xmlObj = etree.parse(particleFile)
                    particles = xmlObj.xpath('Particle')
                    binning = int(particles[0].xpath('InfoTomogram')[0].get('BinningFactor'))
                    refmarkindex = int(particles[0].xpath('InfoTomogram')[0].get('RefMarkIndex'))
                except:
                    print('Default values used for {}:\n\tbin recon = 8\n\t ref mark index = 1'.format(
                        os.path.basename(particleFile)))
                    binning = 8
                    refmarkindex = 1
                # binning = os.popen('cat {} | grep "--referenceMarkerIndex" '.format(a)).read()[:-1]
                # print(binning)
                origin = ['alignment', 'ctf', 'sorted']
                metafiles = [os.path.join(sor, f) for f in os.listdir(sor) if f.endswith('.meta')] + ['']

                import glob
                query = os.path.join(al,'*/alignmentResults.txt')
                alignfiles= [os.path.join(sor, f) for f in glob.glob(query)] + ['']

                values.append([particleFile, True, origin, choices, templates, dimZ, binning, 25, 0, 0, 0, metafiles, alignfiles, fscfiles, ''])
                refmarkindices.append(refmarkindex)

        try:
            self.num_nodes[key].setParent(None)
        except:
            pass

        if values:
            self.valuesBatchParticlePolishing = values
            self.fill_tab(key, headers, types, values, sizes, tooltip=tooltip, nn=True,wname=self.stage + 'BatchParticlePolish_')
            self.tabPolish_widgets = self.tables[key].widgets
            for n, index in enumerate(refmarkindices):
                self.tabPolish_widgets['widget_{}_3'.format(n)].setCurrentIndex(index)

            for row in range(len(values)):
                w = self.tabPolish_widgets['widget_{}_2'.format(row)]
                w.currentIndexChanged.connect(lambda d, r=row, v=values: self.updateChoicesPP(r, v))

            self.particleFilesBatchPolish = particleFiles
            self.pbs[key].clicked.connect(lambda dummy, pid=key, v=values: self.massPolishParticles(pid, v))
        else:
            return

        for i in range(len(values)):
            self.updateChoicesPP(i, values)

    def massPolishParticles(self, pid, values):
        qIDs =[]
        num_nodes = int(self.num_nodes[pid].value())
        num_submitted_jobs = nsj = 0
        values = self.valuesBatchParticlePolishing
        for row in range(self.tables[pid].table.rowCount()):
            if self.tabPolish_widgets['widget_{}_1'.format(row)].isChecked():
                particleXML = values[row][0]  # [row][0]
                tomoindex = particleXML.split('_tomogram_')[-1][:3]
                origin = values[row][2][self.tabPolish_widgets['widget_{}_{}'.format(row, 2)].currentIndex()]
                folder_aligned = values[row][3][self.tabPolish_widgets['widget_{}_{}'.format(row, 3)].currentIndex()]
                template = values[row][4][self.tabPolish_widgets['widget_{}_{}'.format(row, 4)].currentIndex()]

                dimZ = self.tabPolish_widgets['widget_{}_{}'.format(row, 5)].text()
                bin_read = self.tabPolish_widgets['widget_{}_{}'.format(row, 6)].text()
                maxShift = self.tabPolish_widgets['widget_{}_{}'.format(row, 7)].text()
                offx = self.tabPolish_widgets['widget_{}_{}'.format(row, 8)].text()
                offy = self.tabPolish_widgets['widget_{}_{}'.format(row, 9)].text()
                offz = self.tabPolish_widgets['widget_{}_{}'.format(row, 10)].text()
                meta = values[row][11][self.tabPolish_widgets['widget_{}_{}'.format(row, 11)].currentIndex()]
                alignfile = values[row][12][self.tabPolish_widgets['widget_{}_{}'.format(row, 12)].currentIndex()]
                fscfile = values[row][13][self.tabPolish_widgets['widget_{}_{}'.format(row, 13)].currentIndex()]

                fscflag = f'--FSCPath {fscfile}' if fscfile else ''
                metaflag = f'--metaFile {meta}' if meta else ''
                alignflag = f'--alignmentResultsFile {alignfile}' if alignfile else ''

                refid = folder_aligned.split('marker_')[1].split('_')[0]

                # generate output directory for particle polishing results if output directory does not exist.
                destination = os.path.join(f'ParticlePolishing/tomogram_{tomoindex}')
                if not os.path.exists(os.path.join(self.subtomodir, destination)):
                    os.mkdir(os.path.join(self.subtomodir, destination))

                execfilename = os.path.join(self.subtomodir, destination, f'inputParticlePolishing.sh')
                qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchParticlePolish'].values()

                paramsCmd = [self.subtomodir, str(cores * n_nodes), self.pytompath, particleXML, folder_aligned,
                             template, destination, bin_read, maxShift, offx, offy, offz, fscflag, metaflag, alignflag]

                txt = polishParticles.format(d=paramsCmd)
                jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder,
                                                       name='PolishParticles_{}'.format(nsj % num_nodes),
                                                       num_jobs_per_node=cores, time=time, partition=qname,
                                                       modules=modules, num_nodes=n_nodes, cmd=qcmd) + txt

                ID, num = self.submitBatchJob(execfilename, pid, jobtxt)
                nsj += 1
                qIDs.append(ID)


        if nsj:
            self.popup_messagebox('Info', 'Submission Status', f'Submitted {nsj} jobs to the queue.')
            self.addProgressBarToStatusBar(qIDs, key='QJobs', job_description='Particle Polish')

    def addAverageFields(self, mode, title=''):

        tooltip = 'Calculate the average of the particles described in a particle list.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets


        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list if v1 and v2 are not available.', mode='file',
                                    filetype='xml')
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_line(parent, "Outname", mode + 'outName', cstep=-1,
                               tooltip="What name should the output file have (.mrc and .em extentions are allowed)")
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', cstep=0,rstep=1,
                                  tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")





        # Connected Widgets
        self.widgets[mode + 'outFolder'].textChanged.connect(
            lambda d, m=mode: self.updateOutFileName(m))
        self.widgets[mode + 'outName'].textChanged.connect(
            lambda d, m=mode: self.updateOutFileName(m))

        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        # Widgets Updated When Other Widgets Are Updated
        self.widgets[mode + 'outFileName'] = QLineEdit('')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')

        # Parameters for execution

        exefilename = [mode + 'outFolder', 'AverageParticleList.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='APL', folder=self.logfolder,
                                                      id='AverageParticleList')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode+ 'outFileName', mode + 'numberMpiCores',
                     mode + 'gpuString', templateAverageParticleList]
        mandatory_fill = [mode + 'particleList', mode + 'outName', mode + 'outFolder']


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)

        # Run Update With Data From Logfile
        self.updateOutFileName(mode)
        self.updateGpuString(mode)

        setattr(self, mode + 'gb_APL', groupbox)
        return groupbox

    def addFRMAlignmentFields(self, mode=None, title=''):
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList', initdir=self.pickpartdir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask', initdir=self.frmdir,
                                    tooltip='Select the mask file.', mode='file', filetype=['em', 'mrc'], cstep=1,
                                    rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask', self.frmdir])
        self.insert_label_line_push(parent, 'Filename Average', mode + 'filenameAverage', initdir=self.frmdir,
                                    tooltip='Choose a filename for the average of all particles.', mode='file',
                                    filetype=['em', 'mrc'], cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Average', rstep=1, cstep=-3, action=self.gen_average,
                               params=[mode + 'particleList', mode + 'filenameAverage', mode + 'outputDir'])

        self.insert_label_line_push(parent, 'Output Directory', mode + 'outputDir',initdir=self.frmdir,
                                    tooltip='Folder in which the output will be written.')
        self.insert_label(parent, rstep=1, cstep=0)
        self.insert_label_spinbox(parent, mode + 'bwMin', 'Min Order SH Polynomial',
                                  value=8, minimum=0, stepsize=1,
                                  tooltip='The minimal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'bwMax', 'Max Order SH Polynomial',
                                  value=64, minimum=0, stepsize=1,
                                  tooltip='The maximal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'frequency', 'Frequency (px)',
                                  value=8, minimum=2, stepsize=1,
                                  tooltip='The maximum frequency to consider for the first iteration (in fourier \n'
                                          'space pixels from the center). If the reference has a known resolution \n'
                                          'you can provide this as the starting max resolution. If the reference has \n'
                                          'no known resolution set to the minimum value of 2')
        self.insert_label_spinbox(parent, mode + 'maxIterations', 'Maximum Iterations',
                                  value=8, minimum=1, stepsize=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
        self.insert_label_spinbox(parent, mode + 'peakOffset', 'Peak Offset',
                                  value=0, minimum=0, stepsize=1,
                                  tooltip='Sets the peak offset.')
        self.insert_label(parent, rstep=1, cstep=0)
        self.insert_label_spinbox(parent, mode + 'pixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox, minimum=0.1, stepsize=0.1, value=1.75)
        self.insert_label_spinbox(parent, mode + 'particleDiameter', 'Particle Diameter (A)',
                                  minimum=10, stepsize=1, value=300, maximum=10000, width=150)
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=0, stepsize=1,
                                  minimum=1, value=1, tooltip='Perform binning (downscale) of subvolumes by factor. '
                                                              'Default=1.')

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updateFRM(m))

        rscore = 'False'
        weightedAv = 'False'
        weighting = ''
        sphere = 'True'
        ad_res = '0.00'
        fsc = '0.50'

        jobfilename = [mode + 'outputDir', 'job_description.xml']  # os.path.join(self.frmdir, 'job_description.xml')
        exefilename = [mode + 'outputDir', 'frmAlignment.sh']  # os.path.join(self.frmdir, 'frmAlignment.sh')

        paramsSbatch = guiFunctions.createGenericDict(fname='FRMAlign', folder=self.logfolder,
                                                      id='FRMAlignment')  # , modules=['openmpi/2.1.1', 'python/2.7', 'lib64/append', 'pytom/dev/gui'])
        paramsJob = [mode + 'bwMin', mode + 'bwMax', mode + 'frequency', mode + 'maxIterations', mode + 'peakOffset',
                     rscore, weightedAv, mode + 'filenameAverage', weighting, mode + 'filenameMask', mode + 'binning',
                     sphere,
                     mode + 'pixelSize', mode + 'particleDiameter', mode + 'particleList', mode + 'outputDir']
        paramsCmd = [self.subtomodir, self.pytompath, jobfilename, mode + 'numberMpiCores', templateFRMSlurm]
        mandatory_fill = [mode + 'particleList', mode + 'filenameMask', mode+'filenameAverage', mode + 'outputDir']

        self.insert_gen_text_exe(parent, mode, xmlfilename=jobfilename, jobfield=True, exefilename=exefilename,
                                 paramsXML=paramsJob + [templateFRMJob], paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)

        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def addGLocalAlignmentFields(self,mode, title=''):
        tooltip = 'Run pytom GLocal routine.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList', initdir=self.pickpartdir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Initial reference model', mode + 'referenceModel', mode='file',
                                    filetype=['em','mrc'], enabled=True,
                                    tooltip='Reference : the initial reference - if none provided average of particle list')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask',self.glocaldir])
        self.insert_label_line_push(parent, 'Output Directory', mode + 'destination', mode='folder',initdir=self.glocaldir,
                                    tooltip='Select the destination directory.')
        self.insert_label_spinbox(parent, mode + 'numIterations', 'Number of Iterations',
                                  minimum=1, value=4, stepsize=1,
                                  tooltip='Sets the number of iterations.')
        self.insert_label_spinbox(parent, mode + 'pixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox, minimum=0.1, value=1.75, stepsize=0.1,
                                  tooltip='Pixelsize in Angstrom ')
        self.insert_label_spinbox(parent, mode + 'particleDiameter', 'Particle Diameter (A)',
                                  minimum=10, stepsize=1, value=300, maximum=10000,
                                  rstep=1, cstep=-1, tooltip='Particle diameter in Angstrom.')
        self.insert_label_spinbox(parent, mode + 'angleShells', 'Number of angular shells',
                                  tooltip='# Angle shells used for angular refinement.',
                                  minimum=1, stepsize=1, value=3, maximum=100,
                                  rstep=1, cstep=-1)
        self.insert_label_spinbox(parent, mode + 'angleIncrement', 'Angular Increment (degrees)',
                                  minimum=1, stepsize=1, value=3, maximum=359, wtype=QDoubleSpinBox,
                                  rstep=1, cstep=-1, tooltip='Angular increment for refinement.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=-1,
                                  stepsize=1,minimum=1,value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0,
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.widgets[mode+'jobName'] = QLineEdit()
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets['referenceCommand'] = QLineEdit(self)
        self.widgets['referenceCommand'].setVisible(False)

        self.widgets[mode + 'referenceModel'].textChanged.connect(lambda dummy, m=mode: self.updateReference(m))
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updatePixelSize(m))
        self.widgets[mode + 'destination'].textChanged.connect(lambda d, m=mode: self.updateJobname(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        self.updateJobname(mode)
        #self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.glocaldir: self.updateOutFolder(mode, p))
        exefilename = [mode+'destination', 'GLocal_Alignment.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='GLocal', folder=self.logfolder, id='GLocalAlignment')
        paramsCmd = [self.subtomodir, self.pytompath, self.pytompath, mode+'particleList', 'referenceCommand',
                     mode+'filenameMask', mode+'numIterations', mode+'pixelSize', mode+'particleDiameter',
                     mode+'binning', mode+'jobName', mode+'destination', mode + 'angleShells',
                     mode + 'angleIncrement', mode + 'numberMpiCores', mode + 'gpuString',  templateGLocal]

        mandatory_fill=[mode+'particleList', mode+'filenameMask', mode+'destination' ]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)

        self.updateGpuString(mode)
        self.updateReference(mode)
        self.updatePixelSize(mode)

        setattr(self, mode + 'gb_GLocal', groupbox)
        return groupbox

    def addCPCAFields(self,mode='', title=''):
        tooltip = 'The CCC is further used for classification. This script computes \n'+\
                  'the eigenvectors of the CCC and projects the data on the first \n'+\
                  'neig eigenvectors. Subsequently, these multidimensional vectors \n'+\
                  'are clustered into nclass groups using a kmeans algorithm.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',initdir=self.glocaldir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',initdir=self.cpcadir,
                                    tooltip='Select/Create an output folder.')
        self.insert_label_line(parent, 'Filename Output particleList', mode + 'outputFilename',
                               tooltip='Filename for generated XML file that includes the assigned classes for each particle. No full path needed.')
        self.insert_label_line_push(parent, 'CCC File', mode + 'cccFile',initdir=self.cpcadir,
                                    tooltip='Select the constrained correlation matrix file from the previous step.', mode='file', filetype='csv')
        self.insert_label_spinbox(parent, mode + 'numEig', text='Number of Eigenvectors',
                                  value=4, minimum=1, stepsize=1, maximum=100,
                                  tooltip='Sets the number of eigenvectors (corresponding to largest eigenvectors) used for clustering.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes',
                                  value=4, minimum=1,stepsize=1, maximum=1000,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_line(parent, 'Prefix', mode + 'prefix', rstep=1, cstep=0,
                               tooltip='Root for generated averages of the corresponding classes. The files will be called "Prefix"_iclass.em.')

        self.widgets[mode + 'particleList'].textChanged.connect(
            lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode, p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        exefilename = [mode + 'outFolder', 'CPCA_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CPCA', folder=self.logfolder, id='CPCA')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'outputFilename',
                     mode + 'cccFile', mode + 'numEig', mode+'numClasses', mode + 'outFolder', mode+'prefix',  templateCPCA]
        mandatory_fill = [mode + 'particleList', mode+'outputFilename', mode + 'outFolder', mode + 'cccFile']

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)

        setattr(self, mode + 'gb_CPCA', groupbox)
        return groupbox

    def addCCCFields(self,mode='', title=''):
        tooltip = 'Calculate the pairwise constrained cross correlation.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',initdir=self.glocaldir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0, initdir=self.cpcadir)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask', self.cpcadir])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',initdir=self.cpcadir,
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'lowpass', 'Lowpass filter (px)',
                                  minimum=0, maximum=1024, stepsize=1, value=20,
                                  tooltip='The lowpass filter is applied to all subtomograms after binning.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=-1,
                                  minimum=1, stepsize=1, value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0,
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'gpuString'] = QLineEdit('')

        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode,p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        exefilename = [mode + 'outFolder', 'CCC_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CCC_Class', folder=self.logfolder,
                                                      id='PairwiseCrossCorrelation')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'filenameMask',
                     mode + 'lowpass', mode + 'binning', mode + 'outFolder', mode + 'numberMpiCores', mode + 'gpuString',
                     templateCCC]
        mandatory_fill = [mode + 'particleList', mode + 'filenameMask', mode + 'outFolder']


        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)
        self.updateGpuString(mode)

        setattr(self, mode + 'gb_CCC', groupbox)
        return groupbox

    def addAC3DFields(self,mode, title=''):
        tooltip = 'Run autofocussed classification.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',initdir=self.glocaldir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Classification Mask', mode + 'filenameClassificationMask', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,initdir=self.acdir,
                                    tooltip='This mask is used for constraining the calculation of the focused mask. (Optional)', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameClassificationMask', self.acdir])
        self.insert_label_line_push(parent, 'Alignment Mask', mode + 'filenameAlignmentMask', mode='file',
                                    filetype=['mrc', 'em'], enabled=True, initdir=self.acdir,
                                    tooltip='This mask is only used for the alignment purpose. Only specify it if the particle list is not aligned (Optional).', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameAlignmentMask', self.acdir])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder', initdir=self.acdir,
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes', stepsize=1, value=4, minimum=2,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_spinbox(parent, mode + 'maxIterations', 'Number Of Iterations',stepsize=1,value=10, minimum=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
        self.insert_label_spinbox(parent, mode + 'binningFactor', 'Binning Factor',stepsize=1,value=1, minimum=1,
                                  tooltip='Sets the binning factor of the input data.')
        self.insert_label_spinbox(parent, mode + 'bwMax', 'Max Bandwidth Reconstruction (px)',
                                  minimum=1, maximum=1024, stepsize=1, value=20,
                                  tooltip='The maximal frequency used for reconstruction.')
        self.insert_label_spinbox(parent, mode + 'peakOffset', 'Max Peak Offset', stepsize=1,value=10, minimum=1,
                                  tooltip='Sets the peak offset.')
        self.insert_label_spinbox(parent, mode+'noisePercentage', 'Noise Percentage',
                                  wtype=QDoubleSpinBox, value=0.1, stepsize=.1,
                                  tooltip='Noise percentage (between 0 and 1). If you estimate your dataset contains certain amount of noise outliers, specify it here.')
        self.insert_label_spinbox(parent, mode + 'partDensThresh', 'Particle Density Threshold',
                                  wtype=QDoubleSpinBox, value=0., minimum=-6.0, maximum=6.0, stepsize=.1,
                                  tooltip='Particle density threshold for calculating the difference map (optional, by default 0).\n'
                                          'Two other most common choise are -2 and 2. -2 means all the values of the subtomogram \n'
                                          'below the -2 sigma will be used for calculating the difference mask (negative values\n'
                                          'count). 2 means all the values of the subtomogram above the 2 sigma will be used for\n'
                                          'calculating t he difference mask (positive values count). Finally, 0 means all the \n'
                                          'values are used for the calculation.')
        self.insert_label_spinbox(parent, mode + 'stdDiffMap', 'STD Threshold Diff Map', rstep=1, cstep=0,
                                  wtype=QDoubleSpinBox, stepsize=.1, minimum=0, maximum=1, value=0.4,
                                  tooltip='STD threshold for the difference map (optional, by default 0.4). This value should be \n'
                                          'between 0 and 1. 1 means only the place with the peak value will be set to 1 in the \n'
                                          'difference map (too much discriminative ability). 0 means all the places with the value\n'
                                          'above the average of STD will be set to 1 (not enough discriminative ability).')

        # Connected Widgets
        self.widgets[mode + 'filenameAlignmentMask'].textChanged.connect(
            lambda d, m=mode: self.updateAlignmentMaskFlag(m))
        self.widgets[mode + 'filenameClassificationMask'].textChanged.connect(
            lambda d, m=mode: self.updateClassificationMaskFlag(m))
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.acpath: self.updateOutFolder(mode, p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        # Widgets Updated When Other Widgets Are Updated
        self.widgets[mode + 'flagAlignmentMask'] = QLineEdit('')
        self.widgets[mode + 'flagClassificationMask'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')

        # Parameters for execution
        exefilename = [mode + 'outFolder', 'AC_Classification.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='AutoFocus', folder=self.logfolder,
                                                      id='AutoFocusClassification')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'flagAlignmentMask',
                     mode + 'flagClassificationMask', mode + 'numClasses', mode + 'bwMax', mode + 'maxIterations',
                     mode + 'peakOffset', mode + 'noisePercentage', mode + 'partDensThresh', mode + 'stdDiffMap',
                     mode + 'outFolder', mode + 'numberMpiCores', mode + 'binningFactor', templateAC]
        mandatory_fill = [mode + 'particleList', mode + 'outFolder']


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=mandatory_fill)

        # Run Update With Data From Logfile
        self.updateAlignmentMaskFlag(mode)
        self.updateClassificationMaskFlag(mode)
        self.updateOutFolder(mode, self.acpath)

        setattr(self, mode + 'gb_AC', groupbox)
        return groupbox

    def addFSCFields(self, mode, title=''):
        tooltip = 'Run fourier shell correlation function.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets

        self.insert_label_line_push(parent, 'Volume 1', mode + 'volume1', mode='file', initdir=self.glocaldir,
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='The first volume path. (Optional)')
        self.insert_label_line_push(parent, 'Volume 2', mode + 'volume2', mode='file', initdir=self.glocaldir,
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='The second volume path. (Optional)')
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList', initdir=self.glocaldir,
                                    tooltip='Select the particle list if v1 and v2 are not available.', mode='file',
                                    filetype='xml')
        self.insert_label_line_push(parent, 'Mask', mode + 'mask',initdir=self.fscdir,
                                    tooltip='Mask (optional, but recomended).', mode='file',
                                    filetype=['em', 'mrc'],cstep=1,rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_fsc_mask,
                               params=[mode + 'mask', self.fscdir])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder', initdir=self.fscdir,
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'fsc', 'FSC cutoff', stepsize=1, value=0.143, minimum=0, maximum=1.,
                                  wtype=QDoubleSpinBox,decimals=3,
                                  tooltip='The FSC criterion. Value between 0.0 and 1.0. Standard values are 0.5, 0.3. or 0.17')
        self.insert_label_spinbox(parent, mode + 'pixelsize', 'Pixelsize (A)',stepsize=1,value=2.62, minimum=.1, decimals=3,
                                  tooltip='Pixelsize in Angstrom', wtype=QDoubleSpinBox)
        self.insert_label_spinbox(parent, mode + 'randomizePhases', 'randomizePhases', stepsize=1,value=0.8, minimum=0, maximum=1,
                                  tooltip='Check validity of FSC using phases randomization beyond spatial frequency'
                                          ' where the uncorrected FSC curve drops below set threshold. Values between '
                                          '0, 1. No phase randomization check performed when value is set to 0.', decimals=3,
                                  wtype=QDoubleSpinBox)
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', cstep=-1,
                                  tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.insert_label_checkbox(parent, mode + 'CombinedResolution', 'Combined Resolution',
                                   tooltip='Check this box to calculate the resolution of the two volumes together.')

        self.insert_label_checkbox(parent, mode + 'plot', 'Plot Results',
                                   tooltip='Check this box to plot the results.', cstep=0, rstep=1)


        # Connected Widgets
        self.widgets[mode + 'volume1'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 1))
        self.widgets[mode + 'volume2'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 2))
        self.widgets[mode + 'particleList'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 0))
        self.widgets[mode + 'mask'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 3))
        self.widgets[mode + 'plot'].stateChanged.connect(
            lambda d, m=mode: self.updateFSCPlotFlag(m))
        self.widgets[mode + 'CombinedResolution'].stateChanged.connect(
            lambda d, m=mode: self.updateFSCCombResFlag(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        # Widgets Updated When Other Widgets Are Updated
        self.widgets[mode + 'flagVolume1'] = QLineEdit('')
        self.widgets[mode + 'flagVolume2'] = QLineEdit('')
        self.widgets[mode + 'flagParticleList'] = QLineEdit('')
        self.widgets[mode + 'flagMask'] = QLineEdit('')
        self.widgets[mode + 'flagPlot'] = QLineEdit('')
        self.widgets[mode + 'flagCombinedResolution'] = QLineEdit('')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')

        # Parameters for execution

        exefilename = [mode + 'outFolder', 'FSC_Validation.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='FSC', folder=self.logfolder,
                                                      id='FSCValidation')
        paramsCmd = [self.fscdir, self.pytompath, mode + 'flagParticleList', mode + 'flagVolume1',
                     mode + 'flagVolume2', mode + 'flagMask', mode + 'outFolder', mode + 'fsc',
                     mode + 'pixelsize', mode + 'randomizePhases', mode + 'flagPlot', mode + 'flagCombinedResolution',
                     mode + 'gpuString', templateFSC]


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch, mandatory_fill=[mode+'outFolder'])

        self.widgets[mode+'queue'].stateChanged.connect(lambda d, m=mode: self.updateLog(m))

        # Run Update With Data From Logfile
        self.updateFSCFlags(mode)
        self.updateFSCPlotFlag(mode)
        self.updateLog(mode)
        self.updateFSCCombResFlag(mode)
        self.updateGpuString(mode)
        self.redoPixelSize(mode)

        setattr(self, mode + 'gb_FSC', groupbox)
        return groupbox

    # Helper functions

    def updateOutFileName(self, mode):
        try:
            folder = self.widgets[mode + 'outFolder'].text()
            fname = self.widgets[mode + 'outName'].text()
            self.widgets[mode + 'outFileName'].setText(os.path.join(folder, fname))
        except Exception as e:
            print(e)

    def updatePolishFlag(self, mode):
        fname = self.widgets[mode + 'polishFile'].text()
        if not fname:
            self.widgets[mode + 'polishFlag'].setText('')
            return

        self.widgets[mode + 'polishFlag'].setText(f'--particlePolishFile {fname}')

    def updateFSCFlag(self, mode):

        fname = self.widgets[mode + 'FSCPath'].text()
        if not fname:
            self.widgets[mode + 'FSCFlag'].setText('')
            return

        self.widgets[mode + 'FSCFlag'].setText(f'--FSCPath {fname}')

    def updateMetaFileFlag(self, mode):

        fname = self.widgets[mode + 'MetaFile'].text()
        if not fname:
            self.widgets[mode + 'MetaFileFlag'].setText('')
            return

        self.widgets[mode + 'MetaFileFlag'].setText(f'--metaFile {fname}')

    def updateAlignmentResultsFlag(self, mode):

        fname = self.widgets[mode + 'AlignmentResultsFile'].text()
        if not fname:
            self.widgets[mode + 'AlignmentResultsFileFlag'].setText('')
            return

        self.widgets[mode + 'AlignmentResultsFileFlag'].setText(f'--alignmentResultsFile {fname}')

    def updateMeta(self,mode):
        pl = self.widgets[mode + 'particlelist'].text()
        try: 
            tomoID = int(pl.split('_tomogram_')[-1][:3])
            tomo = os.path.join(self.tomogram_folder, 'tomogram_{:03d}/sorted/'.format(tomoID))
            a = glob.glob(tomo+'*.meta')
            if not a: print('No meta file found. auto update stopped.')
            a = a[-1]
            self.widgets[mode+'MetaFile'].setText(a)
            #self.widgets[mode + 'AlignedTiltDir'].setText(os.path.dirname(a))
        except Exception as e:
            print(e)
            pass

    def update_alignment_options(self, mode):
        pl = self.widgets[mode + 'particlelist'].text()
        # to get the tomogram id from the particle list  # TODO could be improved by reading tomogram from list
        tomo_id = int(pl.split('_tomogram_')[-1][:3])  # TODO this formattting i really dont like
        self.tomogram = 'tomogram_{:03d}'.format(tomo_id)

        # clear previous alignment choices
        self.widgets[mode + 'alignment'].disconnect()
        self.widgets[mode + 'alignment'].clear()
        # atttempt to find valid alignments
        alignment_dir = os.path.join(self.tomogram_folder, self.tomogram, 'alignment')
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
        # set the alignment options
        self.widgets[mode + 'alignment'].currentIndexChanged.connect(
            lambda d, m=mode: self.update_alignment_choice(m))
        self.widgets[mode + 'alignment'].addItems(alignment_choices)

        # get ctf options
        self.widgets[mode + 'ctfCorrChoice'].disconnect()
        self.widgets[mode + 'ctfCorrChoice'].clear()
        sorted_dir = os.path.join(self.tomogram_folder, self.tomogram, 'sorted')
        ctf_sorted_dir = os.path.join(self.tomogram_folder, self.tomogram, 'ctf',
                                      'sorted_ctf')
        tilt_choices = []
        if len([f for f in os.listdir(sorted_dir) if f.endswith('.mrc')]) > 0:
            tilt_choices.append('sorted')
        if (os.path.exists(ctf_sorted_dir) and 
			len([f for f in os.listdir(ctf_sorted_dir) if f.endswith('.mrc')]) > 0):  # at least
            tilt_choices.append('sorted_ctf')
        # set choices in the widget
        self.widgets[mode + 'ctfCorrChoice'].currentIndexChanged.connect(
            lambda d, m=mode: self.update_ctf_corr_choice(m))
        self.widgets[mode + 'ctfCorrChoice'].addItems(tilt_choices)

    def update_alignment_choice(self, mode):
        folder = os.path.join(self.tomogram_folder, self.tomogram, 'alignment',
                              self.widgets[mode + 'alignment'].currentText(), 'GlobalAlignment')
        try:
            pointer = [d for d in os.listdir(folder) if 'sorted' in d][0]
            self.widgets[mode + 'alignmentResultsFile'].setText(os.path.join(folder, pointer, 'alignmentResults.txt'))

            # if an alignment file has been selected, we can start looking for angle range
            if os.path.exists(self.widgets[mode + 'alignmentResultsFile'].text()):
                try:
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

        except IndexError:
            print('No alignment result files found in alignment directory')

    def update_alignment_choice_batch(self, row_id, table_id):
        w_tomogram = self.tables[table_id].widgets['widget_{}_1'.format(row_id)]
        w_alignment = self.tables[table_id].widgets['widget_{}_3'.format(row_id)]
        w_first_angle = self.tables[table_id].widgets['widget_{}_4'.format(row_id)]
        w_last_angle = self.tables[table_id].widgets['widget_{}_5'.format(row_id)]

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
            self.widgets[mode + 'projectionDir'].setText(os.path.join(self.tomogram_folder, self.tomogram,
                                                                      'ctf', 'sorted_ctf'))
        else:
            self.widgets[mode + 'projectionDir'].setText(os.path.join(self.tomogram_folder, self.tomogram, 'sorted'))

    def updateHeaderChoices(self, rowID, key):
        print('Update Header', rowID, key)
        header, row = self.tables[key].general_widgets[3], self.tab12_widgets['widget_{}_3'.format(rowID)]
        AllItemsGeneral = [header.itemText(i) for i in range(header.count())]
        indPart = [row.itemText(i) for i in range(row.count())]
        header.clear()

        header.addItems([item for item in indPart if not item in AllItemsGeneral and 'CLOSEST' in item])
        header.addItems(AllItemsGeneral)

    def updateAlignmentType(self, rowID, table_id):
        print('Update alignment type')
        values = self.valuesBatchSubtomoReconstruction
        self.tab12_widgets['widget_{}_4'.format(rowID)].clear()
        current_index_f = self.tab12_widgets['widget_{}_2'.format(rowID)].currentIndex()
        current_index_a = self.tab12_widgets['widget_{}_3'.format(rowID)].currentIndex()
        folder = values[rowID][2][current_index_f]
        aligntype = values[rowID][3][current_index_a]

        if not '_closest_' in os.path.basename(folder.lower()):
            if aligntype:
                temp = os.path.join(folder, aligntype)
                origin = [f'{temp}/{f}' for f in os.listdir(temp) if os.path.isdir(f'{temp}/{f}')]
            else:
                origin = []
        else:
            origin = []
            angles = os.path.basename(folder).split('_')[-1]
            if aligntype:
                for ff in values[rowID][2]:
                    temp_angles = os.path.basename(ff).split('_')[-1]
                    if 'CLOSEST' in ff or temp_angles != angles:
                        continue
                    tempaligntype = [f for f in os.listdir(ff) if os.path.isdir(ff + '/' + f)]

                    if tempaligntype:
                        temporigin = [f for f in os.listdir(f'{ff}/{tempaligntype[0]}') if os.path.isdir(f'{ff}/{tempaligntype[0]}/{f}')]
                    else:
                        temporigin = []

                    for to in temporigin:
                        if not to in origin:
                            origin.append(to)

        self.valuesBatchSubtomoReconstruction[rowID][4] = origin
        for item in origin:
            self.tab12_widgets['widget_{}_4'.format(rowID)].addItem(os.path.basename(item))

    def updateChoices(self, rowID, table_id):
        print(f'Update Choices {rowID}')

        values = self.valuesBatchSubtomoReconstruction
        self.tab12_widgets['widget_{}_3'.format(rowID)].clear()
        self.tab12_widgets['widget_{}_4'.format(rowID)].clear()

        current_index = self.tab12_widgets['widget_{}_2'.format(rowID)].currentIndex()

        try:
           folder = values[rowID][2][current_index]
        except:
            return

        particleFile = values[rowID][0]
        a = 'tomogram_' + particleFile.split('_tomogram_')[1][:3]

        print('update for folder', folder)

        if not '_closest_' in os.path.basename(folder.lower()):
            aligntype = [folder + '/' + f for f in os.listdir(folder) if os.path.isdir(folder + '/' + f)]
            if aligntype:
                origin = [f'{aligntype[0]}/{f}' for f in os.listdir(aligntype[0]) if os.path.isdir(f'{aligntype[0]}/{f}')]
            else:
                origin = []
        else:
            aligntype = []
            origin = []
            angles = os.path.basename(folder).split('_')[-1]
            for ff in values[rowID][2]:
                temp_angles = os.path.basename(ff).split('_')[-1]
                if 'CLOSEST' in ff or temp_angles != angles:
                    continue
                tempaligntype = [f for f in os.listdir(ff) if os.path.isdir(ff + '/' + f)]

                if tempaligntype:
                    temporigin = [f for f in os.listdir(f'{ff}/{tempaligntype[0]}') if os.path.isdir(f'{ff}/{tempaligntype[0]}/{f}')]
                else:
                    temporigin = []

                for ta in tempaligntype:
                    if not ta in aligntype:
                        aligntype.append(ta)
                for to in temporigin:
                    if not to in origin:
                        origin.append(to)

        self.valuesBatchSubtomoReconstruction[rowID][3] = aligntype
        self.valuesBatchSubtomoReconstruction[rowID][4] = origin

        for item in aligntype:
            self.tab12_widgets['widget_{}_3'.format(rowID)].addItem(os.path.basename(item))

        if aligntype:
            self.updateAlignmentType(rowID, table_id)
        self.updateHeaderChoices(rowID, table_id)

    def updateChoicesPP(self, rowID, values):
        values = self.valuesBatchParticlePolishing
        self.tabPolish_widgets['widget_{}_3'.format(rowID)].clear()

        current_index = self.tabPolish_widgets['widget_{}_2'.format(rowID)].currentIndex()

        origin = values[rowID][2][current_index]
        particleFile = values[rowID][0]
        a = 'tomogram_' + particleFile.split('_tomogram_')[1][:3]

        folder = os.path.join(self.tomogram_folder, a, origin)

        choices = [folder + '/' + f for f in os.listdir(folder) if
                    'marker_' in f and os.path.isdir(folder + '/' + f)]

        self.valuesBatchParticlePolishing[rowID][3] = choices

        for item in choices:
            self.tabPolish_widgets['widget_{}_3'.format(rowID)].addItem(os.path.basename(item))

    def updateFRM(self, mode):
        item = self.widgets[mode + 'particleList'].text()
        if not item:
            return
        folder, ext = os.path.splitext( os.path.basename(item))
        outputDir = os.path.join(self.frmdir, folder.replace('particleList_', ''))
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        else:
            for i in range(1000):
                if not os.path.exists(outputDir+'_{:03d}'.format(i)):
                    outputDir += '_{:03d}'.format(i)
                    os.mkdir(outputDir)
                    break


        self.widgets[mode+'outputDir'].setText(outputDir)

    def updateReference(self, mode):
        if self.widgets[mode + 'referenceModel'].text():
            self.widgets['referenceCommand'].setText( "--reference " + self.widgets[mode + 'referenceModel'].text())
        else:
            self.widgets['referenceCommand'].setText("")

    def updateGpuString(self, mode):
        id = self.widgets[mode + 'gpuID'].text()
        try:
            a = list(map(int,[el for el in id.split(',') if el != '']))
        except:
            self.widgets[mode + 'gpuID'].setText('')
            self.popup_messagebox('Warning', 'Invalid value in field', 'Impossible to parse gpu IDs, field has been cleared.')
            return

        if len(id) > 0:
            self.widgets[mode + 'gpuString'].setText(f' --gpuID {id}')
            self.widgets[mode + 'numberMpiCores'].setText(f'{len(a)+1}')
        else:
            self.widgets[mode + 'gpuString'].setText('')

            if 'gLocal'in mode and 'GLocalAlignment' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['GLocalAlignment'].values()
            elif 'CCC' in mode and 'PairwiseCrossCorrelation' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['PairwiseCrossCorrelation'].values()
            elif 'FSC' in mode and 'FSCValidation' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['FSCValidation'].values()
            else:
                cores = 20
            self.widgets[mode + 'numberMpiCores'].setText(f'{cores}')

    def updatePixelSize(self, mode):
        from pytom.basic.structures import ParticleList

        try:
            particleList = self.widgets[mode + 'particleList'].text()
            if not particleList: return
            pl = ParticleList()
            pl.fromXMLFile(particleList)
            tomid = pl[0].getPickPosition().getOriginFilename().split('/tomogram_')[1].split('_')[0]
            metaquery = os.path.join( self.tomoanalysis, f'tomogram_{tomid}/sorted/*.meta')
            a = glob.glob(metaquery)
            if len(a):
                pixelsize = guiFunctions.loadstar(a[0], dtype=guiFunctions.datatype)['PixelSpacing'][0]
                self.widgets[mode + 'pixelSize'].setValue(pixelsize)
        except Exception as e:
            print(e)

    def updateJobname(self, mode):
        dest = self.widgets[mode+'destination'].text()
        if dest:
            self.widgets[mode + 'jobName'].setText( os.path.join(dest, 'glocal_input_params_'+os.path.basename(dest)+'.xml'))

    def updateClassificationMaskFlag(self, mode):
        try:
            mask = self.widgets[mode+'filenameClassificationMask'].text()
            if mask:
                self.widgets[mode + 'flagClassificationMask'].setText('\\\n-c {}'.format(mask) )
            else:
                self.widgets[mode + 'flagClassificationMask'].setText('')
        except:
            pass

    def updateAlignmentMaskFlag(self, mode):
        try:
            mask = self.widgets[mode+'filenameAlignmentMask'].text()
            if mask:
                self.widgets[mode + 'flagAlignmentMask'].setText('-m {}'.format(mask) )
            else:
                self.widgets[mode + 'flagAlignmentMask'].setText('-a')
        except:
            pass

    def updateOutFolder(self,mode, path, name='particleList'):
        pl = self.widgets[mode + name].text()
        if not pl or not os.path.exists(pl):
            self.widgets[mode + name].setText('')
            return

        if not self.doAllParticlesExist(pl):
            self.popup_messagebox('Warning', 'Particles do not exist.',
                                  f'Not all particles in {os.path.basename(pl)} have been created.\nPlease be aware that the job might fail. ')
            # self.widgets[mode + name].setText('')
            return

        folder = os.path.basename(pl).replace('particleList_','')[:-4]
        folder = os.path.join(path,folder)
        if os.path.exists(folder):
            folder = folder
            if not folder: return

        if not os.path.exists(folder):
            os.mkdir(folder)
            #os.system('ln -s {}/Subtomograms {}/Subtomograms'.format(self.subtomodir, folder ) )
        try: self.widgets[mode + 'outFolder'].setText(folder)
        except: self.widgets[mode + 'destination'].setText(folder)

    def updateFSCPlotFlag(self, mode):
        t = self.widgets[mode + 'plot'].isChecked()

        if t:
            self.widgets[mode + 'flagPlot'].setText('--plot ')
            self.widgets[mode + 'queue'].setChecked(False)
        else:
            self.widgets[mode + 'flagPlot'].setText('')

    def updateFSCCombResFlag(self, mode):
        t = self.widgets[mode + 'CombinedResolution'].isChecked()

        if t:
            self.widgets[mode + 'flagCombinedResolution'].setText('--combinedResolution')
        else:
            self.widgets[mode + 'flagCombinedResolution'].setText('')

    def updateLog(self, mode):
        if self.widgets[mode + 'queue'].isChecked():
            self.widgets[mode + 'plot'].setChecked(False)

    def updateFSCFlags(self, mode, id='all'):

        for n, (name, flag) in enumerate((('particleList', '--pl'), ('volume1','--v1'), ('volume2', '--v2'),
                                          ('mask', '--mask'))):
            if n != id and id != 'all': continue
            t = self.widgets[mode + name].text()

            if t: self.widgets[mode + 'flag' + name[:1].capitalize() + name[1:]].setText(f'{flag} {t}')
            else: self.widgets[mode + 'flag' + name[:1].capitalize() + name[1:]].setText('')


        di = {'-Even.':'-Odd.', '-Odd.':'-Even.'}
        if id == 1:
            t = self.widgets[mode + f'volume{id}'].text()
            for i in ('-Even.', '-Odd.'):
                if not i in t: continue
                otherHalfMap = t.replace(i,di[i])
                #print(otherHalfMap == self.widgets[mode + f'volume{3-id}'].text(), otherHalfMap, self.widgets[mode + f'volume{3-id}'].text())
                if os.path.exists(otherHalfMap) and self.widgets[mode + f'volume{3-id}'].text() != otherHalfMap:
                    self.widgets[mode + f'volume{3-id}'].setText(otherHalfMap)
                    break

            self.redoPixelSize(mode)

    def createOutFolder(self, mode):
        folder = self.widgets[mode + 'outFolder'].text()
        if not folder: return
        if not os.path.exists(folder):
            os.mkdir(folder)

    def doAllParticlesExist(self, particleList):
        from pytom.basic.structures import ParticleList
        if not particleList: return
        pl = ParticleList()
        pl.fromXMLFile(particleList)
        allParticlesExist = True
        for particle in pl:
            if not os.path.exists(os.path.join(self.subtomodir ,particle.getFilename())):
                allParticlesExist = False
                break
        return allParticlesExist

    def gen_average(self, params):
        key_particleList, key_filename_average, key_outputDir = params
        particleList = self.widgets[key_particleList].text()
        if not particleList:
            self.popup_messagebox('Warning', 'Averaging Failed', 'Averaging did not succeed. No particle list selected.')
            return
        folder = self.widgets[key_outputDir].text()
        if not os.path.exists(folder): os.mkdir(folder)
        output_name = os.path.join( folder, 'average.em')

        out = os.popen('cd {}; average.py -p {} -a {} -c 4'.format(self.subtomodir, particleList, output_name)).read()
        print(out)
        if not os.path.exists(particleList):
            self.popup_messagebox('Warning', 'Averaging Failed',
                                  'Averaging did not succeed, please try again.')
            return
        self.widgets[key_filename_average].setText(output_name)

    def gen_mask(self,params):
        maskfilename = CreateMaskFile(self, maskfname=params[0],folder=params[1])
        maskfilename.show()

    def selectTemplates(self, key):
        self.selectPLWindow.close()
        try: self.templateList.text()
        except: self.templateList = QLineEdit()

        self.selectTemplatesWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.em', '.mrc'],
                                                 outputline=self.templateList, id=key,
                                                 run_upon_complete=self.selectFSCFilterFiles,
                                                 title='Select templates')

    def selectFSCFilterFiles(self, key):
        self.selectTemplatesWindow.close()
        try: self.fscList.text()
        except: self.fscList = QLineEdit()

        self.selectFSCFilesWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.dat'],
                                                 outputline=self.fscList, id=key,
                                                 run_upon_complete=self.populatePartPolishingBatchTable,
                                                 title='Select FSC File for filtering cross-correlation (optional)')

    def redoPixelSize(self, mode):
        from pytom.alignment.alignmentStructures import GLocalSamplingJob

        t = self.widgets[mode + 'volume2'].text()
        if not t: return
        try:
            folder = os.path.dirname(t)
            files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('xml')]
            job = GLocalSamplingJob()
            for f in files:
                try:
                    job.fromXMLFile(f)
                    break
                except Exception as e:
                    pass
            self.widgets[mode + 'pixelsize'].setValue(
                float(job.samplingParameters.sampleInformation.getPixelSize()))
        except Exception as e:
            print(e)

    def gen_fsc_mask(self,params):
        try:
            maskfilename = CreateFSCMaskFile(self, params[0], params[1])
            maskfilename.show()
        except Exception as e:
            print(e)


