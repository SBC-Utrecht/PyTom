import sys
import os
import random
import glob
import numpy
import time

from multiprocessing import Manager, Event, Process
from ftplib import FTP_TLS, FTP
from os.path import dirname, basename

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets

from pytom.basic.files import read
from pytom.gui.guiStyleSheets import *
from pytom.gui.guiSupportCommands import *
from pytom.gui.guiStructures import *
from pytom.gui.fiducialAssignment import FiducialAssignment
from pytom.gui.guiFunctions import avail_gpu
import pytom.gui.guiFunctions as guiFunctions
from pytom.bin.extractTomoNameFromXML import *


class SubtomoAnalysis(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(SubtomoAnalysis, self).__init__(parent)
        self.stage          = 'v04_'
        self.pytompath      = self.parent().pytompath
        self.projectname    = self.parent().projectname
        self.logfolder      = self.parent().logfolder
        self.subtomodir     = self.parent().subtomo_folder
        self.frmdir         = os.path.join(self.subtomodir,'Alignment/FRM')
        self.glocaldir      = os.path.join(self.subtomodir, 'GLocal/FRM')
        self.cpcadir        = os.path.join(self.subtomodir, 'Classification/CPCA')
        self.acdir          = os.path.join(self.subtomodir, 'Classification/AutoFocus')
        self.pickpartdir    = self.parent().particlepick_folder+'/Picked_Particles'
        self.tomogramfolder = os.path.join(self.parent().particlepick_folder, 'Tomograms')
        self.tomogram_folder = self.parent().tomogram_folder
        self.acpath = os.path.join(self.parent().subtomo_folder, 'Classification/AutoFocus')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        headers = ["Reconstruct Subtomograms","Align Subtomograms","Classify Subtomograms"]
        subheaders = [['Single Reconstruction','Batch Reconstruction'],['FRM Alignment','GLocal'],['CPCA','Auto Focus']]

        self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders)

        self.widgets = {}
        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.num_nodes = {}

        self.tabs = {'tab11': self.tab11, 'tab12': self.tab12,
                     'tab21': self.tab21, 'tab22': self.tab22,
                     'tab31': self.tab31, 'tab32': self.tab32,}

        self.tab_actions = {'tab11': self.tab11UI, 'tab12': self.tab12UI,
                            'tab21': self.tab21UI, 'tab22': self.tab22UI,
                            'tab31': self.tab31UI, 'tab32': self.tab32UI}

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if tt in ('tab11', 'tab21', 'tab22', 'tab31', 'tab32'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()
                button = QPushButton('Refresh Tab')
                button.setSizePolicy(self.sizePolicyC)
                button.clicked.connect(self.tab_actions[tt])

                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)

                if tt in ('tab12'):
                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])

                if not tt in ('tab12'):
                    self.tab_actions[tt]()

                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])

    def tab11UI(self):
        key = 'tab11'
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleReconstruction_'

        items += list(self.create_expandable_group(self.createSubtomograms, self.sizePolicyB, 'Single Reconstruction',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab12UI(self):
        try: self.extractLists.text()
        except: self.extractLists = QLineEdit()

        self.mass_extract = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.xml'],
                                        outputline=self.extractLists, run_upon_complete=self.populate_batch_create,
                                        title='Select particlLists')
        pass

    def tab21UI(self):
        key = 'tab21'
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1, t2 =  self.stage + 'inputFiles_', self.stage + 'frmSetttings_',self.stage + 'sampleInformation_'

        items += list(self.create_expandable_group(self.inputFiles, self.sizePolicyB, 'FRM Alignment',
                                                   mode=t0))
        items[-1].setVisible(False)


        for n, item in enumerate(items):
            grid.addWidget(item, n, 0,1,4)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab22UI(self):
        key = 'tab22'

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        items += list(self.create_expandable_group(self.glocal, self.sizePolicyB, 'GLocal Alignment',
                                                   mode=self.stage + 'gLocal_'))
        items[-1].setVisible(False)
        # items += list( self.create_expandable_group(self.createSubtomograms, self.sizePolicyB, 'Extract Subtomograms',
        #                                            mode=self.stage+'extract_') )
        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab31UI(self):
        key = 'tab31'
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'CCC_', self.stage + 'CPCA_'

        items += list(self.create_expandable_group(self.CCC, self.sizePolicyB, 'Pairwise Cross Correlation',
                                                   mode=t0))
        items[-1].setVisible(False)

        items += list(self.create_expandable_group(self.CPCA, self.sizePolicyB, 'CPCA',
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

    def tab32UI(self):
        key = 'tab32'
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'AC_', self.stage + 'CPCA_'

        items += list(self.create_expandable_group(self.ac, self.sizePolicyB, 'Autofocussed Classification',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def createSubtomograms(self, mode=''):

        title = "Single Reconstruction"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particlelist',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Folder with aligned tilt images', mode + 'AlignedTiltDir',
                                    'Select the folder with the aligned tilt images.')
        self.insert_label_line_push(parent, 'Meta file with tilt angles', mode + 'MetaFile',mode='file',filetype='meta',
                                    tooltip='Select the corresponding metafile.')



        self.insert_label_spinbox(parent, mode + 'BinFactorReconstruction', 'Binning factor used in the reconstruction.',
                                  'Defines the binning factor used in the reconstruction of the tomogram from which'+
                                  'the particles are selected.',
                                  minimum=1,stepsize=1,value=8)

        self.insert_label_spinbox(parent,  mode + 'WeightingFactor', 'Apply Weighting (0/1)',
                                  'Sets the weighting scheme applied to the tilt images.\n'+
                                  '0: no weighting.\n1: ramp filter.', minimum=-5, maximum=5, stepsize=1,value=0)

        self.insert_label_spinbox(parent,mode+'SizeSubtomos', 'Size subtomograms.','Sets the size of the subtomograms.',
                                  minimum=10,maximum=1000,stepsize=1,value=128)

        self.insert_label_spinbox(parent, mode+'BinFactorSubtomos', 'Binning Factor Subtomograms.',
                                  'Sets the binning factor of the subtomograms.',rstep=1,
                                  value=1, stepsize=1,minimum=1)

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

        self.widgets[mode + 'particlelist'].textChanged.connect(lambda d, m=mode: self.updateMeta(m))


        execfilename = os.path.join( self.subtomodir, 'Reconstruction/reconstructSubtomograms.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='subtomoReconstr', folder=self.logfolder,
                                                      id='SingleSubtomoReconstruct')
        paramsCmd = [mode+'particlelist', mode+'AlignedTiltDir', mode + 'BinFactorReconstruction',
                     mode+'SizeSubtomos', mode+'BinFactorSubtomos', mode+'OffsetX', mode+'OffsetY', mode+'OffsetZ',
                     self.subtomodir, mode+'WeightingFactor', mode+'MetaFile', '20', extractParticles]

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename,paramsSbatch=paramsSbatch)
        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def updateMeta(self,mode):
        pl = self.widgets[mode + 'particlelist'].text()
        try: 
            tomoID = int(pl.split('_tomogram_')[-1][:3])
            tomo = os.path.join(self.tomogram_folder, 'tomogram_{:03d}/sorted/'.format(tomoID))
            print(tomo)
            a = glob.glob(tomo+'*.meta')
            if not a: print('No meta file found. auto update stopped.')
            a = a[-1]
            self.widgets[mode+'MetaFile'].setText(a)
        except:
            pass

    def populate_batch_create(self):
        self.mass_extract.close()
        particleFilesStart = sorted( self.extractLists.text().split('\n') )
        particleFiles = []

        for particleFile in particleFilesStart:
            if '_tomogram_' in particleFile:
                particleFiles.append(particleFile)
            else:
                pLs = extractParticleListsByTomoNameFromXML(particleFile, directory=os.path.dirname(particleFile))
                for pl in pLs:
                    particleFiles.append(pl)


        particleFiles = sorted(particleFiles)

        id='tab12'
        headers = ["Filename particleList", "Run", "Origin", "Tilt Images", 'Bin factor recon', 'Weighting', "Size subtomos", "Bin subtomos", "Offset X", "Offset Y", "Offset Y", '']
        types = ['txt', 'checkbox', 'combobox', 'combobox', 'lineedit', 'lineedit', 'lineedit','lineedit', 'lineedit', 'lineedit', 'lineedit','txt']
        a=40
        sizes = [0, 0, 80, 80, a, a, a, a, a, a, a]

        tooltip = ['Names of the particleList files', 
                   'Check this box to run subtomogram reconstruction.',
                   'Which folder contains the tilt-images you want to use for subtomo reconstruction?',
                   'Aligned Images',
                   'Binning factor used for the reconstruction.',
                   'Weighting Type.\n0: No Weighting,\n1: Analytical Weighting.\n-1: Ramp Weighting',
                   'Size Subtomograms','Binning factor for subtomograms (--projBinning)','Offset in X-dimension',
                   'Offset in Y-dimension','Offset in Z-dimension']

        values = []
        refmarkindices = []
        for n, particleFile in enumerate( particleFiles ):
            if not particleFile: continue
            base, ext = os.path.splitext(os.path.basename(particleFile).replace('particleList_', '').replace('coords_','').replace('_flipped',''))
            if '_tomogram_' in base: base = 'tomogram_' + base.split('_tomogram_')[1]

            for t in ('WBP', 'INFR'):
                if t in base: base = base.split(t)[0]+t

            if base+'.mrc' in os.listdir(self.tomogramfolder) or base+'.em' in os.listdir(self.tomogramfolder):
                if os.path.exists(os.path.join(self.tomogramfolder, base+'.mrc')):
                    folder = os.popen('ls -alrt {}.mrc'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                elif os.path.exists(os.path.join(self.tomogramfolder, base+'.em')):
                    folder = os.popen('ls -alrt {}.em'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                else:
                    folder=''
                if not folder: continue
                folder = os.path.dirname(folder.split()[-1])

                markerfile  = os.path.join(folder, 'markerfile.em')
                markerdata = read(markerfile,binning=[1,1,1])

                al  = os.path.join(os.path.dirname(os.path.dirname(folder) ),'alignment')
                ctf = os.path.join(os.path.dirname(os.path.dirname(folder) ),'ctf')
                choices = [al+'/'+f for f in os.listdir(al) if 'unweighted_unbinned' in f and os.path.isdir(al+'/'+f)]
                #choices = list(map(str,range(markerdata.sizeZ()))) # + ['closest']
                #a = sorted(glob.glob('{}/Reconstruction*-*.out'.format(folder)))[-1]

                try:
                    from lxml import etree
                    xmlObj = etree.parse(particleFile)
                    particles = xmlObj.xpath('Particle')
                    binning = int( particles[0].xpath('InfoTomogram')[0].get('BinningFactor') )
                    refmarkindex = int(particles[0].xpath('InfoTomogram')[0].get('RefMarkIndex'))
                except:
                    print('Default values used for {}:\n\tbin recon = 8\n\t ref mark index = 1'.format(os.path.basename(particleFile)))
                    binning = 8
                    refmarkindex = 1
                #binning = os.popen('cat {} | grep "--referenceMarkerIndex" '.format(a)).read()[:-1]
                #print(binning)
                origin = ['alignment', 'ctf', 'sorted']
                values.append([particleFile, True, origin, choices, binning, -1, 128, 1, 0, 0, 0, ''])
                refmarkindices.append(refmarkindex)

        try:
            self.num_nodes[id].setParent(None)
        except:
            pass

        if values:
            self.valuesBatchSubtomoReconstruction = values
            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, nn=True)
            self.tab12_widgets = self.tables[id].widgets
            for n, index in enumerate(refmarkindices):
                self.tab12_widgets['widget_{}_3'.format(n)].setCurrentIndex(index)

            for row in range(len(values)):
                w = self.tab12_widgets['widget_{}_2'.format(row)]
                w.currentIndexChanged.connect(lambda d, r=row, v=values: self.update_choices(r, v))

            self.particleFilesBatchExtract = particleFiles
            self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.mass_extract_particles(pid, v))
        else:
            return

        for i in range(len(values)):
            self.update_choices(i, values)

    def update_choices(self, rowID, values):
        if 1:
            values = self.valuesBatchSubtomoReconstruction
            self.tab12_widgets['widget_{}_3'.format(rowID)].clear()

            current_index = self.tab12_widgets['widget_{}_2'.format(rowID)].currentIndex()

            origin = values[rowID][2][current_index]
            particleFile = values[rowID][0]
            a = 'tomogram_' + particleFile.split('_tomogram_')[1][:3]

            folder = os.path.join(self.tomogram_folder, a, origin)
            print(folder)
            choices = [folder + '/' + f for f in os.listdir(folder) if
                        'unweighted_unbinned' in f and os.path.isdir(folder + '/' + f)]

            closest_choices = {}

            for choice in choices:
                f = os.path.dirname(choice)
                try: a = choice.split('unweighted_unbinned_marker_')[1]
                except: continue
                if 'reduced' in a:
                    print(a)
                    try:
                        aS, aE = a.split('_')[2:4]
                        angleS = '{:4.1f}'.format( float(aS) )
                        angleE = '{:4.1f}'.format( float(aE) )

                        if 'ctf' in a: ctf = '_ctf'
                        else: ctf = ''


                        key = '{}/unweighted_unbinned_marker_CLOSEST_reduced_{}_{}{}'.format(f, angleS, angleE, ctf)
                        value = '{}/unweighted_unbinned_marker_CLOSEST_reduced_{}_{}{}'.format(f, aS, aE, ctf)


                        closest_choices[key] = value
                    except:
                        pass

                elif 'unweighted_unbinned_marker_' in a:
                    closest_choices['{}/unweighted_unbinned_marker_CLOSEST'.format(f)] = 1
                elif 'sorted' in a:
                    closest_choices[choice] = 1

            choices  = list(closest_choices.values()) + choices

            self.valuesBatchSubtomoReconstruction[rowID][3] = choices

            for item in choices:
                self.tab12_widgets['widget_{}_3'.format(rowID)].addItem(os.path.basename(item))

        else:#except Exception as e:
            print('Update choices has failed.')
            print(e)
            return

    def mass_extract_particles(self,pid, values):
        num_nodes = int(self.num_nodes[pid].value())
        num_submitted_jobs = nsj = 0
        values = self.valuesBatchSubtomoReconstruction
        for row in range(self.tables[pid].table.rowCount()):
            if self.tab12_widgets['widget_{}_1'.format(row)].isChecked():
                particleXML = values[row][0] #[row][0]
                tomoindex = particleXML.split('_tomogram_')[-1][:3]
                origin = values[row][2][self.tab12_widgets['widget_{}_{}'.format(row,2)].currentIndex()]
                folder_aligned = values[row][3][self.tab12_widgets['widget_{}_{}'.format(row,3)].currentIndex()]
                metafile = glob.glob('{}/03_Tomographic_Reconstruction/tomogram_{}/sorted/*.meta'.format(self.projectname,tomoindex))
                if not metafile: continue
                metafile = metafile[0]
                #q = '{}/03_Tomographic_Reconstruction/tomogram_{}/alignment/unweighted_unbinned_marker_{}'
                #folder_aligned = q.format(self.projectname,tomoindex,ref_marker)
                bin_read = self.tab12_widgets['widget_{}_{}'.format(row,4)].text()
                weight = self.tab12_widgets['widget_{}_{}'.format(row, 5)].text()
                size = self.tab12_widgets['widget_{}_{}'.format(row, 6)].text()
                bin_subtomo = self.tab12_widgets['widget_{}_{}'.format(row,7)].text()
                offx = self.tab12_widgets['widget_{}_{}'.format(row, 8)].text()
                offy = self.tab12_widgets['widget_{}_{}'.format(row, 9)].text()
                offz = self.tab12_widgets['widget_{}_{}'.format(row, 10)].text()


                refid = folder_aligned.split('unweighted_unbinned_marker_')[1].split('_')[0]

                if refid.lower() != 'closest':

                    outname = 'Reconstruction/reconstruct_subtomograms_{:03d}_refmarker_{}.sh'.format(int(tomoindex),refid)
                    execfilename = os.path.join(self.subtomodir, outname)

                    paramsCmd = [particleXML, folder_aligned, bin_read, size, bin_subtomo, offx, offy, offz,
                                 self.subtomodir, weight, metafile, '20']

                    txt = extractParticles.format(d=paramsCmd)
                    jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder, name=os.path.basename(outname[:-3]),
                                                           num_jobs_per_node=20, time=12) + txt

                else:
                    outname = 'Reconstruction/reconstruct_subtomograms_{:03d}_refmarker_{}.sh'
                    outname = outname.format(int(tomoindex), refid)
                    execfilename = os.path.join(self.subtomodir, outname)

                    tomodir = os.path.dirname(os.path.dirname(metafile))

                    reconAlg = None
                    for alg in ('WBP', 'INFR'):
                        if alg in particleXML:
                            reconAlg = alg

                    if not reconAlg:
                        print('FAIL: subtomogram reconstruction for {} failed. No INFR or WBP in xml path to particle.')
                        continue

                    end = 'reconstruction/{}/marker_locations_tomogram_{}_{}_irefmark_*.txt'
                    end = end.format(reconAlg, tomoindex, reconAlg)

                    logfilequery = os.path.join(tomodir, end)
                    logfile = sorted(glob.glob(logfilequery))[0]
                    qname, n_nodes, cores, time = self.qparams['BatchSubtomoReconstruct'].values()

                    paramsCmd = [particleXML, folder_aligned, bin_read, size, bin_subtomo, offx, offy, offz,
                                 self.subtomodir, weight, metafile, logfile, str(cores), 'sorted_ctf_aligned']


                    txt = extractParticlesClosestMarker.format(d=paramsCmd)

                    jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder,singleton=True, num_nodes=n_nodes,
                                                           name='SubtomoRecon_{}'.format(nsj % num_nodes),
                                                           num_jobs_per_node=20, time=time) + txt
                out = open(execfilename, 'w')
                out.write(jobtxt)
                out.close()
                os.system('{} {}'.format(self.qcommand, execfilename))
                nsj += 1

    def inputFiles(self, mode=None):
        title = "FRM Alignment"
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
                                    tooltip='Select the mask file.', mode='file',filetype=['em','mrc'],cstep=1,rstep=0)
        self.insert_pushbutton(parent,'Create',rstep=1,cstep=-3,action=self.gen_mask,params=[mode+'filenameMask'])
        self.insert_label_line_push(parent, 'Filename Average', mode + 'filenameAverage',initdir=self.frmdir,
                                    tooltip='Choose a filename for the average of all particles.', mode='file',
                                    filetype=['em','mrc'],cstep=1,rstep=0)
        self.insert_pushbutton(parent, 'Average', rstep=1, cstep=-3, action=self.gen_average,
                               params=[mode + 'particleList', mode + 'filenameAverage', mode + 'outputDir'])

        self.insert_label_line_push(parent, 'Output Directory', mode + 'outputDir',
                                    'Folder in which the output will be written.')
        self.insert_label(parent,rstep=1,cstep=0)
        self.insert_label_spinbox(parent, mode + 'bwMin', 'Min Order SH Polynomial',
                                  value=8,minimum=0,stepsize=1,
                                  tooltip='The minimal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'bwMax', 'Max Order SH Polynomial',
                                  value=64, minimum=0, stepsize=1,
                                  tooltip='The maximal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'frequency', 'Frequency (px)',
                                  value=8, minimum=0, stepsize=1,
                                  tooltip='The minimal frequency used for reconstruction.')
        self.insert_label_spinbox(parent,  mode + 'maxIterations', 'Maximum Iterations',
                                  value=8, minimum=1, stepsize=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
        self.insert_label_spinbox(parent, mode + 'peakOffset', 'Peak Offset',
                                  value=0, minimum=0, stepsize=1,
                                   tooltip='Sets the peak offset.')
        self.insert_label(parent, rstep=1, cstep=0)
        self.insert_label_spinbox(parent,mode+'pixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox,minimum=0.1,stepsize=0.1,value=1.75)
        self.insert_label_spinbox(parent,mode+ 'particleDiameter','Particle Diameter (A)',rstep=1,cstep=0,
                                  minimum=10, stepsize=1, value=300, maximum=10000, width=150)

        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updateFRM(m))

        rscore = 'False'
        weightedAv = 'False'
        weighting = ''
        binning_mask = '1'
        sphere = 'True'
        ad_res = '0.00'
        fsc = '0.50'

        jobfilename = [mode + 'outputDir', 'job_description.xml']#os.path.join(self.frmdir, 'job_description.xml')
        exefilename = [mode + 'outputDir', 'frmAlignment.sh'] #os.path.join(self.frmdir, 'frmAlignment.sh')

        paramsSbatch = guiFunctions.createGenericDict(fname='FRMAlign', folder=self.logfolder, id='FRMAlignment') #, modules=['openmpi/2.1.1', 'python/2.7', 'lib64/append', 'pytom/dev/gui'])
        paramsJob = [mode+'bwMin',mode+'bwMax',mode+'frequency',mode+'maxIterations', mode+'peakOffset',
                     rscore, weightedAv, mode+'filenameAverage', weighting, mode+'filenameMask', binning_mask, sphere,
                     mode+'pixelSize', mode+'particleDiameter', mode+'particleList', mode+'outputDir']
        paramsCmd = [self.subtomodir, self.pytompath, jobfilename, templateFRMSlurm]

        self.insert_gen_text_exe(parent, self.stage, xmlfilename=jobfilename, jobfield=True, exefilename=exefilename,
                                 paramsXML=paramsJob + [templateFRMJob], paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def updateFRM(self,mode):
        item = self.widgets[mode + 'particleList'].text()
        if not item:
            return
        folder, ext = os.path.splitext( os.path.basename(item))
        outputDir = os.path.join(self.frmdir, folder.replace('particleList_', ''))
        if not os.path.exists(outputDir): os.mkdir(outputDir)
        else:
            for i in range(1000):
                if not os.path.exists(outputDir+'_{:03d}'.format(i)):
                    outputDir += '_{:03d}'.format(i)
                    os.mkdir(outputDir)
                    break


        self.widgets[mode+'outputDir'].setText(outputDir)

    def gen_average(self, params):
        print(params)
        key_particleList, key_filename_average, key_outputDir = params
        particleList = self.widgets[key_particleList].text()
        if not particleList:
            self.popup_messagebox('Warning', 'Averaging Failed', 'Averaging did not succeed. No particle list selected.')
            return
        folder = self.widgets[key_outputDir].text()
        if not os.path.exists(folder): os.mkdir(folder)
        output_name = os.path.join( folder, 'average.em')
        out = os.popen('cd {}; average.py -p {} -a {} '.format(self.subtomodir, particleList, output_name)).read()
        if not os.path.exists(output_name):
            self.popup_messagebox('Warning', 'Averaging Failed',
                                  'Averaging did not succeed, please try again.')
            return
        self.widgets[key_filename_average].setText(output_name)

    def gen_mask(self,params):
        maskfilename = CreateMaskFile(self, maskfname=params[-1])
        maskfilename.show()

    def referenceUpdate(self, mode):
        if self.widgets[mode + 'referenceModel'].text():
            self.widgets['referenceCommand'].setText( "--reference " + self.widgets[mode + 'referenceModel'].text())
        else:
            self.widgets['referenceCommand'].setText("")

    def glocal(self,mode):
        title = "GLocal alignment"
        tooltip = 'Run pytom GLocal routine.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Initial reference model', mode + 'referenceModel', mode='file',
                                    filetype=['em','mrc'], enabled=True,
                                    tooltip='Reference : the initial reference - if none provided average of particle list')
        self.widgets[mode + 'referenceModel'].textChanged.connect(lambda dummy,mode=mode: self.referenceUpdate(mode))
        self.widgets['referenceCommand'] = QLineEdit(self)
        self.widgets['referenceCommand'].setVisible(False)
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask'])
        self.insert_label_line_push(parent, 'Output Directory', mode + 'destination', mode='folder',
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
                                  minimum=1, stepsize=1, value=3, maximum=359,
                                  rstep=1, cstep=-1, tooltip='Angular increment for refinement.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=0,
                                  stepsize=1,minimum=1,value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')

        self.widgets[mode+'jobName'] = QLineEdit()
        self.widgets[mode + 'destination'].textChanged.connect(lambda d, m=mode: self.update_jobname(m))
        self.update_jobname(mode)
        glocalpath = os.path.join(self.subtomodir, 'Alignment/GLocal')
        exefilename = os.path.join(glocalpath, 'GLocal_Alignment.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='GLocal', folder=self.logfolder, id='GLocalAlignment')
        paramsCmd = [self.subtomodir, self.pytompath, self.pytompath, mode+'particleList', 'referenceCommand',
                     mode+'filenameMask', mode+'numIterations', mode+'pixelSize', mode+'particleDiameter',
                     mode+'binning', mode+'jobName', mode+'destination', mode + 'angleShells',
                     mode + 'angleIncrement', templateGLocal]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_GLocal', groupbox)
        return groupbox

    def update_jobname(self, mode):
        dest = self.widgets[mode+'destination'].text()
        if dest:
            self.widgets[mode + 'jobName'].setText( os.path.join(dest, 'glocal_results_'+os.path.basename(dest)+'.xml'))

    def CPCA(self,mode=''):
        title = "Classify CPCA"
        tooltip = 'The CCC is further used for classification. This script computes \n'+\
                  'the eigenvectors of the CCC and projects the data on the first \n'+\
                  'neig eigenvectors. Subsequently, these multidimensional vectors \n'+\
                  'are clustered into nclass groups using a kmeans algorithm.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_line(parent, 'Filename Output particleList', mode + 'outputFilename',
                               tooltip='Filename for generated XML file that includes the assigned classes for each particle. No full path needed.')
        self.insert_label_line_push(parent, 'CCC File', mode + 'cccFile',
                                    'Select the constrained correlation matrix file from the previous step.', mode='file', filetype='csv')
        self.insert_label_spinbox(parent, mode + 'numEig', text='Number of Eigenvectors',
                                  value=4, minimum=1, stepsize=1,
                                  tooltip='Sets the number of eigenvectors (corresponding to largest eigenvectors) used for clustering.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes',
                                  value=4, minimum=1,stepsize=1,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_line(parent, 'Prefix', mode + 'prefix', rstep=1, cstep=0,
                               tooltip='Root for generated averages of the corresponding classes. The files will be called "Prefix"_iclass.em.')

        self.widgets[mode + 'particleList'].textChanged.connect(
            lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode, p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        exefilename = [mode + 'outFolder', 'CPCA_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CPCA', folder=self.logfolder, id='CPCA')
        paramsCmd = [mode+'outFolder', self.pytompath, mode + 'particleList', mode + 'outputFilename',
                     mode + 'cccFile', mode + 'numEig', mode+'numClasses', mode+'prefix',  templateCPCA]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_CPCA', groupbox)
        return groupbox

    def CCC(self,mode=''):
        title = "Pairwise Constrained Cross Correlation"
        tooltip = 'Calculate the pairwise constrained cross correlation.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask'])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'lowpass', 'Lowpass filter (px)',
                                  minimum=0, maximum=1024, stepsize=1, value=20,
                                  tooltip='The lowpass filter is applied to all subtomograms after binning.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=0,
                                  minimum=1, stepsize=1, value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')

        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode,p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        exefilename = [mode + 'outFolder', 'CCC_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CCC_Class', folder=self.logfolder,
                                                      id='PairwiseCrossCorrelation')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'filenameMask',
                     mode + 'lowpass', mode + 'binning', mode + 'outFolder', templateCCC]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_CCC', groupbox)
        return groupbox

    def ac(self,mode):
        title = "Autofocussed Classification"
        tooltip = 'Run autofocussed classification.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Classification Mask', mode + 'filenameClassificationMask', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='This mask is used for constraining the calculation of the focused mask. (Optional)', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameClassificationMask'])
        self.insert_label_line_push(parent, 'Alignment Mask', mode + 'filenameAlignmentMask', mode='file', filetype=['em', 'mrc'], enabled=True,
                                    tooltip='This mask is only used for the alignment purpose. Only specify it if the particle list is not aligned.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameAlignmentMask'])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes', stepsize=1, value=4, minimum=2,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_spinbox(parent, mode + 'maxIterations', 'Number Of Iterations',stepsize=1,value=10, minimum=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
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
                                  tooltip='Particle density threshold for calculating the difference map (optional, by default 0). Two other most common choise are -2 and 2. -2 means all the values of the subtomogram below the -2 sigma will be used for calculating the difference mask (negative values count). 2 means all the values of the subtomogram above the 2 sigma will be used for calculating the difference mask (positive values count). Finally, 0 means all the values are used for the calculation.')
        self.insert_label_spinbox(parent, mode + 'stdDiffMap', 'STD Threshold Diff Map', rstep=1, cstep=0,
                                  wtype=QDoubleSpinBox, stepsize=.1, minimum=0, maximum=1, value=0.4,
                                  tooltip='STD threshold for the difference map (optional, by default 0.4). This value should be between 0 and 1. 1 means only the place with the peak value will be set to 1 in the difference map (too much discriminative ability). 0 means all the places with the value above the average of STD will be set to 1 (not enough discriminative ability).')

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

        # Parameters for execution
        exefilename = [mode + 'outFolder', 'AC_Classification.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='AutoFocus', folder=self.logfolder,
                                                      id='AutoFocusClassification')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'flagAlignmentMask',
                     mode + 'flagClassificationMask', mode + 'numClasses', mode + 'bwMax', mode + 'maxIterations',
                     mode + 'peakOffset', mode + 'noisePercentage', mode + 'partDensThresh', mode + 'stdDiffMap',
                     mode + 'outFolder', templateAC]


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        # Run Update With Data From Logfile
        self.updateAlignmentMaskFlag(mode)
        self.updateClassificationMaskFlag(mode)
        self.updateOutFolder(mode, self.acpath)

        setattr(self, mode + 'gb_AC', groupbox)
        return groupbox

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
        if not pl: return

        folder = os.path.basename(pl).replace('particleList_','')[:-4]
        folder = os.path.join(path,folder)
        if os.path.exists(folder):
            folder = folder
            if not folder: return

        if not os.path.exists(folder): 
            os.mkdir(folder)
            #os.system('ln -s {}/Subtomograms {}/Subtomograms'.format(self.subtomodir, folder ) )
        self.widgets[mode + 'outFolder'].setText(folder) 

    def createOutFolder(self, mode):
        folder = self.widgets[mode + 'outFolder'].text()
        if not folder: return
        if not os.path.exists(os.path.join(folder, 'Subtomograms')):
            os.system('ln -s {}/Subtomograms {}/Subtomograms'.format(self.subtomodir, folder ) )
