import os
import copy
import pickle
import glob
import pyqtgraph as pg
import multiprocessing
import getpass
import re

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from pytom.basic.files import pdb2mrc
from pytom.gui.mrcOperations import read_mrc
from pytom.gui.guiStyleSheets import MAINC, BARS, DEFAULT_STYLE_PROGRESSBAR, WHITE
from pytom.gui.guiFunctions import initSphere
from pytom.gui import guiFunctions
from pytom.gui.guiSupportCommands import templateConvertData
import numpy as np

from scipy.ndimage.filters import gaussian_filter
from ftplib import FTP_TLS
import lxml.etree as et

from multiprocessing import Manager, Event


colorsHEX = ['00ffff', 'f5f5dc', '0000ff', 'a52a2a', '7fff00', 'd2691e', 'daa520', 'ff7f50', '00ffff',
                     'dc143c', '00008b', '006400', '7fffd4', 'ff00ff', 'ffd700', '008000', '4b0082', 'f0e68c', 'add8e6',
                     '90ee90', 'ff00ff', '800000', '808000', 'ffa500', 'ff4500', 'da70d6', 'ffc0cb', '800080', 'dda0dd',
                     'fa8072', 'fa8072', '008080', 'ffff00', '40e0d0', '9acd32', ] * 15

html = '<div style="text-align: center"><span style="color: #{}; font-size: 20pt;">{}</span></div>'


class BrowseWindowRemote(QMainWindow):
    '''This class creates a new window for remote file browsing'''
    def __init__(self, parent=None, initdir='/',filter=[''],search='file',credentials=['','',''],outputline='',
                 validate=True):
        super(BrowseWindowRemote, self).__init__(parent)
        self.setGeometry(50, 50, 800, 300)
        self.pathdisplay = QLineEdit()
        self.pathdisplay.setStyleSheet(WHITE)
        self.pathdisplay.setEnabled(False)
        self.splitter0 = splitter0 = QSplitter(Qt.Vertical)
        self.topleft = QListWidget()
        self.topleft.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.topleft.setStyleSheet(WHITE)
        self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)
        self.topright = QListWidget()
        self.topright.setStyleSheet(WHITE)
        if search == 'file': self.topright.itemDoubleClicked.connect(self.select_file)
        self.splitter1 = splitter1 = QSplitter(Qt.Horizontal)
        splitter1.addWidget(self.topleft)
        splitter1.addWidget(self.topright)
        splitter1.setSizes([300,500])

        bottom_layout = QHBoxLayout()
        ok_button = QPushButton(text='OK')
        ok_button.clicked.connect(self.select_item)
        cancel_button = QPushButton(text='Cancel')
        cancel_button.clicked.connect(self.cancel)

        bottom = QWidget()
        bottom_layout.addWidget(ok_button)
        bottom_layout.addWidget(cancel_button)
        bottom_layout.addStretch(1)
        bottom.setLayout(bottom_layout)

        splitter0.addWidget(self.pathdisplay)
        splitter0.addWidget(splitter1)
        splitter0.addWidget(bottom)

        self.success = False
        self.ftps=None
        self.initdir = initdir
        self.folderpath = initdir
        self.filters = filter
        self.search = search
        self.outputline = outputline

        try:
            self.servername = str(credentials[0])
            self.username   = str(credentials[1])
            self.password   = str(credentials[2])
        except:
            return


        #self.servername,self.username,self.password = 'emsquare1.science.uu.nl','emuser','#99@3584cg'

        #print(self.servername,self.username,self.password)

        try:
            self.connect_ftp_server(self.servername, self.username, self.password)
            self.setCentralWidget(splitter0)
            self.success = True
            self.add_folders()

        except:
            if validate:
                QMessageBox().critical(self, "Credentials not valid.",
                                       "Please provide a valid combination of servername, username, and password",
                                                 QMessageBox.Ok)
                #self.destroy(True,True)
                self.close()


    def select_item(self):
        if self.outputline: self.outputline.setText(self.folderpath)
        self.close()

    def select_file(self):
        self.folderpath = os.path.join(self.folderpath, self.topright.currentItem().text())
        self.outputline.setText(self.folderpath)
        self.close()

    def cancel(self):
        #self.pathselect = self.initdir
        self.close()

    def connect_ftp_server(self,servername, username, password):
        self.ftps = FTP_TLS(servername,username,password)
        self.ftps.prot_p()

    def append(self,line):
        words = line.split(None, 3)
        datestamp = words[0]
        timestamp = words[1]
        folder = words[2]
        filename = words[3]
        size = None
        if folder == '<DIR>':
            folder = True
        else:
            size = int(folder)
            folder = False
        if folder:
            self.subdirs.append(filename)
        else:
            self.matchingfiles.append(filename)

    def add_folders(self, folders=[]):
        self.data()

        self.topleft.clear()
        self.topright.clear()

        self.topleft.insertItems(0,self.subdirs)
        self.topright.insertItems(0,self.matchingfiles)

        #if self.search == 'file': self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)

    def data(self):
        self.subdirs = [os.pardir]
        self.matchingfiles = []


        self.ftps.cwd(self.folderpath)
        self.ftps.retrlines('LIST', self.append)

        matchingfiles = []
        for pat in self.filters:
            matchingfiles += [line for line in self.matchingfiles if line.endswith(pat.split('.')[-1])]
        self.matchingfiles = matchingfiles

    def repopulate_folder_list(self):

        extra = self.topleft.currentItem().text()
        if extra != '..':
            self.folderpath = os.path.join(self.folderpath, extra)

        elif len(self.folderpath) >1:
            self.folderpath = os.path.dirname(self.folderpath)
        self.pathdisplay.setText(self.folderpath)
        self.add_folders()


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:

    finished
        No data

    error
        `tuple` (exctype, value, traceback.format_exc() )

    result
        `object` data returned from processing, anything

    '''

    finished_queue = pyqtSignal(object)
    finished_mcor = pyqtSignal()
    finished_collect = pyqtSignal()
    error = pyqtSignal(tuple)
    result1 = pyqtSignal(object)
    result2 = pyqtSignal(object)
    success = pyqtSignal()
    results = pyqtSignal(object)
    startMessage = pyqtSignal(object)


class Worker(QRunnable):

    def __init__(self, fn=None, args=[], sig=True, results=False):
        super(Worker, self).__init__()
        # QRunnable.__init__(self)
        #QObject.__init__(self)

        #super(Worker, self).__init__()
        self.signals = WorkerSignals()
        self.fn = fn
        self.args = tuple(list(args))
        self.finishedJob = False
        self.results = results

        if sig: self.args = tuple(list(args)+[self.signals])


    def run(self):
        result = self.fn(*self.args)
        self.finishedJob = True
        if self.results:
            self.signals.results.emit(result)

    def start(self):
        QThreadPool.globalInstance().start(self)


class TimedMessageBox(QMessageBox):
    '''This class creates a Message Box that will disappear after time out is reached'''
    def __init__(self, parent=None, timeout=3, info=None, type=''):
        super(TimedMessageBox, self).__init__(parent)

        from PyQt5.QtCore import QTimer

        self.setWindowTitle(info[1])
        self.message = info[2]
        self.time_to_wait = timeout
        self.setText("{1}\n(closing automatically in {0} seconds.)".format(self.time_to_wait,self.message))
        self.setStandardButtons(info[3])
        self.timer = QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.changeContent)
        # if type == 'TimedInfo':
        #     self.information(*info)
        # elif type == 'TimedWarning':
        #     self.warning(*info)
        self.timer.start()
        self.exec_()

    def changeContent(self):
        self.setText("{1}\n(closing automatically in {0} seconds.)".format(self.time_to_wait,self.message))
        if self.time_to_wait <= 0:
            self.close()
        self.time_to_wait -= 1

    def closeEvent(self, event):
        print('closing window')
        self.timer.stop()

        event.accept()


class CommonFunctions():
    '''Class with commonly used functions'''
    def __init__(self,parent=None):
        pass

    def size_policies(self):

        self.sizePolicyA = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.sizePolicyA.setHorizontalStretch(1)
        self.sizePolicyA.setVerticalStretch(1)

        self.sizePolicyB = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.sizePolicyB.setHorizontalStretch(1)
        self.sizePolicyB.setVerticalStretch(0)

        self.sizePolicyC = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sizePolicyC.setHorizontalStretch(0)
        self.sizePolicyC.setVerticalStretch(0)

    def add_toolbar(self, decider, new=False, open=True, save=True,
                    newText='Create New Project', openText='Load Project', saveText="Save Project"):
        self.setStyleSheet(MAINC)
        bar = self.menuBar()
        tb = QToolBar()
        tb.setStyleSheet(BARS)
        self.addToolBar(tb)

        if new:
            new = QAction(QIcon("{}/gui/Icons/new_project4.png".format(self.pytompath)), newText, self)
            tb.addAction(new)

        if open:
            load = QAction(QIcon("{}/gui/Icons/OpenProject.png".format(self.pytompath)), openText, self)
            tb.addAction(load)

        if save:
            s = QAction(QIcon("{}/gui/Icons/SaveProject.png".format(self.pytompath)), saveText, self)
            tb.addAction(s)
            #tb.actionTriggered[QAction].connect(self.save_particleList)

        tb.actionTriggered[QAction].connect(decider)

    def insert_slider(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0, tickinterval=1,
                        alignment=Qt.AlignCenter, tooltip='', logvar=False, width=0, value=15, minimum=0, maximum=50):

        widget = QSlider(Qt.Horizontal)
        widget.setStyleSheet("QSlider{background:none;}")
        widget.setMinimum(minimum)
        widget.setMaximum(maximum)
        widget.setValue(value)
        widget.setTickPosition(QSlider.TicksBelow)
        widget.setTickInterval(tickinterval)

        if tooltip: widget.setToolTip(tooltip)
        if width: widget.setFixedWidth(width)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan, alignment)
        setattr(self, wname, widget)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_checkbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,
                        alignment=Qt.AlignCenter, tooltip='', logvar=False,width=0):

        widget = QCheckBox(self)
        widget.setStyleSheet("QCheckBox{background:none;}")
        if logvar:
            if wname in self.logbook.keys():
                widget.setChecked(self.logbook[wname])
            else:
                self.logbook[wname] = False
            widget.setObjectName(wname)
            widget.stateChanged.connect(lambda dummy, name=wname, text='checkbox':
                                        self.update_logbook(name=name, text=text))

        if text:
            widget.setText(text)
        else:
            sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
        if tooltip: widget.setToolTip(tooltip)
        if width: widget.setFixedWidth(width)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan, alignment)
        setattr(self, wname, widget)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label(self, parent, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,transparent=True,
                     tooltip='', alignment=Qt.AlignLeft, value=None, width=0, sizepolicy=''):
        widget = QLabel(self)

        if transparent:
            widget.setStyleSheet("QLabel{background:transparent;} \nQToolTip { background-color: white; }")

        if text:
            widget.setText(text)
        elif width == 0:
            sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
        if tooltip: widget.setToolTip(tooltip)
        if sizepolicy:widget.setSizePolicy(sizepolicy)
        if width: widget.setFixedWidth(width)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan, alignment)
        if alignment: widget.setAlignment(alignment)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_combobox(self, parent, wname, labels, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,
                        width=80, tooltip='', logvar=False):
        widget = QComboBox(self)
        widget.addItems(labels)

        if logvar:
            if wname in self.logbook.keys():
                for n, label in enumerate(labels):
                    if label == self.logbook[wname]:
                        widget.setCurrentIndex(n)
            else:
                self.logbook[wname] = labels[0]
            widget.setObjectName(wname)
            widget.currentIndexChanged.connect(lambda dummy, name=wname, text='combobox':
                                               self.update_logbook(name=name, text=text))

        if tooltip: widget.setToolTip(tooltip)
        widget.setStyleSheet("QComboBox{background:white; selection-background-color: #1989ac;}")
        widget.setFixedWidth(width)
        sizePolicy = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
        widget.setSizePolicy(sizePolicy)

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget

        self.row += rstep
        self.column += cstep

    def insert_module(self, parent, wname='', cstep=0, rstep=1, rowspan=1, columnspan=1, options=[], mode=''):
        widget = SelectModules(self, modules=options, mode=mode)
        if wname: self.widgets[wname] = widget
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label_modules(self, parent, wname, text='', rstep=1, cstep=-1, width=0, tooltip='', height=0,
                             logvar=False, options=[], mode=''):
        self.insert_label(parent, text=text, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_module(parent, wname=wname, cstep=cstep, rstep=rstep, options=options, mode=mode)

    def insert_pgSpinbox(self, parent, wname, text='',rowspan=1, columnspan=1, rstep=0, cstep=0, decimals=0,width=0,
                         height=0, tooltip='', logvar=False, value=None, maximum=None, minimum=0, stepsize=1,suffix=''):

        widget = pg.SpinBox(parent=parent)
        self.widgets[wname] = widget

        if maximum: widget.setMaximum(maximum)
        if minimum: widget.setMinimum(minimum)
        if suffix: widget.setSuffix(suffix)
        if stepsize: widget.setSingleStep(stepsize)
        if value: widget.setValue(value)
        if decimals>=0: widget.setDecimals(decimals)
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)

        if logvar:
            try:
                widget.setText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.valueChanged.connect(lambda ignore, name=wname, text='lineedit':
                                       self.update_logbook(name=name, text=text))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_spinbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,width=0, height=0,
                       validator=None, password=False, tooltip='', logvar=False, value=None, att=False, enabled=True,
                       maximum=0, minimum=0, stepsize=0, widgetType=QSpinBox, decimals=0):

        widget = widgetType(self)
        widget.setLocale(QLocale('.'))
        self.widgets[wname] = widget
        widget.setEnabled(enabled)
        widget.setAutoFillBackground(True)
        widget.setStyleSheet("background-color: rgb(255,255,255);")

        widget.setPalette(QPalette(Qt.white ) )
        widget.setKeyboardTracking(False)

        if logvar:
            try:
                widget.setValue(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.valueChanged.connect(lambda ignore, name=wname, text='spinbox':
                                       self.update_logbook(name=name, text=text))

        if password: widget.setEchoMode(QLineEdit.Password)
        if validator: widget.setValidator(validator)
        if tooltip: widget.setToolTip(tooltip)
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)
        if minimum: widget.setMinimum(minimum)
        if maximum: widget.setMaximum(maximum)
        if stepsize: widget.setSingleStep(stepsize)
        if value: widget.setValue(value)
        if att and not logvar: widget.textChanged.connect(lambda ignore, name=wname: self.set_attribute(wname))
        if decimals:
            widget.setDecimals(decimals)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_lineedit(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,width=0, height=0,
                        validator=None, password=False, tooltip='', logvar=False, value=None, att=False, enabled=True):
        widget = QLineEdit(self)

        if logvar:
            try:
                widget.setText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.textChanged.connect(lambda ignore, name=wname, text='lineedit':
                                       self.update_logbook(name=name, text=text))
        if password: widget.setEchoMode(QLineEdit.Password)
        if validator: widget.setValidator(validator)
        if tooltip: widget.setToolTip(tooltip)
        if att and not logvar:
            widget.textChanged.connect(lambda ignore, name=wname:
                                       self.set_attribute(wname))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)


        widget.setStyleSheet("QLineEdit{background:white;}; QLineEdit:disabled{background:white; color:black;} ")
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)
        if value: widget.setText(str(value))
        widget.setEnabled(enabled)
        self.row += rstep
        self.column += cstep

    def set_attribute(self, wname):
        txt = self.widgets[wname].text()

        if wname[-1] == '1':
            wname = wname[:-1]
            setattr(self, wname, txt)
            self.widgets[wname].setText(txt)
        else:
            setattr(self, wname, txt)
            wname = wname + '1'
            self.widgets[wname].setText(txt)

    def insert_pushbutton(self, parent, text='', rowspan=1, columnspan=1, rstep=0, cstep=0, action=None,
                          tooltip='', width=0, params=['file', '', '', False], iconpath='', wname='', state=True):
        widget = QPushButton(self)

        if wname:
            self.widgets[wname] = widget
        widget.setEnabled(state)
        if tooltip: widget.setToolTip(tooltip)
        if text:
            widget.setText(text)
        elif width==0:
            sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)

        if action:
            if params:
                widget.clicked.connect(lambda ignore, p=params: action(p))
            else:
                widget.clicked.connect(lambda ignore: action())
        if iconpath:
            icon = QPixmap(iconpath)
            widget.setIcon(icon)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        if width: widget.setFixedWidth(width)
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_textfield(self, parent, wname, width=500, height=400, logvar=True,
                         rstep=0, cstep=0, rowspan=1, columnspan=1):

        widget = QPlainTextEdit(self)
        widget.setSizePolicy(self.sizePolicyB)
        widget.setStyleSheet("QPlainTextEdit{background:white;}")
        widget.setObjectName(wname)


        if logvar:
            try:
                widget.setPlainText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            #widget.textChanged.connect(lambda ignore, name=wname, text='textfield':
            #                           self.update_logbook(name=name, text=text))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label_spinbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=1, cstep=-1, width=0, height=0,
                             validator=None, tooltip='', value=None, att=False, enabled=True, maximum=0, minimum=0,
                             wtype=QSpinBox, stepsize=0, logvar=True, decimals=0):

        self.insert_label(parent, text=text, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_spinbox(parent, wname, validator=validator, width=width, enabled=enabled,
                            maximum=maximum, minimum=minimum, cstep=cstep, rstep=rstep, value=value,
                            widgetType=wtype, stepsize=stepsize, logvar=logvar, decimals=decimals)

    def insert_label_line_push(self, parent, textlabel, wname, tooltip='', text='', cstep=-2, rstep=1, validator=None,
                               mode='folder', remote=False, pushtext='Browse', width=150, filetype='',action='',
                               initdir='', enabled=False, action2=None, pushtext2='# Particles'):

        if action == '':
            action = self.browse


        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_lineedit(parent, wname, cstep=1, logvar=True, text='',validator=validator,width=width,enabled=enabled)

        if initdir:
            params = [mode, self.widgets[wname], filetype, remote, initdir]
        else:
            params = [mode, self.widgets[wname], filetype, remote]

        if not action2 is None:
            self.insert_pushbutton(parent, cstep=1, rstep=0, text=pushtext,width=100,
                                   action=action, params=params)
            self.insert_pushbutton(parent, cstep=cstep, rstep=rstep, text=pushtext2, width=100,
                                   action=action2, params=wname)
        else:
            self.insert_pushbutton(parent, cstep=cstep, rstep=rstep, text=pushtext,width=100,
                                   action=action, params=params)

    def insert_label_push(self, parent, textlabel, wname, tooltip='', cstep=-1, rstep=1, pushtext='Surprise Me!', action='',
                          width=100, params=[]):

        if action == '':
            return
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_pushbutton(parent, cstep=cstep, rstep=rstep, text=pushtext, width=width, action=action, params=params, wname=wname)

    def insert_label_line(self, parent, textlabel, wname, tooltip='', value='', cstep=-1, rstep=1, validator=None,
                          width=150,logvar=True, enabled=True):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_lineedit(parent, wname, cstep=cstep, rstep=rstep, value=value, logvar=logvar, validator=validator,
                             width=width, enabled=enabled)

    def insert_label_checkbox(self, parent, wname, textlabel, tooltip='', cstep=-1, rstep=1, logvar=True):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_checkbox(parent, wname, logvar=logvar, cstep=cstep, rstep=rstep,alignment=Qt.AlignLeft)

    def insert_checkbox_label(self, parent, wname, textlabel, tooltip='', cstep=-1, rstep=1, logvar=True, width=200):
        self.insert_checkbox(parent, wname, logvar=logvar, alignment=Qt.AlignLeft, cstep=1)
        self.insert_label(parent, text=textlabel, cstep=cstep, rstep=rstep, alignment=Qt.AlignRight, tooltip=tooltip,
                          width=width)

    def insert_checkbox_label_line(self, parent, wname, textlabel, wname2, tooltip='', cstep=-2, rstep=1, logvar=True,
                                   width=150, labelwidth=200, validator=None, value='', enabled=False, sizepolicy=None,
                                   alignment=None):
        self.insert_checkbox(parent, wname, logvar=logvar, alignment=Qt.AlignLeft, cstep=1)
        self.insert_label(parent, text=textlabel, cstep=1, tooltip=tooltip, width=labelwidth,
                          sizepolicy=sizepolicy, alignment=Qt.AlignRight)
        self.insert_lineedit(parent, wname2, cstep=cstep, rstep=rstep, value=value, logvar=logvar, validator=validator,
                             width=width, enabled=enabled)

    def insert_checkbox_label_spinbox(self, parent, wname, textlabel, wname2, tooltip='', cstep=-2, rstep=1, logvar=True,
                                      width=150, validator=None, enabled=True, maximum=100, minimum=0, decimals=0,
                                      wtype=None, value=1, stepsize=1):
        self.insert_checkbox(parent, wname, logvar=logvar, alignment=Qt.AlignLeft, cstep=1)
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_spinbox(parent, wname2, validator=validator, width=width, enabled=enabled,
                            maximum=maximum, minimum=minimum, cstep=cstep, rstep=rstep, value=value,
                            widgetType=wtype, stepsize=stepsize, logvar=logvar, decimals=decimals)

    def insert_label_combobox(self, parent, textlabel, wname, labels, rowspan=1, columnspan=1, rstep=1, cstep=-1,
                        width=150, tooltip='', logvar=False ):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_combobox(parent, wname, labels, rowspan=rowspan, columnspan=columnspan, rstep=rstep, cstep=cstep,
                             width=width, logvar=logvar)

    def insert_label_action_label(self, parent, actiontext, action='', cstep=0, rstep=1, sizepolicy='', params='', wname=''):
        self.insert_label(parent,alignment=Qt.AlignRight,sizepolicy=sizepolicy,rstep=1)
        self.insert_pushbutton(parent, text=actiontext, rstep=rstep, cstep=cstep, action=action, params=params, wname=wname)
        self.insert_label(parent,sizepolicy=sizepolicy,rstep=1)

    def insert_label_action(self, parent, actiontext, action='', cstep=0, rstep=1, sizepolicy='', params='', wname='',
                            width=150):
        self.insert_label(parent,alignment=Qt.AlignRight,sizepolicy=sizepolicy,rstep=1)
        self.insert_pushbutton(parent, text=actiontext, rstep=1, cstep=cstep, action=action, params=params, wname=wname,
                               width=width)

    def insert_gen_text_exe(self, parent, mode, action='', paramsAction=[], paramsXML=[], paramsCmd=[],
                            paramsSbatch={}, xmlfilename='', exefilename='exe.sh', jobfield=False, gpu=False,
                            queue=True, cs=3, mandatory_fill=[]):

        if not queue:
            self.insert_label_action(parent, 'Generate command', cstep=1, rstep=-1, sizepolicy=self.sizePolicyB,
                                           action=self.gen_action, wname=mode + 'GeneratePushButton',
                                           params=[[mode + 'XMLText'] + paramsXML,
                                                   [mode + 'CommandText'] + paramsCmd,
                                                   paramsSbatch])
            self.widgets[mode+'queue'] = QCheckBox()


        if queue:
            self.insert_label_action_label(parent, 'Generate command', cstep=1, rstep=-1, sizepolicy=self.sizePolicyB,
                                           action=self.gen_action, wname=mode + 'GeneratePushButton',
                                           params=[[mode + 'XMLText'] + paramsXML,
                                                   [mode + 'CommandText'] + paramsCmd,
                                                   paramsSbatch, mandatory_fill])
            self.insert_checkbox(parent,mode + 'queue',text='queue',cstep=-cs,rstep=1,logvar=True,alignment=Qt.AlignLeft)
        else:
            self.column -= cs
            self.row += 1

        if jobfield:
            self.insert_textfield(parent, mode + 'XMLText', columnspan=cs, rstep=0, cstep=3, width=600,logvar=False)
            if gpu:
                self.insert_checkbox(parent,mode + 'gpuRun',text='gpu',cstep=-1,rstep=1,logvar=True, alignment=Qt.AlignTop | Qt.AlignLeft)

            self.insert_label(parent, alignment=Qt.AlignRight, rstep=1, cstep=-3, sizepolicy=self.sizePolicyB)
        self.insert_textfield(parent, mode + 'CommandText', columnspan=cs, rstep=1, cstep=cs-1, width=600, logvar=False)

        if queue:
            self.insert_label_action_label(parent, 'Execute command', rstep=1, action=self.exe_action, wname=mode+'ExecutePushButton',
                                       params=[exefilename, mode+'CommandText', xmlfilename, mode+'XMLText', action,
                                               paramsAction, mode])
        else:
            self.insert_label_action(parent, 'Execute command', rstep=1, action=self.exe_action, wname=mode+'ExecutePushButton',
                                           params=[exefilename, mode + 'CommandText', xmlfilename, mode + 'XMLText',
                                                   action, paramsAction, mode])

    def exe_action(self, params):

        mode = params[6]
        if not self.widgets[mode + 'ExecutePushButton'].isEnabled():
            return

        self.widgets[mode + 'ExecutePushButton'].setEnabled(False)
        QApplication.processEvents()

        if params[4]: params[4](params[5])

        if 1:
            # Check if one needs to write pytom related XML file.
            if params[2]:
                if len(params[2][0]) > 1:
                    try:
                        tempfilename = os.path.join(self.widgets[params[2][0]].text(),self.widgets[params[2][1]].text())
                    except:
                        tempfilename = os.path.join(self.widgets[params[2][0]].text(), params[2][1])
                else:
                    tempfilename = params[2]

                jobfile = open(tempfilename,'w')
                jobfile.write(self.widgets[params[3]].toPlainText())
                jobfile.close()

            print(params[0], params[0][0])
            # Write executable file
            if len(params[0][0]) > 1:
                exefilename = os.path.join(self.widgets[params[0][0]].text(), params[0][1])
            else:
                exefilename = params[0]

            exefile = open(exefilename, 'w')

            exefile.write(self.widgets[params[1]].toPlainText())
            exefile.close()

            if self.widgets[mode+'queue'].isChecked():
                dd = os.popen('{} {}'.format(self.qcommand, exefilename))
                print('Submitted')
                text = dd.read()[:-1]
                id = text.split()[-1]
                self.popup_messagebox('Info','Submitted job to the queue', text)

                logcopy = os.path.join(self.projectname, f'LogFiles/{id}_{os.path.basename(exefilename)}')
                self.addProgressBarToStatusBar([id], key='QJobs',
                                               job_description=os.path.basename(exefilename)[:-3])
                try: os.system(f'cp {exefilename} {logcopy}')
                except Exception as e:
                    print(e)
                    
            else:
                import time
                ID = self.getLocalID()
                try:
                    self.localqID[ID] = 0
                except:
                    self.localqID = {}
                    self.localqID[ID] = 0

                self.activeProcesses[ID] = Worker(fn=self.submit_local_job, args=[exefilename, ID], sig=False)
                self.threadPool.start(self.activeProcesses[ID])

                # check = Worker(fn=self.checkLocalJob, args=[ID], sig=False)
                # check.start()
                self.addProgressBarToStatusBar(['Local_'+ID], key='QJobs',
                                               job_description=os.path.basename(exefilename)[:-3])
                # while self.localqID[ID] == 0:
                #     time.sleep(1)
                # self.popup_messagebox('Info', 'Local Job Finished', f'Finished Job {ID}')
                # os.system('sh {}'.format(params[0]))
        # except Exception as e:
        #     print(e)
        #     print ('Please check your input parameters. They might be incomplete.')



        self.widgets[mode + 'ExecutePushButton'].setEnabled(True)

    def countNumberOfParticles(self, wname):
        filename = self.widgets[wname].text()
        if filename and os.path.exists(filename):
            numberParticles = os.popen(f'grep Shift {filename} | wc -l').read()[:-1]
            self.popup_messagebox('Info', 'Number of Particles in ParticleList', f'Number of particles in {os.path.basename(filename)}: {numberParticles}')

    def checkLocalJob(self, ID):
        import time
        while self.localqID[ID] == 0:
            time.sleep(1)
        self.popup_messagebox('Info', 'Local Job Finished', f'Finished Job {ID}')


    def finishedJob(self, keys=None):
        self.popup_messagebox('Info', 'Job Finished', 'Job Finished')

    def executeJob(self, command, signals=None):
        os.system(command)

        print('Motioncor finished')
        signals.finished_mcor.emit()

    def gen_action(self, params):
        try:
            if params[3]:
                for field in params[3]:
                    if self.widgets[field].text() == '':
                        if not self.silent:  # TODO no attribute silent in commonfunctions??
                            self.popup_messagebox('Warning', 'Required field is empty',
                                                  'One of the required paramaters has not been supplied, '
                                                  'please supply it and regenerate script before pressing execute.')
                        return
        except Exception as e:
            print(e)

        mode = params[1][0][:-len('CommandText')]

        id = params[2]['id']
        if id or self.widgets[mode+'queue'].isChecked():  # TODO remove if statement here?
            partition, num_nodes, cores, time, modules, qcmd = self.qparams[id].values()

        for key in params[1][1:-1]:
            if 'numberMpiCores' in key and params[2]['id']:
                try:
                    ss = self.widgets[mode + 'gpuString'].text()
                    if ss == '':
                        update = True
                    else:
                        update = False
                except:
                    update = True
                if update:
                    self.widgets[key].setText(str(num_nodes*cores))

        for i in range(2):
            if params[i][0] in self.widgets.keys():
                text = params[i][-1]
                if params[i][1:-1]:
                    d = []
                    for a in params[i][1:-1]:
                        if type(a) == type([]):
                            d.append(os.path.join(self.widgets[a[0]].text(), a[1]))
                        elif a in self.widgets.keys():
                            datatype = type(self.widgets[a])
                            if datatype in (QSpinBox, QDoubleSpinBox, QLineEdit): d.append(self.widgets[a].text())
                            elif datatype == QComboBox: d.append(self.widgets[a].currentText())
                            elif datatype == QCheckBox: d.append(str(int(self.widgets[a].isChecked())))
                            else: pass
                        elif type(str(a)) == type(''):
                            d.append(a)
                        else: pass
                    self.widgets.keys()
                    text = text.format( d=d )
                if i==0: self.widgets[params[i][0]].setPlainText(text)

        # Check if user wants to submit to queue. If so, add queue header.
        if self.widgets[mode+'queue'].isChecked():
            d = params[2]
            folder = d['folder']
            suffix = d['suffix']

            # not needed it seems
            # num_jobs_per_node = d['num_jobs_per_node']
            # time = d['time']
            # partition = d['partition']
            # num_nodes = d['num_nodes']
            # id = d['id']
            # modules = d['modules']
            # if id:
            #     partition, num_nodes, num_jobs_per_node, time, modules, qcmd = self.qparams[id].values()

            try:
                gpus = self.widgets[mode+'gpuID'].text()
            except:
                gpus=''

            if type(folder) == type([]):

                outfolder = ''

                for extra in folder:
                    try:
                        extra = self.widgets[extra].text()
                    except:
                        pass
                    outfolder = os.path.join( outfolder, extra)
                folder = outfolder

            if self.custom:
                text = self.genSettingsWidgets['v00_QParams_CustomHeaderTextField'].toPlainText() + '\n\n' + text
            else:
                # this generates the queue header based on the parameters in self.qparams
                text = guiFunctions.gen_queue_header(name=d['fname'], folder=folder, cmd=qcmd, modules=modules,
                                                     qtype=self.qtype, partition=partition, time=time,suffix=suffix,
                                                     num_jobs_per_node=cores, num_nodes=num_nodes,
                                                     gpus=gpus) + text

        self.widgets[params[i][0]].setPlainText(text)

    def browse(self, params):
        mode = params[0]
        line = params[1]
        filetype = params[2]

        try:
            remote = params[3].isChecked()
        except:
            remote = False


        if line.text():
            if not remote:
                initdir = os.path.join('/', line.text())
            else:
                initdir = line.text()
        else:
            try:
                initdir = params[4]
            except:
                try:
                    initdir = self.projectname
                except:
                    initdir='./'


        try: initdir = initdir.text()
        except: pass


        if filetype:
            try:
                filters = "Image files (*.{})".format( filetype.currentText() )
            except:
                if type(filetype) == str:
                    filetype=[filetype]
                    filters = "Image files (*.{})".format( filetype )
                if type(filetype) == list:
                    filters = ''
                    for f in filetype:
                        filters += "{} files (*.{});; ".format(f, f)
                    filters = filters[:-3]
        else:
            filters = "Image files (*.*)"

        if not remote:
            if mode == 'folder':
                line.setText(QFileDialog.getExistingDirectory(self, 'Open file', initdir))
            if mode == 'file':

                filename, dummy = QFileDialog.getOpenFileName(self, 'Open file', initdir, filters)
                line.setText(str(filename))

        if remote:
            credentials = [self.servername, self.username, self.password]
            if not initdir: initdir = '/'
            filters = ['*.{}'.format(filetype.currentText())]

            b = BrowseWindowRemote(self, initdir=initdir, search=mode, filter=filters, outputline=line,
                                   credentials=credentials)

    def update_logbook(self,name='',text=''):

        if name:
            if text == 'checkbox':
                self.logbook[name] = self.widgets[name].isChecked()
            if text == 'lineedit':
                self.logbook[name] = self.widgets[name].text()
            if text == 'combobox':
                self.logbook[name] = self.widgets[name].currentText()
            elif text == 'spinbox':
                self.logbook[name] = self.widgets[name].value()
            elif text == 'textfield':
                self.logbook[name] = self.widgets[name].toPlainText()

    def create_expandable_group(self, action, sizeP, text, mode='', setVisible=False):
        a = QCheckBox(text=text)
        a.setSizePolicy(sizeP)

        try:
            b = action(mode=mode, title=text)
        except:
            b = action(mode=mode)

        a.clicked.connect(lambda ignore,  w=a, w2=b: self.toggle_groupbox_visible(w, w2))
        b.clicked.connect(lambda ignore, w2=a,  w=b: self.toggle_groupbox_visible(w, w2))
        b.setVisible(setVisible)
        return a, b

    def toggle_groupbox_visible(self,widget,widget2):
        widget.setVisible(False)
        widget2.setChecked(True==widget.isChecked())
        widget2.setVisible(True)

    def popup_messagebox(self, messagetype, title, message):
        import time
        if messagetype == 'Info':
            QMessageBox().information(self, title, message, QMessageBox.Ok)

        elif messagetype == 'Error':
            QMessageBox().critical(self, title, message, QMessageBox.Ok)

        elif messagetype == 'Warning':
            QMessageBox().warning(self, title, message, QMessageBox.Ok)

        elif messagetype == 'TimedInfo':
            TimedMessageBox(self, timeout=4, info=(self, title, message, QMessageBox.Ok), type=messagetype)

    def update_pb(self, pb, value):
        print( 'update: ', value.text() )
        pb.setValue( int(value.text()) )

    def fill_tab(self, id, headers, types, values, sizes, tooltip=[], wname='v02_batch_aligntable_', connect=0,
                 nn=False,
                 sorting=False, runbutton=True, runtitle='Run', addQCheckBox=True):
        try:
            self.tables[id].setParent(None)
            self.pbs[id].setParent(None)
            self.ends[id].setParent(None)
        except:

            pass

        self.tables[id] = SimpleTable(headers, types, values, sizes, tooltip=tooltip, connect=connect, id=id)
        self.widgets['{}{}'.format(wname, id)] = self.tables[id]

        if nn:
            num_nodes = QSpinBox()
            num_nodes.setValue(2)
            num_nodes.setRange(1, 99)
            num_nodes.setPrefix('Num Nodes: ')

        try:
            self.num_nodes[id] = num_nodes
        except:
            self.num_nodes = {}
            self.num_nodes[id] = 0
        if runbutton:
            print(f'Generating PB {wname}ExecutePushButton')
            self.pbs[id] = QPushButton(runtitle)
            self.pbs[id].setSizePolicy(self.sizePolicyC)
            self.widgets[f'{wname}ExecutePushButton'] = self.pbs[id]
        self.ends[id] = QWidget()
        self.ends[id].setSizePolicy(self.sizePolicyA)

        for n, a in enumerate((self.tables[id], self.num_nodes[id], self.checkbox[id], self.pbs[id], self.ends[id])):
            if n == 1 and nn == False:
                continue

            if n == 2:
                self.widgets[f'{id}_queue'] = a
                a.setEnabled(self.qtype != 'none')

                if addQCheckBox:
                    self.table_layouts[id].addWidget(a)
            else:
                self.table_layouts[id].addWidget(a)

    def create_groupbox(self,title,tooltip,sizepol):
        groupbox = QGroupBox(title)
        groupbox.setToolTip(tooltip)
        groupbox.setCheckable(True)
        groupbox.setChecked(False)
        groupbox.setSizePolicy(sizepol)
        parent = QGridLayout()
        groupbox.setLayout(parent)
        return groupbox, parent

    def submitBatchJob(self, execfilename, id, command, threaded=True):
        outjob = open(execfilename, 'w')
        outjob.write(command)
        outjob.close()

        if self.checkbox[id].isChecked():
            dd = os.popen('{} {}'.format(self.qcommand, execfilename))
            text = dd.read()[:-1]
            try:
                ID = text.split()[-1]
            except:
                return None, -1
            logcopy = os.path.join(self.projectname, f'LogFiles/{ID}_{os.path.basename(execfilename)}')
            os.system(f'cp {execfilename} {logcopy}')
            return ID, 1
        else:
            ID = self.getLocalID()
            self.localqID[ID] = 0

            if threaded:
                print('submitting threaded job')
                self.activeProcesses[ID] = Worker(fn=self.submit_local_job, args=[execfilename, ID], sig=False, results=True)
                self.threadPool.start(self.activeProcesses[ID])
            else:
                self.submit_local_job(execfilename, ID)

            return 'Local_'+ID, 0

    def multiSeq(self, func, params, wID=0, threaded=False):
        import time
        maxtime = 600
        for execfilename, pid, job in params:
            sleeptime = 0
            ID, num = func(execfilename, pid, job, threaded=threaded)
            if num == 0:
                self.localJobs[wID].append(ID)
                # while self.localqID[ID.split('_')[-1]] == 0 and sleeptime < maxtime:
                #     time.sleep(5)
                #     sleeptime += 5

    def getLocalID(self):
        import os
        try:
            idfile = os.path.join(self.logfolder, 'Local/.idcounter.txt')

            if os.path.exists(idfile):
                ID = open(idfile, 'r').readlines()[0].split()[0].zfill(8)
            else:
                ID = '0'.zfill(8)
        except:
            ID = '0'.zfill(8)

        out = open(idfile, 'w')
        out.write(str(int(ID) + 1) + '\n')
        out.close()

        return ID

    def submit_local_job(self, execfilename, ID):
        try:
            logcopy = os.path.join(self.logfolder, f'Local/{ID}-{os.path.basename(execfilename)}')
            os.system(f'cp {execfilename} {logcopy}')
            jobstring = f'sh {execfilename} >> {os.path.splitext(logcopy)[0]}.out'
            self.localJobStrings[ID] = jobstring
            os.system(jobstring)

            self.localqID[ID] = 1
        except:
            self.localqID[ID] = 2

    def submit_multi_job(self, exefilename, command, queue=True):

        try:
            exefile = open(exefilename, 'w')
            exefile.write(command)
            exefile.close()

            if queue:

                dd = os.popen('{} {}'.format(self.qcommand, exefilename))
                text = dd.read()[:-1]
                id = text.split()[-1]
                logcopy = os.path.join(self.projectname, f'LogFiles/{id}_{os.path.basename(exefilename)}')
                os.system(f'cp {exefilename} {logcopy}')
            else:
                os.system('sh {}'.format(exefilename))
        except:
            print ('Please check your input parameters. They might be incomplete.')

    def addProgressBarToStatusBar(self, qids=[], key='', job_description='Queue', num_submitted_jobs=0):
        try:
            if not key or not len(qids): return

            if num_submitted_jobs > 0:
                counters = [0,] * num_submitted_jobs
            else:
                counters = [0, ] * len(qids)

            manager = Manager()

            ID = qids[0]

            self.progressBarCounters[ID] = manager.list(counters)
            self.generateStatusBar(len(counters), ID, job_description)

            proc = Worker(fn=self.checkRun, args=(ID, qids, job_description, num_submitted_jobs))
            proc.signals.result1.connect(self.updateProgressBar)
            proc.signals.finished_queue.connect(self.deleteProgressBar)
            proc.start()
            event = Event()
            self.queueEvents[job_description] = event
        except Exception as e:
            print(e)
            print('''No statusbar created.''')

    def generateStatusBar(self, nrJobs, key, job_description):
        widget = QWidget(self)
        layout = QHBoxLayout()
        widget.setLayout(layout)

        self.progressBars[key] = QProgressBar()

        self.progressBars[key].setSizePolicy(self.sizePolicyA)
        self.progressBars[key].setStyleSheet(DEFAULT_STYLE_PROGRESSBAR)
        self.progressBars[key].setFormat(f'{job_description}: %v/%m')
        self.progressBars[key].setSizePolicy(self.sizePolicyB)
        self.progressBars[key].setMaximum(nrJobs)
        layout.addWidget(self.progressBars[key])
        self.statusBar.addPermanentWidget(self.progressBars[key], stretch=1)

    def whichRunning(self, qids, num_submitted_jobs):
        qids_tot = [self.localJobs[wid] for wid in qids if type(wid) == int]
        qt = []
        for q in qids_tot:
            qt += q

        [qt.append(q) for q in qids if type(q) == str]

        qids = qt

        runningJobs = [rid[:-1] for rid in os.popen(" squeue | awk 'NR > 1 {print $1}' ").readlines()]
        inQueue = [str(int(qid)) for qid in qids if qid in runningJobs]
        inQueue += [qid[6:] for qid in qids if str(qid).startswith('Local_') and self.localqID[qid[6:]] == 0]
        if len(qt) < num_submitted_jobs:
            inQueue += [0,]*(num_submitted_jobs-len(qt))
        return inQueue

    def whichLocalJobsRunning(self, qids):

        completedJobs = [[id for id in self.localJobs[qid] if self.localqID[id] > 0] for qid in qids]

    def checkRun(self, id, qids, job_description, num_submitted_jobs, signals):
        import time
        total_in = num_submitted_jobs if num_submitted_jobs > 0 else len(qids)
        inQueue = self.whichRunning(qids, num_submitted_jobs)
        while len(inQueue):
            total = total_in - len(inQueue)
            signals.result1.emit([id, total])
            time.sleep(1)
            inQueue = self.whichRunning(qids, num_submitted_jobs)
            # print(f'Number of active threads: {self.threadPool.activeThreadCount()}')
            if self.queueEvents[job_description].is_set():
                print(f'exit {job_description} loop')
                return

        qids_tot = [self.localJobs[wid] for wid in qids if type(wid) == int]
        qt = []
        for q in qids_tot: qt += q

        [qt.append(q) for q in qids if type(q) == str]

        qids = qt

        #self.popup_messagebox("Info", "Completion", f'Finished Queue Jobs {job_description}')
        signals.finished_queue.emit([id, job_description, qids])
        #self.popup_messagebox("Info", "Completion", f'Finished Queue Jobs {job_description}')

    def updateProgressBar(self, keys):
        key, total = keys
        self.progressBars[key].setValue(total)

    def deleteProgressBar(self, keys):
        key, job_description, qids = keys
        self.statusBar.removeWidget(self.progressBars[key])
        if len(qids)>1:
            self.popup_messagebox("Info", "Completion", f'Finished {job_description} Jobs {qids[0]}-{qids[-1]}')
        elif len(qids) > 0:
            self.popup_messagebox("Info", "Completion", f'Finished {job_description} Job {qids[0]}')

    def retrieveJobID(self,results):
        ID, num = results
        print(f'finished job {ID}')
        self.popup_messagebox('TimedInfo', 'Finished Job', f'Finished job {ID}')

    def activate_stage(self, stage_id):
        if not self.parent().parent().parent().stage_buttons[stage_id].isEnabled():
            self.parent().parent().parent().stage_buttons[stage_id].setEnabled(True)
            self.parent().parent().parent().logbook['00_framebutton_{}'.format(self.parent().parent().parent().targets[stage_id][0])] = True
            self.parent().parent().parent().save_logfile()

    def update_cores_flag(self, widget_cores, widget_cores_flag, suppress_message=False):
        """
        Set the cores command line flag for an algorithm.
        Pass widget that holds user cores input, and one that holds the flag.
        @param widget_cores: text widget for user gpu ids
        @param widget_cores_flag: text widget with gpu flag
        @param suppress_message: dont show qmessage box if input is wrong
        """
        n_cores = widget_cores.text()
        if n_cores == '':
            widget_cores_flag.setText('')
            return
        try:
            n_cores = int(n_cores)
        except ValueError:
            widget_cores.setText('')
            widget_cores_flag.setText('')

            if suppress_message:
                return
            else:
                self.popup_messagebox('Warning', 'Invalid value in field',
                                      'Can only fill in integer number of cpu cores')
                return

        widget_cores_flag.setText(f'--numProcesses {n_cores} ')

    def update_gpu_flag(self, widget_gpu_ids, widget_gpu_flag, single_gpu=False, widget_mpi_cores=None,
                        suppress_message=False):
        """
        Set the gpu command line flag for an algorithm.
        Pass widget that holds user gpu ids input, and one that holds the flag.

        @param widget_gpu_ids: text widget for user gpu ids
        @param widget_gpu_flag: text widget with gpu flag
        @param single_gpu: restrict to parse only one gpu
        @param widget_mpi_cores: set mpi number of cores widget for algos that use openmpi for gpu parallel
        @param suppress_message: dont show qmessage box if input is wrong
        """
        gpu_ids = widget_gpu_ids.text()

        if gpu_ids == '':
            widget_gpu_flag.setText('')
            return

        if single_gpu:
            try:
                gpu = int(gpu_ids)
            except ValueError:
                widget_gpu_ids.setText('')
                widget_gpu_flag.setText('')

                if suppress_message:
                    return
                else:
                    self.popup_messagebox('Warning', 'Invalid value in field',
                                          'Only single gpu allowed for this algorithm.')
                    return

            widget_gpu_flag.setText(f'--gpuID {gpu} ')

        else:
            try:
                gpu = list(map(int, [g for g in gpu_ids.split(',') if g != '']))
            except ValueError:
                widget_gpu_ids.setText('')
                widget_gpu_flag.setText('')

                if suppress_message:
                    return
                else:
                    self.popup_messagebox('Warning', 'Invalid value in field',
                                          'Could not parse gpu field, should be of form: 0,3,4 ')
                    return

            widget_gpu_flag.setText(f'--gpuID {gpu_ids} ')

            if widget_mpi_cores is not None:
                widget_mpi_cores.setText(f'{len(gpu)+1}')


class SelectFiles(BrowseWindowRemote, CommonFunctions):
    def __init__(self, parent=None, initdir='/', filter=[''], search='file', credentials=['', '', ''], outputline='',
                 validate=False, run_upon_complete=print, title='', id=''):
        self.finished = False
        super(SelectFiles, self).__init__(parent, initdir=initdir, filter=[''], search='file', credentials=['', '', ''],
                                          outputline=outputline, validate=validate)
        self.setWindowTitle(title)
        self.run_upon_complete = run_upon_complete
        self.filters = filter
        self.id = id
        self.matchingfiles = []
        self.selectedfiles = []
        if self.outputline.text(): self.selectedfiles += [ll for ll in self.outputline.text().split('\n')]
        self.setCentralWidget(self.splitter0)
        self.add_folders()
        self.show()

    def select_item(self):
        if not self.selectedfiles:
            self.popup_messagebox('Info', 'No Files Selected', 'You did not select any files. Please double click the file(s) you want to select.\nThe file will move to the right when selected.')
            return
        self.outputline.setText('\n'.join(self.selectedfiles))
        self.finished = True
        self.run_upon_complete(self.id)

    def select_file(self):
        item = self.topright.currentItem().text()

        for n, file in enumerate(self.selectedfiles):
            if file == item:
                self.selectedfiles = self.selectedfiles[:n] + self.selectedfiles[n + 1:]
                break
        self.add_folders()

    def data(self):
        self.subdirs = [os.pardir]
        self.matchingfiles = []
        self.append(os.listdir(self.folderpath))

        matchingfiles = []
        for pat in self.filters:
            matchingfiles += [line for line in self.matchingfiles if line.endswith(pat.split('.')[-1])]
        self.matchingfiles = matchingfiles

    def append(self, line):
        for filename in line:
            folder = os.path.isdir(os.path.join(self.folderpath, filename))

            if folder:
                self.subdirs.append(filename)
            else:
                selected = False
                for file in self.selectedfiles:
                    if os.path.join(self.folderpath, filename) == file:
                        selected = True
                        break
                if not selected: self.matchingfiles.append(filename)

        self.matchingfiles = sorted(self.matchingfiles)
        self.subdirs = sorted(self.subdirs)

    def add_folders(self, folders=[]):
        self.data()

        self.topleft.clear()
        self.topright.clear()

        self.topleft.insertItems(0, self.subdirs + self.matchingfiles)
        self.topright.insertItems(0, self.selectedfiles)

    def repopulate_folder_list(self):

        extra = self.topleft.currentItem().text()
        potential = os.path.join(self.folderpath, extra)
        if extra != '..' and os.path.isdir(potential):
            self.folderpath = potential

        elif len(self.folderpath) > 1 and extra == '..':
            self.folderpath = os.path.dirname(self.folderpath)
        else:
            self.selectedfiles.append(potential)

        self.pathdisplay.setText(self.folderpath)
        self.add_folders()


class CreateMaskFile(QMainWindow, CommonFunctions):
    def __init__(self,parent=None,maskfname='',folder=''):
        super(CreateMaskFile,self).__init__(parent)
        w = QWidget(self)
        l = QGridLayout()

        self.logbook = self.parent().logbook
        w.setLayout(l)
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = l
        self.insert_label_line(parent,'Sixe X (px)',wname='size_template_x',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent,'Size Y (px)',wname='size_template_y',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent,'Size Z (px)',wname='size_template_z',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent, 'Radius Spherical Mask', wname='radius',rstep=1, cstep=-1,
                               validator=QIntValidator())
        self.insert_label_spinbox(parent, 'smooth_factor', text='Smoothing Factor', rstep=1, cstep=0, value=0,
                                  wtype=QDoubleSpinBox,stepsize=0.1)
        self.insert_pushbutton(parent, 'Create', action=self.generate_mask,
                               params=['size_template', 'radius','smooth_factor',maskfname, folder])
        self.setCentralWidget(w)


        self.widgets['size_template_x'].textChanged.connect(self.set_sizes)

        self.show()

    def set_sizes(self):
        size = self.widgets['size_template_x'].text()
        for name in ('size_template_y','size_template_z'):
            self.widgets[name].setText(size)

    def generate_mask(self,params):
        try:
            radius = int(self.widgets['radius'].text())
            smooth = float(self.widgets['smooth_factor'].value())
            sizeX = int(self.widgets['size_template_x'].text())
            sizeY = int(self.widgets['size_template_y'].text())
            sizeZ = int(self.widgets['size_template_z'].text())
            try:
                folder = params[4] if params[4] else self.parent().projectname
                print(folder)
                fname = os.path.join(folder, f'Mask_{sizeX}_{radius}_{smooth:.2f}.mrc')
            except Exception as e:
                print(e)
                tomoname = os.path.basename(self.parent().widgets['v03_TemplateMatch_tomoFname'].text())
                filename, file_extension = os.path.splitext(tomoname)
                if not os.path.exists(os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename)):
                    os.mkdir(os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename))
                fname = os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename, 'TM_mask.mrc')

            maskfilename = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.mrc')[0])
            if not maskfilename:
                return

            if maskfilename and not maskfilename.endswith('.mrc'): maskfilename += '.mrc'
            try:
                success = initSphere(sizeX, sizeY, sizeZ, radius=radius, smooth=smooth, filename=maskfilename)
                if success:
                    self.parent().widgets[params[3]].setText(maskfilename)
                else:
                    self.popup_messagebox('Error', 'Mask Generation Failed',
                                          'Generation of the mask failed. Please select an existing mask, or generate a mask yourself.')

            except:
                self.popup_messagebox('Error','Mask Generation Failed', 'Generation of the mask failed. Please select an existing mask, or generate a mask yourself.')

            self.close()
        except Exception as e:
            print(e)
            print('Mask generation failed.')


class ConvertEM2PDB(QMainWindow, CommonFunctions):
    def __init__(self,parent,emfname='',folder='./'):
        super(ConvertEM2PDB, self).__init__(parent)
        self.folder = folder
        w = QWidget(self)
        l = QGridLayout()
        self.logbook = self.parent().logbook
        w.setLayout(l)
        self.widgets = {}
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = l

        self.insert_label_line(parent, 'PDB ID', wname='pdb_id', cstep=-1, rstep=1,width=150)

        self.insert_label_line(parent, 'Chain', wname='chain', rstep=1, cstep=-1, value='A')

        self.insert_label_spinbox(parent, 'volume_size', text='Volume size', rstep=1, cstep=-1,
                                  stepsize=4,minimum=8,maximum=4096,value=128,
                                  tooltip='Volume length(size) in all dimensions')

        self.insert_label_spinbox(parent, 'pixel_size', text='Pixel Size', rstep=1, cstep=0, value=3.5,
                                  wtype=QDoubleSpinBox, stepsize=0.1, minimum=1.0, maximum=100,
                                  tooltip='Pixel size of output volume (in Angstrom)')

        self.insert_pushbutton(parent, 'Create', action=self.pdb2emr,
                               params=['size_template', 'radius', 'smooth_factor', emfname])

        self.setCentralWidget(w)
        self.show()

    def pdb2emr(self,params):
        pdbid = self.widgets['pdb_id'].text()
        if not pdbid: return
        chain = self.widgets['chain'].text()
        cube_size = int(self.widgets['volume_size'].text())
        pixel_size = float(self.widgets['pixel_size'].value())
        os.system('getpdb {} {}'.format(pdbid, self.folder))
        out_fname = str(QFileDialog.getSaveFileName(self, 'Save EM model.', self.folder, filter='*.mrc')[0])
        if not out_fname: return
        if not out_fname.endswith('.mrc'): out_fname += '.mrc'
        fname_pdb = f'{self.folder}/{pdbid}.pdb'

        try:
            pdb2mrc(fname_pdb, pixel_size, cube_size, chain=chain, fname=out_fname, invertDensity=True)
        except Exception as e:
            print(e)
            return

        self.parent().widgets[params[-1]].setText(out_fname)
        self.close()


class CreateFSCMaskFile(QMainWindow, CommonFunctions):
    def __init__(self, parent, emfname='',folder='./'):
        super(CreateFSCMaskFile, self).__init__(parent)
        self.folder = folder
        w = QWidget(self)
        l = QGridLayout()
        self.logbook = self.parent().logbook
        w.setLayout(l)
        self.widgets = {}
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = l

        self.insert_label_line_push(parent, 'Volume (Filtered)', 'volume', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='Volume path.', initdir=self.folder)
        self.insert_label_line_push(parent, 'Mask (Optional)', 'mask', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='The mask is used to only select a part of the model for resolution '
                                            'determination.', initdir=self.folder)
        self.insert_label_spinbox(parent, 'numstd', text='Threshold: #std below mean', rstep=1, cstep=-1,
                                  minimum=0, maximum=100, value=1, wtype=QDoubleSpinBox,
                                  tooltip='This parameter sets the threshold value for what is a particle.\n '
                                          'Threshold = mean signal - num_stds * std signal. ')
        self.insert_label_spinbox(parent, 'smooth', rstep=1, cstep=-1, wtype=QDoubleSpinBox,
                                  tooltip='std for the gaussian kernel used for smoothing of the edge of the mask.',
                                  minimum=0, maximum=100, value=2,
                                  text='Smoothing Edges')
        self.insert_label_spinbox(parent, 'cycles', text='Number of Dilation Cycles', rstep=1, cstep=0,
                                  stepsize=1,minimum=0,maximum=100,value=3,
                                  tooltip='Number of dilation cycles. Creates a less structured mask')

        self.insert_pushbutton(parent, 'Create', action=self.generate,
                               params=['volume', 'mask', 'numstd', 'smooth', 'cycles', emfname])

        self.setCentralWidget(w)
        self.show()

    def generate(self,params):
        from pytom.bin.gen_mask import gen_mask_fsc
        from pytom.agnostic.io import read



        data = read(self.widgets[params[0]].text())
        if self.widgets[params[1]].text():
            mask = read(self.widgets[params[1]].text())
        else:
            mask= None
        numstd = float(self.widgets[params[2]].value())
        smooth = float(self.widgets[params[3]].value())
        cycles = int(self.widgets[params[4]].value())

        fname = os.path.join(self.folder, f'mask_{numstd:.2f}_{smooth:.2f}_{cycles}.mrc')
        out_fname = str(QFileDialog.getSaveFileName(self, 'Save model as.', fname, filter='*.mrc')[0])
        if not out_fname: return
        if not out_fname.endswith('.mrc'): out_fname += '.mrc'

        gen_mask_fsc(data, cycles, out_fname, numstd, smooth, maskD=mask)

        self.parent().widgets[params[5]].setText(out_fname)
        self.close()


class MyCircleOverlay(pg.EllipseROI):
    def __init__(self, pos, size, label='', **args):
        self.path = None
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = True

        if label:
            self.label = QGraphicsTextItem(label, self)
            # for d in dir(self.label): print (d)
            self.label.setHtml(html.format(colorsHEX[int(label)], int(label)))
            self.label.setPos(
                QPoint(self.boundingRect().center().x() - (self.label.boundingRect().width() / 2) * 0,
                              self.state['size'][1] / 1.5))


def circle(pos, size=2, label='', color=Qt.blue, pensize=16):
    pensize = 0.02*40./pensize
    return MyCircleOverlay(pos=pos, size=size, label=label, pen=QPen(color, pensize), movable=False)


class SimpleTable(QMainWindow, CommonFunctions):

    def __init__(self, headers, types, values, sizes=[], tooltip=[] , connect=0, sorting=False, id=''):
        super(SimpleTable, self).__init__()
        self.size_policies()
        self.setWindowFlags(Qt.FramelessWindowHint)

        #self.scrollarea.setGeometry(0, 0, 700, 1000)
        central_widget = QWidget(self)

        #self.scrollarea.setWidget(central_widget)

        #self.setGeometry(10, 10, 700, 300)  # Create a central widget
        grid_layout = QGridLayout(self)
        grid_layout.setContentsMargins(0,0,0,0) # Create QGridLayout
        central_widget.setLayout(grid_layout)  # Set this layout in central widget
        table = QTableWidget(self)
        table.setColumnCount(len(headers))  # Set three columns
        table.setRowCount(len(values))  # and one row
        header = table.horizontalHeader()
        self.widgets = {}
        # Set the table headers
        table.setHorizontalHeaderLabels(headers)
        table.setVerticalHeaderLabels(['', ] * len(values))
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.headers = headers
        hh = table.horizontalHeader()
        table.verticalHeader().hide()


        options= []
        for i in range(len(headers)):
            options.append([])
        # Set the alignment to the headers
        for i in range(len(headers)):

            if i+1 == len(headers):
                hh.setResizeMode(i, QHeaderView.Stretch)
            elif sizes[i] == 0:
                hh.setSectionResizeMode(i,QHeaderView.ResizeToContents)
            else:
                hh.setSectionResizeMode(i,sizes[i])

            if i+1 != len(headers):
                table.horizontalHeaderItem(i).setTextAlignment(Qt.AlignHCenter)
            else:
                table.horizontalHeaderItem(i).setTextAlignment(Qt.AlignLeft)

            if len(tooltip) > i: table.horizontalHeaderItem(i).setToolTip(tooltip[i])

            for v in range(len(values)):

                # Fill the first line
                if types[i] == 'txt':
                    widget = QWidget()
                    data= "{}".format(values[v][i].split('/')[-1]) if '/' in values[v][i] else values[v][i]
                    cb = QLabel( data )
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(10, 0, 10, 0)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v, i)] = cb
                    widget.setStyleSheet('background: white;')

                elif types[i] == 'txt2':
                    widget = QWidget()

                    try:
                        a = values[v][i].split('/')[-2]
                        if 'import_' in a:
                            data = f"{a}/{values[v][i].split('/')[-1]}"
                        else:
                            data= "{}".format(values[v][i].split('/')[-1]) if '/' in values[v][i] else values[v][i]
                    except:
                        data= "{}".format(values[v][i].split('/')[-1]) if '/' in values[v][i] else values[v][i]

                    cb = QLabel( data )
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(10, 0, 10, 0)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v, i)] = cb
                    widget.setStyleSheet('background: white;')

                elif types[i] == 'sort_txt':
                    table.setSortingEnabled(True)
                    item = QTableWidgetItem()
                    item.setTextAlignment(Qt.AlignCenter)
                    item.setData(Qt.EditRole, values[v][i])
                    table.setItem(v, i, item)
                    self.widgets['widget_{}_{}'.format(v, i)] = item


                elif types[i] == 'checkbox' and values[v][i]:
                    widget = QWidget()
                    cb = QCheckBox()
                    if connect:
                        cb.stateChanged.connect(lambda d, ID=id, rowID=v, columnID=i: connect(ID, rowID, columnID))
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    if values[v][i] ==16:
                        cb.setChecked(True)
                        cb.setEnabled(False)
                    if i+1==len(headers): layoutCheckBox.setAlignment(Qt.AlignLeft)
                    layoutCheckBox.setContentsMargins(10, 0, 10, 0)
                    table.setCellWidget(v,i,widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = cb
                    widget.setStyleSheet('background: white;')

                elif types[i] == 'lineedit':
                    widget = QWidget()
                    layoutCheckBox = QHBoxLayout(widget)
                    le = QLineEdit()
                    le.setFixedWidth(table.columnWidth(i))
                    le.setText(str(values[v][i]))
                    layoutCheckBox.addWidget(le)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = le
                    widget.setStyleSheet('background: white;')

                    if values[v][i] == 'UNCHANGED':
                        widget.setEnabled(False)

                elif types[i] == 'spinbox':
                    widget = QSpinBox()
                    widget.setValue(values[v][i])

                elif types[i] in ('combobox', 'comboboxF'):

                    widget = QWidget()
                    cb = QComboBox()
                    l = QVBoxLayout(widget)
                    l.addWidget(cb)
                    cb.setContentsMargins(0, 0, 0, 0)
                    l.setContentsMargins(0, 0, 0, 0)
                    for value in values[v][i]:
                        if types[i].endswith('F'):
                            try:
                                val = "/".join(value.split('/')[-2:])
                            except:
                                val = value
                        else:
                            val = value.split('/')[-1]
                        if not val in options[i]:
                            options[i] +=[val]

                        cb.addItem(val)

                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = cb
                    widget.setStyleSheet('background: white; selection-background-color: #1989ac;')
                else:
                    widget = QWidget()
                    cb = QLabel('')
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                    table.setCellWidget(v, i, widget)
                    widget.setStyleSheet('background: white;')

        #if not sizes: table.resizeColumnsToContents()
        self.sizePolicyC = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sizePolicyC.setHorizontalStretch(0)
        self.sizePolicyC.setVerticalStretch(0)
        self.sizePolicyB = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.sizePolicyB.setVerticalStretch(0)



        self.types = types
        self.general_widgets = []
        self.table2 = QTableWidget()
        self.table2.setColumnCount(len(headers))
        self.table2.setRowCount(1)
        self.table2.verticalHeader().hide()
        self.table2.horizontalHeader().hide()
        self.table2.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table2.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        #table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        #table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        applyOptionSum = 0

        for n, t in enumerate(types):
            if t == 'checkbox':
                widget = QWidget()
                cb = QCheckBox()
                l = QVBoxLayout(widget)
                l.addWidget(cb)
                #cb.setStyleSheet('QCheckBox::indicator{ width: 10px; height: 10px; }')
                cb.setFixedWidth(table.columnWidth(n))
                widget.setSizePolicy(self.sizePolicyC)
                cb.setToolTip('{} All'.format( headers[n].capitalize() ))

                cb.setContentsMargins(10, 0, 10, 0)
                l.setContentsMargins(10, 0, 10, 0)

                self.general_widgets.append(cb)
                self.table2.setCellWidget(0, n, widget)
                cb.stateChanged.connect(lambda dummy, rowIndex=n, c=t: self.on_changeItem(rowIndex,c))
                applyOptionSum +=1
            elif t in ( 'combobox', 'comboboxF'):
                widget = QWidget()
                cb = QComboBox()
                l = QVBoxLayout(widget)
                l.addWidget(cb)

                cb.setContentsMargins(0, 0, 0, 0)
                l.setContentsMargins(0, 0, 0, 0)
                cb.setFixedWidth(table.columnWidth(n))
                # for value in values[0][n]:
                #     if t.endswith('F'):
                #         try : v = "/".join(value.split('/')[-2:])
                #         except: v =value
                #     else:
                #         v = value.split('/')[-1]
                #
                #     cb.addItem(v)

                for v in options[n]:
                    cb.addItem(v)

                #glayout.addWidget(cb)
                widget.setStyleSheet('selection-background-color: #1989ac;')
                self.general_widgets.append(cb)
                self.table2.setCellWidget(0, n, widget)
                cb.currentTextChanged.connect(lambda dummy, rowIndex=n, c=t.replace('F', ''): self.on_changeItem(rowIndex, c))
                applyOptionSum +=1

            elif t == 'lineedit':
                widget = QWidget()
                layoutCheckBox = QHBoxLayout(widget)
                le = QLineEdit()
                le.setText(str(values[0][n]))
                layoutCheckBox.addWidget(le)
                layoutCheckBox.setAlignment(Qt.AlignCenter)
                layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                self.table2.setCellWidget(0, n, widget)
                self.general_widgets.append(le)
                applyOptionSum += 1
                le.textChanged.connect(lambda dummy, rowIndex=n, c=t: self.on_changeItem(rowIndex, c))
            else:
                widget = QWidget()
                a = QLabel('')
                l = QVBoxLayout(widget)
                l.setAlignment(Qt.AlignCenter)
                if n == 0:
                    a.setText('Apply to all')
                    l.setAlignment(Qt.AlignCenter)
                l.addWidget(a)
                a.setFixedWidth(table.columnWidth(n))
                widget.setSizePolicy(self.sizePolicyC)
                l.setContentsMargins(0, 0, 0, 0)
                self.table2.setCellWidget(0, n, widget)
                self.general_widgets.append(a)

        #glayout.addWidget(self.table2)
        self.table = table
        self.table2.setMaximumHeight(self.table2.rowHeight(0))
        self.table.setMinimumHeight(min(400, self.table2.rowHeight(0)*(len(values)+.8)))

        if applyOptionSum: grid_layout.addWidget(self.table2, 0, 0)
        grid_layout.addWidget(table, 1, 0)
        #self.table2.horizontalHeaderItem(len(headers) - 1).setTextAlignment(Qt.AlignLeft)

        for i in range(self.table.columnCount()):
            self.table2.setColumnWidth(i, self.table.columnWidth(i))

        self.setCentralWidget(central_widget)

        self.sliderBar1 = self.table2.horizontalScrollBar()
        self.sliderBar2 = self.table.horizontalScrollBar()

        self.sliderBar2.valueChanged.connect(lambda d, s1=self.sliderBar2, s2=self.sliderBar1: self.SyncScroll(s1,s2))

    def SyncScroll(self,slider1, slider2):
        sliderValue = slider1.value()
        slider2.setValue(sliderValue)

    def on_changeItem(self, rowIndex, widgetType):
        for i in range(self.table.rowCount()):
            if 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='combobox':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setCurrentText(self.general_widgets[rowIndex].currentText())
            elif 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='checkbox':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setChecked(self.general_widgets[rowIndex].isChecked())
            elif 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='lineedit':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setText(self.general_widgets[rowIndex].text())

        self.table.horizontalHeader().setResizeMode(len(self.headers)-1, QHeaderView.Stretch)

        for i in range(self.table.columnCount()):
            self.table2.setColumnWidth(i, self.table.columnWidth(i))


class GuiTabWidget(QWidget, CommonFunctions):
    def __init__(self, parent=None, headers=[],offx=0,offy=0,dimx=900,dimy=721,logbook=[]):

        super(GuiTabWidget, self).__init__(parent)

        self.addTabs(headers=headers, offx=offx, offy=offy, dimx=dimx, dimy=dimy,soff=50)

    def addTabs(self, headers, widget=QWidget, subheaders=[], offx=0,offy=0,dimx=900,dimy=721,soff=0,
                sizeX=900,sizeY=700, tabUIs=None, tabs=None, tab_actions=None, static_tabs=None):

        self.size_policies()
        self.scrollarea = QScrollArea(self)
        if soff: self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        #self.scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scrollarea.setSizePolicy(self.sizePolicyA)
        self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setContentsMargins(0,0,0,0)
        #self.scrollarea.setGeometry(offx, offy, dimx, dimy)

        self.scrollarea.setGeometry(0,0,sizeX,sizeY)
        self.scrollarea.setFrameShape(QFrame.NoFrame)

        self.tabWidget = QTabWidget(self.scrollarea)
        self.tabWidget.setContentsMargins(0,0,0,0)
        self.tabWidget.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)

        self.scrollarea.setWidget(self.tabWidget)
        self.tabs = []
        self.tabs2 = {}
        self.tab_actions2 = {}

        try:
            self.scrollareas.append(self.scrollarea)
            self.scrolloffset.append(soff)
        except:
            self.scrollareas = [self.scrollarea]
            self.scrolloffset = [soff]

        try:
            self.logbook = self.parent().logbook
            self.widgets = {}
            self.statusBar = self.parent().sbar
            self.scrollarea.setStyleSheet('background: #{};'.format(self.parent().middlec))
        except:
            pass

        for n, header in enumerate(headers):
            if len(subheaders) >= len(headers) and len(subheaders[n]) > 1:
                subtab = widget(headers=subheaders[n], dimx=dimx - 50)
                self.scrollareas += subtab.scrollareas
                self.scrolloffset += subtab.scrolloffset
                tab = QWidget(self.tabWidget)
                tab.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
                layout = QVBoxLayout()
                layout.addWidget(subtab)
                tab.setLayout(layout)
                for m, subheader in enumerate(subheaders[n]):
                    setattr(self, 'tab{}{}'.format(n + 1, m + 1), subtab.tabs[m])

                    self.tabs2[f'tab{n+1}{m+1}'] = getattr(self, f'tab{n+1}{m+1}')
                    tabs[f'tab{n+1}{m+1}'] = getattr(self, f'tab{n+1}{m+1}')
                    if tabUIs[n][m]:
                        tab_actions[f'tab{n+1}{m+1}'] = tabUIs[n][m]
                        self.tab_actions2['tab{}{}'.format(n + 1, m + 1)] = tabUIs[n][m]

            else:
                tab = QWidget()
            self.tabs.append(tab)
            tab.setObjectName(header)
            setattr(self, 'tab{}'.format(n + 1), tab)

            try:
                tabs[f'tab{n+1}'] = getattr(self, f'tab{n+1}')
                self.tabs2[f'tab{n+1}'] = getattr(self, f'tab{n+1}')
                if tabUIs[n]:
                    tab_actions[f'tab{n+1}'] = tabUIs[n][0]
                    self.tab_actions2[f'tab{n+1}'] = tabUIs[n][0]

            except Exception as e:
                pass
            self.tabWidget.addTab(tab, header)

        # self.tabs2 = {'tab1': self.tab1,
        #              'tab2':  self.tab2,
        #              'tab31': self.tab31, 'tab32': self.tab32,
        #              'tab41': self.tab41, 'tab42': self.tab42, 'tab43': self.tab43,
        #              'tab51': self.tab51, 'tab52': self.tab52, 'tab53': self.tab53}
        #

        # self.tab_actions2 = {'tab1':  self.tab1UI,
        #                     'tab2':  self.tab2UI,
        #                     'tab31': self.tab31UI, 'tab32': self.tab32UI,
        #                     'tab41': self.tab41UI, 'tab42': self.tab42UI, 'tab43': self.tab43UI,
        #                     'tab51': self.tab51UI, 'tab52': self.tab52UI, 'tab53': self.tab53UI}

        if not static_tabs is None:

            for i in range(len(headers)):
                t = 'tab{}'.format(i + 1)
                empty = 1 * (len(subheaders[i]) == 0)

                for j in range(len(subheaders[i]) + empty):
                    tt = t + (str(j + 1) * (1 - empty))

                    if static_tabs[i][j]:  # tt in ('tab2', 'tab31', 'tab41', 'tab42', 'tab51', 'tab52'):
                        self.table_layouts[tt] = QGridLayout()
                    else:
                        self.table_layouts[tt] = QVBoxLayout()

                    self.tables[tt] = QWidget()
                    self.pbs[tt] = QWidget()
                    self.ends[tt] = QWidget()
                    self.ends[tt].setSizePolicy(self.sizePolicyA)
                    self.checkbox[tt] = QCheckBox('queue')

                    if not static_tabs[i][j]:  # tt in ('tab1','tab32', 'tab43', 'tab53'):
                        button = QPushButton('Refresh Tab')
                        button.setSizePolicy(self.sizePolicyC)
                        button.clicked.connect(lambda d, k=tt, a=self.tab_actions[tt]: a(k))
                        self.table_layouts[tt].addWidget(button)
                        self.table_layouts[tt].addWidget(self.ends[tt])

                    else:  # if tt in ('tab2','tab31','tab41', 'tab42', 'tab51', 'tab52'):

                        self.tab_actions[tt](tt)

                    tab = self.tabs_dict[tt]
                    tab.setLayout(self.table_layouts[tt])

    def addGeneralVariables(self):

        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.tomogram_folder = self.parent().tomogram_folder
        self.rawnanographs_folder = self.parent().rawnanographs_folder
        self.motioncor_folder = self.parent().motioncor_folder
        self.particlepick_folder = self.parent().particlepick_folder
        self.tomoanalysis = self.parent().tomogram_folder
        self.subtomodir = self.parent().subtomo_folder
        self.subtomofolder = self.parent().subtomo_folder

        self.fscdir = os.path.join(self.subtomodir, 'Validation')
        self.polishfolder = os.path.join(self.subtomodir, 'ParticlePolishing')
        self.frmdir = os.path.join(self.subtomodir, 'Alignment/FRM')
        self.glocaldir = os.path.join(self.subtomodir, 'Alignment/GLocal')
        self.cpcadir = os.path.join(self.subtomodir, 'Classification/CPCA')
        self.acdir = os.path.join(self.subtomodir, 'Classification/AutoFocus')
        self.acpath = os.path.join(self.subtomodir, 'Classification/AutoFocus')

        self.pickpartdir = os.path.join(self.particlepick_folder, 'Picked_Particles')
        self.pickpartfolder = os.path.join(self.particlepick_folder, 'Picked_Particles')
        self.tomogramfolder = os.path.join(self.particlepick_folder, 'Tomograms')
        self.templatematchfolder = os.path.join(self.particlepick_folder, 'Template_Matching')
        self.ccfolder = os.path.join(self.templatematchfolder, 'cross_correlation')

        self.silent = self.parent().silent
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand
        self.widgets = {}
        self.progressBarCounters = {}
        self.progressBars = {}
        self.queueEvents = self.parent().qEvents
        self.localqID = {}
        self.activeProcesses = {}
        self.threadPool = self.parent().threadPool
        self.localJobs = {}
        self.tabs_dict, self.tab_actions = {}, {}
        self.qparams = self.parent().qparams
        self.localJobStrings = {}
        self.modes = {}

        self.widgets['pytomPath'] = QLineEdit()
        self.widgets['pytomPath'].setText(self.parent().pytompath)

        self.binningFactorIMOD = 4

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}

        self.queue_job_names = []

        self.workerID = 0
        self.TMCounter = 0
        self.ECCounter = 0


class KeyPressGraphicsWindow(pg.GraphicsWindow):
    sigKeyPress = pyqtSignal(object)

    def __init__(self, *args, **kwargs):
        super(KeyPressGraphicsWindow,self).__init__(*args, **kwargs)

    def keyPressEvent(self, ev):
        self.scene().keyPressEvent(ev)
        self.sigKeyPress.emit(ev)

    def addViewBox2(self, row=None, col=None, rowspan=1, colspan=1, **kargs):
        """
        Create a ViewBox and place it in the next available cell (or in the cell specified)
        All extra keyword arguments are passed to :func:`ViewBox.__init__ <pyqtgraph.ViewBox.__init__>`
        Returns the created item.
        """
        print('\n\n\ng')
        vb = ViewBoxItems(**kargs)
        self.addItem(vb, row, col, rowspan, colspan)
        return vb


class ViewBoxItems(pg.graphicsItems.ViewBox.ViewBox):

    def __init__(self, *args, **kwargs):
        super(ViewBoxItems, self).__init__(*args, **kwargs)

    def addItems(self, items, ignoreBounds=False):
        """
        Add a QGraphicsItem to this view. The view will include this item when determining how to set its range
        automatically unless *ignoreBounds* is True.
        """

        scene = self.scene()
        for i, item in enumerate(items):
            if item.zValue() < self.zValue():
                item.setZValue(self.zValue() + 1)

            if scene is not None and scene is not item.scene():
                scene.addItem(item)  ## Necessary due to Qt bug: https://bugreports.qt-project.org/browse/QTBUG-18616

            item.setParentItem(self.childGroup)

        if not ignoreBounds:
            self.addedItems += items

        self.updateAutoRange()


class CreateMaskTMOld(QMainWindow, CommonFunctions):
    def __init__(self,parent=None, fname=''):
        super(CreateMaskTM,self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 1000, 1000)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.add_toolbar(self.open_load)
        self.logbook = {}
        self.radius = 50
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0,0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(250, 750), border=True)
        self.leftimage  = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(750, 750), border=True)
        self.centimage  = w.addViewBox(row=0, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(750, 250), border=True)
        self.bottomimage  = w2.addViewBox(row=0, col=0 )
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]



        self.layout.addWidget(w1,0,0)
        self.layout.addWidget(w,0,1)
        self.layout.addWidget(w2,1,1)
        self.layout.addWidget(self.operationbox,1,0)
        self.title = self.parent().parent().widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle( "Create Mask File for: {}".format(os.path.basename(self.title ) ))
        self.centcanvas.wheelEvent = self.wheelEvent
        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        #self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centcanvas.sigMouseReleased.connect(self.empty)
        self.load_image()
        self.leftimage.setXRange(0, self.vol.shape[0])

        self.add_controls(self.layout_operationbox)

        self.subtomo_plots = PlotterSubPlots(self)
        self.subtomo_plots.show()
        QApplication.processEvents()

    def wheelEvent(self, event):

        step = event.angleDelta().y()/120
        increment = int(self.widgets['step_size'].text())*step
        if self.slice+increment < self.vol.shape[0] and self.slice+increment > -1:
            self.update_circles(increment)
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked()==False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = -1*int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.subtomo_plots.close()
            self.close()

    def add_controls(self,prnt):
        vDouble = QDoubleValidator()
        vInt = QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt,'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.',width=100)

        self.insert_label(prnt, text='Size Selection: ',cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ',cstep=-1,rstep=1, value=100, minimum=1, maximum=1000,width=100, stepsize =20)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt,text='Step Size',cstep=1,alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt,wname='step_size',cstep=-1, value=10, rstep=1,width=100,
                            minimum=1, maximum=int(self.vol.shape[0]/4))

        self.insert_label(prnt,text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt,'Number Particles','numSelected',validator=vInt,width=100)

        self.insert_checkbox(prnt, 'show_mask', cstep=1)
        self.insert_label(prnt, text='Show Mask', cstep=0, rstep=1, alignment=Qt.AlignLeft)
        self.widgets['show_mask'].stateChanged.connect(self.show_maskTM)
        #self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def show_maskTM(self):
        w = self.widgets['show_mask']
        width = self.widgets['show_mask'].text()

        if w.isChecked():
            Y,Z,X= np.meshgrid(np.arange(float(self.dim)),np.arange(float(self.dim)),np.arange(float(
                self.dim)))

            if len(self.particleList):
                try:
                    self.mask += 0
                except:
                    self.mask = np.zeros_like(self.vol)

                for x,y,z,r in self.particleList:

                    x = self.dim-x
                    y = self.dim-y
                    z = self.dim-z
                    X -= self.dim - x -0.5
                    Y -= self.dim - y -0.5
                    Z -= self.dim - z -0.5
                    R = np.sqrt(X**2+Y**2+Z**2)
                    self.mask[R<=r] = 1
                    X += self.dim - x -0.5
                    Y += self.dim - y -0.5
                    Z += self.dim - z -0.5

                self.vol = self.volg = self.vol*self.mask

            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

    def open_load(self, q):
        if q.text() == 'Save Project': self.save_particleList()
        elif q.text() == 'Load Project': self.load_particleList()

    def save_particleList(self):
        try: self.mask += 0
        except: 
            print('Please create a mask before saving.')
            return
        import mrcfile
        base, ext = os.path.splitext(self.title)
        fname = os.path.join(self.parent().parent().pickpartfolder, 'mask_' + os.path.basename(self.title.replace(ext, '')))
        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', fname, filter="MRC File (*.mrc);; EM File (*.em)")[0])
        if not fname.endswith('.mrc'): fname += '.mrc'
        if not self.mask.sum(): self.mask +=1
        mrcfile.new(fname, self.mask, overwrite=True)
        self.parent().widgets['filenameMask'].setText(fname)
    
    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "MRC File (*.mrc);; EM File (*.em)")[0])
        if not filename: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol  = read(self.title)
            self.mask = copy.deepcopy( vol2npy(vol) )
            self.mask = self.mask.T

        elif self.title.endswith('mrc'):
            self.mask = read_mrc(self.title)

        if not self.widgets['show_mask'].isChecked():
            self.widgets['show_mask'].setChecked(True)
        else:
            self.show_maskTM()

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        tree = et.parse(self.xmlfile)


        for particle in tree.xpath("Particle"):
            remove = True

            position = []

            for tag in ('X', 'Y', 'Z'):
                position.append(float(particle.xpath('PickPosition')[0].get(tag)))

            x, y, z = position

            for cx, cy, cz, score in self.particleList:
                if abs(x - cx) < 1 or abs(y - cy) < 1 or abs(z - cz) < 1:
                    remove = False
                    break
            if remove:
                self.remove_element(particle)

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
            except:
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def sizeChanged(self):
        a = self.widgets['size_selection: '].text()
        if not a: return


        plist = copy.deepcopy(self.particleList)
        
        self.radius = int(self.widgets['size_selection: '].text()) / 2
        

    def stateScoreChanged(self):
        pass

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

    def replot(self):
        crop = self.vol[int(self.slice), :, :]
        self.img1m.setImage(image=crop.T)

        #self.centcanvas.removeItem(self.hist)
        #self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(np.median(crop)-crop.std()*3,np.median(crop)+crop.std()*3)
        #self.centcanvas.addItem(self.hist)

    def mouseHasMoved(self, evt):

        try:
            pos = self.centimage.mapSceneToView( evt.scenePos() )
        except:
            return

        add = True
        remove = evt.button()==2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x,y,z,s) in enumerate(self.particleList):
            if 1:#np.sqrt( (x-pos.x())**2 + (y-pos.y())**2 + (z-self.slice)**2 ) < self.radius:
                add = False
                if remove:
                    self.remove_point(n-num_deleted_items,z)
                    self.subtomo_plots.delete_subplot([x,y,z])
                    num_deleted_items += 1

        if add and not remove:
            X, Y = pos.x(), pos.y()

            self.add_points(pos, X, Y, self.slice, self.radius, self.radius,add=add,score=self.radius)
            self.subtomo_plots.add_subplot(self.vol, self.particleList[-1])

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def remove_from_coords(self,coords):
        cx,cy,cz = coords[:3]
        for n, (x,y,z,s) in enumerate(self.particleList):
            if np.sqrt( (x-cx)**2 + (y-cy)**2 + (z-cz)**2 ) < self.radius:
                self.remove_point(n, z)
                self.subtomo_plots.delete_subplot([x, y, z])
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0.):

        self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz)), score])
        if radius < 1: return
        pos.setX( cx-radius)
        pos.setY( cy-radius)
        self.circles_cent.append(circle(pos, size=(radius)*2))

        if abs(cz - self.slice) < self.radius:
            self.centimage.addItem(self.circles_cent[-1])

        pos.setX(cx - self.radius)
        pos.setY(cz - self.radius)
        self.circles_bottom.append(circle(pos, size=self.radius * 2))
        if add:
            self.bottomimage.addItem(self.circles_bottom[-1])

        pos.setX(cz - self.radius)
        pos.setY(cy - self.radius)
        self.circles_left.append(circle(pos, size=self.radius * 2))
        if add:
            self.leftimage.addItem(self.circles_left[-1])

    def remove_point(self, n, z, from_particleList=True):
        if from_particleList:
            self.particleList = self.particleList[:n] + self.particleList[n + 1:]

        if abs(z - self.slice) <= self.radius:
            self.centimage.removeItem(self.circles_cent[n])
        self.circles_cent = self.circles_cent[:n] + self.circles_cent[n + 1:]

        self.leftimage.removeItem(self.circles_left[n])
        self.circles_left = self.circles_left[:n] + self.circles_left[n + 1:]

        self.bottomimage.removeItem(self.circles_bottom[n])
        self.circles_bottom = self.circles_bottom[:n] + self.circles_bottom[n + 1:]

    def remove_all(self):
        for image in self.image_list:
            for child in image.allChildren():
                if type(child) == MyCircleOverlay:
                    image.removeItem(child)

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.particleList = []

    def load_image(self):
        if not self.title: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol  = read(self.title)
            self.vol = copy.deepcopy( vol2npy(vol) )
            self.vol = self.vol.T
        elif self.title.endswith('mrc'):
            self.vol = read_mrc(self.title)

        #self.vol[self.vol < -4.] = -4.
        self.backup = self.vol.copy()
        self.origin = self.vol.copy()

        self.vol[ self.vol < self.vol.min()] = self.vol.min()

        self.volg = self.vol
        self.mask = np.zeros_like(self.vol)
        self.dim = self.vol.shape[0]
        self.slice = self.d = int(self.dim / 2)


        self.img1a = pg.ImageItem(self.vol.sum(axis=1))
        self.img1b = pg.ImageItem(self.vol.sum(axis=2))
        self.img1m = pg.ImageItem(self.vol[int(self.slice), :, :])

        self.leftimage.addItem(self.img1a)
        self.centimage.addItem(self.img1m)
        self.bottomimage.addItem(self.img1b)

        self.leftcanvas.setAspectLocked(True)
        self.bottomcanvas.setAspectLocked(True)

        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.centcanvas.addItem(self.hist)

        self.replot_all()

    def update_circles(self, update=0):
        plist = copy.deepcopy(self.particleList)
        self.remove_all()

        for n, (cx, cy, cz, cs) in enumerate(plist):
            
            radius = 0
            if abs(cz - self.slice) <= cs:
                radius = np.sqrt(cs ** 2 - (cz - self.slice)**2)

            if self.xmlfile and len(plist) > 100:
                add = False
            else: add = True
            self.add_points(self.pos, cx, cy, cz, cs, radius,add=add,score=cs)
        self.slice += update


class CreateMaskTM(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, fname=''):
        super(CreateMaskTM, self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 1000, 1000)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.add_toolbar(self.open_load)
        self.logbook = {}
        self.radius = 50
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(250, 750), border=True)
        self.leftimage = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(750, 750), border=True)
        self.centimage = w.addViewBox(row=0, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(750, 250), border=True)
        self.bottomimage = w2.addViewBox(row=0, col=0)
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]

        self.layout.addWidget(w1, 0, 0)
        self.layout.addWidget(w, 0, 1)
        self.layout.addWidget(w2, 1, 1)
        self.layout.addWidget(self.operationbox, 1, 0)
        self.title = self.parent().parent().widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle("Create Mask File for: {}".format(os.path.basename(self.title)))
        self.centcanvas.wheelEvent = self.wheelEvent
        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        # self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centcanvas.sigMouseReleased.connect(self.empty)
        self.load_image()
        self.leftimage.setXRange(0, self.vol.shape[0])

        self.add_controls(self.layout_operationbox)

        QApplication.processEvents()

    def wheelEvent(self, event):

        step = event.angleDelta().y() / 120
        increment = int(self.widgets['step_size'].text()) * step
        if self.slice + increment < self.vol.shape[0] and self.slice + increment > -1:
            self.update_circles(increment)
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked() == False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = -1 * int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.close()

    def add_controls(self, prnt):
        vDouble = QDoubleValidator()
        vInt = QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, text='Size Selection: ', cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ', cstep=-1, rstep=1, value=100, minimum=1, maximum=1000,
                            width=100, stepsize=20)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt, text='Step Size', cstep=1, alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt, wname='step_size', cstep=-1, value=10, rstep=1, width=100,
                            minimum=1, maximum=int(self.vol.shape[0] / 4))

        self.insert_label(prnt, text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt, 'Number Particles', 'numSelected', validator=vInt, width=100)

        self.insert_checkbox(prnt, 'show_mask', cstep=1)
        self.insert_label(prnt, text='Show Mask', cstep=0, rstep=1, alignment=Qt.AlignLeft)
        self.widgets['show_mask'].stateChanged.connect(self.show_maskTM)
        # self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def show_maskTM(self):
        w = self.widgets['show_mask']
        width = self.widgets['show_mask'].text()

        if w.isChecked():
            # Y,Z,X= np.meshgrid(np.arange(float(self.dim)),np.arange(float(self.dim)),np.arange(float(self.dim)))

            if len(self.particleList):
                try:
                    self.mask += 0
                except:
                    self.mask = np.zeros_like(self.vol)

                self.mask *= 0

                for x, y, z, r in self.particleList:
                    dim = r * 2. + 2.
                    Y, Z, X = np.meshgrid(np.arange(dim), np.arange(dim), np.arange(dim))

                    X -= dim / 2 - 0.5
                    Y -= dim / 2 - 0.5
                    Z -= dim / 2 - 0.5
                    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

                    startx, endx = int(max(0, x - dim // 2)), int(min(self.dim, x + dim // 2))
                    starty, endy = int(max(0, y - dim // 2)), int(min(self.dim, y + dim // 2))
                    startz, endz = int(max(0, z - dim // 2)), int(min(self.dim, z + dim // 2))

                    rsx, rex, rsy, rey, rsz, rez = 0, dim, 0, dim, 0, dim

                    if endx - startx < dim:
                        if x - dim // 2 < 0:
                            rsx = dim // 2 - x
                        else:
                            rex = dim // 2 - x + self.dim
                    if endy - starty < dim:
                        if y - dim // 2 < 0:
                            rsy = dim // 2 - y
                        else:
                            rey = dim // 2 - y + self.dim
                    if endz - startz < dim:
                        if z - dim // 2 < 0:
                            rsz = dim // 2 - z
                        else:
                            rez = dim // 2 - z + self.dim

                    rsx, rex, rsy, rey, rsz, rez = map(int, (rsx, rex, rsy, rey, rsz, rez))
                    R = R[rsz:rez, rsy:rey, rsx:rex]

                    self.mask[startz:endz, starty:endy, startx:endx][R <= r] = 1

                self.vol = self.volg = self.vol * (0.25 + self.mask)

            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

    def open_load(self, q):
        if q.text() == 'Save Project':
            self.save_particleList()
        elif q.text() == 'Load Project':
            self.load_particleList()

    def save_particleList(self):
        try:
            self.mask += 0
        except:
            print('Please create a mask before saving.')
            return
        import mrcfile
        base, ext = os.path.splitext(self.title)
        fname = os.path.join(self.parent().parent().pickpartfolder,
                             'mask_' + os.path.basename(self.title.replace(ext, '')))
        fname = str(
            QFileDialog.getSaveFileName(self, 'Save particle list.', fname, filter="MRC File (*.mrc);; EM File (*.em)")[
                0])
        if not fname.endswith('.mrc'): fname += '.mrc'
        if not self.mask.sum(): self.mask += 1
        mrcfile.new(fname, self.mask, overwrite=True)
        self.parent().widgets['filenameMask'].setText(fname)

    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "MRC File (*.mrc);; EM File (*.em)")[0])
        if not filename: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol = read(self.title)
            self.mask = copy.deepcopy(vol2npy(vol))
            self.mask = self.mask.T

        elif self.title.endswith('mrc'):
            self.mask = read_mrc(self.title)

        if not self.widgets['show_mask'].isChecked():
            self.widgets['show_mask'].setChecked(True)
        else:
            self.show_maskTM()

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        tree = et.parse(self.xmlfile)

        for particle in tree.xpath("Particle"):
            remove = True

            position = []

            for tag in ('X', 'Y', 'Z'):
                position.append(float(particle.xpath('PickPosition')[0].get(tag)))

            x, y, z = position

            for cx, cy, cz, score in self.particleList:
                if abs(x - cx) < 1 or abs(y - cy) < 1 or abs(z - cz) < 1:
                    remove = False
                    break
            if remove:
                self.remove_element(particle)

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
            except:
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

        #self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def sizeChanged(self):
        a = self.widgets['size_selection: '].text()
        if not a: return

        plist = copy.deepcopy(self.particleList)

        self.radius = int(self.widgets['size_selection: '].text()) / 2

    def stateScoreChanged(self):
        pass

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

    def replot(self):
        crop = self.vol[int(self.slice), :, :]
        self.img1m.setImage(image=crop.T)

        # self.centcanvas.removeItem(self.hist)
        # self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(np.median(crop) - crop.std() * 3, np.median(crop) + crop.std() * 3)
        # self.centcanvas.addItem(self.hist)
        del crop

    def mouseHasMoved(self, evt):

        try:
            pos = self.centimage.mapSceneToView(evt.scenePos())
        except:
            return

        add = True
        remove = evt.button() == 2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x, y, z, s) in enumerate(self.particleList):
            if np.sqrt((x - pos.x()) ** 2 + (y - pos.y()) ** 2 + (z - self.slice) ** 2) < self.radius:
                add = False
                if remove:
                    self.remove_point(n - num_deleted_items, z)

                    num_deleted_items += 1
            add=True

        if add and not remove:
            X, Y = pos.x(), pos.y()
            self.particleList.append([int(round(X)), int(round(Y)), int(round(self.slice)), self.radius])

            self.add_points(pos, X, Y, self.slice, self.radius, self.radius, add=add, score=self.radius)

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def remove_from_coords(self, coords):
        cx, cy, cz = coords[:3]
        for n, (x, y, z, s) in enumerate(self.particleList):
            if np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) < self.radius:
                self.remove_point(n, z)
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0., new=True):

        if new:
            self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz)), score])

        if radius < 1: return

        pos.setX(cx - radius)
        pos.setY(cy - radius)

        if abs(cz - self.slice) < self.radius:
            self.circles_cent.append(circle(pos, size=(radius) * 2))
            self.centimage.addItem(self.circles_cent[-1])

        if new:
            pos.setX(cx - self.radius)
            pos.setY(cz - self.radius)
            self.circles_bottom.append(circle(pos, size=self.radius * 2))
            self.bottomimage.addItem(self.circles_bottom[-1])

            pos.setX(cz - self.radius)
            pos.setY(cy - self.radius)
            self.circles_left.append(circle(pos, size=self.radius * 2))
            self.leftimage.addItem(self.circles_left[-1])

    def remove_point(self, n, z, from_particleList=True):

        if n < len(self.circles_cent):
            if abs(z - self.slice) <= self.radius:
                self.centimage.removeItem(self.circles_cent[n])
                self.circles_cent.pop(n)  # = self.circles_cent[:n] + self.circles_cent[n + 1:]

                self.leftimage.removeItem(self.circles_left[n])
                self.circles_left.pop(n)  # = self.circles_left[:n] + self.circles_left[n + 1:]

                self.bottomimage.removeItem(self.circles_bottom[n])
                self.circles_bottom.pop(n)  # = self.circles_bottom[:n] + self.circles_bottom[n + 1:]

                if from_particleList:
                    self.particleList.pop(n)
                    # self.particleList = self.particleList[:n] + self.particleList[n + 1:]

    def remove_all(self):
        for image in self.image_list:
            for child in image.allChildren():
                if type(child) == MyCircleOverlay:
                    image.removeItem(child)

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.particleList = []

    def remove_cent(self):
        for child in self.centimage.allChildren():
            if type(child) == MyCircleOverlay:
                self.centimage.removeItem(child)

        for i in range(len(self.circles_cent)):
            self.circles_cent.pop(0)

    def load_image(self):
        if not self.title: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol = read(self.title)
            self.vol = copy.deepcopy(vol2npy(vol))
            self.vol = self.vol.T
        elif self.title.endswith('mrc'):
            self.vol = read_mrc(self.title)

        # self.vol[self.vol < -4.] = -4.
        self.backup = self.vol.copy()
        self.origin = self.vol.copy()

        self.vol[self.vol < self.vol.min()] = self.vol.min()

        self.volg = self.vol
        self.mask = np.zeros_like(self.vol)
        self.dim = self.vol.shape[0]
        self.slice = self.d = int(self.dim / 2)

        self.img1a = pg.ImageItem(self.vol.sum(axis=1))
        self.img1b = pg.ImageItem(self.vol.sum(axis=2))
        self.img1m = pg.ImageItem(self.vol[int(self.slice), :, :])

        self.leftimage.addItem(self.img1a)
        self.centimage.addItem(self.img1m)
        self.bottomimage.addItem(self.img1b)

        self.leftcanvas.setAspectLocked(True)
        self.bottomcanvas.setAspectLocked(True)

        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.centcanvas.addItem(self.hist)

        self.replot_all()

    def update_circles(self, update=0):
        self.slice += update
        self.remove_cent()

        for n, (cx, cy, cz, cs) in enumerate(self.particleList):
            radius = 0
            if abs(cz - self.slice) <= cs:
                radius = np.sqrt(cs ** 2 - (cz - self.slice) ** 2)

            add = False
            self.add_points(self.pos, cx, cy, cz, cs, radius, add=add, score=cs, new=False)

        self.centimage.addItems(self.circles_cent)


class InfiniteLinePlotter(pg.InfiniteLine):
    pass


class LinearRegionItem(pg.LinearRegionItem):
    def lineMoveFinished(self):

        try:
            if self.saveZLimits:
                outfile = open(self.filename, 'w')
                a, b = int(self.lines[0].value()), int(self.lines[1].value())
                outfile.write(f'{min(a,b)} {max(a,b)}')
                outfile.close()
                print(f'written {self.filename}')
        except Exception as e:
            print(e)
            print('file NOT written')

        self.other.lines[0].setValue(self.lines[0].value())
        self.other.lines[1].setValue(self.lines[1].value())
        self.sigRegionChangeFinished.emit(self)


class ParticlePicker(QMainWindow, CommonFunctions):
    def __init__(self,parent=None, fname=''):
        super(ParticlePicker,self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 800, 800)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.add_toolbar(self.open_load, saveText='Save Particle List', openText='Open Particle List')
        self.logbook = {}
        self.radius = 8
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0,0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.txtfile = ''
        self.filetype = 'txt'
        self.activate = 0
        self.activateLine = 0
        self.angleLine = 0
        self.red_circle = 0
        self.mask = None

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(200, 600), border=True)
        self.leftimage  = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(600, 600), border=True)
        self.datalabel = pg.LabelItem(justify='right')
        self.centcanvas.addItem(self.datalabel)
        self.centimage  = w.addViewBox2(row=1, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)

        self.target = w3 = pg.ImageView()
        self.title = parent.widgets['v03_manualpp_tomogramFname'].text()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(600, 200), border=True)
        self.bottomimage  = w2.addViewBox(row=0, col=0 )
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]

        self.layout.addWidget(w1,0,0)
        self.layout.addWidget(w,0,1)
        self.layout.addWidget(w2,1,1)
        self.layout.addWidget(self.operationbox,1,0)
        self.title = parent.widgets['v03_manualpp_tomogramFname'].text()
        self.TOMNAME = self.title
        if not self.title: self.title = 'Dummy Data'

        self.folder = os.path.dirname(self.title)
        self.setWindowTitle("Manual Particle Selection From: {}".format( os.path.basename(self.title)))
        self.centcanvas.wheelEvent = self.wheelEvent

        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        self.leftimage.scene().sigMouseClicked.connect(self.mouseHasMovedLeft)
        self.bottomimage.scene().sigMouseClicked.connect(self.mouseHasMovedBottom)

        self.centcanvas.sigMouseReleased.connect(self.empty)

        self.slice =0
        self.pen = pg.mkPen(color='#117a42', width=2)
        self.lineB = InfiniteLinePlotter([self.slice, self.slice], angle=0, movable=False, pen=self.pen)
        self.lineL = InfiniteLinePlotter([self.slice, self.slice], angle=90, movable=False, pen=self.pen)


        self.load_image()
        self.leftimage.setXRange(0, self.vol.shape[0])

        self.add_controls(self.layout_operationbox)

        self.subtomo_plots = PlotterSubPlots(self, size_subtomo=self.radius*2)
        self.subtomo_plots.show()
        self.subtomo_plots.hide()

        self.centimage.scene().sigMouseMoved.connect(self.updateLabel)

        self.bottomimage.addItem(self.lineB)
        self.leftimage.addItem(self.lineL)

        self.replot()
        QApplication.processEvents()

    def closeEvent(self, event):

        close = QMessageBox()
        close.setText("Are you sure you want to close the particle picking windows?")
        close.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
        close = close.exec()

        if close == QMessageBox.Yes:
            self.subtomo_plots.destroy()
            event.accept()

        else:
            event.ignore()

    def wheelEvent(self, event):
        try:
            step = event.angleDelta().y()/120
            increment = int(self.widgets['step_size'].text())*step
            if self.activateLine:
                self.angleLine = (self.angleLine + step)
                if (self.angleLine) > 90: self.angleLine -= 180
                if (self.angleLine) < -90: self.angleLine += 180

                self.line.setAngle(self.angleLine)

            elif self.slice+increment < self.vol.shape[0] and self.slice+increment > -1:
                self.slice += increment
                self.update_circles()
                self.replot()
        except Exception as e:
            print(e)

    def keyPressEvent(self, evt):
        if Qt.Key_P == evt.key():
            self.subtomo_plots.close()
            self.subtomo_plots.show()

        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked()==False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.slice += update
                self.replot()
                self.update_circles()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = int(self.widgets['step_size'].text())
                self.slice -= update
                self.replot()
                self.update_circles()

        if evt.key() == Qt.Key_Escape:
            self.close()

        if evt.key() == Qt.Key_A:
            self.activate = 1 - self.activate
            if self.activate:
                self.adjustZLimits()
            else:
                self.resetZLimits()
                self.replot_all()

        if evt.key() == Qt.Key_L:
            try:
                self.activateLine = 1 - self.activateLine
                if self.activateLine:
                    self.drawLine()
                else:
                    self.drawLine(insert=False)
            except Exception as e:
                print(e)

        # if evt.key() == Qt.Key_K:
        #     try:
        #         self.activateLine = 1 - self.activateLine
        #         if self.activateLine:
        #             self.drawLine(id=1)
        #         else:
        #             self.drawLine(insert=False)
        #     except Exception as e:
        #         print(e)

    def add_controls(self,prnt):
        vDouble = QDoubleValidator()
        vInt = QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt,'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.',width=100)

        self.insert_label(prnt, text='Size Selection: ',cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ',cstep=-1,rstep=1, value=16, minimum=1, maximum=1000,width=100)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt,text='Step Size',cstep=1,alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt,wname='step_size',cstep=-1, value=1, rstep=1,width=100,
                            minimum=1, maximum=int(self.vol.shape[0]/4))

        self.insert_label(prnt,text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt,'# Selected Particles','numSelected',validator=vInt,width=100, enabled=False)
        self.insert_label_spinbox(prnt,'minScore','Minimal Score', value=0.01,wtype=QDoubleSpinBox, minimum=0, maximum=1,
                                  stepsize=0.05, width=100, decimals=4)
        self.insert_label_spinbox(prnt, 'maxScore', 'Maximal Score', value=1., wtype=QDoubleSpinBox, minimum=0,
                                  maximum=10., stepsize=0.05, width=100, decimals=4, cstep=-2)

        self.widgets['minScore'].valueChanged.connect(self.stateScoreChanged)
        self.widgets['maxScore'].valueChanged.connect(self.stateScoreChanged)
        
        self.insert_checkbox(prnt, 'apply_mask', cstep=1)
        self.insert_label(prnt, text='Apply mask', cstep=0, rstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt, 'filenameMask', cstep=1, rstep=-1, value='')
        params = ['file', self.widgets['filenameMask'], ['mrc', 'em'], False, self.parent().pickpartfolder]
        self.insert_pushbutton(prnt, 'Browse', action=self.browse, params=params, width=100)
        # self.insert_pushbutton(prnt, 'Create', rstep=1, cstep=-1, action=self.gen_mask)
        self.widgets['apply_mask'].stateChanged.connect(self.stateMaskChanged)

    def gen_mask(self):
        self.genmask = CreateMaskTM(self)
        self.genmask.show()

    def stateMaskChanged(self):
        w = self.widgets['apply_mask']
        fname = self.widgets['filenameMask'].text()

        if w.isChecked():

            if os.path.exists(fname) and os.path.splitext(fname)[1] in ['.mrc', '.em']:

                # read as em or mrc file
                if fname.endswith('.em'):
                    from pytom.basic.files import read
                    from pytom_numpy import vol2npy
                    vol = read(fname)
                    self.mask = copy.deepcopy(vol2npy(vol))
                    self.mask = self.mask.T  # why transpose ?

                elif self.title.endswith('mrc'):
                    self.mask = read_mrc(fname)

                if self.mask is not None:
                    # check that the shape is the same as tomogram
                    if not all([s1 == s2 for s1, s2 in zip(self.mask.shape, self.vol.shape)]):
                        self.mask = None

                        # give user warning box
                        msg = QMessageBox()
                        msg.setText("Mask dimensions are not matching with tomogram")
                        msg.setStandardButtons(QMessageBox.Ok)
                        msg.exec()

                    else:  # make mask binary if not already binary
                        self.mask = (self.mask - self.mask.min()) / (self.mask.max() - self.mask.min())
                        self.mask = ((self.mask >= 0.5) * 1).astype(np.int8)

            else:
                # give user warning box
                msg = QMessageBox()
                msg.setText("Invalid tomogram mask file")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec()

        else:
            self.mask = None

        if self.filetype == 'txt':
            if os.path.exists(self.txtfile):
                self.load_txtFile(self.txtfile)
        elif self.filetype == 'xml':
            if os.path.exists(self.xmlfile):
                self.load_xmlFile(self.xmlfile)

    def open_load(self, q):
        if q.text() == 'Save Particle List': self.save_particleList()
        elif q.text() == 'Open Particle List': self.load_particleList()

    def save_particleList(self):

        if self.activateLine == 1:
            self.save_angle()
            return

        if self.filetype == 'xml':
            self.save_particleList_xml()
        else:
            self.save_particleList_txt()

    def save_particleList_xml(self):

        try:
            self.remove_deselected_particles_from_XML()
        except Exception as e:
            print('particleList not saved.')
            print(e)

    def save_particleList_txt(self):
        ext = '.'+os.path.basename(self.title).split('.')[-1]
        fname = os.path.join(self.parent().pickpartfolder,'coords_' + os.path.basename(self.title.replace(ext,'')))

        fname = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.txt')[0])
        if not fname.endswith('.txt'): fname += '.txt'

        if fname:

            folder = os.path.dirname(self.title)
            if 'INFR' in os.path.basename(self.title):
                key = 'INFR'
            else:
                key ='WBP'

            inputJobName = "cat {}/{}_reconstruction.sh | grep '{}' | {}"


            from pytom.bin.templateMatchingCandidateExtractionSingleGPU import getBinningFactorAndReferenceMarker

            binning, refmark = getBinningFactorAndReferenceMarker(self.title)


            out = open(fname, 'w')
            out.write('#\n')
            out.write('#TOMONAME \t{}\n'.format(self.TOMNAME))
            out.write(f'#MARKERINDEX \t{binning}\n')
            out.write(f'#BINNINGFACTOR \t{refmark}\n')
            out.write('#\n')

            for x, y, z, s in self.particleList:
                outtext =  '{:8.0f} {:8.0f} {:8.0f}\n'.format(x,y,z)
                out.write(outtext)

    def load_particleList(self):
        initdir = self.parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir,
                                                   "PyTom ParticleList(*.xml);; Coordinate File (*.txt)")[0])
        if not filename: return

        if filename.endswith('.txt'):
            self.filetype = 'txt'
            self.load_txtFile(filename)
            self.txtfile = filename

        elif filename.endswith('.xml'):
            self.filetype = 'xml'
            self.load_xmlFile(filename)
            self.xmlfile = filename

        self.subtomo_plots.close()
        self.subtomo_plots.show()

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        from pytom.basic.structures import ParticleList, Particle, PickPosition
        from pytom.basic.score import FLCFScore
        from pytom.agnostic.io import read
        import numpy
        import random

        asked = False

        cc = 180./ np.pi

        plOld = ParticleList()
        plOld.fromXMLFile(self.xmlfile)

        folder = os.path.dirname(plOld[0].getFilename())
        wedge = plOld[0].getWedge()
        origin = plOld[0].getPickPosition().getOriginFilename()

        plNew = ParticleList()

        tempPL = []

        prefix = 'man_slct_'
        id = 0

        anglelist = read(f'{self.pytompath}/angles/angleLists/angles_11_15192.em')

        for p in plOld:
            fname = os.path.basename(p.getFilename())
            if prefix in fname:
                tid = int(fname.split('.')[0][len(prefix):])
                id = id if tid < id else tid

            x,y,z = p.getPickPosition().toVector()
            print(x, x.__class__)
            s = p.getScore().getValue()
            tempPL.append([x,y,z,s])

        for x,y,z,s in self.particleList:
            found =False

            for n,[cx,cy,cz,s0] in enumerate(tempPL):
                if abs(x - cx) < 1 and abs(y - cy) < 1 and abs(z - cz) < 1:
                    print(n)
                    plNew.append(plOld[n])
                    found = True
                    break


            if not asked and not found:
                close = QMessageBox()
                close.setText("Do you want to save the manually added particles to your particleList?")
                close.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
                close = close.exec()

                if close == QMessageBox.Yes:
                    add = True
                else:
                    add = False

                asked = True

            if not found and add:
                print(id)
                rot = list(random.choice(anglelist) * cc)
                particle = Particle(filename=f'{folder}/{prefix}{id}.em', rotation=rot, shift=[0,0,0], wedge=wedge,
                                    className = 0, pickPosition=PickPosition(x=x,y=y,z=z, originFilename=origin),
                                    score=FLCFScore(1E-6))
                plNew.append(particle)
                id += 1

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                plNew.toXMLFile(fname)
            except Exception as e:
                print(e)
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

    def load_xmlFile(self, filename):

        self.remove_all()

        xmlObj = et.parse(filename)
        particles = xmlObj.xpath('Particle')

        self.slice = int(self.vol.shape[0] / 2)
        dz, dy, dx = self.vol.shape
        for p in particles:
            try:
                score = float( p.xpath('Score')[0].get('Value') )
            except:
                score = 0.5

            if self.max_score >= score >= self.min_score:
                
                x,y,z = map(float, (p.xpath('PickPosition')[0].get('X'),
                                    p.xpath('PickPosition')[0].get('Y'),
                                    p.xpath('PickPosition')[0].get('Z')))

                if self.mask is not None:
                    if not self.mask[int(x), int(y), int(z)]:
                        continue

                if abs(dx/2-x) > dx/2 or abs(dy/2-y) > dy/2 or abs(dz/2-z) > dz/2:
                    print(f'particle at ({x},{y},{z}) not added')
                    continue

                self.particleList.append([int(round(x)), int(round(y)), int(round(z)), score])

        z_sorted_ids = sorted(range(len(self.particleList)), key=lambda k: self.particleList[k][2])

        # initialize the circles
        if len(self.particleList):
            self.update_circles()

        self.widgets['numSelected'].setText(str(len(self.particleList)))
        self.replot()
        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def load_txtFile(self, filename):

        self.remove_all()

        particlePos = [map(float, line.split()) for line in open(filename).readlines()
                       if line != '' and not '#' in line]

        self.particleList = []

        self.slice = int(self.vol.shape[0] / 2)
        dz, dy, dx = self.vol.shape

        for x,y,z in particlePos:

            if self.mask is not None:
                if not self.mask[int(x), int(y), int(z)]:
                    continue

            if abs(dx / 2 - x) > dx / 2 or abs(dy / 2 - y) > dy / 2 or abs(dz / 2 - z) > dz / 2:
                print(f'particle at ({x},{y},{z}) not added')
                continue

            self.particleList.append([int(round(x)), int(round(y)), int(round(z)), self.radius])
        #     self.z_sorted.append(int(round(z)))
        #
        # z_sorted = sorted(self.sorted)
        # z_sorted_ids = sorted(range(len(self.particleList)), key=lambda k: self.particleList[k][2])

        if len(self.particleList):
            self.update_circles()

        self.widgets['numSelected'].setText(str(len(self.particleList)))
        self.replot()
        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def sizeChanged(self):
        a = self.widgets['size_selection: '].text()

        if not a:
            return

        self.radius = int(self.widgets['size_selection: '].text()) / 2

        self.update_circles()

        self.subtomo_plots.size_subtomo = self.radius*2
        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def stateScoreChanged(self):
        self.min_score = float(self.widgets['minScore'].value())
        self.max_score = float(self.widgets['maxScore'].value())
        if self.filetype == 'txt':
            if os.path.exists(self.txtfile):
                self.load_txtFile(self.txtfile)
        elif self.filetype == 'xml':
            if os.path.exists(self.xmlfile):
                self.load_xmlFile(self.xmlfile)

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

    def replot(self):
        crop = self.vol[int(self.slice), :, :]
        self.img1m.clear()
        self.img1m.setImage(image=crop.T, autoDownsample=True)

        #self.centcanvas.removeItem(self.hist)
        #self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(np.median(crop)-crop.std()*3,np.median(crop)+crop.std()*3)
        self.writeLabelText(z=self.slice)

        if self.activateLine == 0 and self.activate == 0:
            self.lineB.setValue(self.slice)
            self.lineL.setValue(self.slice)

        #self.centcanvas.addItem(self.hist)

    def mouseHasMoved(self, evt):

        try:
            pos = self.centimage.mapSceneToView( evt.scenePos() )
        except:
            return

        add = True
        remove = evt.button()==2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[2] or pos.y() >= self.vol.shape[1]:
            #print(pos.x(), pos.y())
            return

        num_deleted_items = 0
        for n, (x,y,z,s) in enumerate(self.particleList):
            if np.sqrt( (x-pos.x())**2 + (y-pos.y())**2 + (z-self.slice)**2 ) < self.radius:
                add = False
                if remove:
                    self.remove_point(n-num_deleted_items,z)
                    self.subtomo_plots.delete_subplot([x,y,z])
                    self.update_circles()

                    num_deleted_items += 1

        if add and not remove:
            X, Y = pos.x(), pos.y()
            self.particleList.append([int(round(X)), int(round(Y)), int(round(self.slice)), self.radius])
            self.update_circles()
            #self.add_points(pos, X, Y, self.slice, 0, self.radius,add=add)
            self.subtomo_plots.add_subplot(self.vol, self.particleList[-1])

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def updateLabel(self, pos):

          ## using signal proxy turns original arguments into a tuple
        pos = self.centimage.mapSceneToView(pos)
        if pos.x() >= 0 and pos.y() >= 0 and pos.x() < self.vol.shape[2] and pos.y() < self.vol.shape[1]:

            self.datalabel.setText(f"value = {self.vol[int(self.slice)][int(pos.y())][int(pos.x())]:7.3f} x = {pos.x():6.0f} y={pos.y():6.0f} z = {self.slice:4.0f}")

    def writeLabelText(self,v=0, x=0,y=0,z=0):
        self.datalabel.setText(f"value = {v:7.3f} x = {x:6.0f} y={y:6.0f} z = {z:4.0f}")

    def mouseHasMovedBottom(self, evt):
        pos = self.bottomimage.mapSceneToView( evt.scenePos() )
        if pos.y() < 0 or pos.y() >= self.vol.shape[2]:
            return
        step = pos.y() - self.slice
        self.slice += step
        self.update_circles()
        self.replot()

    def mouseHasMovedLeft(self, evt):
        pos = self.leftimage.mapSceneToView( evt.scenePos() )
        if pos.x() < 0 or pos.x() >= self.vol.shape[1]: return

        step = pos.x()-self.slice
        self.slice += step
        self.update_circles()
        self.replot()

    def remove_from_coords(self, coords):
        cx,cy,cz = coords[:3]
        for n, (x,y,z,s) in enumerate(self.particleList):
            if np.sqrt( (x-cx)**2 + (y-cy)**2 + (z-cz)**2 ) < self.radius:
                self.remove_point(n, z)
                self.subtomo_plots.delete_subplot([x, y, z])
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0.):

        pos.setX(cx - radius)
        pos.setY(cy - radius)
        self.circles_cent += [circle(pos, size=(radius) * 2)]

        if add:
            pos.setX(cx - self.radius)
            pos.setY(cz - self.radius)
            self.circles_bottom.append(circle(pos, size=self.radius * 2))
            self.bottomimage.addItem(self.circles_bottom[-1])

            pos.setX(cz - self.radius)
            pos.setY(cy - self.radius)
            self.circles_left.append(circle(pos, size=self.radius * 2))
            self.leftimage.addItem(self.circles_left[-1])

    def remove_point(self, n, z, from_particleList=True):
        if from_particleList:
            self.particleList.pop(n)

    def remove_all(self, keepPL=False):
        [self.centimage.removeItem(child) for child in self.circles_cent]
        [self.bottomimage.removeItem(child) for child in self.circles_bottom]
        [self.leftimage.removeItem(child) for child in self.circles_left]
        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        if not keepPL:
            self.particleList = []

    def load_image(self):
        if not self.title: return
        self.tomogram_name = self.title

        try:
            if self.title.endswith('em'):
                from pytom.basic.files import read
                from pytom_numpy import vol2npy
                vol  = read(self.title)
                self.vol = copy.deepcopy( vol2npy(vol) )
                self.vol = self.vol.T
            elif self.title.split('.')[-1] in ('mrc', 'mrcs', 'rec', 'st', 'map'):
                self.vol = read_mrc(self.title)

            # not really necessary...
            #self.vol[self.vol < -4.] = -4.
            self.backup = self.vol.copy()
            self.origin = self.vol.copy()

            self.vol[ self.vol < self.vol.min()] = self.vol.min()

            self.dim = self.vol.shape[0]
            self.slice = self.d = int(self.dim / 2)


            self.img1a = pg.ImageItem(self.vol.sum(axis=2))
            self.img1b = pg.ImageItem(self.vol.sum(axis=1).T)
            self.img1m = pg.ImageItem(self.vol[int(self.slice), :, :])

            self.resetZLimits(True)
            self.leftimage.addItem(self.img1a)
            self.centimage.addItem(self.img1m)
            self.bottomimage.addItem(self.img1b)

            self.leftcanvas.setAspectLocked(True)
            self.bottomcanvas.setAspectLocked(True)

            self.hist = pg.HistogramLUTItem()
            self.hist.setImageItem(self.img1m)
            self.centcanvas.addItem(self.hist, row=1, col=1)

            self.setAngle()

            self.replot_all()

        except Exception as e:
            print(f'ERROR: {e}')
            self.close()

    def setAngle(self):
        try:
            folder = os.path.dirname(os.popen(f'ls -alrt {self.tomogram_name}').read().split()[-1])

            specimen_angle_file = os.path.join(folder, 'specimen_rotation_angle.txt')
            if os.path.exists(specimen_angle_file):
                angle = float(open(specimen_angle_file,'r').read()[:-1])
            else:
                angle=0
        except Exception as e:
            print(e)
            angle=0
        folder = os.path.dirname(os.popen(f'ls -alrt {self.tomogram_name}').read().split()[-1])
        if os.path.exists(os.path.join(folder, 'WBP_Reconstruction.sh')):
            d = open(os.path.join(folder, 'WBP_Reconstruction.sh'), 'r').read()[:-1]
            try:
                used_angle = float(d.split('specimenAngle ')[-1].split()[0])
            except Exception as e:
                print('Error: ', e)
                used_angle = 0
        else:
            used_angle = 0

        self.angleLine = (angle-used_angle)*-1

    def update_circles(self):

        if type(self.red_circle) == MyCircleOverlay:
            self.centimage.removeItem(self.red_circle)
            self.red_circle = 0

        nn = 0
        added = 0
        num_circles = len(self.circles_cent)

        # if the points were ordered by z height, the loop could be faster
        # also coordinates could be stored in numpy array to speed up.
        # might also be that the ROI init from pyqtgraph is the bottleneck, because the
        # loop itself should be relatively fast for ~500 particles.
        for n, (cx, cy, cz, cs) in enumerate(self.particleList):

            radius = self.radius

            # if n % 100 == 0:
            #     print(n)

            # this check if particle n has a circle that should be plotted in this slice
            if abs(cz - self.slice) >= self.radius:
                continue

            # if it should be plotted, calculate radius of this slice of the sphere
            if abs(cz - self.slice) < self.radius:
                radius = np.sqrt(self.radius ** 2 - (cz - self.slice)**2)

            # if we have a particle list larger than 100 we should not add to side views to save time
            # but why do we need to update this while scrolling, could be set immediately
            if self.xmlfile and len(self.particleList) > 100:
                add = False
            else:
                add = True

            # if nn is smaller than number of circles, update one of the circles
            if nn < num_circles:
                self.circles_cent[nn].setPos(cx-radius, cy-radius, update=False)
                self.circles_cent[nn].setSize(radius*2)
            else:  # else add a new circle
                self.add_points(self.pos, cx, cy, cz, cs, radius, add=add)
                added += 1

            # nn counts if the particle sphere should be displayed this slice
            nn += 1

        # only add the new circles to the image
        self.centimage.addItems(self.circles_cent[-added:])

        # remove = 0
        if nn < num_circles:
            # is it faster to remove the circles ?
            for circle in self.circles_cent[nn:]:
                circle.setSize(0.001)

    def adjustZLimits(self):
        self.replot_all()
        self.resetZLimits(True)
        for r in [self.rgnleft, self.rgnbott]:
            r.saveZLimits = True
            r.filename = os.path.join(self.folder, 'z_limits.txt')

    def resetZLimits(self, insert=False):

        for image in [self.leftimage, self.bottomimage]:
            for child in image.allChildren():
                if type(child) == LinearRegionItem or type(child) == InfiniteLinePlotter:
                    image.removeItem(child)
        if not insert:
            self.lineB = InfiniteLinePlotter([self.slice, self.slice], angle=0, movable=False, pen=self.pen)
            self.lineL = InfiniteLinePlotter([self.slice, self.slice], angle=90, movable=False, pen=self.pen)
            self.bottomimage.addItem(self.lineB)
            self.leftimage.addItem(self.lineL)
            return
        try:
            tt = os.path.dirname(os.popen(f'ls -alrt {self.tomogram_name}').read().split()[-1])
            self.folder = tt
            zlimitfile = os.path.join(tt, 'z_limits.txt')

            if os.path.exists(zlimitfile):
                z_start, z_end = list(map(int, open(zlimitfile).readlines()[0].split()[:2]))
            else:
                z_start, z_end = 0, self.vol.shape[0]
                zlimitfile = ''
            self.rgnleft = LinearRegionItem([z_start, z_end])
            self.rgnbott = LinearRegionItem(values=(z_start, z_end), orientation=1)
            self.rgnbott.filename = zlimitfile
            self.rgnleft.filename = zlimitfile
            self.rgnleft.other = self.rgnbott
            self.rgnbott.other = self.rgnleft
            self.leftimage.addItem(self.rgnleft)
            self.bottomimage.addItem(self.rgnbott)

        except Exception as e:
            print('\n\n\n', e)

    def drawLine(self, id=0, insert=True):
        for image in [self.leftimage, self.bottomimage]:
            for child in image.allChildren():
                if type(child) == InfiniteLinePlotter:
                    image.removeItem(child)
        if not insert:
            self.lineB = InfiniteLinePlotter([self.slice, self.slice], angle=0, movable=False, pen=pg.mkPen(color='#117a42', width=2))
            self.lineL = InfiniteLinePlotter([self.slice, self.slice], angle=90, movable=False,
                                             pen=pg.mkPen(color='#117a42', width=2))

            self.bottomimage.addItem(self.lineB)
            self.leftimage.addItem(self.lineL)

            return
        try:

            self.line =  InfiniteLinePlotter([int(self.dim / 2),int(self.dim / 2)], angle=self.angleLine, movable=True )

            if 0 == id:
                self.bottomimage.addItem(self.line)

            else:
                self.leftimage.addItem(self.line)

        except Exception as e:
            print('\n\n\n', e)

    def save_angle(self):
        tt = os.path.dirname(os.popen(f'ls -alrt {self.tomogram_name}').read().split()[-1])
        self.folder = tt
        sra_file = os.path.join(tt, 'specimen_rotation_angle.txt')
        if os.path.exists(os.path.join(os.path.dirname(sra_file), 'WBP_Reconstruction.sh')):
            d = open(os.path.join(os.path.dirname(sra_file), 'WBP_Reconstruction.sh'), 'r').read()[:-1]
            try:
                used_angle = float(d.split('specimenAngle ')[-1].split()[0])
            except Exception as e:
                print(e)
                used_angle = 0
        else:
            used_angle = 0

        f = open(sra_file, 'w')
        f.write(f'{(used_angle + self.angleLine * -1):.1f}\n')

        f.close()

        self.popup_messagebox('Info', 'Angle file written', f'Specimen Angle {self.angleLine*-1 + used_angle} '
                                                            f'written in {os.path.basename(sra_file)}.\n'
                                                            f'This includes the specimen angle {used_angle} used for the reconstruction.' )


class QParams():
    """
    Each update is structured in the same way:
    - If a value is provided that will be used to directly update that attribute of the structure.
    - If no value is provided, parent and mode will be used to update the attribute.
    - If parent is provided all the settings will be updated. This means a pickle will be stored, and the qparams for
    pytomGUI head structure will be updated.
    - If current_job is All, this will update the specified parameter for all jobs.
    """
    # TODO custom header needs to be stored as well
    def __init__(self, time=12, queue='defq', nodes=1, cores=20, modules=[], command=''):
        self.time = time
        self.queue = queue
        self.nodes = nodes
        self.cores = cores
        self.modules = modules
        self.command = command

    def update_time(self, parent=None, mode=None, value=None, current_job=''):
        if value is None:
            self.time = parent.widgets[mode + 'maxTime'].value()
        else:
            self.time = value

        if current_job == 'All':
            for key_other in parent.qparams.keys():
                parent.qparams[key_other].update_time(value=self.time)
        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_queue(self, parent=None, mode=None, value=None, current_job=''):
        if value is None:
            self.queue = parent.widgets[mode + 'queueName'].text()
        else:
            self.queue = value

        if current_job == 'All':
            for key_other in parent.qparams.keys():
                parent.qparams[key_other].update_queue(value=self.queue)
        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_nodes(self, parent=None, mode=None, value=None, current_job=''):
        if value is None:
            self.nodes = parent.widgets[mode + 'numberOfNodes'].value()
        else:
            self.nodes = value

        if current_job == 'All':
            for key_other in parent.qparams.keys():
                parent.qparams[key_other].update_nodes(value=self.nodes)
        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_cores(self, parent=None, mode=None, value=None, current_job=''):
        if value is None:
            self.cores = parent.widgets[mode + 'numberOfCores'].value()
        else:
            self.cores = value

        if current_job == 'All':
            for key_other in parent.qparams.keys():
                parent.qparams[key_other].update_cores(value=self.cores)
        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_command(self, parent=None, mode=None, value=None, current_job=''):
        if value is None:
            self.command = parent.widgets[mode + 'command'].text()
        else:
            self.command = value

        if current_job == 'All':
            for key_other in parent.qparams.keys():
                parent.qparams[key_other].update_command(value=self.command)
        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_modules(self, parent=None, mode=None, value=None):
        if value is None:
            self.modules = parent.widgets[mode + 'modules'].getModules()
        else:
            self.modules = value

        if parent is not None:
            # update everything
            self.update_settings(parent)

    def update_settings(self, parent, store=True):
        if store:
            with open(os.path.join(parent.projectname, '.qparams.pickle'), 'wb') as handle:
                pickle.dump(parent.qparams, handle, protocol=pickle.HIGHEST_PROTOCOL)

        parent.parent().qparams = parent.qparams  # updated qparams in the PyTomGUI class structure

        for tab in (parent.parent().CD, parent.parent().TR, parent.parent().PP, parent.parent().SA):
            tab.qparams = parent.qparams

    def values(self):
        return [self.queue, self.nodes, self.cores, self.time, self.modules, self.command]


class DisplayText(QMainWindow, CommonFunctions):
    def __init__(self, parent, type='read'):
        super(DisplayText, self).__init__(parent)
        self.setGeometry(100,100,700,400)
        self.widget = QPlainTextEdit()
        self.widget.setStyleSheet("QPlainTextEdit{background:white;}")
        #self.widget.setEnabled(False)

        layout = QVBoxLayout()
        layout.addWidget(self.widget)
        self.setLayout(layout)
        self.setCentralWidget(self.widget)
        self.projectname = self.parent().projectname
        self.pytompath = self.parent().pytompath
        self.add_toolbar(self.processtrigger, open=True, save=True, openText='Open File', saveText='Save File')

        if type == 'edit':

            button = QPushButton('Save')
            button.clicked.connect(self.saveText)
            layout.addWidget(button)

    def processtrigger(self, q):
        if   q.text() == 'New':          self.readText(folder=self.projectname)
        elif q.text() == 'Open File':         self.readText(folder=self.projectname)
        elif q.text() == 'Quit':         self.close()
        elif q.text() == 'Save File':         self.saveText(folder=self.projectname)

    def setText(self, text, title=''):
        self.setWindowTitle(title)
        self.widget.clear()
        self.widget.insertPlainText(text)

    def saveText(self, folder='./'):
        filename = str(QFileDialog.getSaveFileName(self, 'Save particle list.', folder)[0])
        if filename:
            txt = self.widget.toPlainText()
            outfile = open(filename, 'w')
            outfile.write(txt)
            outfile.close()

    def readText(self, folder='./'):
        filename = QFileDialog.getOpenFileName(self, 'Open file', folder, "Marker files (*.*)")[0]
        if filename:
            outfile = open(filename, 'r')
            text = outfile.read()
            outfile.close()
            self.setText(text, title=os.path.basename(filename))
            self.show()


class SelectModules(QWidget):
    def __init__(self, parent=None, modules=[], mode=''):
        super(SelectModules, self).__init__(parent)

        myBoxLayout = QVBoxLayout()
        self.setLayout(myBoxLayout)
        #self.setCentralWidget(myQWidget)
        self.generalize = True
        self.toolbutton = QToolButton(self)
        self.toolbutton.setText('Select Modules')
        self.toolmenu = QMenu(self)
        self.toolmenu.setStyleSheet("selection-background-color: #1989ac;")
        self.toolbutton.setMinimumWidth(150)
        myBoxLayout.setContentsMargins(0, 0, 0, 0)
        self.setContentsMargins(0, 0, 0, 0)
        self.actions = []
        self.modules = []
        self.mode = mode
        self.p = parent

        q = "module avail --long 2>&1 | awk 'NR >2 {print $1}'"
        avail = [line for line in os.popen(q).readlines() if not line.startswith('/')
                 and not line.startswith('shared') and not 'intel' in line]
        # only add pytom, imod, and motioncor, rest is contained now in the conda environment
        self.grouped = [mod.strip("\n") for mod in avail if 'motioncor' in mod
                        or 'imod' in mod or 'pytom' in mod or 'python' in mod
                        or 'openmpi' in mod or 'lib64' in mod]
        self.update = True
        for i, name in enumerate(self.grouped):
            action = self.toolmenu.addAction(name)
            action.setCheckable(True)
            self.actions.append(action)
            # this links the action of the widget I think...
            action.toggled.connect(lambda d, index=i: self.updateModules(index))
            self.toolbutton.setMenu(self.toolmenu)

        self.toolbutton.setPopupMode(QToolButton.InstantPopup)
        myBoxLayout.addWidget(self.toolbutton)

        # activate modules for all qparams ?
        self.activateModules(modules)

    def updateModules(self, index):
        # this makes sure this function is not running multiple times
        if self.update == False:
            return

        name = self.actions[index].text()
        origin = name.split('/')[0]

        self.update = False  # start update step
        # this probably ensures that only one version of the module is activated
        for action in self.actions:
            tempName = action.text()
            if name == tempName:
                continue
            tempOrigin = tempName.split('/')[0]
            if origin == tempOrigin and action.isChecked() == True:
                action.setChecked(False)
        self.update = True  # finish update step

        # start settings activated modules based on qparams file
        self.modules = self.getActivatedModules()

        # get current job name
        text = self.p.widgets[self.mode + 'jobName'].currentText()
        try:  # try to update the modules, which should be the widget selected modules
            self.p.qparams[text].update_modules(parent=self.p, mode=self.mode)
        except KeyError:  # if KeyError the widget has not been intialized, pytom is initializing
            pass

        # ensure that module in all is added and removed from the other jobs
        removed = name not in self.modules
        if (text == 'All') and (self.mode + 'modules' in self.p.widgets.keys()):  # remove from all if the module
            # was activated/deactivated for all jobs
            for jobname in self.p.jobnames:
                if jobname != 'All':
                    if removed:
                        self.p.qparams[jobname].modules = [mod for mod in
                                                           self.p.qparams[jobname].modules if mod != name]
                    else:
                        self.p.qparams[jobname].modules = list(set(self.p.qparams[jobname].modules + [name]))
                        # How to discriminate init step?
                    self.p.qparams[jobname].update_settings(parent=self.p)

    def getLatestVersion(self, module_name):
        from distutils.version import StrictVersion
        # remove module name from available
        available = [m for m in self.grouped if module_name in m]
        if len(available) == 0:
            return module_name
        versions = []
        for avail in [m.replace(module_name, '') for m in self.grouped if module_name in m]:
            matches = [x.group() for x in re.finditer('\\d+(\\.\\d+)*', avail)]
            if len(matches) != 0:
                versions.append(max(matches, key=len))
        if len(versions) != 0:
            latest = sorted(versions, key=StrictVersion)[-1]
            return [a for a in available if latest in a][0]
        else:  # no versions numbers found, just return the first encountered module
            return available[0]

    def activateModules(self, modules, block=False):
        # if module name contains exactly pytom, imod, or motioncor2, find the latest versions and activate
        # this can happen uppon processing directory init
        for i, mod in enumerate(modules):
            if mod not in self.grouped and mod in ['pytom', 'imod', 'motioncor2']:
                modules[i] = self.getLatestVersion(mod)

        if block:
            for action in self.actions:
                action.blockSignals(True)

        for action in self.actions:
            if action.text() in modules:
                action.setChecked(True)
            else:
                action.setChecked(False)

        for action in self.actions:
            action.blockSignals(False)

    def getModules(self):
        return self.modules

    def getActivatedModules(self):
        return [action.text() for action in self.actions if action.isChecked()]


class PlotterSubPlots(QMainWindow,CommonFunctions):
    def __init__(self, parent=None, width=800, size_subplot=80, size_subtomo=40, height=1000, offset_x=0, offset_y=0):
        super(PlotterSubPlots,self).__init__(parent)
        self.width = width
        self.height = height
        self.size_subplot = size_subplot
        self.size_subtomo = size_subtomo
        self.dx = offset_x
        self.dy = offset_y

        self.num_subtomo_per_row = round(self.width/self.size_subplot)

        self.size_policies()
        self.superscroll = QScrollArea(self)
        self.superscroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.superscroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.superscroll.setGeometry(self.dx,self.dy,self.width,self.height)

        self.scrollarea = QScrollArea()
        self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setGeometry(self.dx,self.dy,self.width,self.height)
        self.scrollarea.setSizePolicy(self.sizePolicyA)
        self.canvas = KeyPressGraphicsWindow(size=(self.width, self.height), border=None)
        self.scrollarea.setWidget(self.canvas)
        self.superscroll.setWidget(self.scrollarea)
        #self.setCentralWidget(self.scrollarea)
        self.setLayout(QHBoxLayout())
        self.setGeometry(self.dx,self.dy,self.width,self.height)
        self.setWindowTitle('Selected Particles')
        self.init_variables()

    def init_variables(self):
        #layout.addWidget(self.canvas)
        self.vBoxList = []
        self.iItemList = []
        self.coordinates = []
        self.assigned = []
        self.index = [0,0]
        self.num_assigned = 0

    def keyPressEvent(self, evt):
        if evt.key() == Qt.Key_Escape:
            self.close()
            self.parent().close()

    def add_subplot(self, tomo, position):
        [x,y,z,s]  = position
        xmin = max(0, x - int(self.size_subtomo / 2))
        xmax = min(tomo.shape[2], x + int(self.size_subtomo / 2))
        ymin = max(0, y - int(self.size_subtomo / 2))
        ymax = min(tomo.shape[1], y + int(self.size_subtomo / 2))
        subtomo = np.zeros((int(self.size_subtomo),int(self.size_subtomo)),dtype=float)
        subtomo[:xmax-xmin,:ymax-ymin] = (tomo[position[2]-4:position[2]+4].sum(axis=0).T)[xmin:xmax, ymin:ymax]

        self.position = position
        if self.index[1] == 0 or self.index[0]+1 >= len(self.vBoxList):

            self.newViewBox()
            self.iItemList.append( pg.ImageItem(subtomo))
            self.vBoxList[-1].addItem(self.iItemList[-1])
            self.vBoxList[-1].setRange(xRange=[0, self.size_subtomo], yRange=[0, self.size_subtomo], padding=0)
            self.vBoxList[-1].setAspectLocked(True)

            self.index = [self.index[0] + 1, 0]
        else:

            self.iItemList[self.index[0]].setImage(image=subtomo)
            self.assigned[self.index[0]] = 1
            self.num_assigned += 1
            self.coordinates[self.index[0]] = self.position + [self.index[0]]
            self.index = [self.num_assigned, 0]
            self.vBoxList[self.index[0]].setRange(xRange=[0, self.size_subtomo], yRange=[0, self.size_subtomo], padding=0)
            for n, a in enumerate(self.assigned):
                if a == 0:
                    self.index = [n,1]
                    break

    def reset_display_subtomograms(self, particleList, volume):

        self.index = [0,0]
        self.num_assigned = 0
        for n, (x,y,z,s) in enumerate(particleList):
            position = [x,y,z,s]
            xmin = max(0, x - int(self.size_subtomo / 2))
            xmax = min(volume.shape[2], x + int(self.size_subtomo / 2))
            ymin = max(0, y - int(self.size_subtomo / 2))
            ymax = min(volume.shape[1], y + int(self.size_subtomo / 2))

            subtomo = np.zeros((int(self.size_subtomo),int(self.size_subtomo)),dtype=float)
            try:subtomo[:xmax-xmin,:ymax-ymin] = (volume[z-4:z+4].sum(axis=0).T)[xmin:xmax, ymin:ymax]
            except: pass
            if n < len(self.iItemList):
                self.position = [x,y,z,s]
                self.iItemList[self.index[0]].setImage(image=subtomo)
                self.assigned[self.index[0]] = 1
                self.num_assigned += 1
                self.coordinates[self.index[0]] = self.position + [self.index[0]]
                self.vBoxList[self.index[0]].setRange(xRange=[0, self.size_subtomo], yRange=[0, self.size_subtomo], padding=0)
                self.index = [self.num_assigned, 0]
            else:
                self.index = [n,0]
                self.add_subplot(volume, [x,y,z,s])

        for nn in range(len(particleList), len(self.iItemList)):
            self.delete_subplot(self.coordinates[nn], start=len(particleList))

    def delete_subplot(self, position, start=0):
        for x,y,z,s, index in self.coordinates[start:]:
            if x == position[0] and y == position[1] and z == position[2]:
                blank = np.zeros((int(self.size_subtomo), int(self.size_subtomo)), dtype=float)
                try:
                    self.iItemList[index].setImage(image=blank)
                except:
                    continue
                if index < self.index[0]:
                    self.index = [index,1]
                self.assigned[index] = 0
                self.num_assigned -= 1
                self.coordinates[index][:3] = [-1, -1, -1]
                break

    def newViewBox(self):
        row = int(self.index[0] / self.num_subtomo_per_row)
        col = int(self.index[0] % self.num_subtomo_per_row)
        if row*self.size_subplot > self.height:
            self.scrollarea.setGeometry(self.dx,self.dy,self.width,row*self.size_subplot)
        self.vBoxList.append( self.canvas.addViewBox(row=row, col=col) )
        self.vBoxList[-1].setMenuEnabled(False)
        self.vBoxList[-1].setMouseEnabled(False, False)
        if self.index[0] == 0: self.vBoxList[-1].scene().sigMouseClicked.connect(self.mouseClicked)

        self.vBoxList[-1].setGeometry(0, 0, self.size_subplot, self.size_subplot)

        #
        self.assigned.append(1)
        self.coordinates.append(self.position+[self.index[0]])
        self.num_assigned += 1

    def mouseClicked(self,event):
        ID = -2
        try:
            event.pos()

            for id, vb in enumerate(self.vBoxList):
                P = vb.mapSceneToView(event.scenePos())
                if P.x() > 0.0001 and P.y() > 0.0001 and P.x() < self.size_subtomo and P.y() < self.size_subtomo:
                    ID = id
                    break
        except:
            return

        if ID < 0: return

        if event.button() == 2:
            self.parent().remove_from_coords(self.coordinates[ID])

        elif self.coordinates[ID][2] > -1:

            #self.parent().pos.setX(self.coordinates[ID][0] - self.parent().radius)
            #self.parent().pos.setY(self.coordinates[ID][1] - self.parent().radius)
            #self.parent().centimage.addItem(circle(self.parent().pos, size=(self.parent().radius)*1.5, color=Qt.yellow))
            self.setWindowTitle("Error Score: {:6.3f}".format( self.coordinates[ID][3]))
            self.parent().slice = self.coordinates[ID][2]
            self.parent().replot()
            self.parent().update_circles()
            #self.parent().remove_all(keepPL=True)

            self.parent().pos.setX(self.coordinates[ID][0] - self.parent().radius)
            self.parent().pos.setY(self.coordinates[ID][1] - self.parent().radius)
            try:
                self.parent().red_circle.setPos(self.parent().pos)
            except:
                self.parent().red_circle = circle(self.parent().pos, size=(self.parent().radius) * 2, color=Qt.red)
                self.parent().centimage.addItem(self.parent().red_circle)

    def reset_display_subtomograms_old(self, particleList, volume ):
        import time
        t = time.time()

        for child in self.vBoxList:
           self.canvas.removeItem(child)

        print(f'Removal time: {(time.time()-t)*1000:.1f"}')

        self.init_variables()
        import atexit
        from multiprocessing import Process
        from pytom.gui.guiFunctions import loadstar, savestar, kill_proc

        self.num_subtomo_per_row = int(self.width / self.size_subplot)
        #self.index = [0,self.index[1]]

        nr_procs = 0
        procs=[]
        for proc_id in range(nr_procs):
            manager = Manager()
            vb = manager.list([[],]*nr_procs)
            ass = manager.list([[],]*nr_procs)
            coords =  manager.list([[],]*nr_procs)
            ilist =  manager.list([[],]*nr_procs)
            num_ass = manager.list([0,]*nr_procs)
            proc = Process(target=self.update_do,
                           args=(volume, particleList[proc_id::nr_procs]))
            procs.append(proc)
            proc.start()
            atexit.register(kill_proc, proc)

        while procs:
            procs = [proc for proc in procs if proc.is_alive()]

        if nr_procs == 0:
            for n, (x,y,z,s) in enumerate(particleList):
                self.add_subplot( volume, [x,y,z,s] )

            #self.update_subplot(tomo, particleList)
            #self.update_do(volume, particleList,procid, self.size_subtomo, vb, self.canvas, ass, coords, ilist, self.num_subtomo_per_row, num_ass )

        self.show()

    def update_do(self, volume, particleList, procid, size_subtomo, vb, canvas, ass, coords, ilist, num_subtomo_per_row, num_ass ):
        for n, (x,y,z,s) in enumerate(particleList):
            self.add_subplot( volume, [x,y,z,s], procid, size_subtomo, vb, canvas, ass, coords, ilist, num_subtomo_per_row, num_ass )

    def add_subplot_parallel(self, tomo, position, procid, size_subtomo, vBoxList, canvas, assigned, coordinates, iItemList, num_subtomo_per_row, num_assigned):
        [x,y,z,s]  = position
        xmin = max(0, x - int(size_subtomo / 2))
        xmax = min(tomo.shape[1], x + int(size_subtomo / 2))
        ymin = max(0, y - int(size_subtomo / 2))
        ymax = min(tomo.shape[2], y + int(size_subtomo / 2))

        subtomo = np.zeros((int(size_subtomo),int(size_subtomo)),dtype=float)
        subtomo[:xmax-xmin,:ymax-ymin] = (tomo[z-size_subtomo//2:z+size_subtomo//2].sum(axis=0).T)[xmin:xmax, ymin:ymax]


        if self.index[1] == 0 or self.index[0]+1 >= len(vBoxList):

            row = int(self.index[0] / num_subtomo_per_row)
            col = int(self.index[0] % num_subtomo_per_row)
            # if row*self.size_subplot > self.height:
            #     self.scrollarea.setGeometry(self.dx,self.dy,self.width,row*self.size_subplot)
            vBoxList[procid].append( canvas.addViewBox(row=row, col=col) )
            vBoxList[procid][-1].setMenuEnabled(False)
            vBoxList[procid][-1].setMouseEnabled(False, False)
            #    self.vBoxList[-1].scene().sigMouseClicked.connect(self.mouseClicked)

            vBoxList[procid][-1].setGeometry(0, 0, self.size_subplot, self.size_subplot)

            #
            assigned[procid].append(1)
            coordinates[procid].append(position+[self.index[0]])
            num_assigned += 1

            iItemList.append( pg.ImageItem(subtomo))
            vBoxList[procid][-1].addItem(iItemList[-1])
            vBoxList[procid][-1].setRange(xRange=[0, size_subtomo], yRange=[0, size_subtomo], padding=0)
            vBoxList[procid][-1].setAspectLocked(True)

            self.index = [self.index[0] + 1, 0]
        else:

            iItemList[procid][self.index[0]].setImage(image=subtomo)
            assigned[procid][self.index[0]] = 1
            num_assigned[procid] += 1
            coordinates[self.index[0]] = position + [self.index[0]]
            self.index = [num_assigned[procid], 0]
            vBoxList[self.index[0]].setRange(xRange=[0, size_subtomo], yRange=[0, size_subtomo], padding=0)
            for n, a in enumerate(assigned[procid]):
                if a == 0:
                    self.index = [n,1]
                    break


# Select folder with tomograms to load them all at once for batch template matching
class SelectTomogramDir(QMainWindow, CommonFunctions):
    '''This class lets you browse to a directory to load all .mrc or .em tomos in it.'''
    def __init__(self, parent):
        super(SelectTomogramDir, self).__init__(parent)
        self.setGeometry(50, 50, 300, 100)
        self.cwidget = QWidget()
        self.gridLayout = QGridLayout()
        self.setWindowModality(Qt.ApplicationModal)

        #self.gridLayout.setContentrsMargins(10, 10, 10, 10)
        self.setStyleSheet('background: #{};'.format(self.parent().middlec) )
        self.cwidget.setLayout(self.gridLayout)
        self.setCentralWidget(self.cwidget)
        self.fill()

    def fill(self):
        columns, rows = 5, 5

        self.items, self.widgets = [['', ] * columns, ] * rows, {}
        parent = self.gridLayout

        self.row, self.column = 0, 1
        self.insert_label(parent, text='Select tomogram directory', rstep=1, alignment=Qt.AlignHCenter,
                          tooltip='Provide the foldername where a bunch of tomograms are located.')
        self.insert_lineedit(parent, 'tomogramdir', cstep=1)
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Browse',
                               action=self.browse, params=['folder', self.items[self.row][self.column - 1], ''])
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Select',
                               action=self.return_value)

    def return_value(self):
        path = self.widgets['tomogramdir'].text()

        files = [os.path.join(path, fname) for fname in os.listdir(path) if fname.endswith('.em') or fname.endswith(
                '.mrc')]

        self.parent().tomogramlist = files

        if len(files):
            QMessageBox().warning(self, "No tomograms in folder",
                                  "Folder needs to contain tomograms in .mrc or .em format.", QMessageBox.Ok)

        self.close()


# Windows Connected to Icons in header bar.
class NewProject(QMainWindow, CommonFunctions):
    '''This class creates a new windows for browsing'''
    def __init__(self,parent,label):
        super(NewProject, self).__init__(parent)
        self.setGeometry(50, 50, 300, 100)
        self.cwidget = QWidget()
        self.gridLayout = QGridLayout()
        self.setWindowModality(Qt.ApplicationModal)

        #self.gridLayout.setContentrsMargins(10, 10, 10, 10)
        self.label = label
        self.setStyleSheet('background: #{};'.format(self.parent().middlec) )
        self.cwidget.setLayout(self.gridLayout)
        self.setCentralWidget(self.cwidget)
        self.fill()


    def fill(self):
        columns, rows = 5, 5

        self.items, self.widgets = [['', ] * columns, ] * rows, {}
        parent = self.gridLayout

        self.row, self.column = 0, 1
        self.insert_label(parent, text='Project name', rstep=1, alignment=Qt.AlignHCenter,
                          tooltip='Provide the foldername of a new project. The foldername must not exist.')
        self.insert_lineedit(parent, 'projectname', cstep=1)
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Browse',
                               action=self.browse, params=['folder', self.items[self.row][self.column - 1], ''])
        self.insert_pushbutton(parent, cstep=self.column * -1+1, rstep=1, text='Create',
                               action=self.return_value)

    def return_value(self,params):
        path = self.widgets['projectname'].text()

        if os.path.exists(os.path.join(path,'logfile.js')):
            QMessageBox().warning(self, "Folder exists",
                                   "Please provide a non existing folder name.", QMessageBox.Ok)
        else:

            self.label.setText(path)
            self.close()


class GeneralSettings(QMainWindow, GuiTabWidget, CommonFunctions):
    resized = pyqtSignal()

    def __init__(self, parent):
        super(GeneralSettings, self).__init__(parent)  # this inits the inherited class, so QMainWindow?
        self.stage='generalSettings_'
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logbook = self.parent().logbook
        self.setGeometry(0, 0, 900, 500)
        self.qcommanddict = {'slurm': 'sbatch', 'sge': 'qsub', 'torque': 'qsub', 'none': 'none'}

        try:
            queue_system = self.logbook[self.stage + 'queue_system']
            num_cores = int(self.logbook[self.stage + 'num_cores'])
            time = int(self.logbook[self.stage + 'time'])
            queue = self.logbook[self.stage + 'default_queue']
        except:
            queue_system, done4 = QInputDialog.getItem(self, 'Name Queuing system',
                                                                 'Queuing system you are using:',
                                                                 self.qcommanddict.keys())
            if queue_system != 'none':
                queue, done1 = QInputDialog.getText(self, 'Name default queue', 'Enter the name of default queue:', text='defq')
                time, done2 = QInputDialog.getInt(self, 'Default Time-out time', 'After how many hours does a queued job time out:', 24)
                num_cores, done3 = QInputDialog.getInt(self, 'Default number of cores per node', 'How many cores does one node have:',max(1,multiprocessing.cpu_count()-4))
            else:  # if nothing provided set some default values
                queue_system, queue, time, num_cores = 'slurm', 'defq', '24', '20'
            self.logbook[self.stage + 'queue_system'] = queue_system
            self.logbook[self.stage + 'num_cores'] = num_cores
            self.logbook[self.stage + 'time'] = time
            self.logbook[self.stage + 'default_queue'] = queue
            self.parent().save_logfile()

        self.queue_system = queue_system
        self.qname = queue
        self.expTime = time
        self.num_cores = num_cores

        headers = ['Queuing Parameters', 'Data Transfer', 'Tomographic Reconstruction', 'Particle Picking', 'Subtomogram Analysis']
        subheaders = [[], ] * len(headers)

        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, sizeX=800, sizeY=500)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}
        self.widgets = {}
        self.subprocesses = 10
        self.modules = {}

        self.tabs = {'tab1': self.tab1,
                     'tab2': self.tab2,
                     'tab3': self.tab3,
                     'tab4': self.tab4,
                     'tab5': self.tab5,
                     }

        self.tab_actions = {'tab1': self.tab1UI,
                            'tab2': self.tab2UI,
                            'tab3': self.tab3UI,
                            'tab4': self.tab4UI,
                            'tab5': self.tab5UI,

                            }

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if tt in ('tab1', 'tab2', 'tab3', 'tab4', 'tab5'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()

                if tt in ('tab1', 'tab2', 'tab3', 'tab4', 'tab5'):
                    self.tab_actions[tt]()

                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])

        self.resized.connect(self.sizetest)
        self.sizetest()

    def resizeEvent(self, event):
        self.resized.emit()
        return super(GeneralSettings, self).resizeEvent(event)

    def sizetest(self):
        w = self.frameGeometry().width()
        h  = self.frameGeometry().height()

        for scrollarea in self.scrollareas:
            scrollarea.resize(w,h)

    def setQNames(self):
        self.qnames = ['defq', 'fastq']

    def updateJobName(self, mode='v00_QParams_'):
        self.currentJobName = self.widgets[mode + 'jobName'].currentText()
        self.widgets[mode + 'queueName'].setText(self.qparams[self.currentJobName].queue)
        self.widgets[mode + 'command'].setText(self.qparams[self.currentJobName].command)
        self.widgets[mode + 'maxTime'].setValue(self.qparams[self.currentJobName].time)
        self.widgets[mode + 'numberOfNodes'].setValue(self.qparams[self.currentJobName].nodes)
        self.widgets[mode + 'numberOfCores'].setValue(self.qparams[self.currentJobName].cores)
        self.widgets[mode + 'modules'].activateModules(self.qparams[self.currentJobName].modules,
                                                       block=(self.currentJobName == 'All'))

    def tab1UI(self):
        self.jobnames = ['All',
                         'DataCollection', 'MotionCorrection',
                         'SingleAlignment', 'BatchAlignment',
                         'ReconstructWBP', 'ReconstructINFR', 'BatchReconstruct',
                         'CTFDetermination', 'SingleCTFCorrection', 'BatchCTFCorrection',
                         'SingleTemplateMatch','SingleExtractCandidates','BatchTemplateMatch','BatchExtractCandidates',
                         'SingleSubtomoReconstruct', 'BatchSubtomoReconstruct',
                         'SingleParticlePolish', 'BatchParticlePolish',
                         'AverageParticleList',
                         'FRMAlignment','GLocalAlignment',
                         'PairwiseCrossCorrelation', 'CPCA', 'AutoFocusClassification', 'FSCValidation']
        self.currentJobName = self.jobnames[0]

        if os.path.exists(os.path.join(self.projectname, '.qparams.pickle')):
            # load the previous queue settings
            with open(os.path.join(self.projectname, '.qparams.pickle'), 'rb') as handle:
                self.qparams = pickle.load(handle)

            # allow compatability of command option with old pytom projects
            # old pickle will be stored as backup
            copy_flag = False
            for val in self.qparams.values():
                if not hasattr(val, 'command'):
                    val.command = ''
                    copy_flag = True
            if copy_flag:
                os.system(f"cp {os.path.join(self.projectname, '.qparams.pickle')} "
                          f"{os.path.join(self.projectname, '.qparams.pickle.bak')}")

            # update the settings to the parent but dont store the pickle
            [p.update_settings(self, store=False) for p in self.qparams.values()]
        else:
            self.qparams = {}

        for jobname in self.jobnames:
            if not jobname in self.qparams.keys():
                try:
                    self.qparams[jobname] = QParams(queue=self.qname, time=self.expTime,
                                                    cores=self.num_cores, modules=self.parent().modules)
                except:
                    self.qparams[jobname] = QParams(modules=self.parent().modules)
                self.qparams[jobname].update_settings(self, store=True)  # write to pickle

        id = 'tab1'
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.queue_system_options = ['slurm', 'torque', 'sge', 'none']
        parent = self.table_layouts[id]
        mode = 'v00_QParams_'

        w = 150
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_label_combobox(parent, 'Queuing System ', mode + 'qType', self.queue_system_options,
                                    tooltip='Select the queuing system for you cluster.\n', logvar=True)
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_label_combobox(parent, 'Job Submission Parameters', mode + 'jobName', self.jobnames, logvar=True,
                                   tooltip='Select the job for which you want to adjust the queuing parameters.\n')
        self.insert_label_line(parent, 'Submit to: ', mode+'queueName', value=self.qparams[self.currentJobName].queue)
        self.insert_label_spinbox(parent,mode+'numberOfNodes', 'Number of Nodes',
                                  value=self.qparams[self.currentJobName].nodes, minimum=1, stepsize=1,
                                  wtype=QSpinBox)
        self.insert_label_spinbox(parent,mode+'numberOfCores', 'Number of Cores',
                                  value=self.qparams[self.currentJobName].cores, minimum=1, maximum=10000,
                                  stepsize=1, wtype=QSpinBox)
        self.insert_label_spinbox(parent,mode+'maxTime', 'Maximum Time (hours)',
                                  value=self.qparams[self.currentJobName].time, minimum=1, stepsize=1,
                                  wtype=QSpinBox)
        # options should be changed to qparams().modules,  self.parent().modules
        self.insert_label_modules(parent, mode + 'modules', text='Select Modules',
                                  options=self.qparams[self.currentJobName].modules, rstep=1, cstep=-1, mode=mode)
        self.insert_label_line(parent, 'Custom Command', mode + 'command',
                               value=self.qparams[self.currentJobName].command)

        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_label_push(parent, 'Re-activate Pushbuttons', mode + 'activatePushButtons',rstep=1,cstep=-1,width=w,
                               tooltip='Re-activate the pushbuttons of specific or all stages.\nStage is set above.',
                               pushtext='Reactivate!', action=self.reactivatePushButtons, params=mode+'jobName')
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)

        self.insert_checkbox(parent, mode + 'CustomHeader', 'use custom header for queue', cstep=0, rstep=1,
                            alignment=Qt.AlignLeft, logvar=True, columnspan=2)
        self.insert_textfield(parent, mode + 'CustomHeaderTextField', columnspan=5, rstep=1, cstep=0, logvar=True)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

        self.widgets[mode + 'qType'].currentTextChanged.connect(lambda d, m=mode: self.updateQType(m))
        self.widgets[mode + 'CustomHeader'].stateChanged.connect(lambda d, m=mode: self.updateCustomHeader(m))
        self.widgets[mode + 'jobName'].currentTextChanged.connect(lambda d, m=mode: self.updateJobName(m))

        self.widgets[mode + 'queueName'].textChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update_queue(
                                                                 parent=self, mode=m,
                                                                 current_job=self.currentJobName))
        self.widgets[mode + 'command'].textChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update_command(
                                                                 parent=self, mode=m,
                                                                 current_job=self.currentJobName))
        self.widgets[mode + 'numberOfNodes'].valueChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update_nodes(
                                                                 parent=self, mode=m,
                                                                 current_job=self.currentJobName))
        self.widgets[mode + 'numberOfCores'].valueChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update_cores(
                                                                 parent=self, mode=m,
                                                                 current_job=self.currentJobName))
        self.widgets[mode + 'maxTime'].valueChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update_time(
                                                                 parent=self, mode=m,
                                                                 current_job=self.currentJobName))

        self.updateCustomHeader(mode)  # TODO custom header not fully working
        # self.qparams[self.currentJobName].update(mode, self)

        # Some automatic queue selection in case the selected queue is not available. Some robustness for users
        # inexperienced with queueing systems.
        existing_queues = []
        for n, value in enumerate(self.qcommanddict.values()):
            if value != 'none' and os.path.exists(os.popen('which {}'.format(value)).read()[:-1]):
                existing_queues.append(n)
        queue_index = self.queue_system_options.index(self.queue_system)
        if queue_index in existing_queues:
            self.widgets[mode + 'qType'].setCurrentIndex(queue_index)
        elif len(existing_queues) != 0:
            if len(existing_queues) == 1:
                i = 0
            else:
                i = self.queue_system_options.index('slurm') if any(['slurm' == self.queue_system_options[q]
                                                                     for q in existing_queues]) else 0
            other_queue = self.queue_system_options[existing_queues[i]]
            print(f'WARNING! selected queue system ({self.queue_system}) command ' \
                  f'{self.qcommanddict[self.queue_system]} is not found, ' \
                  f'switching to available queue {other_queue} with command ' \
                  f'{self.qcommanddict[other_queue]}.')
            self.widgets[mode + 'qType'].setCurrentIndex(existing_queues[i])
        else:
            print(f'WARNING! your selected queuing system ({self.queue_system}) command '
                  f'{self.qcommanddict[self.queue_system]} is not found on your system.')
            self.widgets[mode + 'qType'].setCurrentIndex(queue_index)

    def updateCustomHeader(self, mode):
        try:
            for tab in (self.parent().CD, self.parent().TR, self.parent().PP, self.parent().SA):
                tab.custom = self.widgets[mode + 'CustomHeader'].isChecked()
                tab.genSettingsWidgets = self.widgets
        except:
            pass

    def updateQType(self, mode):

        qtype = self.widgets[mode + 'qType'].currentText().lower()
        self.queue_system = qtype
        self.logbook[self.stage + 'queue_system'] = self.queue_system
        self.parent().save_logfile()
        qcommand = self.qcommanddict[qtype]

        self.parent().qtype = qtype
        self.parent().qcommand = qcommand
        active = (qtype != 'none')
        try:
            for tab in (self.parent().CD, self.parent().TR, self.parent().PP, self.parent().SA):
                tab.qtype = qtype
                tab.qcommand = qcommand
                for key in tab.widgets.keys():
                    if key.endswith('queue'):
                        if not active:
                            tab.widgets[key].setChecked(False)
                        tab.widgets[key].setEnabled(active)

        except:
            pass

    def showTMPlot(self, mode):
        from pytom.plotting.plottingFunctions import plotTMResults

        normal = self.widgets[mode + 'particleListNormal'].text()
        mirrored = self.widgets[mode + 'particleListMirrored'].text()

        plotTMResults([normal, mirrored], labels=['Normal', 'Mirrored'])

    def tab2UI(self):
        pass

    def showFSCPlot(self, mode):
        from pytom.plotting.plotFSC import plot_FSC
        filename = self.widgets[mode + 'FSCFilename'].text()
        pixel_size = self.widgets[mode + 'PixelSize'].value()
        box_size = self.widgets[mode + 'BoxSize'].value()
        cut_off = self.widgets[mode + 'CutOff'].value()
        show_image = True
        outFname = 'temp.png'
        if filename and outFname:
            plot_FSC(filename, pixel_size, boxsize=box_size, show_image=show_image, c=cut_off)

    def tab3UI(self):
        id = 'tab3'
        self.row, self.column = 0, 2
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_DT'
        w = 150

        self.insert_label(parent, cstep=-2, rstep=2, sizepolicy=self.sizePolicyB, width=w, columnspan=2, rowspan=2)
        self.insert_label_spinbox(parent, mode + 'fidBinning', tooltip='Binning factor used for creating fiducial model in IMOD',
                                  text='Binning Factor IMOD.',cstep=-1,
                                  minimum=1, value=4, maximum=20, width=w)
        self.insert_label_spinbox(parent, mode+'numberOfCores', text='Number of cores',
                                  tooltip='Number of cores used for creation of tomogram directories.',
                                  minimum=1, value=5, maximum=20, width=w, cstep=0, rstep=1)
        self.insert_label(parent, cstep=-1, rstep=1, sizepolicy=self.sizePolicyA, width=w, columnspan=2)

        self.widgets[mode + 'numberOfCores'].valueChanged.connect(lambda d, m=mode: self.updateNumberOfCores(mode))
        self.widgets[mode + 'fidBinning'].valueChanged.connect(lambda d, m=mode: self.updateBinningFactorIMOD(mode))

        self.updateNumberOfCores(mode)

    def updateNumberOfCores(self, mode):
        self.parent().TR.num_parallel_procs = int(self.widgets[mode + 'numberOfCores'].value())

    def updateBinningFactorIMOD(self, mode):
        self.parent().TR.binningFactorIMOD = int(self.widgets[mode + 'fidBinning'].value())

    def tab4UI(self):
        pass

    def tab5UI(self):
        pass

    def reactivatePushButtons(self, id):
        name = self.jobnames[self.widgets[id].currentIndex()]
        if name == 'All': name = ''
        frameID = ['CD', 'TR', 'PP', 'SA']

        for n , frame in enumerate([self.parent().CD, self.parent().TR, self.parent().PP, self.parent().SA]):
            for key in frame.widgets.keys():
                if key.endswith(f'{name}_ExecutePushButton'):
                    frame.widgets[key].setEnabled(True)
                    print(f'{frameID[n]}: {key} has been reactivated')


class PlotWindow(QMainWindow, GuiTabWidget, CommonFunctions):
    resized = pyqtSignal()

    def __init__(self, parent):
        super(PlotWindow, self).__init__(parent)
        self.stage = 'generalSettings_'
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.qtype = None
        self.setGeometry(0, 0, 1050, 490)

        headers = ['Alignment Errors', 'Template Matching Results', 'FSC Curve']
        subheaders = [['Alignment', 'Reconstruction'], [], []] * len(headers)
        static_tabs = [[False, False], [True], [True]]

        tabUIs = [[self.tab32UI, self.tab31UI],
                  [self.tab1UI],
                  [self.tab2UI]]
        self.tabs_dict, self.tab_actions = {}, {}

        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, tabUIs=tabUIs, tabs=self.tabs_dict,
                     tab_actions=self.tab_actions, sizeX=700, sizeY=300)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}

        self.queue_job_names = []

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)

            for j in range(len(subheaders[i]) + empty):
                tt = t + (str(j + 1) * (1 - empty))

                if static_tabs[i][j]:  # tt in ('tab2', 'tab31', 'tab41', 'tab42', 'tab51', 'tab52'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()

                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)
                self.checkbox[tt] = QCheckBox('queue')

                if not static_tabs[i][j]:  # tt in ('tab1','tab32', 'tab43', 'tab53'):
                    button = QPushButton('Refresh Tab')
                    button.setSizePolicy(self.sizePolicyC)
                    button.clicked.connect(lambda d, k=tt, a=self.tab_actions[tt]: a(k))
                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])

                else:  # if tt in ('tab2','tab31','tab41', 'tab42', 'tab51', 'tab52'):
                    self.tab_actions[tt](tt)

                tab = self.tabs_dict[tt]
                tab.setLayout(self.table_layouts[tt])

        self.resized.connect(self.sizetest)
        self.sizetest()

    def resizeEvent(self, event):
        self.resized.emit()
        return super(PlotWindow, self).resizeEvent(event)

    def sizetest(self):
        w = self.frameGeometry().width()
        h = self.frameGeometry().height()

        for scrollarea in self.scrollareas:
            scrollarea.resize(w, h)

    def tab1UI(self, id=''):

        self.row, self.column = 0, 4
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_PlotTM'
        w = 150
        last, reftilt = 10, 5
        self.insert_label(parent, cstep=-4, rstep=1, sizepolicy=self.sizePolicyB, width=w)
        self.insert_label_line_push(parent, 'particleList 1', mode + 'particleListNormal', width=w, mode='file',
                                    filetype='xml', enabled=True,
                                    tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_line_push(parent, 'particleList 2 (mirrored)', mode + 'particleListMirrored', width=w,
                                    mode='file', filetype='xml', enabled=True,
                                    tooltip='Select a particleList which you want to plot.\n', cstep=-2)
        self.insert_label_line(parent, 'Label 1', mode + 'label1', width=w,value='Normal',
                                    tooltip='Label for the first particle list.\n', cstep=-1)
        self.insert_label_line(parent, 'Label 2', mode + 'label2', width=w, value='Mirrored',
                               tooltip='Label for the first particle list.\n', cstep=-1)
        self.insert_label_checkbox(parent, mode + 'plotCrossings', 'Plot crossings of curves',
                               tooltip='Do you want to display the locations where two curves cross?\n', cstep=-0)

        self.insert_pushbutton(parent, 'Plot', action=self.showTMPlot, params=mode, rstep=1, cstep=0)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

    def showTMPlot(self, mode):
        from pytom.plotting.plottingFunctions import plotTMResults

        normal = self.widgets[mode + 'particleListNormal'].text()
        mirrored = self.widgets[mode + 'particleListMirrored'].text()

        plotTMResults([normal, mirrored], plot_cross=self.widgets[mode + 'plotCrossings'].isChecked(),
                      labels=[self.widgets[mode + 'label1'].text(), self.widgets[mode + 'label2'].text()])

    def tab2UI(self, id=''):
        self.row, self.column = 0, 4
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_PlotFSC'
        w = 150
        last, reftilt = 10, 5
        self.insert_label(parent, rstep=1, cstep=-4, sizepolicy=self.sizePolicyB, width=w)
        self.insert_label_line_push(parent, 'FSC File (ascii)', mode + 'FSCFilename', mode='file', width=w,
                                    initdir=self.projectname,
                                    filetype='dat', tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_line_push(parent, 'Folder with FSC Files', mode + 'FSCFolder', mode='folder', width=w,
                                    initdir=self.projectname,
                                    filetype='dat', tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_spinbox(parent, mode + 'BoxSize', text='Dimension of Subtomogram',
                                  tooltip='Box size of 3D object',
                                  value=64, minimum=1, maximum=4000, stepsize=1, width=w)
        self.insert_label_spinbox(parent, mode + 'PixelSize', text='Pixel Size',
                                  tooltip='Pixel size of a voxel in teh object.',
                                  value=2.62, minimum=1, stepsize=1, wtype=QDoubleSpinBox, decimals=2, width=w)
        self.insert_label_spinbox(parent, mode + 'CutOff', text='Resolution Cutoff', value=0, minimum=0, stepsize=0.1,
                                  wtype=QDoubleSpinBox, decimals=3, width=w, cstep=0,
                                  tooltip='Cut-off used to determine the resolution of your object from the FSC curve. \nTypical values are 0.5 or 0.143')


        self.widgets[mode + 'FSCFilename'].textChanged.connect(lambda d, m=mode: self.updateBoxsizeAndPixelsize(m))
        self.widgets[mode + 'FSCFolder'].textChanged.connect(lambda d, m=mode: self.updateBoxsizeAndPixelsize(m))

        self.insert_pushbutton(parent, 'Plot!', action=self.showFSCPlot, params=mode, rstep=1, cstep=0)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

    def showFSCPlot(self, mode):
        from pytom.plotting.plotFSC import plot_FSC
        filename = self.widgets[mode + 'FSCFilename'].text()
        pixel_size = self.widgets[mode + 'PixelSize'].value()
        box_size = self.widgets[mode + 'BoxSize'].value()
        cut_off = self.widgets[mode + 'CutOff'].value()
        directory = self.widgets[mode + 'FSCFolder'].text()
        show_image = True
        outFname = os.path.join(self.projectname, 'Images/temp.png')
        if (filename or directory) and outFname:
            plot_FSC(filename, pixel_size, boxsize=box_size, show_image=show_image, c=cut_off, directory=directory)

    def updateBoxsizeAndPixelsize(self, mode):
        from pytom.alignment.alignmentStructures import GLocalSamplingJob
        from pytom.agnostic.io import read_size
        fscfile = self.widgets[mode+'FSCFilename'].text()
        folder =  self.widgets[mode+'FSCFolder'].text()
        if not fscfile and not folder: return

        try:
            if not folder:
                folder = os.path.dirname(fscfile)
            extensions = ['em', 'mrc']

            files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('em') or f.endswith('mrc')]
            if files:
                x,y,z = read_size(files[0])
                self.widgets[mode + 'BoxSize'].setValue(float(x))
        except Exception as e:
            print(e)

        try:
            if not folder:
                folder = os.path.dirname(fscfile)
            files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('xml') ]
            job = GLocalSamplingJob()
            for f in files:
                try:
                    job.fromXMLFile(f)
                    break
                except Exception as e:
                    pass
            self.widgets[mode + 'PixelSize'].setValue(float(job.samplingParameters.sampleInformation.getPixelSize()))
        except Exception as e:
            print(e)


    def tab31UI(self, id=''):
        try:
            import glob
            headers = ["Name Tomogram", 'Score', 'First Angle', "Last Angle", 'Ref. Image', 'Ref. Marker',
                       'Enp. Rot. Angle', 'Det. Rot. Angle', '']
            types = ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt']
            sizes = [0, 80, 0, 0, 0, 0, 0, 0, 0]

            tooltip = ['Names of existing tomogram folders.',
                       'Alignment Score.',
                       'First angle of tiltimages.',
                       'Last angle of tiltimages.',
                       'Reference image number.',
                       'Reference Marker',
                       'Expected Rotation Angle',
                       'Retrieved Rotation Angle']

            logfiles = sorted(glob.glob('{}/LogFiles/*construction*.out'.format(self.projectname)))

            values = []
            tomograms = {}
            for logfile in logfiles[::-1]:
                # tom = os.popen(f'cat {logfile} | grep "Name of Reconstruction Volume:" | awk "{print $5} " ').read()[:-1]
                dir = os.path.dirname(logfile)
                ids = os.path.basename(logfile).split('-')[0]
                infile = glob.glob(f"{dir}/{ids}_*.sh")
                if not infile:
                    continue
                logdata = open(logfile, 'r').read()
                indata = open(infile[0], 'r').read()

                tomogram = os.path.basename(indata.split('cd ')[1].split('\n')[0])
                if tomogram in tomograms: continue
                tomograms[tomogram] = 1
                alignmentscore = str(np.around(float(logdata.split('Score after optimization: ')[1].split('\n')[0]), 3))
                firstangle = logdata.split('tiltAngle=')[1].split(')')[0]
                lastangle = logdata.split('tiltAngle=')[-1].split(')')[0]
                refindex = indata.split('--referenceIndex ')[1].split(' ')[0]
                refmarker = indata.split('--referenceMarkerIndex ')[1].split(' ')[0]
                expected = indata.split('--expectedRotationAngle ')[1].split(' ')[0]

                angles = [float(i.split(',')[0]) for i in logdata.split('rot=')[1:]]
                det_angle = str(int(round(sum(angles) / len(angles))) % 360)
                values = [[tomogram, alignmentscore, firstangle, lastangle, refindex, refmarker, expected, det_angle,
                           '']] + values
            if not values:
                return

            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, runtitle='Save', addQCheckBox=False)
            self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.save2file(pid, v))

            for i in range(len(values)):
                if float(values[i][1]) < 3:
                    color = 'green'
                elif float(values[i][1]) < 4.5:
                    color = 'orange'
                else:
                    color = 'red'
                self.tables[id].widgets['widget_{}_{}'.format(i, 1)].setStyleSheet("QLabel { color : " + color + "}")
        except:
            pass

    def save2file(self, id, values):
        outname = str(QFileDialog.getSaveFileName(self, 'Save alignment scores.', self.projectname, filter='*.txt')[0])
        if outname and not outname.endswith('.txt'):
            outname += '.txt'

        if not outname:
            return

        outfile = open(outname, 'w')
        try:
            for i in range(len(values)):
                for j in range(1, 4):
                    values[i][j] = float(values[i][j])
                outfile.write('{} {:10.3f} {:10.1f} {:10.1f}    {:4s} {:4s} {:3s}   {:3s}\n'.format(*(values[i][:-1])))
        except:
            outfile.write('# Name Tomogram\n# Alignment Method\n# Origin Folder\n# Alignment Score\n# First Angle\n'
                          '# Last Angle\n# Referce Image ID\n# Reference Marker ID\n# Expected Rotation Angle\n# Determined Rotation Angle\n')
            for tomofolder in self.alignmentResulsDict.keys():
                for refmarker in self.alignmentResulsDict[tomofolder].keys():
                    for first in  self.alignmentResulsDict[tomofolder][refmarker].keys():
                        for last in  self.alignmentResulsDict[tomofolder][refmarker][first].keys():
                            for alignType in self.alignmentResulsDict[tomofolder][refmarker][first][last].keys():
                                for origin in self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType].keys():
                                    d = self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType][origin]

                #
                # results, refmarkers = self.alignmentResulsDict[key].values()
                # for i in results.keys():
                #     for j in range(len(results[i])):
                #         for k in range(1, 4):
                #             results[i][j][k] = float(results[i][j][k])
                                    outfile.write(f'{d[0]:15s} {d[4]:22s} {d[5]:11s} {d[1]:>10s} {d[2]:>8s} {d[3]:>8s} {d[6]:>7s} {d[7]:>7s} {d[8]:>7s} {d[9]:>7s}\n')

        outfile.close()

    def tab32UI(self, id=''):
        try:
            import glob, numpy
            from pytom.gui.guiFunctions import loadstar, datatype
            dname = os.path.dirname
            bname = os.path.basename
            headers = ["name tomogram", 'Score', 'First Angle', "Last Angle", 'Alignment Type', 'Origin', 'Ref. Image', 'Ref. Marker',
                       'Exp. Rot. Angle', 'Det. Rot. Angle', '']
            types = ['txt', 'txt', 'combobox', 'combobox', 'combobox', 'combobox',  'txt', 'combobox', 'txt', 'txt', 'txt']
            sizes = [0, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0]

            tooltip = ['Names of existing tomogram folders.',
                       'Alignment Score.',
                       'First angle of tiltimages.',
                       'Last angle of tiltimages.',
                       'Alignment type: Global or based on fixed Local Marker',
                       'Have you aligned sorted or sorted_ctf corrected images?',
                       'Reference image number.',
                       'Reference Marker',
                       'Expected Rotation Angle', 'Determined Rotation Angle']

            tomofolders = sorted(
                [f for f in os.listdir(f'{self.projectname}/03_Tomographic_Reconstruction/') if f.startswith('tomogram_')])

            values = []
            tomograms = {}
            self.alignmentResulsDict = {}
            for tomofolder in tomofolders:

                first = True
                try:
                    metafile = \
                    [f for f in glob.glob(f'{self.projectname}/03_Tomographic_Reconstruction/{tomofolder}/sorted/*.meta')][
                        0]
                    metadata = loadstar(metafile, dtype=datatype)
                except Exception as e:
                    print(e)
                    continue
                for logfile in sorted(glob.glob(
                        f'{self.projectname}/03_Tomographic_Reconstruction/{tomofolder}/alignment/marker*/*/*/logfile*.txt')):
                    # tom = os.popen(f'cat {logfile} | grep "Name of Reconstruction Volume:" | awk "{print $5} " ').read()[:-1]
                    logdata = open(logfile, 'r').read()
                    try:
                        ctffolder = f'{self.projectname}/03_Tomographic_Reconstruction/{tomofolder}/ctf/sorted_ctf/'
                        if os.path.exists(ctffolder):
                            ctfcorr = [f for f in os.listdir(ctffolder) if '_ctf_' in f and f.endswith('.mrc')]
                        else:
                            ctfcorr = []
                        if len(ctfcorr) > 4:
                            ctf = 16
                        else:
                            ctf = False

                        d = eval(logdata.split("Spawned job")[1].split('\n')[1])
                        first, last = os.path.basename(dname(dname(dname(logfile)))).split('_')[-1].split(',')
                        if not d: continue
                        alignmentscore = str(
                            np.around(float(logdata.split('Score after optimization: ')[1].split('\n')[0]), 3))
                        firstangle = str(np.around(metadata['TiltAngle'][d['firstProj']], 1))
                        lastangle = str(np.around(metadata['TiltAngle'][d['lastProj']], 1))
                        refindex = str(d['ireftilt'])
                        refmarker = str(d['irefmark'])
                        expected = str(int(np.around(180 * float(d['handflip'] / np.pi))))
                        logfal = logdata.split('Alignment successful. See ')[1].split(' ')[0]
                        path = os.path.join(self.projectname, '03_Tomographic_Reconstruction', tomofolder, logfal)
                        angles = guiFunctions.loadstar(path, dtype=guiFunctions.datatypeAR)['InPlaneRotation'].mean()
                        detangle = str(int(round(angles)) % 360)
                        markerPath = dname(dname(dname(logfile)))
                        alignType  = bname(dname(dname(logfile)))
                        origin     = bname(dname(logfile))
                        if not (tomofolder in self.alignmentResulsDict.keys()):
                            self.alignmentResulsDict[tomofolder] = {}
                        if not (refmarker in self.alignmentResulsDict[tomofolder].keys()):
                            self.alignmentResulsDict[tomofolder][refmarker] = {}
                        if not (first in self.alignmentResulsDict[tomofolder][refmarker].keys()):
                            self.alignmentResulsDict[tomofolder][refmarker][first] = {}
                        if not (last in self.alignmentResulsDict[tomofolder][refmarker][first].keys()):
                            self.alignmentResulsDict[tomofolder][refmarker][first][last] = {}
                        if not (alignType in self.alignmentResulsDict[tomofolder][refmarker][first][last].keys()):
                            self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType] = {}
                        if not (origin in self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType].keys()):
                            self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType][origin] = {}

                        results = [tomofolder, alignmentscore, firstangle, lastangle, alignType, origin, refindex, refmarker, expected,
                                   detangle]

                        self.alignmentResulsDict[tomofolder][refmarker][first][last][alignType][origin] = results


                    except Exception as e:
                        print('Error in alignment table: ', e)
                        continue

                try:
                    rrr = list(self.alignmentResulsDict[tomofolder].keys())
                    fff = list(self.alignmentResulsDict[tomofolder][rrr[0]].keys())
                    lll = list(self.alignmentResulsDict[tomofolder][rrr[0]][fff[0]].keys())
                    aaa = list(self.alignmentResulsDict[tomofolder][rrr[0]][fff[0]][lll[0]].keys())
                    ooo = list(self.alignmentResulsDict[tomofolder][rrr[0]][fff[0]][lll[0]][aaa[0]].keys())

                    tt, ss, ff, ll, aa, oo, rr, mm, ee, dd = self.alignmentResulsDict[tomofolder][rrr[0]][fff[0]][lll[0]][aaa[0]][ooo[0]]
                    values.append([tt, ss, fff, lll, aaa, ooo, rr, rrr, ee, dd, ''])
                except Exception as e:
                    print(f'No alignment done for {e}')

            if not values:
                return

            self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, runtitle='Save', addQCheckBox=False)
            self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.save2file(pid, v))

            for i in range(len(values)):
                if float(values[i][1]) < 3:
                    color = 'green'
                elif float(values[i][1]) < 4.5:
                    color = 'orange'
                else:
                    color = 'red'
                self.tables[id].widgets['widget_{}_{}'.format(i, 1)].setStyleSheet("QLabel { color : " + color + "}")
                tom = values[i][0]
                self.tables[id].widgets[f'widget_{i}_2'].currentIndexChanged.connect(
                    lambda d, index=i, ID=id, t=tom: self.update(index, ID, t, 2))
                self.tables[id].widgets[f'widget_{i}_3'].currentIndexChanged.connect(
                    lambda d, index=i, ID=id, t=tom: self.update(index, ID, t, 3))
                self.tables[id].widgets[f'widget_{i}_4'].currentIndexChanged.connect(
                    lambda d, index=i, ID=id, t=tom: self.update(index, ID, t, 4))
                self.tables[id].widgets[f'widget_{i}_5'].currentIndexChanged.connect(
                    lambda d, index=i, ID=id, t=tom: self.update(index, ID, t, 5))
                self.tables[id].widgets[f'widget_{i}_7'].currentIndexChanged.connect(
                    lambda d, index=i, ID=id, t=tom: self.update(index, ID, t, 1))
        except:
            pass

    def update(self, row, ID, tomofolder, priority):
        try:
            v = []
            for p in range(priority):
                if p == 0:
                    v.append(self.tables[ID].widgets[f'widget_{row}_7'].currentText())
                else:
                    v.append(self.tables[ID].widgets[f'widget_{row}_{p+1}'].currentText())
            if '' in v:
                return

            if priority == 1 and v[-1]:
                current_option = self.tables[ID].widgets[f'widget_{row}_2'].currentText()
                options = list(self.alignmentResulsDict[tomofolder][v[0]].keys())

                if current_option in options:
                    v.append(current_option)
                    #self.update(row, ID, tomofolder, priority+1)
                else:
                    v.append(options[0])

                self.tables[ID].widgets[f'widget_{row}_2'].clear()
                [self.tables[ID].widgets[f'widget_{row}_2'].addItem(option) for option in options]


            if priority == 2 and v[-1]:
                current_option = self.tables[ID].widgets[f'widget_{row}_3'].currentText()
                options = list(self.alignmentResulsDict[tomofolder][v[0]][v[1]].keys())

                if current_option in options:
                    v.append(current_option)
                    #self.update(row, ID, tomofolder, priority + 1)
                else:
                    v.append(options[0])
                self.tables[ID].widgets[f'widget_{row}_3'].clear()
                [self.tables[ID].widgets[f'widget_{row}_3'].addItem(option) for option in options]

            if priority == 3 and v[-1]:
                current_option = self.tables[ID].widgets[f'widget_{row}_4'].currentText()
                options = list(self.alignmentResulsDict[tomofolder][v[0]][v[1]][v[2]].keys())

                if current_option in options:
                    v.append(current_option)
                    # self.update(row, ID, tomofolder, priority + 1)
                else:
                    v.append(options[0])
                self.tables[ID].widgets[f'widget_{row}_4'].clear()
                [self.tables[ID].widgets[f'widget_{row}_4'].addItem(option) for option in options]

            if priority == 4 and v[-1]:
                current_option = self.tables[ID].widgets[f'widget_{row}_5'].currentText()
                options = list(self.alignmentResulsDict[tomofolder][v[0]][v[1]][v[2]][v[3]].keys())

                if current_option in options:
                    v.append(current_option)
                    # self.update(row, ID, tomofolder, priority + 1)
                else:
                    v.append(options[0])
                self.tables[ID].widgets[f'widget_{row}_5'].clear()
                [self.tables[ID].widgets[f'widget_{row}_5'].addItem(option) for option in options]

            if priority == 5 and not '' in v:
                tt, ss, ff, ll, aa, oo, rr, mm, ee, dd = self.alignmentResulsDict[tomofolder][v[0]][v[1]][v[2]][v[3]][v[4]]

                for column, text in zip([1,6,8,9],[ss,rr,ee,dd]):
                    self.tables[ID].widgets[f'widget_{row}_{column}'].setText(text)

                self.updateScoreColor(ID, row)
                self.updateTopRow(ID)

        except Exception as e:
            print('ERROR in updating Alignment Table', e)


    def updateTopRow(self, ID):


        columns = [widget for widget in list(self.tables[ID].widgets.keys()) if '0' in widget.split('_')[1]]
        num_columns = len(columns)
        row = [widget for widget in list(self.tables[ID].widgets.keys()) if '0' == widget.split('_')[-1]]
        num_rows = len(row)
        for column in range(num_columns):
            gw = self.tables[ID].general_widgets[column]
            if column not in (2,3,4,5,7):
                continue
            current_items = [gw.itemText(index) for index in range(gw.count())]

            for row in range(num_rows):
                w = self.tables[ID].widgets[f'widget_{row}_{column}']
                temp_items = [w.itemText(index) for index in range(w.count())]
                for t_item in temp_items:
                    if not t_item in current_items:
                        gw.addItem(t_item)



    def updateScoreColor(self, ID, row):
        score = float(self.tables[ID].widgets[f'widget_{row}_1'].text())
        if score < 3:
            color = 'green'
        elif score < 4.5:
            color = 'orange'
        else:
            color = 'red'

        self.tables[ID].widgets[f'widget_{row}_1'].setStyleSheet("QLabel { color : " + color + "}")

    def updateFirst(self, row, ID, tomofolder):
        refmarker = self.tables[ID].widgets[f'widget_{row}_5'].currentText()
        index = self.tables[ID].widgets[f'widget_{row}_3'].currentIndex()
        self.tables[ID].widgets[f'widget_{row}_2'].setCurrentIndex(index)
        for column in (1, 5, 7, 8):
            text = self.alignmentResulsDict[tomofolder]['results'][refmarker][index][column]
            self.tables[ID].widgets[f'widget_{row}_{column}'].setText(text)
        self.updateScoreColor(ID, row)

    def updateLast(self, row, ID, tomofolder):
        refmarker = self.tables[ID].widgets[f'widget_{row}_5'].currentText()
        index = self.tables[ID].widgets[f'widget_{row}_2'].currentIndex()
        self.tables[ID].widgets[f'widget_{row}_3'].setCurrentIndex(index)
        for column in (1, 5, 7, 8):
            text = self.alignmentResulsDict[tomofolder]['results'][refmarker][index][column]
            self.tables[ID].widgets[f'widget_{row}_{column}'].setText(text)
        self.updateScoreColor(ID, row)


class ExecutedJobs(QMainWindow, GuiTabWidget, CommonFunctions):
    resized = pyqtSignal()

    def __init__(self, parent):
        super(ExecutedJobs, self).__init__(parent)
        self.stage = 'generalSettings_'
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logbook = self.parent().logbook
        self.logfolder = os.path.join(self.projectname, 'LogFiles')

        self.qcommand = self.parent().qcommand
        self.setGeometry(0, 0, 900, 550)
        self.size_policies()
        self.progressBarCounters = {}
        self.progressBars = {}
        self.queueEvents = self.parent().qEvents
        self.qtype = None

        headers = ['Local Jobs', 'Queued Jobs']
        subheaders = [[], ] * len(headers)

        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, sizeX=900, sizeY=550)
        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}
        self.widgets = {}
        self.buttons = {}
        self.subprocesses = 10
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        self.tabs = {'tab1': self.tab1,
                     'tab2': self.tab2,
                     }

        self.tab_actions = {'tab1': self.tab1UI,
                            'tab2': self.tab2UI,
                            }

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if tt in ('tab1', 'tab2'):
                    self.table_layouts[tt] = QVBoxLayout()

                button = QPushButton('Refresh Tab')
                button.setSizePolicy(self.sizePolicyC)
                button.clicked.connect(self.tab_actions[tt])

                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)
                self.checkbox[tt] = QCheckBox('queue')

                if tt in ('tab1', 'tab2'):
                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])
                    self.buttons[tt] = button

                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])

        self.resized.connect(self.sizetest)
        self.sizetest()

    def resizeEvent(self, event):
        self.resized.emit()
        return super(ExecutedJobs, self).resizeEvent(event)

    def sizetest(self):
        w = self.frameGeometry().width()
        h = self.frameGeometry().height()

        for scrollarea in self.scrollareas:
            scrollarea.resize(w, h)

    def tab1UI(self):
        self.buttons['tab1'].setEnabled(False)
        jobfilesLocal = [line for line in sorted(os.listdir(self.logfolder + '/Local')) if line.endswith('.out')]
        self.jobFilesLocal = [os.path.join(self.logfolder, 'Local', job) for job in jobfilesLocal]
        self.populate_local()
        self.buttons['tab1'].setEnabled(True)

    def populate_local(self):

        if len(self.jobFilesLocal) == 0:
            return

        id = 'tab1'
        headers = ["Type", "QueueId", "Open Job", "Open Log", "Teminate", "Filename Jobfile Queue", 'Filename Logfile',
                   '']
        types = ['sort_txt', 'sort_txt', 'checkbox', 'checkbox', 'checkbox', 'txt', 'txt', 'txt']
        sizes = [0, 0, 0, 0, 0, 0, 0, 0]

        tooltip = []
        values = []
        added_jobs = []
        processes = [pid for pid in os.popen(f'''ps -u {getpass.getuser()} -f ''').read().split('\n') if pid]

        self.localJobStrings = {}
        self.pids = []
        for frame in self.parent().frames:
            self.localJobStrings.update(frame.localJobStrings)

        for n, logfile in enumerate(reversed(self.jobFilesLocal)):
            queueId = os.path.basename(logfile).split('-')[0]
            jobname = glob.glob(os.path.join(self.logfolder, 'Local', f'{queueId}*.sh'))
            if len(jobname) < 1:
                continue
            jobname = jobname[0]
            added_jobs.append(queueId)
            name = os.path.splitext(os.path.basename(logfile).split('-')[1])[0]

            try:
                pids = [l.split()[1] for l in processes if
                        queueId in self.localJobStrings and self.localJobStrings[queueId] in l]
                if len(pids) == 1:
                    terminate = 1
                    self.pids.append([pids[0]])
                    for i in range(4):
                        [self.pids[-1].append(l.split()[1]) for l in processes if l.split()[2] in self.pids[-1]]

                else:
                    terminate = 0
                    self.pids.append([])
            except Exception as e:
                terminate = 0
                self.pids.append([])

            values.append([name, queueId, 1, 1, terminate, jobname, logfile, ''])

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, sorting=True, connect=self.checkboxUpdate,
                      addQCheckBox=False)

        self.tab1_widgets = self.tables[id].widgets

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.do_something(pid, v))

    def checkboxUpdate(self, id, rowID=0, columnID=2):

        try:
            widgets = self.tables[id].widgets
            self.tab2_widgets
        except:
            return
        other = {2: 3, 3: 2}
        if not columnID in (2, 3):
            return
        status = widgets[f'widget_{rowID}_{columnID}'].isChecked()
        if not status: return
        for i in range(self.tables[id].table.rowCount()):
            if i != rowID:
                try:
                    widgets[f'widget_{i}_{columnID}'].setChecked(False)
                    widgets[f'widget_{i}_{other[columnID]}'].setChecked(False)
                except:
                    pass
            else:
                widgets[f'widget_{i}_{other[columnID]}'].setChecked(False)
        # self.tab2_widgets[f'widget_{i}_2'].setChecked(status)

    def tab2UI(self):
        self.buttons['tab2'].setEnabled(False)
        jobfiles = [line for line in sorted(os.listdir(self.logfolder)) if line.endswith('.out')]
        self.jobFilesQueue = [os.path.join(self.logfolder, job) for job in jobfiles if not job.startswith('local_')]
        self.populate_queue()
        self.buttons['tab2'].setEnabled(True)

    def nthElem(self, elem, n=1):
        return elem[n]

    def populate_queue(self):
        if len(self.jobFilesQueue) == 0:
            return

        id = 'tab2'
        headers = ["Type", "QueueId", "Open Job", "Open Log", "Running", 'Terminate', "Filename Jobfile Queue",
                   'Filename Logfile', '']
        types = ['sort_txt', 'sort_txt', 'checkbox', 'checkbox', 'checkbox', 'checkbox', 'txt', 'txt', 'txt']
        sizes = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        tooltip = []
        values = []

        import getpass
        whoami = getpass.getuser()
        qjobs = [int(line.split()[0]) for line in os.popen(f'squeue -u {whoami} | grep -v JOBID').readlines()]
        added_jobs = []
        for n, logfile in enumerate(self.jobFilesQueue):
            queueId = int(os.path.basename(logfile).split('-')[0])
            jobname = glob.glob(os.path.join(self.logfolder, f'{queueId}*.sh'))
            if len(jobname) < 1:
                continue
            jobname = jobname[0]
            added_jobs.append(queueId)
            running = 1 * (queueId in qjobs)
            name = os.path.splitext(os.path.basename(logfile).split('-')[1])[0]
            values.append([name, queueId, 1, 1, 16 * running, running, jobname, logfile, ''])

        for running in reversed(qjobs):
            if not running in added_jobs:
                queueId = int(running)
                if len(glob.glob(os.path.join(self.logfolder, f'{queueId}*.sh'))) < 1:
                    continue
                values.append(['', int(running), 0, 0, 16, 1, '', '', ""])

        values = sorted(values, key=self.nthElem, reverse=True)

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, sorting=True, connect=self.checkboxUpdate,
                      addQCheckBox=False)

        self.tab2_widgets = self.tables[id].widgets

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.do_something(pid, v))

    def do_something(self, pid, values):
        term = 0
        if pid == 'tab1':
            for row in range(self.tables[pid].table.rowCount()):
                logfile = values[row][6]
                exefile = values[row][5]
                if self.tab1_widgets[f'widget_{row}_3'].isChecked():
                    self.open_resultfile(logfile)
                if self.tab1_widgets[f'widget_{row}_2'].isChecked():
                    self.open_resultfile(exefile)

                if f'widget_{row}_4' in self.tab1_widgets.keys() and self.tab1_widgets[f'widget_{row}_4'].isChecked():
                    if self.pids[row]:
                        for pid in reversed(self.pids[row]):
                            term += 1
                            os.system(f'kill -9 {pid} >& /dev/null')

            if term: self.tab1UI()

        if pid == 'tab2':
            for row in range(self.tables[pid].table.rowCount()):
                logfile = values[row][7]
                jobfile = values[row][6]
                qId = self.tab2_widgets[f'widget_{row}_1'].text()

                try:
                    if self.tab2_widgets[f'widget_{row}_2'].isChecked():
                        self.open_resultfile(jobfile)
                except:
                    pass
                try:
                    if self.tab2_widgets[f'widget_{row}_3'].isChecked():
                        self.open_resultfile(logfile)
                except:
                    pass

                try:
                    if self.tab2_widgets[f'widget_{row}_5'].isChecked():
                        try:
                            os.system('scancel {}'.format(qId))
                            self.tab2_widgets[f'widget_{row}_4'].setChecked(False)
                            self.tab2_widgets[f'widget_{row}_5'].setChecked(False)
                        except:
                            print('Failed to cancel job {}'.format(qId))
                except:
                    pass

    def open_resultfile(self, logfile):
        with open(logfile, 'r') as f:
            txt = f.read()
            try:
                self.d.close()
                self.d.setText(txt, os.path.basename(logfile))
                self.d.show()
            except:
                self.d = DisplayText(self)
                self.d.setText(txt, os.path.basename(logfile))
                self.d.show()


class ConvertData(QMainWindow, GuiTabWidget, CommonFunctions):
    resized = pyqtSignal()
    def __init__(self,parent):
        super(ConvertData, self).__init__(parent)



        self.stage='convertData_'
        self.stage = 'v02_'
        self.addGeneralVariables()

        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logbook = self.parent().logbook
        self.setGeometry(0, 0, 900, 800)
        headers = ['Convert Data']
        subheaders = [[]]
        tabUIs = [[self.tab1UI]]
        static_tabs = [[True]]
        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, tabUIs=tabUIs, tabs=self.tabs_dict,
                     tab_actions=self.tab_actions, static_tabs=static_tabs, sizeY=800)

    def tab1UI(self, key, title='DataConversion'):
        tooltip = ''
        sizepol = self.sizePolicyB
        parent = self.table_layouts[key]
        mode = key
        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB, width=400 )
        self.insert_label_line_push(parent, 'Input file', wname=mode+'InputFile',width=w, initdir=self.tomogramfolder,
                                    tooltip='Select the tomogram file used for template matching.',
                                    filetype=['mrc', 'em', 'rec', 'st', 'log', 'txt', 'star', 'meta', 'xml'], mode='file')
        self.insert_label_line_push(parent,'Folder With Files', mode +'InputFolder',width=w,
                                    tooltip='Select the file where the sorted tiltimages are located.\n')
        self.insert_label_line_push(parent,'Output/Target Folder', mode +'TargetFolder',
                                    'Select the folder where theoutput file(s) are saved.\n',width=w)
        self.insert_label_line(parent, 'Prefix Query (Optional)', mode + 'PrefixQuery',
                                    'With what symbols do the file names begin.\n',width=w)
        self.insert_label_line(parent, 'Suffix Query (Optional)', mode + 'SuffixQuery',
                                    'What do the filenames end with.\n',width=w)
        self.insert_label_combobox(parent,'Output file type',mode+'OutputType',
                                   labels=['em', 'mrc', 'rec', 'st', 'txt', 'meta', 'log', 'star', 'xml'],
                                   tooltip='Select the file type of the output file',
                                   width=w,cstep=-1)
        self.insert_label_line(parent, 'Output Name (Optional)', mode + 'OutputName',
                               'What is the file name off your output file (Optional).\n',width=w)
        self.insert_label_line(parent, 'Wedge Angle(s) (Optional)', mode + 'WedgeAngles',
                               'What are the wedgeAngles. Example: 30 or 30,30.\n',width=w)
        self.insert_label_spinbox(parent, mode + 'PixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox, minimum=0.1, stepsize=0.1, value=1.0)
        self.insert_label_line_push(parent, 'Alignment result file IMOD (.xf)', wname=mode+'AlignXF',width=w, initdir=self.tomogramfolder,
                                    tooltip='Select the tomogram file used for template matching.',
                                    filetype=['xf'], mode='file')
        self.insert_label_line_push(parent, 'Folder Sorted Images', wname=mode+'SortedFolder',width=w, initdir=self.tomogramfolder,
                                    tooltip='Select the fodler with sorted images used for reconstruction.', mode='folder', cstep=-1)

        self.flags_dict = {}
        for name, flag in (('InputFile', '-f'), ('InputFolder', '-d'), ('PrefixQuery', '--prefixQuery'),
                           ('SuffixQuery', '--suffixQuery'), ('OutputName', '--outname'),
                           ('WedgeAngles', '--wedgeAngles'), ('AlignXF', '--alignxf'),
                           ('SortedFolder', '--sortedFolder')):
            self.widgets[mode + 'flag' + name] = QLineEdit()
            self.widgets[mode + name].textChanged.connect(lambda dummy, m=mode, ff=flag, nn=name: self.updateFlag(m, ff, nn))
            self.updateFlag(mode, flag, name)

        execfilename = [mode + 'TargetFolder', 'convert.sh']

        paramsSbatch = guiFunctions.createGenericDict(fname='ConvertData', folder=self.logfolder)
        paramsCmd    = [mode + 'TargetFolder', mode + 'flagInputFile', mode + 'flagInputFolder',
                        mode + 'flagPrefixQuery', mode + 'flagSuffixQuery', mode + 'OutputType',
                        mode + 'flagOutputName', mode + 'PixelSize', mode + 'flagAlignXF',
                        mode + 'flagSortedFolder', mode + 'flagWedgeAngles', templateConvertData]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=execfilename, paramsSbatch = paramsSbatch,
                                 paramsCmd=paramsCmd, mandatory_fill = [mode + 'TargetFolder'])
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[mode].addWidget(label)


    def updateFlag(self, mode, flag, name):
        a = self.widgets[mode + name].text()
        if a:
            self.widgets[mode+'flag'+name].setText(f'{flag} {a} ')
        else:
            self.widgets[mode+'flag'+name].setText('')

    def resizeEvent(self, event):
        self.resized.emit()
        return super(ConvertData, self).resizeEvent(event)

    def sizetest(self):
        w = self.frameGeometry().width()
        h  = self.frameGeometry().height()

        for scrollarea in self.scrollareas:
            scrollarea.resize(w,h)


class View2d(QMainWindow, CommonFunctions):
    def __init__(self, parent):
        super(View2d, self).__init__(parent)
        self.setGeometry(50, 50, 200, 200)
        self.size_policies()
        self.cwidget = QWidget()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.templatematchfolder = os.path.join(self.projectname, '04_Particle_Picking/Template_Matching')
        self.pickpartfolder = os.path.join(self.projectname, '04_Particle_Picking/Picked_Particles')
        self.subtomofolder = os.path.join(self.projectname, '05_Subtomogram_Analysis')
        self.tomogramfolder = os.path.join(self.projectname, '04_Particle_Picking/Tomograms')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        initdir = self.parent().projectname
        self.logbook = self.parent().logbook
        parent = QGridLayout()
        parent.setAlignment(self, Qt.AlignTop)
        mode = 'Viewer2D_'
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=0, rstep=1)
        self.insert_label_line_push(parent, 'Image', mode + 'Filename', initdir=initdir, filetype=['mrc', 'em'],
                                    tooltip='Select an image (2D and 3D are possible).', mode='file', cstep=-2, rstep=1)
        self.insert_label_line_push(parent, 'Folder with 2D/3D Images', mode + 'Foldername', initdir=initdir,
                                    tooltip='Select a folder with 2D or 3D images.', mode='folder', cstep=-2, rstep=1)
        self.insert_label(parent, cstep=0, rstep=1)
        self.insert_label_spinbox(parent, mode + 'binningFactor', 'Binning Factor', value=1, stepsize=1, minimum=1, maximum=10)
        self.insert_label_line(parent, 'Prefix', mode + 'prefix', value='sorted_')
        self.insert_label_combobox(parent, 'File Type', mode + 'filetype', ['mrc', 'em'], cstep=0)

        self.widgets[mode + 'Filename'].textChanged.connect(lambda d, m=mode, f='Filename': self.clearChoice(m, f))
        self.widgets[mode + 'Foldername'].textChanged.connect(lambda d, m=mode, f='Foldername': self.clearChoice(m, f))
        self.insert_pushbutton(parent, 'View!', tooltip='Select 3D mrc or em file.', wname=mode+'pushButtonView2D',
                               rstep=1, cstep=2, action=self.insert_image, params=[parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA, rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB, width=10)

        self.cwidget.setLayout(parent)
        self.setWindowTitle('Select either a file, or a folder with 2D images.')
        self.setCentralWidget(self.cwidget)

    def insert_image(self, will):
        self.viewer = Viewer2D(self)
        if not self.viewer.failed:
            self.viewer.show()

    def clearChoice(self, mode, choice):
        txt = self.widgets[mode + choice].text()
        if choice == 'Filename':
            self.widgets[mode + 'Foldername'].setText('')
        if choice == 'Foldername':
            self.widgets[mode + 'Filename'].setText('')
        self.widgets[mode + choice].setText(txt)


class Viewer2D(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, fname=''):
        super(Viewer2D, self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 800, 800)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.logbook = {}
        self.radius = 8
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'
        self.step_size=1
        self.redraw = True
        self.failed = False

        self.centcanvas = w = KeyPressGraphicsWindow(size=(600, 600), border=True)
        self.centimage = w.addViewBox(row=1, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()
        self.datalabel = pg.LabelItem(justify='right')
        self.centcanvas.addItem(self.datalabel, row=0, col=0)

        self.image_list = [self.centimage]

        self.layout.addWidget(w, 0, 1)
        self.layout.addWidget(self.operationbox, 1, 1)
        self.title = parent.widgets['Viewer2D_Filename'].text()

        if not self.title: self.title = '2D Viewer'
        self.setWindowTitle("2D Image Viewer: {}".format(os.path.basename(self.title)))
        self.centcanvas.wheelEvent = self.wheelEvent

        # self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centimage.scene().sigMouseMoved.connect(self.updateLabel)


        self.add_controls(self.layout_operationbox)

        QApplication.processEvents()
        self.load_image()

    def wheelEvent(self, event):

        step = event.angleDelta().y() / 120
        increment = int(self.step_size * step)
        if self.slice + increment < self.vol.shape[0] and self.slice + increment > -1:
            self.slice += increment
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked() == False)

        if Qt.Key_Right == evt.key():
            if self.slice + self.step_size < self.dim:
                self.slice += self.step_size
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > self.step_size - 1:
                self.slice -= self.step_size
                self.setWindowTitle(f"Slice: {self.slice}")
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.close()

    def updateLabel(self, pos):

          ## using signal proxy turns original arguments into a tuple
        pos = self.centimage.mapSceneToView(pos)
        x, y = pos.x(), pos.y()
        if pos.x() >= 0 and pos.y() >= 0 and pos.x() < self.vol.shape[2] and pos.y() < self.vol.shape[1]:
            self.writeLabelText( self.vol[int(self.slice)][int(y)][int(x)], x, y, self.slice)

    def writeLabelText(self, v=0, x=0, y=0, z=0):
        self.datalabel.setText(f"value = {v:7.3f} x = {x:6.0f} y={y:6.0f} z = {z:4.0f}")

    def add_controls(self, prnt):
        vDouble = QDoubleValidator()
        vInt = QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-2, value='1.', width=100)

        self.insert_checkbox(prnt, 'show_power_spectrum', cstep=1)
        self.insert_label(prnt, text='Power Spectrum', cstep=2, rstep=1, alignment=Qt.AlignLeft)

        #self.insert_checkbox(prnt, 'show_phases', cstep=1)
        #self.insert_label(prnt, text='Phases', cstep=2, alignment=Qt.AlignLeft)
        self.insert_label(prnt,sizepolicy=self.sizePolicyB)

        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)
        self.widgets['show_power_spectrum'].stateChanged.connect(self.showPowerSpectrum)
        #self.widgets['show_phases'].stateChanged.connect(self.showPhases)

        # self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if not self.redraw: return
        self.redraw = False
        self.widgets['show_power_spectrum'].setChecked(False)
        #self.widgets['apply_gaussian_filter'].setChecked(False)
        self.redraw = True

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot()

    def showPowerSpectrum(self):
        w = self.widgets['show_power_spectrum']
        if not self.redraw: return
        self.redraw = False
        #self.widgets['show_phases'].setChecked(False)
        self.widgets['apply_gaussian_filter'].setChecked(False)
        self.redraw = True

        if w.isChecked():
            try:
                self.vol = self.ps.copy()
            except:
                self.ps = self.backup.copy()
                for i in range(self.dim):
                    self.ps[i,:,:] = np.log10(np.abs(np.fft.fftshift(np.fft.fftn(self.ps[i,:,:]))))

                self.vol = self.ps.copy()
        else:
            self.vol = self.backup.copy()
            self.widgets['show_power_spectrum'].setChecked(False)
            self.widgets['apply_gaussian_filter'].setChecked(False)

        self.replot()

    def showPhases(self):
        w = self.widgets['show_phases']

        if not self.redraw: return
        self.redraw = False
        self.widgets['show_power_spectrum'].setChecked(False)
        self.widgets['apply_gaussian_filter'].setChecked(False)
        self.redraw = True

        if w.isChecked():
            try:
                self.vol = self.ph.copy()
            except:
                self.ph = self.vol.copy()
                for i in range(self.dim):
                    self.ph[i,:,:] = np.angle(np.fft.fftshift(np.fft.fftn(self.ph[i,:,:])))
                self.vol = self.ph.copy()
        else:
            self.vol = self.backup.copy()

        self.replot()

    def replot_all(self):
        self.replot()

    def replot(self):
        self.setWindowTitle(os.path.basename(self.fnames[self.slice]))
        crop = self.vol[int(self.slice), :, :]
        self.img1m.setImage(image=crop.T)

        # self.centcanvas.removeItem(self.hist)
        # self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(np.median(crop) - crop.std() * 3, np.median(crop) + crop.std() * 3)
        # self.centcanvas.addItem(self.hist)
        self.writeLabelText(z=self.slice)

    def mouseHasMoved(self, evt):
        pass

    def load_image(self):

        from pytom.agnostic.io import read, read_size
        from pytom.gui.mrcOperations import downsample
        mode = 'Viewer2D_'

        if not self.parent().widgets[mode + 'pushButtonView2D'].isEnabled():
            return

        self.parent().widgets[mode + 'pushButtonView2D'].setEnabled(False)
        QApplication.processEvents()

        folder   = self.parent().widgets[mode + 'Foldername'].text()
        file     = self.parent().widgets[mode + 'Filename'].text()
        prefix   = self.parent().widgets[mode + 'prefix'].text()
        filetype = self.parent().widgets[mode + 'filetype'].currentText()
        bin      = int(self.parent().widgets[mode + 'binningFactor'].value())

        if folder:
            files = sorted([os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(filetype) and f.startswith(prefix)])
        elif file:
            files = [file]

        self.fnames = files

        try:
            dx, dy, dz = read_size(files[0])
            self.vol = np.zeros((len(files), dx//bin, dy//bin))
            for n, fname in enumerate(files):
                data = read(fname).squeeze()
                if len(data.shape) >2:
                    data = data.sum(axis=2)

                if bin> 1:
                    data = downsample(data,bin)
                self.vol[n, :, :] = data

            self.backup = self.vol.copy()

            self.dim = self.vol.shape[0]
            self.slice = self.d = int(self.dim // 2)
            self.img1m = pg.ImageItem(self.vol[int(self.slice), :, :])

            self.centimage.addItem(self.img1m)

            self.hist = pg.HistogramLUTItem()
            self.hist.setImageItem(self.img1m)
            self.centcanvas.addItem(self.hist, row=1, col=1)

            self.replot()
        except Exception as e:
            print('ERROR: ', e)
            self.popup_messagebox('Error', 'Reading has failed', 'The reading of the file(s) has failed.')
            self.close()
            self.failed = True

        self.parent().widgets[mode + 'pushButtonView2D'].setEnabled(True)


class View3d(QMainWindow, CommonFunctions):
    def __init__(self, parent):
        super(View3d, self).__init__(parent)
        self.setGeometry(50, 50, 200, 100)
        self.size_policies()
        self.cwidget = QWidget()
        self.pytompath=self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.templatematchfolder = os.path.join(self.projectname, '04_Particle_Picking/Template_Matching')
        self.pickpartfolder = os.path.join(self.projectname, '04_Particle_Picking/Picked_Particles')
        self.subtomofolder = os.path.join(self.projectname, '05_Subtomogram_Analysis')
        self.tomogramfolder = os.path.join(self.projectname, '04_Particle_Picking/Tomograms')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand

        initdir = self.parent().projectname
        self.logbook = self.parent().logbook
        parent = QGridLayout()
        parent.setAlignment(self, Qt.AlignTop)
        mode = 'Viewer3D_'
        self.widgets= {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent,cstep=0, rstep=1)
        self.insert_label_line_push(parent, '3D object', mode + 'tomogramFname', initdir=initdir,
                                    tooltip='Select the particle list.', mode='file', filetype=['mrc', 'em'], rstep=1)
        self.insert_label_combobox(parent, '3D object', mode + 'sliceDirection', ['x', 'y', 'z'], cstep=0, rstep=1)

        self.insert_pushbutton(parent, 'View!', tooltip='Select 3D mrc or em file.',
                               rstep=1, cstep=2, action=self.insert_image, params=[parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA, rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB, width=10)

        self.widgets[mode + 'sliceDirection'].setCurrentIndex(2)

        self.cwidget.setLayout(parent)
        self.setWindowTitle('Select a file describing a 3D model')
        self.setCentralWidget(self.cwidget)


    def insert_image(self, will=''):
        self.view3d = Viewer3D(self)
        if not self.view3d.failed:
            self.view3d.show()


class Viewer3D(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, fname=''):
        super(Viewer3D, self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 800, 800)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.logbook = {}
        self.radius = 8
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'
        self.failed = False
        self.mode = 'Viewer3D_'
        self.dirId = self.parent().widgets[self.mode + 'sliceDirection'].currentIndex()
        self.slicedir = self.parent().widgets[self.mode + 'sliceDirection'].currentText()

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(200, 600), border=True)
        self.leftimage = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(600, 600), border=True)
        self.centimage = w.addViewBox(row=1, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()
        self.datalabel = pg.LabelItem(justify='right')
        self.centcanvas.addItem(self.datalabel, row=0, col=0,)

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(600, 200), border=True)
        self.bottomimage = w2.addViewBox(row=0, col=0)
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]

        self.layout.addWidget(w1, 0, 0)
        self.layout.addWidget(w, 0, 1)
        self.layout.addWidget(w2, 1, 1)
        self.layout.addWidget(self.operationbox, 1, 0)
        self.title = parent.widgets[self.mode + 'tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle("Manual Particle Selection From: {}".format(os.path.basename(self.title)))
        self.centcanvas.wheelEvent = self.wheelEvent

        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        self.leftimage.scene().sigMouseClicked.connect(self.mouseHasMovedLeft)
        self.bottomimage.scene().sigMouseClicked.connect(self.mouseHasMovedBottom)

        # self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centcanvas.sigMouseReleased.connect(self.empty)
        self.centimage.scene().sigMouseMoved.connect(self.updateLabel)

        self.load_image()
        if not self.failed:
            self.leftimage.setXRange(0, self.vol.shape[0])
            self.add_controls(self.layout_operationbox)

        QApplication.processEvents()

    def wheelEvent(self, event):

        step = event.angleDelta().y() / 120
        increment = int(int(self.widgets['step_size'].value()) * step)
        if self.slice + increment < self.vol.shape[0] and self.slice + increment > -1:
            self.slice += increment
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked() == False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.slice += update
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = -1 * int(self.widgets['step_size'].text())
                self.slice += update
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.subtomo_plots.close()
            self.close()

    def updateLabel(self, pos):

          ## using signal proxy turns original arguments into a tuple
        pos = self.centimage.mapSceneToView(pos)
        x, y = pos.x(), pos.y()
        if pos.x() >= 0 and pos.y() >= 0 and pos.x() < self.vol.shape[(self.currID+2)%3] and pos.y() < self.vol.shape[(self.currID+1)%3]:

            if self.currID == 0:
                self.writeLabelText( self.vol[int(self.slice)][int(y)][int(x)], x, y, self.slice)

            if self.currID == 1:
                self.writeLabelText( self.vol[int(y)][int(self.slice)][int(x)], self.slice, x, y)

            if self.currID == 2:
                self.writeLabelText( self.vol[int(y)][int(x)][int(self.slice)], y, self.slice, x)


    def writeLabelText(self, v=0, x=0, y=0, z=0):
        self.datalabel.setText(f"value = {v:7.3f} x = {x:6.0f} y={y:6.0f} z = {z:4.0f}")

    def add_controls(self, prnt):
        vDouble = QDoubleValidator()
        vInt = QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, text='Step Size', cstep=1, alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt, wname='step_size', cstep=-1, value=1, rstep=1, width=100,
                            minimum=1, maximum=int(self.vol.shape[0] / 4))

        self.insert_label(prnt, text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        # self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

    def replot_all(self):
        self.replot()
        volA, volB = self.getSideWindowsIndices()
        self.img1a.setImage(image=volA)
        self.img1b.setImage(image=volB)

    def sliceVol(self):
        if self.slicedir == 'x':
            crop = self.vol[:, :, int(self.slice)]
        elif self.slicedir == 'y':
            crop = self.vol[:, int(self.slice), :]
        else:
            crop = self.vol[int(self.slice), :, :]
        return crop

    def getSideWindowsIndices(self):
        if self.slicedir == 'x':
            return self.vol.sum(axis=1).T, self.vol.sum(axis=0)
        elif self.slicedir == 'y':
            return self.vol.sum(axis=0), self.vol.sum(axis=2)
        else:
            return self.vol.sum(axis=2), self.vol.sum(axis=1).T

    def replot(self):
        crop = self.sliceVol()
        self.img1m.setImage(image=crop.T)
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(np.median(crop) - crop.std() * 3, np.median(crop) + crop.std() * 3)
        self.writeLabelText(z=self.slice)

    def mouseHasMoved(self, evt):
        pass

    def mouseHasMovedBottom(self, evt):
        pos = self.bottomimage.mapSceneToView(evt.scenePos())
        if pos.y() < 0 or pos.y() >= self.vol.shape[self.currID]: return
        step = pos.y() - self.slice
        self.slice += step
        self.replot()

    def mouseHasMovedLeft(self, evt):
        pos = self.leftimage.mapSceneToView(evt.scenePos())
        if pos.x() < 0 or pos.x() >= self.vol.shape[self.currID]: return

        step = pos.x() - self.slice
        self.slice += step
        self.replot()

    def load_image(self):
        if not self.title: return

        if not os.path.exists(self.title):
            self.popup_messagebox('Error', 'File does not exist', 'File does not exist. Please provide a valid filename.')
            self.failed = True
            return

        try:
            if self.title.endswith('em'):
                from pytom.agnostic.io import read
                from pytom_numpy import vol2npy
                self.vol = read(self.title)
                self.vol = self.vol.T
                #self.vol = np.fft.fftshift(np.abs(np.fft.fftn(self.vol))**2)

            elif self.title.split('.')[-1] in ('mrc', 'mrcs', 'rec', 'st', 'map'):
                self.vol = read_mrc(self.title)

            # self.mask = np.ones_like(self.vol)
            # self.vol[self.vol < -4.] = -4.
            self.backup = self.vol.copy()


            self.vol[self.vol < self.vol.min()] = self.vol.min()

            id = 2 - self.dirId

            self.currID = id
            self.dim = self.vol.shape[id]
            self.slice = self.d = int(self.dim / 2)

            volA, volB = self.getSideWindowsIndices()
            self.img1a = pg.ImageItem(volA)
            self.img1b = pg.ImageItem(volB)
            self.img1m = pg.ImageItem(self.sliceVol().T)

            self.leftimage.addItem(self.img1a)
            self.centimage.addItem(self.img1m)
            self.bottomimage.addItem(self.img1b)

            self.leftcanvas.setAspectLocked(True)
            self.bottomcanvas.setAspectLocked(True)

            self.hist = pg.HistogramLUTItem()
            self.hist.setImageItem(self.img1m)
            self.centcanvas.addItem(self.hist, row=1,col=1)

            self.replot_all()
        except Exception as e:
            print(f'Error: {e}')

