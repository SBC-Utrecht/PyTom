
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
import os
import sys

from pytom.gui.pytomGUI import PyTomGui as pytomGUI

fname = 'E2ETests/FullPipeline'

if not os.path.exists(os.path.join(fname, 'logfile.pickle')):
    os.system(f"touch {os.path.join(fname, 'logfile.pickle')}")



def create_project(qtbot,name,wtype=''):
    widget = pytomGUI()
    qtbot.addWidget(widget)
    widget.open_project(name)
    #widget.label.setText(name)
    wtypes = {'CD': widget.CD, 'TR': widget.TR, 'SA': widget.SA, 'PP': widget.PP}
    window = wtypes[wtype]
    return widget, window

def set_values(window, dd, mode):
    for key, value in dd:
        if type(window.widgets[mode + key]) in (QSpinBox, QDoubleSpinBox):
            window.widgets[mode + key].setValue(value)
        elif type(window.widgets[mode + key]) == QLineEdit:
            window.widgets[mode + key].setText(str(value))

def TR_starting_values(window):
    for key in window.widgets.keys():
        tempWidget = window.widgets[key]
        if type(tempWidget) == QLineEdit:
            assert tempWidget.text().replace(' ', '') == ''
            continue

def test_01_start_gui(qtbot):
    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    assert widget.qcommand == "sbatch"

def test_03_Alignment_form(qtbot):
    mode = 'v02_SingleAlignment_'
    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR
    dd = (("LastAngle", 60), ("LastIndex", 36), ("RefMarkerIndex", 1), ("RefTiltIndex", 18), ("RotationTiltAxis", 0),
          ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", f"{fname}/03_Tomographic_Reconstruction/tomogram_000/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    text = window.widgets[mode + 'CommandText'].toPlainText()
    print(text)

    assert text != ''

def test_03_Alignment_form_pbs(qtbot):
    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR
    qtbot.mouseClick(window.pbs['tab32'], QtCore.Qt.LeftButton)

def test_04_INFR_form_mandatory_fill(qtbot):
    mode = 'v02_ReconstructINFR_'

    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR

    dd = (("LastAngle", 60), ("LastIndex", 37), ("RefMarkerIndex", 1), ("RefTiltIndex", 19), ("RotationTiltAxis", 0),
          ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", ""))#""UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)
    text = window.widgets[mode + 'CommandText'].toPlainText()
    print('INFR', text)
    assert window.widgets[mode + 'CommandText'].toPlainText() == ''

def test_04_INFR_form(qtbot):
    mode = 'v02_ReconstructINFR_'

    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR

    dd = (("LastAngle", 60), ("LastIndex", 37), ("RefMarkerIndex", 1), ("RefTiltIndex", 19), ("RotationTiltAxis", 0),
          ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", f"{fname}/03_Tomographic_Reconstruction/tomogram_000/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    text = window.widgets[mode + 'CommandText'].toPlainText()
    print('INFR', text)
    assert window.widgets[mode + 'CommandText'].toPlainText() != ''

def test_05_WBP_form(qtbot):

    mode = 'v02_ReconstructWBP_'

    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR
    dd = (("LastAngle", 60), ("LastIndex", 37), ("RefMarkerIndex", 1), ("RefTiltIndex", 18), ("RotationTiltAxis", 0),
          ("WeightingType", 1), ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", f"{fname}/03_Tomographic_Reconstruction/tomogram_000/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    assert window.widgets[mode + 'CommandText'].toPlainText() != ''

def test_05_WBP_batch(qtbot):
    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR
    qtbot.mouseClick(window.pbs['tab53'], QtCore.Qt.LeftButton)

def test_batch_recon(qtbot):
    widget = pytomGUI(warn_closing=False)
    qtbot.addWidget(widget)
    widget.open_project(fname)
    window = widget.TR
    qtbot.mouseClick(window.pbs['tab43'], QtCore.Qt.LeftButton)


