
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets


class Dance(QMainWindow):
    def __init__(self, parent=None):
        super(Dance,self).__init__(parent)
        self.button = QPushButton('Run )')
        self.label = QLabel('Dummy')
        self.button.clicked.connect(self.update)

    def update(self):
        self.label.setText('Updated')


def test_hello(qtbot):
    widget = Dance()
    qtbot.addWidget(widget)

    # click in the Greet button and make sure it updates the appropriate label
    #qtbot.mouseClick(widget.button, QtCore.Qt.LeftButton)
    
    assert widget.label.text() == "Dummy"

def test_updated(qtbot):
    widget = Dance()
    qtbot.mouseClick(widget.button, QtCore.Qt.LeftButton)
    assert widget.label.text() == "Updated"

