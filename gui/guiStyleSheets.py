DEFAULT_STYLE_PROGRESSBAR = """
QProgressBar{
    border: 2px solid grey;
    border-radius: 5px;
    text-align: center
}

QProgressBar::chunk {
    background-color: lightblue;
    width: 10px;
    margin: 1px;
}
"""

DEFAULT_STYLE_TB = """
QToolBar {
    background-image: url("background10.jpg");
}
QToolBar::item:hover {
background: transparent;
}
"""

MAINC = "background: #fcfaf1;"

BARS = "background: #1989ac;"

WHITE = 'background: white;'

QLINE = 'QLineEdit:disabled {color:#54d777; background-color:white;}'
