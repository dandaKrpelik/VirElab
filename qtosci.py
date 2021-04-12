from osci import *
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *#QApplication, QPushButton
from PyQt5.QtGui import QIcon, QPixmap

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

app = QApplication([])


def on_button_clicked():
	alert = QMessageBox()
	alert.setText('Clicked button!')
	alert.exec()

class MplCanvas(FigureCanvasQTAgg):

	def __init__(self, parent=None, width=5, height=4, dpi=100):
		fig = plt.Figure(figsize=(width, height), dpi=dpi)
		self.axes = fig.add_subplot(111)
		super(MplCanvas, self).__init__(fig)

class MainWindow(QMainWindow):
	def __init__(self, *args, **kwargs):
		super(MainWindow, self).__init__(*args, **kwargs)
		
		bk = QPixmap('lab_meta/OSCI_unit_bkg.png')
		bkg = QLabel('')
		bkg.setPixmap(bk)
		
		self.sc = MplCanvas(self, width=5, height=4, dpi=100)
		self.sc.axes.plot([0,1,2],[2,4,1])
		
		layout = QVBoxLayout()
		layout.addWidget(bkg)
		layout.addWidget(self.sc)
		
		#self. setCentralWidget(layout)
		self.foo = QWidget()
		self.foo.setLayout(layout)
		
		self. setCentralWidget(self.foo)
		
		self.show()

	def plot(self, ax):
		pass
		
	
w = MainWindow()

#button = QPushButton('Hello')
#button.clicked.connect(on_button_clicked)
#button.show()
#app.exec()
