import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton
from setup_model import s, ref, prm, gps_born, likelihood
from obj import get_objective
from pyDOE2 import lhs
from mcmc2 import MCMC_adaptive
from tqdm import tqdm
from scipy.optimize import minimize
from scipy.stats.qmc import LatinHypercube
from scipy.stats import qmc


class MCMCPlot(QWidget):
    def __init__(self):
        super().__init__()

        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)

        self.setLayout(layout)

        self.x_data = []
        self.y_data = []

    def obj(self, x):
        result = get_objective(x, ref, prm, gps_born, likelihood)
        return result

    def MCMC_live_plot(self, x0, x0_length):
        xsto = np.zeros((50000, x0_length))
        outsto = np.zeros(50000)
        for i in tqdm(range(50000)):
            xsto[i] = x0
            outsto[i] = -self.obj(x0)

            # Update live plot
            self.x_data.append(i)
            self.y_data.append(outsto[i])
            self.ax.clear()
            self.ax.plot(self.x_data, self.y_data)
            self.canvas.draw()

            # Perform MCMC step
            x0 = xsto[i] + np.random.normal(0, 0.1, x0_length)
            new_out = -self.obj(x0)
            if new_out > outsto[i]:
                outsto[i] = new_out
            else:
                x0 = xsto[i]

app = QApplication(sys.argv)
window = MCMCPlot()
window.setWindowTitle('MCMC Live Plot')
window.show()
sys.exit(app.exec_())
