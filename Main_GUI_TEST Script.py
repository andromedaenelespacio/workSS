# Title: Main GUI Script
# Organisation : Super-Sharp Space Systems
# Author: Valentin Rummel
# Date 23/08/2025
# Testing a GUI Execution Script
# -------------------------------------------------------------------------------------------------------------
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtWidgets, QtGui, QtCore
from GUI.MainWindow import Ui_MainWindow
from GUI.NewWindow import Ui_NewWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavToolbar
from matplotlib.figure import Figure
import numpy as np
import orekit
from org.orekit.time import AbsoluteDate, TimeScalesFactory, TimeComponents
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions

from orekit.pyhelpers import setup_orekit_curdir
vm = orekit.initVM()
setup_orekit_curdir()

from attitude.attitude_model_test import Att_Model
from Force.Force_Model import Force_Model
from orbit_config_GUI import OrbitConfigGUI
from Orbit.orbit_propagator_ephemeris_gen import generate_ephemeris
from Config.celestialbody_config import InitialSunPos, SunPos, EarthDef
from Config.spacecraft_config import SpacecraftConfig
from Config.groundstation_config import GSConfig
from Utils.AdvPlotTest_GUI import OrbitScrollTracker, AttitudeScrollTracker, GroundPlot
from Utils.GroundStationVis import GSVis
from czml_generator import czml_gen
from Utils.LTAN_to_RAAN import ltan_to_raan


# Generate .py files from .ui files:
#
# 'pyuic5 mainwindow.ui -o MainWindow.py' must be entered into the Git Bash terminal (cd to GUI folder first) to generate the MainWindow.py file after a Qt Creator update
# 'pyuic5 newwindow.ui -o NewWindow.py' must be entered into the Git Bash terminal (cd to GUI folder first) to generate the NewWindow.py file after a Qt Creator update

# Syntax edits in the MainWindow.py script must also be made: 
#
# From bottom most to top most error-
# self.label.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop)
# self.label.setLayoutDirection(QtCore.Qt.LeftToRight)
# self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
# self.centralwidget.setLayoutDirection(QtCore.Qt.LeftToRight)
# MainWindow.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)

# Syntax edits in the NewWindow.py script must also be made: 
#
# From bottom most to top most error-
# self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)



############## 
# GUI

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        fig = Figure(dpi=100)
        self.ax = fig.add_subplot(111, projection='3d')
        super().__init__(fig)

class Mpl2Canvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        fig = Figure(dpi=100)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)

class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    paramsChanged = QtCore.pyqtSignal(float, int, QtCore.QDateTime, float, float, float, QtCore.QTime) # duration, stepsize, datetime, alt_a, alt_p, inc, ltan

    def __init__(self, *args, obj=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.setupUi(self)

        # Window title
        self.setWindowTitle("SuperSharp Systems Toolkit")
        # Window icon
        home_dir = Path.cwd()
        icon_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Icon-INV.png"
        self.setWindowIcon(QtGui.QIcon(str(icon_path)))

        # Initialise sim configuration
        utc = TimeScalesFactory.getUTC()
        initialDate = AbsoluteDate(2000, 1, 1, 0, 0, 00.000, utc)
        self.timeEdit.setTime(QtCore.QTime(10, 30, 0))
        # Frame
        self.inertialFrame = FramesFactory.getEME2000()
        self.ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        # Simulation duration
        self.duration = 6e+3  # (s)
        # Simulation stepsize
        self.stepSize = 10
        # canonical datetime
        self.qdt = self.dateTimeEdit.dateTime()
        # orbit parameters
        self.alt_a = 5e+5  # apogee altitude in meters
        self.alt_p = 5e+5  # perigee altitude in meters
        self.inc = 97.4    # inclination in deg
        self.ltan = self.timeEdit.time()   # local time of ascending node
        # Initial sun position/body
        sun = InitialSunPos(initialDate, self.inertialFrame)
        self.initial_sun_pos = sun.initial_sun_pos
        self.sun_body = sun.sun_body
        # Earth definition
        E = EarthDef(self.ITRF)
        self.mu = E.mu
        self.radius_E = E.radius_E
        self.earth = E.earth
        # Compute initial RAAN from current DateTimeEdit + LTAN timeEdit
        py = self.qdt.toPyDateTime()
        initialDate = AbsoluteDate(py.year, py.month, py.day, py.hour, py.minute,
                                py.second + py.microsecond/1e6, utc)
        self.raan = ltan_to_raan(initialDate, self.timeEdit.time(), self.inertialFrame)
        #
        # --> add configurable orbit parameters...
        #

        # Initial structured class imports
        ## Configs
        self.sc = SpacecraftConfig()
        self.gc = GSConfig()
        self.config1 = OrbitConfigGUI(initialDate, self.duration, self.stepSize, self.alt_a, self.alt_p, self.inc, self.raan, self.mu, self.radius_E, self.earth)
        ## Functional modules
        self.att = Att_Model(self.config1.inertialFrame, self.config1.earth)
        self.F =  Force_Model(self.sc.cross_section_drag, self.sc.drag_coeff, self.sc.cross_section_srp, self.sc.reflection_coeff, self.config1.earth, self.sun_body, self.initial_sun_pos)
        ## Ephemeris generation
        self.ephemeris1 = generate_ephemeris(
        self.config1.a, self.config1.e, self.config1.i, self.config1.omega,
        self.raan, self.config1.av, self.config1.PAT, self.config1.inertialFrame,
        self.config1.initialDate, self.config1.mu, self.config1.radius_E,
        self.config1.duration, self.config1.stepSize, self.config1.ITRF, self.config1.earth, self.F.SRP_force, self.F.drag_force, self.att.att_provider)

        # Validators to keep inputs resonable
        self.lineEdit.setValidator(QtGui.QDoubleValidator(0.0, 1e12, 6, self)) # duration
        self.lineEdit_2.setValidator(QtGui.QIntValidator(1, 10**9, self))      # stepsize
        self.lineEdit_5.setValidator(QtGui.QIntValidator(0, 180, self))        # inclination

        # Seed the edits with current values
        self.lineEdit.setText(str(self.duration))
        self.lineEdit_2.setText(str(self.stepSize))
        self.lineEdit_3.setText(str(self.alt_a))
        self.lineEdit_4.setText(str(self.alt_p))
        self.lineEdit_5.setText(str(self.inc))

        # Connect line edits
        self.lineEdit.editingFinished.connect(self._on_params_edited)
        self.lineEdit_2.editingFinished.connect(self._on_params_edited)
        self.lineEdit_3.editingFinished.connect(self._on_params_edited)
        self.lineEdit_4.editingFinished.connect(self._on_params_edited)
        self.lineEdit_5.editingFinished.connect(self._on_params_edited)

        # Connect dateTime settings
        self.dateTimeEdit.dateTimeChanged.connect(self._on_main_datetime_changed)
        self.timeEdit.timeChanged.connect(self._on_ltan_changed)

        # Check box toggling/pushing
        self.checkBox.clicked.connect(self.the_plot_button_was_toggled)
        self.checkBox_1.clicked.connect(self.the_plot1_button_was_toggled)
        self.checkBox_2.clicked.connect(self.the_plot2_button_was_toggled)
        self.checkBox_3.clicked.connect(self.czml_button_was_toggled)
    

        # Window Logo display
        # robust image path (relative to this script folder)
        img_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Logo-INV.png"
        pm = QtGui.QPixmap(str(img_path))
        self.label.setPixmap(pm)
        self.label.setScaledContents(True)
        # make the label stay centered in the layout
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        # prevent the label from auto-filling the whole row
        self.label.setSizePolicy(QtWidgets.QSizePolicy.Maximum,
                                QtWidgets.QSizePolicy.Maximum)
        # and wrap it in a layout that centers it
        container = QtWidgets.QHBoxLayout()
        container.addStretch(1)         # spacer on the left
        container.addWidget(self.label) # your logo
        container.addStretch(1)         # spacer on the right
        self.verticalLayout.insertLayout(0, container)


    def _qtime_to_timecomponents(self, qt: QtCore.QTime):
        sec = qt.second() + qt.msec()/1000.0
        return TimeComponents(qt.hour(), qt.minute(), sec)

    @QtCore.pyqtSlot(QtCore.QDateTime)
    def _on_main_datetime_changed(self, qdt):
        self.qdt = qdt
        self._recompute_main(qdt)
        self.paramsChanged.emit(self.duration, self.stepSize, self.qdt, self.alt_a, self.alt_p, self.inc, self.ltan)
    
    @QtCore.pyqtSlot(QtCore.QTime)
    def _on_ltan_changed(self, qtime: QtCore.QTime):
        self.ltan = qtime
        self._recompute_main(self.qdt)
        self.paramsChanged.emit(self.duration, self.stepSize, self.qdt,
                                self.alt_a, self.alt_p, self.inc,
                                self.ltan)

    @QtCore.pyqtSlot()
    def _on_params_edited(self):
        try:
            self.duration = float(self.lineEdit.text())
            self.stepSize = int(self.lineEdit_2.text())
            self.alt_a = float(self.lineEdit_3.text())
            self.alt_p = float(self.lineEdit_4.text())
            self.inc = float(self.lineEdit_5.text())
        except ValueError:
            return
        self._recompute_main(self.qdt)
        self.paramsChanged.emit(self.duration, self.stepSize, self.qdt, self.alt_a, self.alt_p, self.inc, self.ltan)

    def _recompute_main(self, qdt: QtCore.QDateTime):
        # keep self.config1 & self.ephemeris1 updated if you need them elsewhere
        self.config1.duration = self.duration
        self.config1.stepSize = self.stepSize
        self.config1.alt_a = self.alt_a
        self.config1.alt_p = self.alt_p
        utc = TimeScalesFactory.getUTC()
        py = qdt.toPyDateTime()
        initialDate = AbsoluteDate(py.year, py.month, py.day, py.hour, py.minute,
                            py.second + py.microsecond/1e6, utc)
        self.raan = ltan_to_raan(initialDate, self.ltan, self.inertialFrame)
        self.config1 = OrbitConfigGUI(initialDate, self.duration, self.stepSize, self.alt_a, self.alt_p, self.inc, self.raan, self.mu, self.radius_E, self.earth)
        self.ephemeris1 = generate_ephemeris(
            self.config1.a, self.config1.e, self.config1.i, self.config1.omega,
            self.raan, self.config1.av, self.config1.PAT, self.config1.inertialFrame,
            self.config1.initialDate, self.config1.mu, self.config1.radius_E,
            self.config1.duration, self.config1.stepSize, self.config1.ITRF, self.config1.earth, self.F.SRP_force, self.F.drag_force, self.att.att_provider)
        
        ##################################
        # Ground station visibility
        Gvis = GSVis()
        ground_station_visibility = Gvis.is_visible(self.gc.g_station_frame, self.gc.min_elev, self.gc.max_range, self.config1.inertialFrame, self.ephemeris1.dat, self.ephemeris1.pos)
        #print(ground_station_visibility)          
        if any(ground_station_visibility) is True:
            print(f"{self.gc.name} ground station is visible")
        else:
            print(f"No ground stations visible")
        ##################################


    # Sat Propogation Window
    class SatPropWindow(QtWidgets.QMainWindow, Ui_NewWindow):
        def __init__(self, outerclass, *args, obj=None, **kwargs):
            super().__init__(outerclass, *args, **kwargs)
            self.setupUi(self)
            self.outerclass = outerclass

            self.canvas = MplCanvas(self)
            self.plotWidget.layout().addWidget(self.canvas)

            # Window title
            self.setWindowTitle("Sat Propagator Plot")
            # Window icon
            home_dir = Path.cwd()
            icon_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Icon-INV.png"
            icon = QtGui.QIcon(str(icon_path))
            self.setWindowIcon(icon)

            # Slider setting link
            self.horizontalSlider.valueChanged.connect(self.on_slider_change)
            self.horizontalSlider.setSingleStep(10)
            # Outer class param change signal link
            self.outerclass.paramsChanged.connect(self.on_params_changed)
            self.on_params_changed(self.outerclass.duration, self.outerclass.stepSize, self.outerclass.qdt, self.outerclass.alt_a, self.outerclass.alt_p, self.outerclass.inc, self.outerclass.ltan)

        @QtCore.pyqtSlot(float, int, QtCore.QDateTime, float, float, float, QtCore.QTime)
        def on_params_changed(self, duration, stepSize, qdt, alt_a, alt_p, inc, ltan):

            # Use current DateTimeEdit as the initial date
            utc = TimeScalesFactory.getUTC()
            py = qdt.toPyDateTime()
            initialDate = AbsoluteDate(py.year, py.month, py.day, py.hour, py.minute,
                                    py.second + py.microsecond/1e6, utc)

            sc = SpacecraftConfig()
            config1 = OrbitConfigGUI(initialDate, duration, stepSize, alt_a, alt_p, inc, self.outerclass.raan, self.outerclass.mu, self.outerclass.radius_E, self.outerclass.earth)
            att = Att_Model(config1.inertialFrame, config1.earth)
            F =  Force_Model(sc.cross_section_drag, sc.drag_coeff, sc.cross_section_srp, sc.reflection_coeff, config1.earth, self.outerclass.sun_body, self.outerclass.initial_sun_pos)
            eph = generate_ephemeris(config1.a, config1.e, config1.i, config1.omega,
                                    self.outerclass.raan, config1.av, config1.PAT, config1.inertialFrame,
                                    config1.initialDate, config1.mu, config1.radius_E,
                                    config1.duration, config1.stepSize, config1.ITRF, config1.earth, F.SRP_force, F.drag_force, att.att_provider)
            sun = SunPos(eph.dat, config1.inertialFrame)

            # Update tracker
            self.canvas.ax.clear()
            self.tracker = OrbitScrollTracker(self.canvas.ax, eph.x_vals, eph.y_vals, eph.z_vals,
                                            eph.attitude_x, eph.attitude_y, eph.attitude_z,
                                            eph.dat, config1.radius_E, eph.pos, sun.sun_pos)
            self.horizontalSlider.setRange(0, self.tracker.max_index)
            self.on_slider_change(min(self.tracker.index, self.tracker.max_index))
            self.toolbar = NavToolbar(self.canvas, self)            # define interactive toolbar
            self.plotWidget.layout().addWidget(self.toolbar)        # add toolbar widget
            self.canvas.draw()
        
        def on_slider_change(self, val):
            # clamp & assign
            self.tracker.index = np.clip(val, 0, self.tracker.max_index)
            # update
            self.tracker.update()

    # Sat Attitude Window
    class SatAttWindow(QtWidgets.QMainWindow, Ui_NewWindow):
        def __init__(self, outerclass, *args, obj=None, **kwargs):
            super().__init__(outerclass, *args, **kwargs)
            self.setupUi(self)
            self.outerclass = outerclass

            self.canvas = MplCanvas(self)
            self.plotWidget.layout().addWidget(self.canvas)

            # Window title
            self.setWindowTitle("Sat Attitude Plot")
            # Window icon
            home_dir = Path.cwd()
            icon_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Icon-INV.png"
            icon = QtGui.QIcon(str(icon_path))
            self.setWindowIcon(icon)

            # Slider setting link
            self.horizontalSlider.valueChanged.connect(self.on_slider_change)
            self.horizontalSlider.setSingleStep(10)
            # Outer class param change signal link
            self.outerclass.paramsChanged.connect(self.on_params_changed)
            self.on_params_changed(self.outerclass.duration, self.outerclass.stepSize, self.outerclass.qdt, self.outerclass.alt_a, self.outerclass.alt_p, self.outerclass.inc, self.outerclass.ltan)      

        @QtCore.pyqtSlot(float, int, QtCore.QDateTime, float, float, float, QtCore.QTime)
        def on_params_changed(self, duration, stepSize, qdt, alt_a, alt_p, inc, ltan):

            # Use current DateTimeEdit as the initial date
            utc = TimeScalesFactory.getUTC()
            py = qdt.toPyDateTime()
            initialDate = AbsoluteDate(py.year, py.month, py.day, py.hour, py.minute,
                                    py.second + py.microsecond/1e6, utc)

            sc = SpacecraftConfig()
            config1 = OrbitConfigGUI(initialDate, duration, stepSize, alt_a, alt_p, inc, self.outerclass.raan, self.outerclass.mu, self.outerclass.radius_E, self.outerclass.earth)
            att = Att_Model(config1.inertialFrame, config1.earth)
            F =  Force_Model(sc.cross_section_drag, sc.drag_coeff, sc.cross_section_srp, sc.reflection_coeff, config1.earth, self.outerclass.sun_body, self.outerclass.initial_sun_pos)
            eph = generate_ephemeris(config1.a, config1.e, config1.i, config1.omega,
                                    self.outerclass.raan, config1.av, config1.PAT, config1.inertialFrame,
                                    config1.initialDate, config1.mu, config1.radius_E,
                                    config1.duration, config1.stepSize, config1.ITRF, config1.earth, F.SRP_force, F.drag_force, att.att_provider)
            sun = SunPos(eph.dat, config1.inertialFrame)

            # Update tracker
            self.canvas.ax.clear()
            self.atttracker = AttitudeScrollTracker(self.canvas.ax, eph.attitude_x, eph.attitude_y, eph.attitude_z, eph.dat, config1.radius_E, eph.pos, sun.sun_pos)
            self.horizontalSlider.setRange(0, self.atttracker.max_index)
            self.on_slider_change(min(self.atttracker.index, self.atttracker.max_index))
            self.toolbar = NavToolbar(self.canvas, self)            # define interactive toolbar
            self.plotWidget.layout().addWidget(self.toolbar)        # add toolbar widget
            self.canvas.draw()
        
        def on_slider_change(self, val):
            # clamp & assign
            self.atttracker.index = np.clip(val, 0, self.atttracker.max_index)
            # update
            self.atttracker.update()

    # Sat Ground Plot Window
    class SatGroundWindow(QtWidgets.QMainWindow, Ui_NewWindow):
        def __init__(self, outerclass, *args, obj=None, **kwargs):
            super().__init__(outerclass, *args, **kwargs)
            self.setupUi(self)
            self.outerclass = outerclass

            # Window title
            self.setWindowTitle("Sat Ground Plot")
            # Window icon
            home_dir = Path.cwd()
            icon_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Icon-INV.png"
            icon = QtGui.QIcon(str(icon_path))
            self.setWindowIcon(icon)

            self.canvas = Mpl2Canvas(self)
            self.plotWidget.layout().addWidget(self.canvas)
            self.groundplot = GroundPlot(self.canvas.ax, self.outerclass.config1.inertialFrame, self.outerclass.ephemeris1.dat, self.outerclass.ephemeris1.pos)
            self.toolbar = NavToolbar(self.canvas, self)            # define interactive toolbar
            self.plotWidget.layout().addWidget(self.toolbar)        # add toolbar widget         
            self.canvas.draw()

    # Sat 2D Plot Window
    '''
    class IIDPlotWindow(QtWidgets.QMainWindow, Ui_NewWindow):
        def __init__(self, outerclass, *args, obj=None, **kwargs):
            super().__init__(outerclass, *args, **kwargs)
            self.setupUi(self)
            self.outerclass = outerclass
        
            # Window title
            self.setWindowTitle("Sat Ground Plot")
            # Window icon
            home_dir = Path.cwd()
            icon_path = home_dir / "GUI" / "Supersharp-A-Satlantis-Company-Icon-INV.png"
            icon = QtGui.QIcon(str(icon_path))
            self.setWindowIcon(icon)

            self.canvas = Mpl2Canvas(self)
            self.plotWidget.layout().addWidget(self.canvas)

            # Plotting Test
            fig3, ax3 = plt.subplots()
            x = [absolutedate_to_datetime(date) for date in self.dat]
            # Orbital velcoity plot
            ax3.plot(x, self.p_vals, color='red', linewidth=1)
            # Current Position
            current_p = self.p_vals[-1]
            current_t = x[-1]
            ax3.scatter(current_t, current_p, color='red', s=10, label='Current Velocity', marker='o')
            # Axes
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
            ax3.xaxis.set_major_locator(mdates.AutoDateLocator())
            ax3.set_xlabel('Date & Time (UTC)')
            ax3.set_ylabel('Distance from Earth Center [m]')
            ax3.set_xlim(np.min(x))
            ax3.set_title('Orbit Radius')
            ax3.grid(True)

            self.canvas.draw()
    '''


    # CZM Visualisation Web Window
    class CZM:
        def __init__(self, config1, sc, att, F):
            self.config1 = config1
            self.sc = sc
            self.att = att
            self.F = F

        def run(self):
            model_path = "bluemoon.glb"
            output_path = "satellite.czml"
            czml_generate = czml_gen(self.config1, self.sc, self.att, self.F, model_path, output_path)
            czml_generate.run()


    # Window toggling functions
    def the_plot_button_was_toggled(self, checked):
        if checked:
            # New Window
            self.spwindow = self.SatPropWindow(self)
            self.spwindow.show()
        else:
            self.spwindow.hide()

    def the_plot1_button_was_toggled(self, checked):
        if checked:
            # New Window
            self.sawindow = self.SatAttWindow(self)
            self.sawindow.show()
        else:
            self.sawindow.hide()

    def the_plot2_button_was_toggled(self, checked):
        if checked:
            # New Window
            self.gpwindow = self.SatGroundWindow(self)
            self.gpwindow.show()
        else:
            self.gpwindow.hide()
            
    def czml_button_was_toggled(self, checked):
        if checked:
            # Web
            print("New CZM Visualisation Generated")
            self.czm_class = self.CZM(self.config1, self.sc, self.att, self.F)
            self.czm_class.run()
        else:
            print("CZM Visualisation Deactivated")


# App execution
app = QtWidgets.QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()