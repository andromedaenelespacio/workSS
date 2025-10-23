import math
import orekit
from orekit.pyhelpers import setup_orekit_curdir
from org.orekit.frames import FramesFactory
from org.orekit.orbits import PositionAngleType
from org.orekit.utils import IERSConventions

vm = orekit.initVM()
setup_orekit_curdir()

class OrbitConfigGUI:
    def __init__(self, initialDate, duration, stepSize, alt_a, alt_p, inc, raan, mu, radius_E, earth):
        # Frames
        self.inertialFrame = FramesFactory.getEME2000()
        self.ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        # AbsoluteDate 
        self.initialDate = initialDate
        # Simulation duration
        self.duration = duration
        # Simulation stepsize
        self.stepSize = stepSize
        # Initial orbit parameters
        self.alt_a = alt_a           # apogee altitude in meters
        self.alt_p = alt_p           # perigee altitude in meters
        self.inc = inc               # inclination in deg
        self.raan = raan             # right ascension of ascending node
        self.mu = mu                 # grav constant
        self.radius_E = radius_E     # earth equitorial radius
        self.earth = earth           # earth body definition

        # Initial orbital parameter comps/calcs
        r_a = self.radius_E + self.alt_a            # apogee radius in meters
        r_p = self.radius_E + self.alt_p            # perigee radius in meters
        self.a = ((r_a + r_p) /2)                   # semi major axis in meters
        self.e = ((r_a - r_p)/(r_a + r_p))          # eccentricity
        self.i = math.radians(self.inc)             # inclination
        self.omega = math.radians(0)                # perigee argument
        self.av = 0                                 # true anomaly
        self.PAT = PositionAngleType.TRUE