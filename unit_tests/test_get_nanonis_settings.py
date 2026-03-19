# go one folder up
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.pyNanonisMeasurements.measurementClasses.MeasurementBase import MeasurementBase

import time

## Initialize communication with Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                      # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface

NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=14000)  # This is how you establish a TCP connection. NTCP is the connection handle.
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

meas = MeasurementBase(NMod)                        # Initialize base measurement class with nanonis modules

settings = meas.nanonisSettingsGet(False)
print("Retrieved settings:")
print(settings)