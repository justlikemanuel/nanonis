# scripts to test the transfer_finder unit
#from libs.AWG_M8195A_interface import M8195A
from transfer_finder import transferFinder
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules

# Establish TCP connection to Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number


# connect to nanonis and load modules
NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

# test functions
def get_session_path():
    return NMod.Util.SessionPathGet()

print("Session path:", get_session_path())
