# go one folder up
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules

# Establish TCP connection to Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number


# connect to nanonis and load modules
NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

# measure communication time for a bias set command

import time
def measure_communication_time(num_measurements=100):
    start_time = time.perf_counter()
    for i in range(num_measurements):
        NMod.Bias.Set(0.1) # set bias to 0.1 V
    end_time = time.perf_counter()
    total_time = end_time - start_time
    average_time = total_time / num_measurements
    print(f"Average communication time for BiasSet command: {average_time:.6f} seconds")    

if __name__ == "__main__":
    measure_communication_time(num_measurements=1000)
    measure_communication_time(num_measurements=1)