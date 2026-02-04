from pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from pyNanonisMeasurements.nanonisTCP import NanonisModules

import time

TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)  # This is how you establish a TCP connection. NTCP is the connection handle.
                                                    # Check your nanonis version in Nanonis > help > info and enter the RT Engine number
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

# check if controller is on
state = NMod.ZCtl.OnOffGet()
print("current state: ", state)

# turn on
NMod.ZCtl.OnOffSet(1)
time.sleep(0.1)
# check if controller is on now
state = NMod.ZCtl.OnOffGet()
print("current state: ", state)

wait_time = 5
for i in range(wait_time):
    time.sleep(1)
    print(f"waited for {i+1} seconds.")

# turn off
NMod.ZCtl.OnOffSet(0)

# check if controller is off now
state = NMod.ZCtl.OnOffGet()
print("current state: ", state)

NTCP.close_connection()                             # Close the connection.

