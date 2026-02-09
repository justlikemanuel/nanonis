from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
#from transfer_finder import transfer_finder
import time

TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)  # This is how you establish a TCP connection. NTCP is the connection handle.
                                                    # Check your nanonis version in Nanonis > help > info and enter the RT Engine number
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

# get current xy position
x,y = NMod.FolMe.XYPosGet(Wait_for_newest_data=True)
print("current xy position: ", x, y)

# move to xy position
wait_end_move = False
NMod.FolMe.XYPosSet(8e-9, -5e-9, wait_end_move)



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

print("Z position: ", NMod.ZCtl.ZPosGet())

NTCP.close_connection()                             # Close the connection.

