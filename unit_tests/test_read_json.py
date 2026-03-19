import json
# read file of format

"""

{
"type": "dummy_data",
"version": [0, 0, 1],
"header": "dummy header",
"start_time": "2026-03-19_12-54-40",
"end_time": "2026-03-19_12-55-25",
"Data": {
	"Channel names": ["frequency (Hz)", "compensation_amplitude (V)", "Current (A)", "Bias (V)", "z-controller setpoint (m)", "Lockin", "Dummy2"],
	"values": [
		[1000000.0, 0.2, 4.9322681355823406e-11, 2.0, 4.105975293100528e-08, 1.7501804828643799, 2.5591986307915704e-09],
		[2000000.0, 0.2, 4.964970448662065e-11, 2.0, 4.105975293100528e-08, 1.7501804828643799, -7.913727451125396e-10],
		[3000000.0, 0.2, 4.947997914173108e-11, 2.0, 4.105975293100528e-08, 1.7501804828643799, -1.8488382380610346e-08],
		[4000000.0, 0.2, 4.9645364208483755e-11, 2.0, 4.105975293100528e-08, 1.7501804828643799, 8.934092576851071e-09],
		[5000000.0, 0.2, 4.951729304369934e-11, 2.0, 4.105975293100528e-08, 1.7501804828643799, 6.734954816067784e-09]
		]
	}
}


"""
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

# test functions
def get_session_path():
    return NMod.Util.SessionPathGet()

#####################################
######### Read json file ############
#####################################

folder = get_session_path()
filename = "transfer_function_measurement_2026-03-19_12-55-25"
filepath = os.path.join(folder, filename + ".json")

with open(filepath, "r") as f:
	data = json.load(f)

print("Data type:", data["type"])
print("Version:", data["version"])
print("Header:", data["header"])
print("Start time:", data["start_time"])
print("End time:", data["end_time"])
print("Channel names:", data["data"]["channel names"])
print("Values:")
for value_set in data["data"]["values"]:
	print(value_set)

frequencies = []
amplitudes = []

for value_set in data["data"]["values"]:
	frequencies.append(value_set[0])
	amplitudes.append(value_set[1]) 

# plot frequencies vs amplitudes
import matplotlib.pyplot as plt
figsize = (8,6)
plt.figure(figsize=figsize)

title = "Frequencies vs Amplitudes"
fontsize = 14
plt.title(title, fontsize=fontsize+2, fontweight="bold")
plt.xlabel("Frequency (Hz)", fontsize=fontsize, fontweight="bold")
plt.ylabel("Compensation Amplitude (V)", fontsize=fontsize, fontweight="bold")
plt.plot(frequencies, amplitudes, marker="x", linestyle="-")
plt.grid()
plt.tight_layout()
plt.show()