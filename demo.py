# test the transfer finder function
from transfer_finder import transferFinder
import numpy as np
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.AWG_M8195A_interface.M8195A_transfer import M8195A_transfer

# Establish TCP connection to Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

# Establish connection to AWG
awg_id = "TCPIP0::localhost::inst0::INSTR"
awg01 = M8195A_transfer(awg_id)


# connect to nanonis and load modules
NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules


# set up all parameters 
"""

class transferFinder:
    def __init__(self, nanonis_module,
                i_rec_integration_time_s,
                amplitude_guess_mode,
                old_transfer_function,
                atom_tracking_parameters,
                atom_tracking_time_s, atom_tracking_interval,
                awg_reference, awg_settling_time,
                sweep_frequencies, reference_frequency,
                use_active_state,
                reference_amplitudes,
                measurement_voltage,
                filename, header,
                 ):
"""

print("Lets hope this works...")

# set up some parameters
atom_tracking_parameters = {
        "Igain": 570e-12,
        "Frequency": 10.0,
        "Amplitude": 100e-12,
        "Phase": 0.0,
        "SwitchOffDelay": 0.5
}

tf_finder = transferFinder(
    nanonis_module=NMod,
    i_rec_integration_time_s=0.1,
    amplitude_guess_mode="reference_irec",
    old_transfer_function=None,
    atom_tracking_settings=atom_tracking_parameters,
    atom_tracking_time_s=4, atom_tracking_interval=5,
    awg_reference=awg01, awg_settling_time=0.1,
    sweep_frequencies=[10, 100, 1000], reference_frequency=10,
    use_active_state=True,
    reference_amplitudes=[0.01, 0.05, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
    measurement_voltage=0.5,
    filename="transfer_function_measurement",
    header={"comment": "This is a test measurement for the transfer function optimization."}
)

print("Starting transfer function optimization...")

# prepare the system for transfer function measurement
tf_finder.prepare_measurement()

print("Preparation complete.")

# record Irec for references
tf_finder.record_irec_for_references()

print("Reference recording complete.")

# find default transfer function
tf_finder.find_reference_parameters()

print("Reference parameter finding complete.")

# measure transfer function
tf_finder.measure_transfer_function_for_all_frequencies()

# save the transfer function to a file
folder_path = "measurements"
# create the folder if it does not exist
import os
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
tf_finder.save_data(folder_path=folder_path)

# plot the reference Irec values for the reference amplitudes
iref_folder = "irec_reference_values"
if not os.path.exists(iref_folder):
    os.makedirs(iref_folder)
tf_finder.save_reference_irec_values(folder_path=iref_folder)
tf_finder.plot_reference_irec_values(folder_path=iref_folder)
print("Transfer function measurement complete.")