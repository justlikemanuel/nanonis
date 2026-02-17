# test the transfer finder function
from transfer_finder import transferFinder
import numpy as np
from nlibs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from nlibs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from nlibs.AWG_M8195A_interface.M8195A import M8195A

# Establish TCP connection to Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

# Establish connection to AWG
instrument_id = "TCPIP0::localhost::inst0::INSTR"
awg01 = M8195A(instrument_id)



NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)  # This is how you establish a TCP connection. NTCP is the connection handle.
                                                    # Check your nanonis version in Nanonis > help > info and enter the RT Engine number
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules


# set up all parameters 
transfer_finder_params = {
    "nanonis_module": NMod,
    "safe_voltage_V": 2.067, 
    "safe_current_A": 1e-9,
    "safe_x_position_m": 3e-9,
    "safe_y_position_m": 3e-9,
    "old_transfer_function": None, # if you have a previous transfer function, you can provide it here as a dict with frequencies as keys and amplitudes as values. This will be used as a starting point for the optimization.
    "atom_tracking_parameters": {
        "Igain": 570e-12,
        "Frequency": 10.0,
        "Amplitude": 100e-12,
        "Phase": 0.0,
        "SwitchOffDelay": 0.5
    },  
    "i_rec_integration_time_s": 1,
    "max_allowed_amplitude_uV": 500_000,
    "amplitude_guess_mode": "fixed",
    "atom_tracking_time_s": 5,
    "atom_tracking_interval": 10,
    "z_off_delay_s": 2,
    "awg_reference": awg01,
    "sweep_frequencies": [1000, 10000, 100000],
    "default_frequency": 10000,
    "reference_amplitudes_uV": [100_000, 150_000, 200_000, 250_000, 300_000, 350_000, 
                                400_000, 450_000, 500_000],
    "default_amplitude_uV": 200_000,
    "threshold_voltage_V": 1.5,
    "version": "1.0",
    "filename": "transfer_function.json",
    "header": "Transfer function measurement"
}

print("Lets hope this works...")
tf_finder = transferFinder(**transfer_finder_params)

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