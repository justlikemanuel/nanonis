# test the transfer finder function
from transfer_finder import transferFinder
import numpy as np
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.AWG_M8195A_interface.M8195A_transfer import M8195A_transfer

import logging
logging.basicConfig(
    #filename='experiment.log', 
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("experiment")


# Establish TCP connection to Nanonis
TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

# Establish connection to AWG
logger.info("Connecting to AWG...")
awg_id = "TCPIP0::localhost::inst0::INSTR"
awg01 = M8195A_transfer(awg_id)

logger.info("Connected to AWG. Now connecting to Nanonis...")
# connect to nanonis and load modules
NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)
NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

logger.info("Connected to Nanonis and loaded modules.")

# set up parameters 
atom_tracking_parameters = {
        "Igain": 570e-12,
        "Frequency": 10.0,
        "Amplitude": 100e-12,
        "Phase": 0.0,
        "SwitchOffDelay": 0.5
}

tf_finder = transferFinder(
    nanonis_module=NMod,
    atom_tracking_settings=atom_tracking_parameters,
    sweep_frequencies=[1e6, 2e6, 3e6, 4e6, 5e6],
    reference_frequency=1e4,
    awg_reference=awg01,
    data_channels=[
        ("Input 2 (V)", "Lockin"),
        ("Input 3 (V)", "Dummy2"), 
    ],
    measurement_voltage=0.5,
)

logger.info("Starting transfer function optimization...")

# prepare the system for transfer function measurement
tf_finder.prepare_measurement()

logger.info("Preparation complete.")

# record Irec for references
tf_finder.record_irec_for_references()

logger.info("Reference recording complete.")

# measure transfer function
tf_finder.measure_transfer_function_for_all_frequencies()

# save the transfer function to a file
folder_path = "measurements"
# create the folder if it does not exist
import os
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
tf_finder.save_data()

# plot the reference Irec values for the reference amplitudes
iref_folder = "irec_reference_values"
if not os.path.exists(iref_folder):
    os.makedirs(iref_folder)
tf_finder.save_reference_irec_values()
tf_finder.plot_reference_irec_values()
logger.info("Transfer function measurement complete.")