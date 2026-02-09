# module to run a measurement that allows computing the transfer function for some frequencies
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.AWG_M8195A_interface.M8195A import M8195A
import time
import numpy as np 
import json

class transfer_finder:
    def __init__(self, nanonis_module, safe_voltage_V, save_height_m, save_current_A,
                 save_x_position, save_y_position, 
                 sweep_frequencies, default_frequency,
                 reference_amplitudes_uV,
                 i_rec_integration_time, version, filename, header,
                 atom_tracking_period_s, atom_tracking_interval, z_off_delay_s, awg_reference):
        
        self.nanonis_module = nanonis_module
        self.safe_voltage = safe_voltage_V
        self.save_height_m = save_height_m
        self.save_current_A = save_current_A
        self.buffer_time_s = 1e-3
        self.awg = awg_reference
        self.i_rec_integration_time = i_rec_integration_time

        # save recovery parameters
        self.save_x_position = save_x_position
        self.save_y_position = save_y_position
        
        # Recording parameters:
        self.sweep_frequencies = sweep_frequencies
        self.default_frequency = default_frequency
        self.default_amplitude_uV = reference_amplitudes_uV[-1]
        self.reference_amplitudes_uV = reference_amplitudes_uV

        self.atom_tracking_period_s = atom_tracking_period_s
        self.atom_tracking_interval = atom_tracking_interval
        self.height_averaging_period_s = z_off_delay_s

        # logging parameters
        self.irec_vs_reference_amplitudes = []
        self.filename = filename
        self.version = version
        self.header = header

        self.default_i_rec = None # current value at the default amplitude
        self.max_amplitude_uV = None # this will be set during the recording routine, as it depends on the recorded Irec values

        self.recorded_data = None # format: [channel_name, values]

    # function to get the current Irec value
    def get_irec(self):
        """
        Function to get the current Irec value. The integration time can be set in the constructor.

        Returns
        Irec (float): The recorded Irec value in Amperes.
        """
        # set the integration time
        # get the Irec value
        irec = self.nanonis_module.Cur.Get()

        return irec

    # record Irec vs the reference amplitudes
    def record_irec_for_references(self):

        # configure AWG to output the reference signal
        # channel_number, frequency_Hz, num_periods = 1, starting_amplitude_uV = 75_000, approximate_frequency = False):
        self.awg.configure_continuous_sine_wave(channel_number=1, frequency=self.default_frequency, 
                                                num_periods=1, starting_amplitude_uV=self.reference_amplitudes_uV[0], approximate_frequency=False)
        
        for index, amplitude_uV in enumerate(self.reference_amplitudes_uV):
            self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=amplitude_uV)
            time.sleep(self.i_rec_integration_time)
            irec = self.get_irec()

            # add to list
            self.irec_vs_reference_amplitudes.append((amplitude_uV, irec))
            #print(f"Reference amplitude: {amplitude_uV} uV, recorded Irec: {irec} A")

            # for the last value: this is a a separate variable
            if index == len(self.reference_amplitudes_uV)-1:
                self.default_i_rec = irec

        # increase the amplitude until irec is 1.5 times the default irec
        # TODO: find better spacing parameter, maybe mean difference between reference amplitudes?
        amplitude_spacing_uV = self.reference_amplitudes_uV[-1] - self.reference_amplitudes_uV[-2]
        max_amplitude_uV = self.reference_amplitudes_uV[-1]

        # consider max iteration number to avoid infinite loop
        max_iterations = 10000 # random number as of now, TODO: find better parameter for this
        iteration = 0
        while self.get_irec() < 1.5 * self.default_i_rec and iteration < max_iterations:
            max_amplitude_uV += amplitude_spacing_uV
            self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=max_amplitude_uV)
            time.sleep(self.i_rec_integration_time)
            iteration += 1

        self.max_amplitude_uV = max_amplitude_uV
        print(f"Max amplitude found: {max_amplitude_uV} uV, recorded Irec: {self.get_irec()} A")


    # function to tune awg amplitude for a specific frequency to match reference irec
    def tune_awg_amplitude_for_frequency(self, frequency_Hz, starting_amplitude_uV = 100_000,
                                         tolerance=0.01, max_iterations=100):
        """
        Function to tune the AWG amplitude for a specific frequency to match the reference Irec value (self.default_i_rec). The function iteratively adjusts the amplitude until the recorded Irec value is within the specified tolerance of the reference Irec.
        Parameters        frequency_Hz (float): The frequency for which to tune the AWG amplitude.
                        starting_amplitude_uV (float): The starting amplitude for the tuning process in microvolts.
                        tolerance (float): The acceptable relative difference between the recorded Irec and the reference Irec
                        max_iterations (int): The maximum number of iterations to perform to avoid infinite loops.
        Returns
        tuned_amplitude_uV (float): The tuned amplitude in microvolts that achieves the desired Irec within the specified tolerance.
        """
        # TODO: Find proper startig and reference amplitude for the tuning process.
        # configure AWG to output the reference signal
        self.awg.configure_continuous_sine_wave(channel_number=1, frequency=frequency_Hz, 
                                                num_periods=1, starting_amplitude_uV=starting_amplitude_uV, approximate_frequency=False)

        tuned_amplitude_uV = starting_amplitude_uV
        iteration = 0
        upper_bound_irec = self.default_i_rec * (1 + tolerance)
        lower_bound_irec = self.default_i_rec * (1 - tolerance)
        while (self.get_irec() > upper_bound_irec or self.get_irec() < lower_bound_irec) and iteration < max_iterations:
            
            # TODO: find better tuning strategy, e.g. proportional control based on the difference between recorded Irec and default Irec, instead of just increasing or decreasing by a fixed percentage
            # or linear sweep
            # Nicolaj mentioned something as the momentum method
            if self.get_irec() > upper_bound_irec:
                tuned_amplitude_uV *= 0.9 # decrease amplitude by 10%
            else:
                tuned_amplitude_uV *= 1.1 # increase amplitude by 10%

            self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=tuned_amplitude_uV)
            time.sleep(self.i_rec_integration_time)
            iteration += 1

        # log the result
        print(f"Tuned amplitude for frequency {frequency_Hz} Hz: {tuned_amplitude_uV} uV, recorded Irec: {self.get_irec()} A, iterations: {iteration}")
        return tuned_amplitude_uV
    # if an error occurs, execute this command
    def escape_routine(self):
        """
        Function which is executed if an error occurs. 
        Sets the z-controller value to its default position and activates the controller.
    
        """

        print("Error occured!")
        print("Recovering to default state.")

        # TODO: check with Nicolaj for best sequence of operations
        # move to xy position
        wait_end_move = False
        self.nanonis_module.FolMe.XYPosSet(self.save_x_position, self.save_y_position, wait_end_move)

        # set the bias mode to the desired value
        self.nanonis_module.Bias.Set(self.safe_voltage)
            
        # deactivate controller
        if self.nanonis_module.ZCtl.OnOffGet()==1:
            self.nanonis_module.ZCtl.OnOffSet(0)

        # not in if-block, as in the worst case, the controller has been switched off just before the error occurred
        time.sleep(self.height_averaging_period_s + self.buffer_time_s)

        # set the desired current
        self.nanonis_module.ZCtl.SetpntSet(self.save_current_A)
        #time.sleep(199)
        time.sleep(5)

        # activate controller
        if self.nanonis_module.ZCtl.OnOffGet()==0:
            self.nanonis_module.ZCtl.OnOffSet(1)

    # function to save the recorded data
    def save_data(self):
        
        # save as json
        # save as a json file
        filename = self.filename

        # TODO: also save settings
        data = {
            "type": "dummy_data",
            "version": self.version,
            "Header": self.header,
            "Data": {
                "channel names": [ channel_name for channel_name, _ in self.recorded_data],
                "values": [ values for _, values in self.recorded_data]
            }
        }
        
        with open(filename, 'w') as f:
            json.dump(
                data,
                f,
                indent=4
            )


if __name__ == '__main__':

    print("This is a unit test")

    TCP_IP  = '127.0.0.1'                               # Local host
    TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
    version = 14000                                     # Nanonis RT Engine version number

    NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)  # This is how you establish a TCP connection. NTCP is the connection handle.
                                                        # Check your nanonis version in Nanonis > help > info and enter the RT Engine number
    NMod = NanonisModules.NanonisModules(NTCP)          # Load all nanonis modules

    print("Connected to the device")

    # test escape_routine
    save_voltage_V = 1.8
    save_height_m = 50e-9
    save_current_A = 14e-12
    off_delay_s = 0.1
    experiment = transfer_finder(nanonis_module=NMod, safe_voltage_V=save_voltage_V, 
                                 save_height_m=save_height_m, 
                                 save_current_A=save_current_A,
                                 off_delay_s=off_delay_s
                                 )

    experiment.escape_routine()

