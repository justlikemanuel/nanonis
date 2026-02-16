# module to run a measurement that allows computing the transfer function for some frequencies
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.AWG_M8195A_interface.M8195A import M8195A
import time
import numpy as np 
import json

class transfer_finder:
    def __init__(self, nanonis_module, safe_voltage_V, save_current_A,
                save_x_position, save_y_position, 
                i_rec_integration_time,
                max_allowed_amplitude_uV,
                amplitude_guess_mode,
                old_transfer_function,
                atom_tracking_parameters,
                atom_tracking_time_s, atom_tracking_interval,
                z_off_delay_s, awg_reference,
                sweep_frequencies, default_frequency,
                reference_amplitudes_uV,
                version, filename, header,
                 ):
        
        """
        Class to measure the transfer function for specified frequencies.

        Args:
            - nanonis_module: The NanonisModules object to interact with the Nanonis system.
            - safe_voltage_V: The voltage to set in case of an error to protect the sample and tip.
            - save_current_A: The current to set in case of an error to protect the sample and tip.
            - save_x_position: The x position to set in case of an error to protect the sample and tip.
            - save_y_position: The y position to set in case of an error to protect the sample and tip.
            - i_rec_integration_time: The integration time to use when recording the Irec value.
            - max_allowed_amplitude_uV: The maximum amplitude to apply with the AWG during the tuning process, in microvolts.
            - amplitude_guess_mode: The strategy to use for estimating the starting amplitude for the tuning process. Options are "known" and "half". "known" uses the recorded Irec values for the reference amplitudes to find the two reference frequencies that are closest to the desired frequency, and performs a linear interpolation to estimate the Irec value at the desired frequency, then finds the corresponding amplitude. "half" assumes 0.5 transmission to estimate the starting amplitude.
            - old_transfer_function: The old transfer function to use for estimating the starting amplitude for the tuning process. 
                Must be an array of frequency and transfer function value.
            - atom_tracking_parameters: Dictionary containing the parameters to use for atom tracking, e.g. a dictionary
            - atom_tracking_time_s: The time to track the atom for after each measurement step.
            - atom_tracking_interval: The inerval after how many measurement steps should be executed again.
            - z_off_delay_s: The time to wait after switching off the z-controller before setting
                            the current and switching the controller back on, to allow for height averaging.
            - awg_reference: The reference to the AWG to use for outputting the reference
            - sweep_frequencies: The frequencies for which to measure the transfer function.
            - default_frequency: The default frequency to use for measuring the reference Irec value.
            - reference_amplitudes_uV: The reference amplitudes to use for recording the Irec values, in microvolts.
            - version: The version of the experiment, to be saved in the data file.
            - filename: The name of the file to save the data to.
            - header: The header to save in the data file, e.g. a description of the experiment and the settings used.
        """

        
        self.nanonis_module = nanonis_module
        self.safe_voltage = safe_voltage_V
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
        self.max_allowed_amplitude_uV = max_allowed_amplitude_uV

        self.atom_tracking_parameters = atom_tracking_parameters
        # apply the parameters to the atom tracking module
        self.nanonis_module.ATrack.PropsSet(**self.atom_tracking_parameters)

        self.atom_tracking_period_s = atom_tracking_time_s
        self.atom_tracking_interval = atom_tracking_interval
        self.height_averaging_period_s = z_off_delay_s

        # logging parameters
        self.irec_vs_reference_amplitudes = []
        self.filename = filename
        self.version = version
        self.header = header


        # compensation parameters
        self.amplitude_guess_mode = amplitude_guess_mode
        self.old_transfer_function = old_transfer_function
        self.current_transfer_function = [] # list of tuples (frequency_Hz, transfer_function)
        self.default_i_rec = None # current value at the default amplitude
        self.default_transfer_function = None # transfer function at the default frequency and default amplitude
        self.max_amplitude_uV = None # voltage to achieve 1.5*Irec. this will be set during the recording routine, as it depends on the recorded Irec values

        self.recorded_data_headers = ["frequency_Hz", "transfer_function"] # TODO: add all parameters
        self.recorded_data_values = [] # list of tuples (frequency_Hz, transfer_function)

    ###################################################
    ############## Positioning functions ##############
    ###################################################
    
    # return to default state
    def return_to_default_state(self):

        # TODO: check with Nicolaj for best sequence of operations
        # move to xy position
        wait_end_move = False
        self.nanonis_module.FolMe.XYPosSet(self.save_x_position, self.save_y_position, wait_end_move)

        # TODO: what is better, set off-delay to 0 and then deactivate, or deactivate immediately and wait for some time? 

        # set off delay to 0
        self.nanonis_module.ZCtl.OffDelaySet(0)

        # deactivate controller
        if self.nanonis_module.ZCtl.OnOffGet()==1:
            self.nanonis_module.ZCtl.OnOffSet(0)

        
        # TODO: find best waiting time
        time.sleep(self.buffer_time_s)

        
        # set the bias mode to the desired value
        self.nanonis_module.Bias.Set(self.safe_voltage)
        
        # set the desired current
        self.nanonis_module.ZCtl.SetpntSet(self.save_current_A)
        
        # TODO: find proper waiting time before and after switching the controller on
        # activate z-controller
        self.nanonis_module.ZCtl.OnOffSet(1)             
 
    # if an error occurs, execute this command
    def escape_routine(self):
        """
        Function which is executed if an error occurs. 
        Sets the z-controller value to its default position and activates the controller.
    
        """
        # TODO: change print statements to logging messages
        print("Error occured!")
        print("Recovering to default state.")

        self.return_to_default_state()

    # function to go to measurement position and height
    def prepare_measurement(self):
        """
        Function to prepare the measurement position by moving to the specified xy position and setting the z-controller to the specified current setpoint.     
        """

        # move to default setpoint
        self.return_to_default_state()

               
        # ensure the z_off_delay is set to the desired value
        self.nanonis_module.ZCtl.OffDelaySet(self.height_averaging_period_s)

        # turn off the z-controller to allow for height averaging
        self.nanonis_module.ZCtl.OnOffSet(0)


    # function to track the atom
    def track_atom(self, tracking_time_s=None, new_I_gain=None, new_Frequency=None, new_Amplitude=None, new_Phase=None, new_SwitchOffDelay=None):
        """
        Function to track the atom for a fixed interval.
        """
        """# ensure modulation and controller are off
        # TODO: Raise error if any is on?
        if self.nanonis_module.ATrack.StatusGet('Controller')=='on':
            self.nanonis_module.ATrack.CtrlSet('Controller','off')
        if self.nanonis_module.ATrack.StatusGet('Modulation')=='on':
            self.nanonis_module.ATrack.CtrlSet('Modulation','off')
        """

        # apply new parameters if specified, otherwise keep old parameters
        # could be made more elegant
        if any(param is not None for param in [new_I_gain, new_Frequency, new_Amplitude, new_Phase, new_SwitchOffDelay]):
            self.nanonis_module.ATrack.PropsSet(Igain=new_I_gain, Frequency=new_Frequency, Amplitude=new_Amplitude, Phase=new_Phase, SwitchOffDelay=new_SwitchOffDelay)
        
        if tracking_time_s is None:
            tracking_time_s = self.atom_tracking_period_s

        # turn modulation and controller on
        self.nanonis_module.ATrack.CtrlSet('Modulation','on') # actually not needed, as activating the controller automatically starts the modulation
        self.nanonis_module.ATrack.CtrlSet('Controller','on')

        # track for the specified time
        time.sleep(tracking_time_s)

        # turn modulation and controller off
        self.nanonis_module.ATrack.CtrlSet('Modulation','off')
        self.nanonis_module.ATrack.CtrlSet('Controller','off') # actually not needed, as switching off the modulation automatically deactivates the controller


    #####################################################
    ############### Measurement functions ###############
    #####################################################
    

    # function to get the current Irec value
    def get_irec(self):
        """
        Function to get the current Irec value. The integration time can be set in the constructor.

        Returns
        Irec (float): The recorded Irec value in Amperes.
        """
        # TODO: Check if cur.get or signal.valget get is better


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
        max_current = self.default_i_rec

        while max_current < 1.5 * self.default_i_rec and iteration < max_iterations:
            max_amplitude_uV += amplitude_spacing_uV
            # check if the amplitude exceeds the maximum allowed amplitude to protect the sample and tip
            if max_amplitude_uV > self.max_allowed_amplitude_uV:
                print(f"Maximum allowed amplitude of {self.max_allowed_amplitude_uV} uV exceeded. Stopping the tuning process to protect the sample and tip.")
                print(f"Achieved Irec at maximum allowed amplitude: {max_current} A, which is {max_current/self.default_i_rec:.2f} times the default Irec.")
                # TODO: Log the achieved Irec value at the maximum allowed amplitude
                break
            self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=max_amplitude_uV)
            time.sleep(self.i_rec_integration_time)
            max_current = self.get_irec()
            iteration += 1

        self.max_amplitude_uV = max_amplitude_uV
        print(f"Max amplitude found: {max_amplitude_uV} uV, recorded Irec: {self.get_irec()} A")


    # function to tune awg amplitude for a specific frequency to match reference irec
    def tune_awg_amplitude_for_frequency(self, frequency_Hz, starting_amplitude_uV = 100_000,
                                         tolerance=0.01, max_iterations=100):
        """
        Function to tune the AWG amplitude for a specific frequency to match the reference Irec value (self.default_i_rec). 
        The function iteratively adjusts the amplitude until the recorded Irec value is within the specified tolerance of the reference Irec.
        Args:
            - frequency_Hz (float): The frequency for which to tune the AWG amplitude.
            - starting_amplitude_uV (float): The starting amplitude for the tuning process in microvolts.
            - tolerance (float): The acceptable relative difference between the recorded Irec and the reference Irec
            - max_iterations (int): The maximum number of iterations to perform to avoid infinite loops.

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
            # CAUTION: too much current rips the sample/tip
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
    
    # function to estimate the starting amplitude for the tuning process
    def estimate_starting_amplitude_for_frequency(self, frequency_Hz, mode):
        """
        Function to estimate the starting amplitude for the tuning process based on the recorded Irec values for the reference amplitudes. The function can use different strategies for the estimation, e.g. linear interpolation between the two reference amplitudes that are closest to the desired frequency, or using the transfer function at the default frequency to extrapolate the expected Irec value at the desired frequency and then find the corresponding amplitude.

        Args:
            - frequency_Hz (float): The frequency for which to estimate the starting amplitude.
            - mode (str): The strategy to use for the estimation. 
                Options are:
                - "known": use the recorded Irec values for the reference amplitudes to find the two reference frequencies that are closest to the desired frequency, and perform a linear interpolation to estimate the Irec value at the desired frequency, then find the corresponding amplitude.
                - "half": assume 0.5 transmission
        """
        starting_amplitude_uV = 0
        if mode == "known":
            # check if frequency is within the old transfer function range (1.st column)
            old_frequencies = self.old_transfer_function[:,0]
            if frequency_Hz in old_frequencies:
                # if the frequency is already in the old transfer function, use the corresponding transfer function value to estimate the starting amplitude
                idx = np.where(old_frequencies == frequency_Hz)[0][0]
                transfer_function_value = self.old_transfer_function[idx,1]
                starting_amplitude_uV = self.default_amplitude_uV * self.default_transfer_function / transfer_function_value

        if mode == "half":
            # assume 0.5 transmission
            starting_amplitude_uV = self.default_amplitude_uV * self.default_transfer_function / 0.5
        
        if mode == "closest":
            # find the transfer function value for the closest frequency of previously evaluated frequencies
            if len(self.current_transfer_function) > 0:
                current_frequencies = np.array([freq for freq, _ in self.current_transfer_function])
                closest_idx = np.argmin(np.abs(current_frequencies - frequency_Hz))
                _, closest_transfer_function = self.current_transfer_function[closest_idx]
                starting_amplitude_uV = self.default_amplitude_uV * self.default_transfer_function / closest_transfer_function

            print(f"Estimated starting amplitude for frequency {frequency_Hz} Hz: {starting_amplitude_uV} uV using mode {mode}")

        return starting_amplitude_uV
    
        
    # function to measure the transfer function for a single frequency
    def measure_transfer_function_for_frequency(self, frequency_Hz):
        """
        Function to measure the transfer function for a single frequency. The function tunes the AWG amplitude to match the reference Irec value, then records the Irec value for the tuned amplitude and the default amplitude, and computes the transfer function as the ratio of these two Irec values.

        Parameters
        frequency_Hz (float): The frequency for which to measure the transfer function.

        Returns
        transfer_function (float): The measure transfer function for the specified frequency.
        """
        amplitude_guess_mode = self.amplitude_guess_mode
        starting_amplitude_uV = self.estimate_starting_amplitude_for_frequency(frequency_Hz, mode=amplitude_guess_mode)
        tuned_amplitude_uV = self.tune_awg_amplitude_for_frequency(frequency_Hz=frequency_Hz, starting_amplitude_uV=starting_amplitude_uV)

        # compute transfer function based on the ratio of the applied and default amplitude, and the transfer function at the default amplitude.
        transfer_function = self.default_transfer_function * (self.default_amplitude_uV / tuned_amplitude_uV)
        
        # append to current transfer function list
        self.current_transfer_function.append((frequency_Hz, transfer_function))
        
        return transfer_function
    

    # function to actually measure the transfer function for the specified frequencies
    def measure_transfer_function_for_all_frequencies(self):
        """
        Function to measure the transfer function for all specified frequencies. The function iterates over the list of frequencies, measures the transfer function for each frequency using the measure_transfer_function_for_frequency function, and stores the results in a list.
        Adds all recorded values to the recorded_data attribute.
        
        transfer_functions (list of float): The list of measured transfer functions for the specified frequencies.
        """

        # iterate over all frequencies and measure the transfer function for each frequency
        for index, frequency_Hz in enumerate(self.sweep_frequencies):
            print(f"Measuring transfer function for frequency {frequency_Hz} Hz ({index+1}/{len(self.sweep_frequencies)})")
            transfer_function = self.measure_transfer_function_for_frequency(frequency_Hz)
            self.recorded_data_values.append((frequency_Hz, transfer_function))

            # execute atom tracking after a specified number of measurement steps
            if (index+1) % self.atom_tracking_interval == 0:
                print(f"Executing atom tracking after {index+1} measurement steps.")
                self.track_atom()

        # return to default state after the measurement is done
        self.return_to_default_state()
        # TODO: implement logging
        # return 0

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
                "channel names": [ channel_name for channel_name in self.recorded_data_headers],
                "values": [ values for _, values in self.recorded_data_values]
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
    save_current_A = 14e-12
    off_delay_s = 0.1
    experiment = transfer_finder(nanonis_module=NMod, safe_voltage_V=save_voltage_V, 

                                 save_current_A=save_current_A,
                                 off_delay_s=off_delay_s
                                 )

    experiment.escape_routine()

