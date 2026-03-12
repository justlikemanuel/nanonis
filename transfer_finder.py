# module to run a measurement that allows computing the transfer function for some frequencies
import matplotlib
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules

import time
import numpy as np 
import json
import datetime

class transferFinder:
    def __init__(self, nanonis_module,
                i_rec_integration_time_s,
                amplitude_guess_mode,
                old_transfer_function,
                atom_tracking_settings,
                atom_tracking_time_s, atom_tracking_interval,
                awg_reference, awg_settling_time,
                sweep_frequencies, reference_frequency,
                use_active_state,
                reference_amplitudes,
                measurement_voltage,
                filename, header,
                 ):
        
        """
        Class to measure the transfer function for specified frequencies.

        Args:
            - nanonis_module: The NanonisModules object to interact with the Nanonis system.
            - i_rec_integration_time_s: The integration time to use for recording the Irec value.
            - amplitude_guess_mode: The strategy to use for estimating the starting amplitude for the tuning process. Options are "known" and "half". "known" uses the recorded Irec values for the reference amplitudes to find the two reference frequencies that are closest to the desired frequency, and performs a linear interpolation to estimate the Irec value at the desired frequency, then finds the corresponding amplitude. "half" assumes 0.5 transmission to estimate the starting amplitude.
            - old_transfer_function: The old transfer function to use for estimating the starting amplitude for the tuning process. 
                Must be an array of frequency and transfer function value.
            - atom_tracking_settings: Dictionary containing the settings to use for atom tracking, e.g. a dictionary of parameters.
            - atom_tracking_time_s: The time to track the atom for after each measurement step.
            - atom_tracking_interval: The inerval after how many measurement steps should be executed again.
                                        If atom tracking shall not be used, set this to a value larger than the number of measurement steps.
            - awg_reference: The reference to the AWG to use for outputting the reference
            - awg_settling_time: The time to wait after changing the AWG settings before recording the Irec value, to allow the system to stabilize.
            - sweep_frequencies: The frequencies for which to measure the transfer function.
            - reference_frequency: The frequency to use for measuring the reference Irec value.
            - reference_amplitudes: The reference amplitudes to use for recording the Irec values, in Volts.
            - measurement_voltage: The voltage used for which the measurement shall be run.
            - filename: The name of the file to save the data to.
            - header: The header to save in the data file, e.g. a description of the experiment and the settings used.
        """

        # dummy parameters (should be used with the constructor)
        self.height_averaging_period_s = 0.5 # time to wait for the height to stabilize after switching off the z-controller, before recording the Irec value for the transfer function measurement
        self.default_frequency = reference_frequency
        
        # get current x_pos, y_pos, voltage and current to be able to return to this state
        self.initial_x_position_m = self.nanonis_module.FolMe.XPosGet()
        self.initial_y_position_m = self.nanonis_module.FolMe.YPosGet()
        self.initial_voltage = self.nanonis_module.Bias.Get()
        self.initial_current_A = self.nanonis_module.ZCtl.SetpntGet()
        
        # get current z-controller settings
        self.zcontroler_switch_off_delay_s = self.nanonis_module.ZCtl.SwitchOffDelayGet() # p_gain, time_constant, i_gain
        gain = self.nanonis_module.ZCtl.GainGet() # p_gain, time_constant, i_gain
        self.p_gain = gain[0]
        self.time_constant = gain[1]

        # AWG parameters
        self.awg = awg_reference
        self.awg_settling_time = awg_settling_time

        # nanonis        
        self.nanonis_module = nanonis_module
        self.i_rec_integration_time_s = i_rec_integration_time_s

        # escape routine
        self.voltage_tolerance = 1e-5
        self.current_tolerance = 0.1e-12
        self.is_in_error_state = False # flag to indicate if the system is in an error state, e.g. due to tip crash or excessive current
        self.escape_routine_voltage_step = 1e-3 # voltage step to apply in the escape routine
        self.escape_routine_current_step = 1e-12 # current step to apply in the escape routine       
        self.escape_routine_time_constant = 0.5  # time step between the update of the bias voltage


        # Sweep parameters:
        self.sweep_frequencies = sweep_frequencies
        self.measurement_voltage = measurement_voltage
    
   

        self.atom_tracking_settings = atom_tracking_settings
        # apply the settings to the atom tracking module
        self.nanonis_module.ATrack.PropsSet(**self.atom_tracking_settings)

        self.atom_tracking_period_s = atom_tracking_time_s
        self.atom_tracking_interval = atom_tracking_interval

        # logging parameters
        self.start_time = time.strftime("%Y-%m-%d_%H-%M-%S")
        self.session_path = self.get_session_path()
        self.filename = filename
        self.version = [1, 0, 0]
        self.header = header
        self.irec_vs_reference_amplitudes = []

        # compensation parameters
        self.amplitude_guess_mode = amplitude_guess_mode
        self.old_transfer_function = old_transfer_function
        self.current_transfer_function = [] # list of tuples (frequency_Hz, transfer_function)
        self.default_i_rec = None # current value at the default amplitude
        self.default_transfer_function = None # transfer function at the default frequency and default amplitude
        self.max_amplitude_uV = None # voltage to achieve 1.5*Irec. this will be set during the recording routine, as it depends on the recorded Irec values

        self.recorded_data_headers = ["frequency in Hz", "transfer_function", "tuned_amplitude in uV"] # TODO: add all parameters
        self.recorded_data_values = [] # list of tuples (frequency_Hz, transfer_function)

    ###################################################
    ############## Positioning functions ##############
    ###################################################
    
    # return to default state
    def return_to_default_state(self):

        # TODO: check with Nicolaj for best sequence of operations
        # move to xy position
        wait_end_move = False
        self.nanonis_module.FolMe.XYPosSet(self.safe_x_position_m, self.safe_y_position_m, wait_end_move)

        # TODO: what is better, set off-delay to 0 and then deactivate, or deactivate immediately and wait for some time? 

        # set off delay to 0
        self.nanonis_module.ZCtl.SwitchOffDelaySet(0)

        # deactivate controller
        if self.nanonis_module.ZCtl.OnOffGet()==1:
            self.nanonis_module.ZCtl.OnOffSet(0)

        
        # TODO: find best waiting time
        time.sleep(self.buffer_time_s)

        
        # set the bias mode to the desired value
        self.nanonis_module.Bias.Set(self.safe_voltage)
        
        # set the desired current
        self.nanonis_module.ZCtl.SetpntSet(self.safe_current_A)
        
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

       
        # Method: recursive try
        """
        1. Activate z-controller    
        2. Find current voltage and current values
        3. Set the bias voltage to the 
        4. if another error occurs, get again the current voltage and current values
        """
        #self.return_to_default_state()

        if not self.is_in_error_state:
            self.is_in_error_state = True

        try:
            # set parameters on z-controller and activate
            self.nanonis_module.ZCtl.SwitchOffDelaySet(0)            
            self.nanonis_module.ZCtl.GainSet(self.p_gain, self.time_constant) # nochmal ausschalten davor?
            self.nanonis_module.ZCtl.OnOffSet(1)

            # find current voltage value and change to default voltage
            measured_voltage = self.nanonis_module.Bias.Get()

            diff_voltage = measured_voltage - self.initial_voltage

            # iterate the voltage in the direction of the default voltage until we are back in the pre-measurement state
            while np.abs(diff_voltage) > self.voltage_tolerance:
                if abs(diff_voltage) < self.escape_routine_voltage_step:
                    step_voltage = diff_voltage
                else:
                    step_voltage = self.escape_routine_voltage_step * np.sign(diff_voltage)

                new_voltage = measured_voltage - step_voltage
                self.nanonis_module.Bias.Set(new_voltage)
                time.sleep(self.escape_routine_time_constant)

                measured_voltage = self.nanonis_module.Bias.Get()
                diff_voltage = measured_voltage - self.initial_voltage   

            # iterate the current in the direction of the default current until we are back in the pre-measurement state
            measured_current = self.nanonis_module.ZCtl.SetpntGet()
            diff_current = measured_current - self.initial_current_A

            while np.abs(diff_current) > self.current_tolerance:
                if abs(diff_current) < self.escape_routine_current_step:
                    step_current = diff_current
                else:
                    step_current = self.escape_routine_current_step * np.sign(diff_current)

                new_current = measured_current - step_current
                self.nanonis_module.ZCtl.SetpntSet(new_current)
                time.sleep(self.escape_routine_time_constant)

                measured_current = self.nanonis_module.ZCtl.SetpntGet()
                diff_current = measured_current - self.initial_current_A

        
            
        except Exception as e:
            print(f"Error during escape routine: {e}. Retrying...")
            self.escape_routine()

    # function to go to measurement position and height
    def prepare_measurement(self):
        """
        Function to prepare the measurement position by moving to the specified xy position and setting the z-controller to the specified current setpoint.     
        """

        # move to default setpoint
        #self.return_to_default_state()

        # run atom tracking
        print("Starting atom tracking.")
        self.track_atom(tracking_time_s=self.atom_tracking_period_s)

        # ensure the z_off_delay is set to the desired value
        self.nanonis_module.ZCtl.SwitchOffDelaySet(self.height_averaging_period_s)

        # turn off the z-controller to allow for height averaging
        self.nanonis_module.ZCtl.OnOffSet(0)
        print("Z-controller turned off for height averaging.")


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
    
    # function to get the session path
    def get_session_path(self):
        """
        Function to get the session path of the current measurement. This can be used for saving the data in the same folder as the measurement data.
        Returns
            - session_path (str): The path to the current measurement session.
        """
        return self.nanonis_module.Util.SessionPathGet()

    # function to get the current Irec value
    def get_irec(self):
        """
        Function to get the current Irec value. The integration time can be set in the constructor.

        Returns
        Irec (float): The recorded Irec value in Amperes.
        """
        # TODO: Check if cur.get or signal.valget get is better

        ####### signal.valget method:
        current_values_A = []
        start_time = time.perf_counter()
        while time.perf_counter() - start_time < self.i_rec_integration_time_s:
            i_tun = self.nanonis_module.Sig.ValGet(0, True) # channel 0 is current, we wait for next update
            current_values_A.append(i_tun)

        # compute the average current value over the integration time
        irec = np.mean(current_values_A)

        ######## cur.get method:
        # set the integration time
        # get the Irec value
        #irec = self.nanonis_module.Cur.Get()

        return irec

    # function to find a threshold value in a curve
    def find_threshold(self, x, y):
        """
        Function to compute the value x at the threshold in y(x)
        Params: 
            - x: list of x values
            - y: list of y values
        """
        # find x where y has the highest change
        dy = np.diff(y)
        dx = np.diff(x)
        derivative = dy/dx
        max_derivative_index = np.argmax(derivative)
        threshold_x = x[max_derivative_index]

        return threshold_x
    
    # function to measure the reference Irec value
    def find_reference_parameters(self):
        """
        Function to measure the reference Irec value at the default frequency and default amplitude. This value will be used as a reference for tuning the AWG amplitude for the other frequencies.
        """

        # configure AWG to output the reference signal
        self.awg.configure_continuous_sine_wave
        # turn on the AWG output
        self.awg.start_playing()
        # wait for the system to stabilize
        time.sleep(0.1)
        
        # record the Irec value for the default frequency and default amplitude
        self.default_i_rec = self.get_irec()

        # turn off the AWG output
        self.awg.stop_playing()
        print("Found Irec. AWG OFF.")

        # TODO: find how from this, we can compute the transfer function at the default frequency and default amplitude, which will be used as a reference for estimating the starting amplitude for the tuning process for the other frequencies.
        # compute the transfer function at the default frequency and default amplitude, which will be used as a reference for estimating the starting amplitude for the tuning process for the other frequencies.
        amplitudes, irec_values = zip(*self.irec_vs_reference_amplitudes)
        # CAUTION: Not all provided reference amplitudes can be played based on limited resolution. Hence, use the actually played amplitudes
        observed_threshold_voltage_uV = self.find_threshold(x=amplitudes, y=irec_values)
        self.default_transfer_function = self.measurement_voltage / observed_threshold_voltage_uV * 1e-6
        print(f"Measured reference Irec: {self.default_i_rec} A at default frequency {self.default_frequency} Hz and default amplitude {self.default_amplitude_uV} uV. Computed transfer function at default frequency: {self.default_transfer_function}")


    # record Irec vs the reference amplitudes
    def record_irec_for_references(self):

        # configure AWG to output the reference signal
        # channel_number, frequency_Hz, num_periods = 1, starting_amplitude_uV = 75_000, approximate_frequency = False):
        self.awg.configure_continuous_sine_wave(channel_number=1, frequency_Hz=self.default_frequency, 
                                                num_periods=1, starting_amplitude_uV=self.reference_amplitudes_uV[0], approximate_frequency=False)
        
        # activate the output of the AWG
        self.awg.start_playing()
        print("AWG ON")

        for index, amplitude_uV in enumerate(self.reference_amplitudes_uV):
            matched_amplitude_uV = self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=amplitude_uV)
            time.sleep(self.i_rec_integration_time_s)
            irec = self.get_irec()

            # add to list
            self.irec_vs_reference_amplitudes.append((matched_amplitude_uV, irec))
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

        
        # record Irec values while increasing the amplitude until we reach 1.5 times the default Irec, which will be used as a reference for the maximum amplitude to apply during the tuning process for the other frequencies, to protect the sample and tip from excessive current.
        while max_current < 1.5 * self.default_i_rec and iteration < max_iterations:
            max_amplitude_uV += amplitude_spacing_uV
            # check if the amplitude exceeds the maximum allowed amplitude to protect the sample and tip
            if max_amplitude_uV > self.max_allowed_amplitude_uV:
                print(f"Maximum allowed amplitude of {self.max_allowed_amplitude_uV} uV exceeded. Stopping the tuning process to protect the sample and tip.")
                print(f"Achieved Irec at maximum allowed amplitude: {max_current} A, which is {max_current/self.default_i_rec:.2f} times the default Irec.")
                # TODO: Log the achieved Irec value at the maximum allowed amplitude
                break
            matched_amplitude_uV = self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=max_amplitude_uV)
            time.sleep(self.i_rec_integration_time_s)
            max_current = self.get_irec() # should be monotonically increasing, so we can just check the last value to see if we have reached the threshold

             # add data point to list
            self.irec_vs_reference_amplitudes.append((matched_amplitude_uV, max_current))
            iteration += 1

        self.max_amplitude_uV = max_amplitude_uV
        print(f"Max amplitude found: {max_amplitude_uV} uV, recorded Irec: {self.get_irec()} A")

        # stop the AWG output
        self.awg.stop_playing()
        print("AWG OFF")
   
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

        # TODO: Find proper starting and reference amplitude for the tuning process.
        # configure AWG to output the reference signal
        # verify that the starting amplitude is integer
        starting_amplitude_uV = int(starting_amplitude_uV)
        self.awg.configure_continuous_sine_wave(channel_number=1, frequency_Hz=frequency_Hz, 
                                                num_periods=1, starting_amplitude_uV=starting_amplitude_uV, approximate_frequency=False)

        # turn on the AWG output
        self.awg.start_playing()
        print(f"AWG ON for frequency {frequency_Hz} Hz, starting amplitude {starting_amplitude_uV} uV")
        time.sleep(0.1) # wait for the system to stabilize

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

            # clip to int. the awg function will clip to the resolution and return the applied value
            tuned_amplitude_uV = int(tuned_amplitude_uV)
            tuned_amplitude_uV = self.awg.update_continuous_sine_wave_amplitude(new_amplitude_uV=tuned_amplitude_uV)
            time.sleep(self.i_rec_integration_time)
            iteration += 1

        # turn off the AWG output
        self.awg.stop_playing()

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

        if mode == "fixed":
            starting_amplitude_uV = 200_000
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
        print(f"Estimated starting amplitude for frequency {frequency_Hz} Hz: {starting_amplitude_uV} uV using mode {amplitude_guess_mode}")
        tuned_amplitude_uV = self.tune_awg_amplitude_for_frequency(frequency_Hz=frequency_Hz, starting_amplitude_uV=starting_amplitude_uV)

        # compute transfer function based on the ratio of the applied and default amplitude, and the transfer function at the default amplitude.
        transfer_function = self.default_transfer_function * (self.default_amplitude_uV / tuned_amplitude_uV)
        
        # append to current transfer function list
        self.current_transfer_function.append((frequency_Hz, transfer_function))
        
        return transfer_function, tuned_amplitude_uV
    

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
            transfer_function, tuned_amplitude_uV = self.measure_transfer_function_for_frequency(frequency_Hz)
            self.recorded_data_values.append((frequency_Hz, transfer_function, tuned_amplitude_uV))

            # execute atom tracking after a specified number of measurement steps
            if (index+1) % self.atom_tracking_interval == 0:
                print(f"Executing atom tracking after {index+1} measurement steps.")

                # turn AWG off
                self.awg.stop_playing()
                print("AWG OFF for atom tracking")

                # tracking
                self.track_atom()

                # turn AWG on again
                self.awg.start_playing()
                print("AWG ON after atom tracking")

        # return to default state after the measurement is done
        self.return_to_default_state()
        # TODO: implement logging
        # return 0


    # function to save the reference Irec values for the reference amplitudes
    def save_reference_irec_values(self):
        # save the recorded Irec values for the reference amplitudes as a json file
        filename = f"{self.session_path}/reference_irec_values_{self.default_frequency}Hz_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        
        data = {
            "default_frequency_Hz": self.default_frequency,
            "reference_amplitudes_uV": [ amplitude for amplitude, _ in self.irec_vs_reference_amplitudes],
            "irec_values_A": [ irec for _, irec in self.irec_vs_reference_amplitudes]
        }
        
        with open(filename, 'w') as f:
            json.dump(
                data,
                #f,
                indent=4
            )
        print(f"Saved reference Irec values for the reference amplitudes to {filename}")

    
    # function to plot the reference Irec values for the reference amplitudes
    def plot_reference_irec_values(self):

        import matplotlib.pyplot as plt
        matplotlib.use('Agg') # use non-interactive backend to avoid issues on headless systems
        # save the recorded Irec values for the reference amplitudes as a json file
        fixsize = (12,6)
        plt.figure(figsize=fixsize)
        font_size = 16
        title_size = 18
        plt.rcParams.update({'font.size': font_size})
        title = "Reference Irec values for the reference amplitudes at " + str(self.default_frequency) + " Hz"
        plt.title(title, fontsize=title_size, fontweight='bold')

        amplitudes, irec_values = zip(*self.irec_vs_reference_amplitudes)
        irec_values = np.array(irec_values)
        amplitudes = 1e-6 * np.array(amplitudes) # convert to volts for plotting
        print("Shape of amplitudes:", np.shape(amplitudes))
        print("Shape of irec_values:", np.shape(irec_values))
        plt.plot(amplitudes, irec_values, marker='x', linestyle='-', color='blue')
        plt.xlabel("Amplitude (uV)", fontweight='bold')
        plt.ylabel("Irec (A)", fontweight='bold')
        plt.grid()
        plt.tight_layout()

        datetime_string = time.strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{self.session_path}/reference_irec_values_{self.default_frequency}Hz_{datetime_string}.png"
        dpi = 300 # save the plot with high resolution

        plt.savefig(filename, dpi=dpi)
        print(f"Saved reference Irec values for the reference amplitudes plot to {filename}")



    # function to save the recorded data
    def save_data(self):
        # get current date and time
        current_time = time.strftime("%Y%m%d-%H%M%S")
        filename = f"{self.session_path}/{self.filename}_{current_time}"

        print(f"Saving data to {filename}...")
        print("Data to be saved:")
        print(self.recorded_data_values)


        # TODO: also save settings
        data = {
            "type": "dummy_data",
            "version": self.version,
            "start_time": self.start_time,
            "end_time": current_time,
            "Header": self.header,
            "Data": {
                "channel names": [ channel_name for channel_name in self.recorded_data_headers],
                "values": [ values for values in self.recorded_data_values]
                                }
        }
        
        with open(filename, 'w') as f:
            json.dump(
                data,
                #f,
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
    safe_voltage_V = 1.8
    safe_current_A = 14e-12
    off_delay_s = 0.1
    experiment = transferFinder(nanonis_module=NMod, safe_voltage_V=safe_voltage_V, 

                                 safe_current_A=safe_current_A,
                                 off_delay_s=off_delay_s
                                 )

    experiment.escape_routine()

