# module to run a measurement that allows computing the transfer function for some frequencies
import matplotlib
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
from libs.pyNanonisMeasurements.nanonisTCP import NanonisModules
from libs.pyNanonisMeasurements.measurementClasses.MeasurementBase import MeasurementBase
from libs.regulator.pi_controller import PIController

import time
import numpy as np 
import json
import datetime

import logging
logger = logging.getLogger("transfer_finder")

class transferFinder:
    def __init__(self, nanonis_module,
                i_rec_integration_time_s = 0.1,
                amplitude_guess_mode = "closest",
                old_compensation_amplitudes = None,
                atom_tracking_settings = None,
                atom_tracking_time = 0.4,
                atom_tracking_interval = 5,
                awg_reference = None, 
                awg_settling_time = 0.1,
                lockin_frequency = 1e3,
                tuning_pgain = 0.5,
                tuning_integration_time_constant = 1.0,
                sweep_frequencies = None,
                reference_frequency = None,
                reference_amplitude = 0.5,
                max_sweep_amplitude = 0.5,
                num_sweep_amplitudes = 5,
                data_channels = None,
                use_active_state = False,
                measurement_voltage = 0.5,
                irec_tolerance = 0.1e-13,
                header = "dummy header",
                filename = "transfer_function_measurement",
                 ):
        
        """
        Class to measure the transfer function for specified frequencies.

        Args:
            - nanonis_module: The NanonisModules object to interact with the Nanonis system.
            - i_rec_integration_time_s: The integration time to use for recording the Irec value.
            - amplitude_guess_mode: The strategy to use for estimating the starting amplitude for the tuning process. Options are "known" and "half". "known" uses the recorded Irec values for the reference amplitudes to find the two reference frequencies that are closest to the desired frequency, and performs a linear interpolation to estimate the Irec value at the desired frequency, then finds the corresponding amplitude. "half" assumes 0.5 transmission to estimate the starting amplitude.
            - old compensation_amplitudes (list of tuples): The list of previously evaluated frequencies and their corresponding compensation amplitudes.
            - atom_tracking_settings: Dictionary containing the settings to use for atom tracking, e.g. a dictionary of parameters.
            - atom_tracking_time: The time to track the atom for after each measurement step.
            - atom_tracking_interval: The inerval after how many measurement steps should be executed again.
                                        If atom tracking shall not be used, set this to a value larger than the number of measurement steps.
            - awg_reference: The reference to the AWG to use for outputting the reference
            - awg_settling_time: The time to wait after changing the AWG settings before recording the Irec value, to allow the system to stabilize.
            - lockin_frequency: The frequency of the lock-in amplifier.
            - tuning_pgain: The proportional gain to use for the tuning process.
            - tuning_integration_time_constant: The time constant for the integral action in the tuning process.
            - tolerance: The tolerance for the tuning process (relative to the current measured at reference amplitude and reference frequency).
            - sweep_frequencies: The frequencies for which the transfer function shall be measured.
            - reference_frequency: The frequency at which the reference Irec value shall be recorded.
            - reference_amplitude: The amplitude at which the reference Irec value shall be recorded.
            - max_sweep_amplitude: The maximum amplitude to use for the sweep, to protect the tip and sample. The actual amplitudes used for the sweep will be generated based on this value and the number of amplitudes.
            - num_sweep_amplitudes: The number of amplitudes to use for the sweep between 0 and the maximum sweep amplitude.
            - data_channels: A tupel of nanonis channel label and desired header for the recorded data.
            - use_active_state: TODO: check with Nicolaj again.
            - measurement_voltage: The voltage used for which the measurement shall be run.
            - irec_tolerance: The tolerance for the Irec value when comparing to the reference Irec value for the compensation amplitude tuning.
            - filename: The name of the file to save the data to.
            - header: The header to save in the data file, e.g. a description of the experiment and the settings used.
        """

        
        # dummy parameters (TODO: should be used with the constructor)
        self.height_averaging_period_s = 0.5 # time to wait for the height to stabilize after switching off the z-controller, before recording the Irec value for the transfer function measurement
        self.reference_frequency = reference_frequency
        self.max_allowed_amplitude = 1 # maximum allowed amplitude in Volts to protect the sample and tip, TODO: find better parameter for this, maybe based on the recorded Irec values for the reference amplitudes
        
        
        # AWG parameters
        self.awg = awg_reference
        self.awg_settling_time = awg_settling_time
        self.lockin_frequency = lockin_frequency

        # create integrator
        self.tuning_controller = PIController(Kp=tuning_pgain, Ti=tuning_integration_time_constant, 
                                              dt=0.1, V_min=0.0, V_max=self.max_allowed_amplitude)
        self.irec_tolerance = irec_tolerance


        # nanonis        
        self.nanonis_module = nanonis_module
        self.i_rec_integration_time_s = i_rec_integration_time_s

        # get current x_pos, y_pos, voltage and current to be able to return to this state
        x_pos, y_pos = self.nanonis_module.FolMe.XYPosGet(Wait_for_newest_data=True)
        self.initial_x_position_m = x_pos
        self.initial_y_position_m = y_pos
        self.initial_voltage = self.nanonis_module.Bias.Get()
        self.initial_current_A = self.nanonis_module.ZCtl.SetpntGet()
        
        # get current nanonis settings
        self.initial_z_controler_switch_off_delay_s = self.nanonis_module.ZCtl.SwitchOffDelayGet() # p_gain, time_constant, i_gain
        gain = self.nanonis_module.ZCtl.GainGet() # p_gain, time_constant, i_gain
        self.initial_z_p_gain = gain[0]
        self.initial_z_time_constant = gain[1]
        self.initial_tracking_settings = self.nanonis_module.ATrack.PropsGet() # Igain, Frequency, Amplitude, Phase, SwitchOffDelay


        # escape routine
        self.voltage_tolerance = 1e-5
        self.current_tolerance = 0.1e-12
        self.is_in_error_state = False # flag to indicate if the system is in an error state, e.g. due to tip crash or excessive current
        self.escape_routine_voltage_step = 1e-3 # voltage step to apply in the escape routine
        self.escape_routine_current_step = 1e-12 # current step to apply in the escape routine       
        self.escape_routine_time_constant = 0.5  # time step between the update of the bias voltage


        # Sweep parameters:
        self.sweep_frequencies = sweep_frequencies
        self.max_sweep_amplitude = max_sweep_amplitude
        self.num_amplitudes = num_sweep_amplitudes

        # generate reference amplitudes
        step = self.max_sweep_amplitude / self.num_amplitudes
        self.sweep_amplitudes = [step * (i+1) for i in range(self.num_amplitudes)]
        logger.info(f"Generated reference amplitudes: {self.sweep_amplitudes} V")

        self.measurement_voltage = measurement_voltage
        self.reference_amplitude = reference_amplitude
        self.old_compensation_amplitudes = old_compensation_amplitudes
    
   
        # atom tracking
        self.atom_tracking_settings = atom_tracking_settings
        # apply the settings to the atom tracking module
        self.nanonis_module.ATrack.PropsSet(**self.atom_tracking_settings)
        self.atom_tracking_time = atom_tracking_time
        self.atom_tracking_interval = atom_tracking_interval

        # logging parameters
        self.start_time = time.strftime("%Y-%m-%d_%H-%M-%S")
        self.session_path = self.get_session_path()
        self.filename = filename
        self.version = [0, 0, 1]
        self.header = header
        self.irec_vs_sweep_amplitudes = []

        # TODO: use indices or names?
        self.data_indices = [
                            0, # current
                            24, # Bias
                            30, # z-controller setpoint
        ]



        self.recorded_data_headers = [
            "frequency (Hz)",
            "compensation_amplitude (V)",
            "Current (A)",
            "Bias (V)",
            "z-controller setpoint (m)",
        ]

        #self.nanonis_channel_names = []

        # when using channel indices for logging
        channel_names = self.nanonis_module.Sig.NamesGet()

        # add the headers for the recorded data channels
        for channel_label, channel_header in data_channels:
            self.recorded_data_headers.append(channel_header)

            # when using channel labels for logging
            # self.nanonis_channel_names.append(channel_label)

            # when using data indices for logging
            # if channel name is not in the list of channel names, raise error
            if channel_label not in channel_names:
                raise ValueError(f"Channel label {channel_label} not found in the list of channel names: {channel_names}")
            channel_index = channel_names.index(channel_label)
            self.data_indices.append(channel_index)

        self.recorded_data_values = [] # list of tuples (frequency, tuned_amplitude, current, bias, z_controller_setpoint,...)

        # compensation parameters
        self.amplitude_guess_mode = amplitude_guess_mode
        self.reference_i_rec = None # current value at the reference amplitude


        # store all parameters in a dictionary for logging
        self.nanonis_settings = {
                    "i_rec_integration_time_s": i_rec_integration_time_s,
                    "atom_tracking_settings": atom_tracking_settings,
                    "atom_tracking_time": atom_tracking_time,
                    "atom_tracking_interval": atom_tracking_interval,
                    "data_channels": data_channels,
                    "measurement_voltage": measurement_voltage,
                    "use_active_state": use_active_state,
                    "initial_x_position_m": self.initial_x_position_m,
                    "initial_y_position_m": self.initial_y_position_m,
                    "initial_bias_voltage": self.initial_voltage,
                    "initial_current_A": self.initial_current_A,
                    "z-controller_settings": {
                                    "switch_off_delay_s": self.initial_z_controler_switch_off_delay_s,
                                    "p_gain": self.initial_z_p_gain,
                                    "time_constant": self.initial_z_time_constant,
                                   }
                    }
        
        self.awg_settings = {
                    "awg_settling_time": awg_settling_time,
                    "lockin_frequency": lockin_frequency,
                    "sweep_frequencies": sweep_frequencies,
                    "reference_frequency": reference_frequency,
                    "reference_amplitude": reference_amplitude,
                    "max_sweep_amplitude": max_sweep_amplitude,
                    "num_sweep_amplitudes": num_sweep_amplitudes,
                }
        
        self.tuning_settings = {
                    "tuning_pgain": tuning_pgain,
                    "tuning_integration_time_constant": tuning_integration_time_constant,
                    "irec_tolerance": irec_tolerance,
                }

        print(f"Session path: {self.session_path}")        
        
        # get the nanonis parameters (to put into the logged file)
        #meas = MeasurementBase(self.nanonis_module)
        #self.nanonis_parameters = meas.nanonisSettingsGet(False)
        self.nanonis_parameters = {"Dummy_parameter": 0} # the nanonis function above just broke...

    ###################################################
    ############## Positioning functions ##############
    ###################################################
 
    # return to default state
    def return_to_starting_state(self):

        """
        Sets all parameters as before the measurement.
        TODO: check what exactly is needed
        """

        try:
            # TODO: check with Nicolaj for best sequence of operations
            # move to xy position
            wait_end_move = False
            self.nanonis_module.FolMe.XYPosSet(self.initial_x_position_m, self.initial_y_position_m, wait_end_move)

            # TODO: what is better, set off-delay to 0 and then deactivate, or deactivate immediately and wait for some time? 

            # set off delay to 0
            self.nanonis_module.ZCtl.SwitchOffDelaySet(0)

            # deactivate controller
            if self.nanonis_module.ZCtl.OnOffGet()==1:
                self.nanonis_module.ZCtl.OnOffSet(0)

            
            # TODO: find best waiting time

            # set the bias mode to the desired value
            self.nanonis_module.Bias.Set(self.initial_voltage)
            
            # set the desired current
            self.nanonis_module.ZCtl.SetpntSet(self.initial_current_A)
            
            # TODO: find proper waiting time before and after switching the controller on
            # activate z-controller
            self.nanonis_module.ZCtl.OnOffSet(1)    

        except Exception as e:
            print(f"Error while returning to starting state: {e}. Executing escape routine.")
            self.escape_routine()         
 
    # if an error occurs, execute this command
    def escape_routine(self):
        """
        Function which is executed if an error occurs. 
        Sets the z-controller value to its default position and activates the controller.
    
        """
        # TODO: change print statements to logging messages
        print("Error occured!")
        print("Recovering to default state.")

        # turn off AWG
        self.awg.stop_playing()

       
        # Method: recursive try
        """
        1. Activate z-controller    
        2. Find current voltage and current values
        3. Set the bias voltage to the pre-measurement voltage
        4. if another error occurs -> recursive call of the escape routine
        """

        if not self.is_in_error_state:
            self.is_in_error_state = True

        try:
            # set parameters on z-controller and activate
            self.nanonis_module.ZCtl.SwitchOffDelaySet(0)            
            self.nanonis_module.ZCtl.GainSet(self.initial_z_p_gain, self.initial_z_time_constant) # nochmal ausschalten davor?
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

        # TODO: what shall happen after recovering?
        exit(1)


    # function to go to measurement position and height
    def prepare_measurement(self):
        """
        Function to prepare the measurement.
        """
        try:
            # TODO: implementation of active state
        
            if self.atom_tracking_interval < len(self.sweep_frequencies):
                # run atom tracking
                print("Starting atom tracking.")
                self.track_atom()
            # ensure the z_off_delay is set to the desired value
            self.nanonis_module.ZCtl.SwitchOffDelaySet(self.height_averaging_period_s)

            # turn off the z-controller to allow for height averaging
            self.nanonis_module.ZCtl.OnOffSet(0)
            print("Z-controller turned off for height averaging.")
        
        except Exception as e:
            print(f"Error while preparing measurement: {e}. Executing escape routine.")
            self.escape_routine()

    # function to track the atom
    def track_atom(self):
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
        try:
            # turn modulation and controller on
            self.nanonis_module.ATrack.CtrlSet('Controller','on')

            # track for the specified time
            time.sleep(self.atom_tracking_time)
            # turn modulation and controller off
            self.nanonis_module.ATrack.CtrlSet('Modulation','off')
            
        except Exception as e:
            print(f"Error while tracking atom: {e}. Executing escape routine.")
            self.escape_routine()


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

    # record Irec vs the reference amplitudes
    def record_irec_for_references(self):
        
        """
        Function to record the Irec value for the reference amplitude, and all sweep amplitudes at the reference frequency.
        """
        try:
            # configure AWG to output the reference signal
            self.awg.configure_continuous_sine_wave(frequency=self.reference_frequency, lockin_frequency=self.lockin_frequency, starting_amplitude=self.reference_amplitude)
            
            # activate the output of the AWG and measure at reference amplitude
            self.awg.start_playing()
            time.sleep(self.awg_settling_time)
            self.reference_i_rec = self.get_irec()
            print("AWG ON")

            for amplitude in self.sweep_amplitudes:

                matched_amplitude = self.awg.update_continuous_sine_wave_amplitude(new_amplitude=amplitude)
                time.sleep(self.awg_settling_time)
                irec = self.get_irec()

                # add data to list
                self.irec_vs_sweep_amplitudes.append((matched_amplitude, irec))

            # stop the AWG output
            self.awg.stop_playing()
            print("AWG OFF")

        except Exception as e:
            print(f"Error while recording Irec for reference amplitudes: {e}. Executing escape routine.")
            self.escape_routine()
   
    # function to tune awg amplitude for a specific frequency to match reference irec
    def tune_awg_amplitude_for_frequency(self, frequency, starting_amplitude = 0.1,
                                         tolerance=0.01, max_iterations=100):
        """
        Function to tune the AWG amplitude for a specific frequency to match the reference Irec value (self.reference_i_rec). 
        The function iteratively adjusts the amplitude until the recorded Irec value is within the specified tolerance of the reference Irec.
        Args:
            - frequency (float): The frequency for which to tune the AWG amplitude.
            - starting_amplitude (float): The starting amplitude for the tuning process in Volts.
            - tolerance (float): The acceptable relative difference between the recorded Irec and the reference Irec
            - max_iterations (int): The maximum number of iterations to perform to avoid infinite loops.

        Returns
        tuned_amplitude (float): The tuned amplitude in microvolts that achieves the desired Irec within the specified tolerance.
        """
        try:
            # TODO: Find proper starting and reference amplitude for the tuning process.
            # configure AWG to output the reference signal
            self.awg.configure_continuous_sine_wave(frequency=frequency,
                                                    lockin_frequency=self.lockin_frequency,
                                                    starting_amplitude=starting_amplitude,
                                                )

            # turn on the AWG output
            self.awg.start_playing()
            print(f"AWG ON for frequency {frequency} Hz, starting amplitude {starting_amplitude} V")
            time.sleep(self.awg_settling_time)

            tuned_amplitude = starting_amplitude
            print("----------------------------------------")
            print(f"Starting tuning for frequency {frequency} Hz. Starting amplitude: {tuned_amplitude} V, reference Irec: {self.reference_i_rec} A")
            iteration = 0
            upper_bound_irec = self.reference_i_rec * (1 + tolerance)
            lower_bound_irec = self.reference_i_rec * (1 - tolerance)

            i_rec = self.get_irec()

            while ((i_rec > upper_bound_irec or i_rec < lower_bound_irec) 
                    and iteration < max_iterations):
                
                self.tuning_controller.update(I_ref=self.reference_i_rec, I_meas=i_rec)
                tuned_amplitude = self.tuning_controller.V_out
                
                self.awg.update_continuous_sine_wave_amplitude(new_amplitude=tuned_amplitude)
                time.sleep(self.i_rec_integration_time_s)
                i_rec = self.get_irec()
                iteration += 1
                
                """                
                # TODO: find better tuning strategy, e.g. proportional control based on the difference between recorded Irec and default Irec, instead of just increasing or decreasing by a fixed percentage
                # or linear sweep
                # CAUTION: too much current rips the sample/tip
                # Nicolaj mentioned something as the momentum method
                if self.get_irec() > upper_bound_irec:
                    tuned_amplitude *= 0.9 # decrease amplitude by 10%
                else:
                    tuned_amplitude *= 1.1 # increase amplitude by 10%

                # the awg function will clip to the resolution and return the applied value
                tuned_amplitude = self.awg.update_continuous_sine_wave_amplitude(new_amplitude=tuned_amplitude)
                time.sleep(self.i_rec_integration_time_s)
                iteration +=1 
                """

            # turn off the AWG output
            self.awg.stop_playing()

            # log the result
            print(f"Tuned amplitude for frequency {frequency} Hz: {tuned_amplitude} V, recorded Irec: {self.get_irec()} A, iterations: {iteration}")
            
            return tuned_amplitude
    
        except Exception as e:
            print(f"Error while tuning AWG amplitude for frequency {frequency} Hz: {e}. Executing escape routine.")
            self.escape_routine()

    
    # function to estimate the starting amplitude for the tuning process
    def estimate_starting_amplitude_for_frequency(self, frequency, mode):
        """
        Function to estimate the starting amplitude for the tuning process based on the recorded Irec values for the reference amplitudes. The function can use different strategies for the estimation, e.g. linear interpolation between the two reference amplitudes that are closest to the desired frequency, or using the transfer function at the default frequency to extrapolate the expected Irec value at the desired frequency and then find the corresponding amplitude.

        Args:
            - frequency (float): The frequency for which to estimate the starting amplitude.
            - mode (str): The strategy to use for the estimation. 
                Options are:
                - "known": 
                - "closest": use the recorded Irec values for the reference amplitudes to find the reference frequency that is closest to the desired frequency, and perform a linear interpolation to estimate the Irec value at the desired frequency, then find the corresponding amplitude.
                - "half": assume 0.5 transmission 
                - "fixed": use a fixed value for the starting amplitude, currently 0.2 V
                - idea: 50% transmission -> problem: what is the desired amplitude at the output of the channel?
        """
        if mode not in ["known", "closest", 
                        #"half", 
                        "fixed"]:
            raise ValueError(f"Invalid mode for estimating starting amplitude: {mode}. Valid options are 'known', 'closest', 'half', 'fixed'.")
        
        starting_amplitude = 0
        if mode == "known":
            # could be modified to not only use the same but a close old frequency
            # check if frequency is within the old transfer function range (1.st column)
            old_frequencies = self.old_compensation_amplitudes[:, 0]
            old_amplitudes = self.old_compensation_amplitudes[:, 1]

            if frequency in old_frequencies:
                # if the frequency is already in the old transfer function, use the corresponding transfer function value to estimate the starting amplitude
                idx = np.where(old_frequencies == frequency)[0][0]
                starting_amplitude = old_amplitudes[idx]
    
        if mode == "closest":
            # find the transfer function value for the closest frequency of previously evaluated frequencies
            if len(self.recorded_data_values) > 0:
                recorded_frequencies = [ data[0] for data in self.recorded_data_values]
                recorded_amplitudes = [ data[1] for data in self.recorded_data_values]

                closest_index = np.argmin(np.abs(np.array(recorded_frequencies) - frequency))
                starting_amplitude = recorded_amplitudes[closest_index]
            
            else:
                # if there are no previously recorded frequencies, use the default amplitude as a starting point
                starting_amplitude = 0.2 # As for fixed value, TODO: find better value
    
        if mode == "fixed":
            starting_amplitude = 0.2 # fixed value for testing, TODO: find better parameter for this

        print(f"Estimated starting amplitude for frequency {frequency} Hz: {starting_amplitude} V using mode {mode}")

        return starting_amplitude
    
        
    # function to measure the transfer function for a single frequency
    def measure_transfer_function_for_frequency(self, frequency):
        """
        Loggs all desired values and adds the row to the recorded data list.
        """
        try:
            starting_amplitude = self.estimate_starting_amplitude_for_frequency(frequency=frequency, mode=self.amplitude_guess_mode)
            #print(f"Estimated starting amplitude for frequency {frequency} Hz: {starting_amplitude} V using mode {self.amplitude_guess_mode}")
            tuned_amplitude = self.tune_awg_amplitude_for_frequency(frequency=frequency, starting_amplitude=starting_amplitude)

            # get data for all elements in the data_indices list and add the values to the recorded data list
            values  = self.nanonis_module.Sig.ValsGet(signal_index_list=self.data_indices, wait_for_newest_data=True)
            data_list = [frequency, tuned_amplitude, *values]
            self.recorded_data_values.append(data_list)
            return 0

        except Exception as e:
            print(f"Error while measuring transfer function for frequency {frequency} Hz: {e}. Executing escape routine.")
            self.escape_routine()


    # function to actually measure the transfer function for the specified frequencies
    def measure_transfer_function_for_all_frequencies(self):
        """
        Function to measure the transfer function for all specified frequencies. The function iterates over the list of frequencies, measures the transfer function for each frequency using the measure_transfer_function_for_frequency function, and stores the results in a list.
        Adds all recorded values to the recorded_data attribute.
        
        transfer_functions (list of float): The list of measured transfer functions for the specified frequencies.
        """

        try:
            # iterate over all frequencies and measure the transfer function for each frequency
            for index, frequency in enumerate(self.sweep_frequencies):
                #print(f"Measuring transfer function for frequency {frequency} Hz ({index+1}/{len(self.sweep_frequencies)})")
                
                # perform the measurement
                self.measure_transfer_function_for_frequency(frequency)

                # execute atom tracking after a specified number of measurement steps
                if (index+1) % self.atom_tracking_interval == 0:
                    print(f"Executing atom tracking after {index+1} measurement steps.")

                    # TODO: Do/Check anything on AWG?

                    # tracking
                    self.track_atom()

            # return to default state after the measurement is done
            self.return_to_starting_state()
            return 0

        except Exception as e:
            print(f"Error while measuring transfer function for all frequencies: {e}. Executing escape routine.")
            self.escape_routine()

    # function to save the reference Irec values for the reference amplitudes
    def save_reference_irec_values(self):
        # save the recorded Irec values for the reference amplitudes as a json file
        filename = f"{self.session_path}/reference_irec_values_{self.reference_frequency}Hz_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        
        data = {
            "reference_frequency": self.reference_frequency,
            "sweep_amplitudes": [ amplitude for amplitude, _ in self.irec_vs_sweep_amplitudes],
            "irec_values_A": [ irec for _, irec in self.irec_vs_sweep_amplitudes]
        }
        
        with open(filename, 'w') as f:
            json.dump(
                data,
                f,
                indent=4
            )
        print(f"Saved reference Irec values for the reference amplitudes to {filename}")
        
        return 0

    
    # function to plot the reference Irec values for the reference amplitudes
    def plot_reference_irec_values(self):

        import matplotlib.pyplot as plt
        matplotlib.use('Agg') # use non-interactive backend to avoid issues on headless systems
        # save the recorded Irec values for the reference amplitudes as a json file
        figsize = (12,6)
        plt.figure(figsize=figsize)
        font_size = 16
        title_size = 18
        plt.rcParams.update({'font.size': font_size})
        title = "Reference Irec values for the reference amplitudes at " + str(self.reference_frequency) + " Hz"
        plt.title(title, fontsize=title_size, fontweight='bold')

        amplitudes, irec_values = zip(*self.irec_vs_sweep_amplitudes)
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
        filename = f"{self.session_path}/reference_irec_values_{self.reference_frequency}Hz_{datetime_string}.png"
        dpi = 300 # save the plot with high resolution

        plt.savefig(filename, dpi=dpi)
        print(f"Saved reference Irec values for the reference amplitudes plot to {filename}")
        return 0


    # helper function to dump a dictionary to a json file
    def dump_dict_to_json(self, file, dict_data,  dict_name = None):
        
        if dict_name:
            file.write(f'"{dict_name}": {{\n')

            for index, (key, value) in enumerate(dict_data.items()):
                file.write(f'\t"{key}": {json.dumps(value)}')
                if index < len(dict_data)-1:
                    file.write(",\n")
                else:
                    file.write("\n")

            file.write("\t},\n")

        else:
            for key, value in dict_data.items():
                file.write(f'"{key}": {json.dumps(value)},\n')

    # function to save the recorded data
    def save_data(self):
        # get current date and time
        current_time = time.strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{self.session_path}/{self.filename}_{current_time}.json"

        print(f"Saving data to {filename}...")
        print("Data headers:")
        print(self.recorded_data_headers)   

        print("Data to be saved:")
        print(self.recorded_data_values)

        preliminaries = {
            "type": "dummy_data",
            "version": self.version,
            "header": f"{self.header}",
            "start_time": f"{self.start_time}",
            "end_time": f"{current_time}",
        }

        # TODO: also save settings
        with open(filename, 'w') as f:
            
            f.write("{\n") # start of json file

            # preliminaries
            self.dump_dict_to_json(file=f, dict_data=preliminaries)

            # nanonis settings
            self.dump_dict_to_json(file=f, dict_data=self.nanonis_settings, dict_name="nanonis_settings")
            # nanonis parameters
            self.dump_dict_to_json(file=f, dict_data=self.nanonis_parameters, dict_name="nanonis_parameters")
  
            # awg settings
            self.dump_dict_to_json(file=f, dict_data=self.awg_settings, dict_name="awg_settings")

            # data section
            f.write('"data": {\n')

            # data headers
            f.write('\t"channel names": ' + json.dumps(self.recorded_data_headers) + ",\n")
            
            # data
            f.write('\t"values": [\n') # start of values list
            for index in range(len(self.recorded_data_values)):
                data_values = self.recorded_data_values[index]
                line = "\t\t" + json.dumps(data_values) # convert list of values to json string and add indentation for better readability
               
                # add comma after each line except the last one
                if index < len(self.recorded_data_values)-1:
                    line += "," 

                f.write(line + "\n")

            # end values list and data section
            f.write("\t\t]\n")
            f.write("\t}\n")

            # end of json file
            f.write("}\n")


        logger.info(f"Data saved to {filename}.")

        return 0