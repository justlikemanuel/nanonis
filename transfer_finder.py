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
                 open_nanonis_settings_gui = False,
                 height_averaging_time = 0.2,
                integration_time = 0.1,
                amplitude_guess_mode = "closest",
                old_measurement_file = None,
                old_compensation_amplitudes = None,
                atom_tracking_settings = None,
                atom_tracking_time = 0.4,
                atom_tracking_interval = 5,
                awg_reference = None, 
                awg_settling_time = 0.1,
                granularity_frequency = 1e5,
                lockin_frequency = 1e3,
                tuning_pgain = 0.5,
                tuning_integration_time_constant = 1.0,
                max_tune_iterations = 10,
                sweep_frequencies = None,
                reference_frequency = None,
                reference_STM_amplitude = 0.5,
                reference_transmission = 0.5,
                #max_sweep_amplitude = 0.5,
                #num_sweep_amplitudes = 5,
                data_channels = None,
                active_state_current = None,
                active_state_voltage = None,
                measurement_voltage = 0.5,
                irec_tolerance = 0.1e-13,
                header = "dummy header",
                filename = "transfer_function_measurement",
                communication_time = 1e-4, # TODO: find value!
                slew_rate = 0.1, # V/s, TODO: find value!
                 ):
        
        """
        Class to measure the transfer function for specified frequencies.

        Args:
            - nanonis_module: The NanonisModules object to interact with the Nanonis system.
            - open_nanonis_settings_gui: Flag which opens the parameter selection GUI if true.
            - height_averaging_time: The time to wait for the height to stabilize after switching off the z-controller.
            - integration_time: The integration time to use for recording the Irec value.
            - amplitude_guess_mode: The strategy to use for estimating the starting amplitude for the tuning process. Options are "known" and "half". "known" uses the recorded Irec values for the reference amplitudes to find the two reference frequencies that are closest to the desired frequency, and performs a linear interpolation to estimate the Irec value at the desired frequency, then finds the corresponding amplitude. "half" assumes 0.5 transmission to estimate the starting amplitude.
            - old compensation_amplitudes (list of tuples): The list of previously evaluated frequencies and their corresponding compensation amplitudes.
            - old_measurement_file (str): The path to the old measurement file from which to read parameters.
            - atom_tracking_settings: Dictionary containing the settings to use for atom tracking, e.g. a dictionary of parameters.
            - atom_tracking_time: The time to track the atom for after each measurement step.
            - atom_tracking_interval: The inerval after how many measurement steps should be executed again.
                                        If atom tracking shall not be used, set this to a value larger than the number of measurement steps.
            - awg_reference: The reference to the AWG to use for outputting the reference
            - awg_settling_time: The time to wait after changing the AWG settings before recording the Irec value, to allow the system to stabilize.
            - granularity_frequency: The granularity frequency to use for the AWG waveform generation. This should be chosen based on the desired frequencies to ensure that the generated waveforms meet the granularity requirements of the AWG.
            - lockin_frequency: The frequency of the lock-in amplifier.
            - tuning_pgain: The proportional gain to use for the tuning process.
            - tuning_integration_time_constant: The time constant for the integral action in the tuning process.
            - tolerance: The tolerance for the tuning process (relative to the current measured at reference amplitude and reference frequency).
            - sweep_frequencies: The frequencies for which the transfer function shall be measured.
            - reference_frequency: The frequency at which the reference Irec value shall be recorded.
            - reference_STM_amplitude: The amplitude at the STM (!!) for which the reference Irec value shall be recorded.
            - max_sweep_amplitude: The maximum amplitude to use for the sweep, to protect the tip and sample. The actual amplitudes used for the sweep will be generated based on this value and the number of amplitudes.
            - num_sweep_amplitudes: The number of amplitudes to use for the sweep between 0 and the maximum sweep amplitude.
            - data_channels: Labels of the channels in Nanonis that shall be logged.
            - use_active_state: TODO: check with Nicolaj again.
            - measurement_voltage: The voltage used for which the measurement shall be run.
            - irec_tolerance: The tolerance for the Irec value when comparing to the reference Irec value for the compensation amplitude tuning.
            - filename: The name of the file to save the data to.
            - header: The header to save in the data file, e.g. a description of the experiment and the settings used.
            - communication_time: The time to wait after each communication with the Nanonis system, to ensure that the system has time to process the command and update the values. This can help to prevent errors due to too fast communication. TODO: find value!
            - slew_rate: The maximum slew rate to use for the voltage changes, to protect the tip and sample. This can be used in the ramping functions to ensure that the voltage is changed in a way that does not exceed this slew rate.       
        """
                
        # dummy parameters (TODO: should be used with the constructor)
        self.max_allowed_amplitude = 1 # maximum allowed amplitude in Volts to protect the sample and tip, TODO: find better parameter for this, maybe based on the recorded Irec values for the reference amplitudes
        self.communication_time = communication_time # time to wait after each communication with the Nanonis system, to ensure that the system has time to process the command and update the values. This can help to prevent errors due to too fast communication. TODO: find value!
        
        # AWG parameters
        self.awg = awg_reference
        self.awg_settling_time = awg_settling_time
        self.granularity_frequency = granularity_frequency
        self.lockin_frequency = lockin_frequency

        # create integrator
        self.tuning_controller = PIController(Kp=tuning_pgain, Ti=tuning_integration_time_constant, 
                                              dt=0.1, V_min=0.1, V_max=self.max_allowed_amplitude)
        self.irec_tolerance = irec_tolerance


        # nanonis        
        self.nanonis_module = nanonis_module
        self.height_averaging_time = height_averaging_time
        self.integration_time = integration_time
        self.current_index = 1 # TODO: find from nanonis

        # get current nanonis settings
        x_pos, y_pos = self.nanonis_module.FolMe.XYPosGet(Wait_for_newest_data=True)
        self.initial_x_position_m = x_pos
        self.initial_y_position_m = y_pos
        self.initial_voltage = self.nanonis_module.Bias.Get()
        self.initial_current_A = self.nanonis_module.ZCtl.SetpntGet()        
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
        self.slew_rate = slew_rate # maximum slew rate to use for the voltage changes in the escape routine, to protect the tip and sample

        # Sweep parameters:
        self.sweep_frequencies = sweep_frequencies
        #self.max_sweep_amplitude = max_sweep_amplitude
        #self.num_amplitudes = num_sweep_amplitudes
        #self.irec_vs_sweep_amplitudes = []

        # generate reference amplitudes
        #step = self.max_sweep_amplitude / self.num_amplitudes
        #self.sweep_amplitudes = [step * (i+1) for i in range(self.num_amplitudes)]
        #logger.info(f"Generated reference amplitudes: {self.sweep_amplitudes} V")

        self.measurement_voltage = measurement_voltage
        self.reference_transmission = reference_transmission
        self.reference_STM_amplitude = reference_STM_amplitude
        self.reference_amplitude = self.reference_STM_amplitude / self.reference_transmission
        self.reference_frequency = reference_frequency
        self.old_compensation_amplitudes = old_compensation_amplitudes
    
   
        # atom tracking
        self.atom_tracking_settings = atom_tracking_settings
        self.nanonis_module.ATrack.PropsSet(**self.atom_tracking_settings)
        self.atom_tracking_time = atom_tracking_time
        self.atom_tracking_interval = atom_tracking_interval

        # logging parameters
        self.start_time = time.strftime("%Y-%m-%d_%H-%M-%S")
        self.session_path = self.get_session_path()
        self.filename = filename
        self.version = [0, 0, 1]
        self.header = header

        # active state parameters
        self.active_state_current = active_state_current
        self.active_state_voltage = active_state_voltage

        self.recorded_data_headers = [
            "frequency (Hz)",
            "compensation_amplitude (V)",
        ]
        self.nanonis_channels = data_channels

        self.recorded_data_headers.extend(self.nanonis_channels)
        self.recorded_data_values = [] # list of tuples (frequency, tuned_amplitude, current, bias, z_controller_setpoint,...)

        # compensation parameters
        self.amplitude_guess_mode = amplitude_guess_mode
        self.reference_i_rec = None # current value at the reference amplitude

        # keep track of current desired paramters for the escape routine
        self.current_desired_voltage = self.initial_voltage
        self.current_desired_current = self.initial_current_A

        ################# For Logging ##################
        # store all parameters in a dictionary for logging
        self.nanonis_settings = {
                    "integration_time": integration_time,
                    "atom_tracking_settings": atom_tracking_settings,
                    "atom_tracking_time": atom_tracking_time,
                    "atom_tracking_interval": atom_tracking_interval,
                    "data_channels": data_channels,
                    "measurement_voltage": measurement_voltage,
                    "active_state_current": active_state_current,
                    "active_state_voltage": active_state_voltage,   
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
                    "granularity_frequency": granularity_frequency,
                    "lockin_frequency": lockin_frequency,
                    "sweep_frequencies": sweep_frequencies,
                    "reference_frequency": reference_frequency,
                    "reference_STM_amplitude": reference_STM_amplitude,
                    "reference_transmission": reference_transmission,
                    #"max_sweep_amplitude": max_sweep_amplitude,
                    #"num_sweep_amplitudes": num_sweep_amplitudes,
                }
        
        self.tuning_settings = {
                    "tuning_pgain": tuning_pgain,
                    "tuning_integration_time_constant": tuning_integration_time_constant,
                    "irec_tolerance": irec_tolerance,
                    "max_tune_iterations": max_tune_iterations,
                }

        print(f"Session path: {self.session_path}")        
        
        # get the nanonis parameters (to put into the logged file)
        meas = MeasurementBase(self.nanonis_module)
        #self.nanonis_parameters = meas.nanonisSettingsGet(open_nanonis_settings_gui)
        self.nanonis_parameters = {"Dummy_parameter": 0} # the nanonis function above just broke...


    # Function to check validity of the settings and parameters
    def check_settings_validity(self):
        """
        Function to check the validity of the settings and parameters. This can include checks such as:
            - Are there any nanonis settings that may cause failure?
            - Are the specified frequencies within the range of the AWG and lock-in amplifier?
            - Is the specified integration time reasonable for the measurement?
            - Are the specified amplitudes within the limits of the AWG and safe for the sample and tip?
            - Are the specified Nanonis channels correctly configured in Nanonis?
            - etc.
        If any check fails, an error should be raised and the measurement should not be started.
        """
        # TODO: implement function :)
        pass


    ###################################################
    ############## Positioning functions ##############
    ###################################################

    # function to (safely) maneeuver to a desired state
    def maneeuver_to_state(self, desired_voltage, desired_current):
        """
        Method to achieve a state given by a desired voltage and desired current.

        Args:
            - desired_voltage: The desired bias voltage in Volts.
            - desired_current: The desired current setpoint in Amperes.

        
        """
        try:
            # check if z-controller is on
            if self.nanonis_module.ZCtl.OnOffGet()==1:
                # Z-controller is ON

                if self.check_bias_polarity(desired_voltage):

                    #polarity matches
                    # TODO: Think about "current" (aktueller) vs "current" (Strom) -> better naming :/
                    current_setpoint = self.nanonis_module.ZCtl.SetpntGet()

                    if current_setpoint < desired_current:
                        self.nanonis_module.ZCtl.SetpntSet(desired_current)
                else:
                    # polarity does not match, turn off z-controller and ramp to desired current
                    self.turn_off_z_controller_and_wait()
                    self.ramp_current_z_off(desired_voltage, desired_current)        
            else:
                # case Z-controller is OFF and/or polarity does not match
                self.ramp_current_z_off(desired_voltage, desired_current)

            self.ramp_bias(desired_voltage)

            # set target current
            self.nanonis_module.ZCtl.SetpntSet(self.current_desired_current)

            return 0
        
        except Exception as e:
            print(f"Error while maneeuvering to state with voltage {desired_voltage} V and current {desired_current} A: {e}. Executing escape routine.")
            print(f"Error in line {e.__traceback__.tb_lineno}")
            self.escape_routine()

    def ramp_current_z_off(self, desired_voltage, desired_current):
        
        """
        Function to handle the case when the z-controller is off.
        First, the voltage is ramped to achieve the current. 
        Then, the z-controller is re-activated.

        Args:
            - desired_voltage: The desired bias voltage in Volts.
            - desired_current: The desired current setpoint in Amperes.
        """

        try:
            if not self.check_bias_polarity(desired_voltage):
                self.bias_to_zero()

            # ensure current is set to desired current and tune voltage
            self.nanonis_module.ZCtl.SetpntSet(desired_current)
            self.tune_voltage_to_current_setpoint(desired_current)

            self.nanonis_module.ZCtl.OnOffSet(1) # turn on z-controller

            return 0
        
        except Exception as e:
            print(f"Error while handling left branch: {e}. Executing escape routine.")
            print(f"Error in line {e.__traceback__.tb_lineno}")
            self.escape_routine()


    # helper function to iteratively set the bias to a desired voltage
    def ramp_bias(self, new_voltage, total_time = 0.05):
        """
        Function that tunes the bias iteratively to a desired voltage
        Args:
            - new_voltage: The desired voltage setpoint in Volts.
            - time: The total duration of the function in seconds.

        """
        current_voltage = self.nanonis_module.Bias.Get()
        diff_voltage = new_voltage - current_voltage

        # ensure slew rate is not exceeded
        time_from_slew_rate = abs(diff_voltage) / self.slew_rate
        total_time = max(total_time, time_from_slew_rate)

        # compute voltage increment
        num_steps = int(np.ceil(diff_voltage / (self.slew_rate * total_time)))       
        if num_steps == 0:
            self.nanonis_module.Bias.Set(new_voltage)
            return 0
        
        step_voltage = diff_voltage / num_steps
        time_per_step = total_time / num_steps
        additional_waiting_time = time_per_step - self.communication_time

        try:
            for _ in range(num_steps):
                current_voltage += step_voltage
                new_voltage = current_voltage - step_voltage
                self.nanonis_module.Bias.Set(new_voltage)

                if (additional_waiting_time > 0):
                    time.sleep(additional_waiting_time)

            return 0
        
        except Exception as e:
            print(f"Error while ramping bias: {e}. Executing escape routine.")
            print(f"Error in line {e.__traceback__.tb_lineno}")
            self.escape_routine()

    # helper function to achieve a desired current while controler is off
    def ramp_to_current(self, desired_current, desired_voltage):
        """
        Function to achieve a desired current setpoint while the z-controller is off, by iteratively adjusting the bias voltage and measuring the current until the desired current is reached within a specified tolerance.

        Args:
            - desired_current: The desired current setpoint in Amperes.
            - desired_voltage: The desired voltage setpoint in Volts.
        """

        # check polarity and set voltage to 0 if polarity does not match
        if not self.check_bias_polarity(desired_voltage):
            self.ramp_bias(0, 0.05) # TODO: find better parameters for the time and voltage step in this function
        # iteratively adjust voltage until achieve current is above the desired current setpoint
        voltage_step = self.slew_rate * self.communication_time

        while self.get_irec() < desired_current:            
            new_voltage = self.nanonis_module.Bias.Get() + voltage_step # TODO: find better strategy for voltage adjustment
            self.nanonis_module.Bias.Set(new_voltage)

        return 0
    
    # helper function to check polarity of measured and desired voltage
    def check_bias_polarity(self, desired_voltage):
        measured_voltage = self.nanonis_module.Bias.Get()
        polarity_measured = np.sign(measured_voltage)
        polarity_desired = np.sign(desired_voltage)

        if polarity_measured != polarity_desired:
            return False
        else:
            return True
        
    
    # helper function to tune the voltage to a desired current setpoint
    def tune_voltage_to_current_setpoint(self, current_setpoint):
        """
        Function to tune the bias voltage to reach a desired current setpoint using a PI controller.

        Args:
            - current_setpoint: The desired current setpoint in Amperes.
        """
        # find polarity of the desired voltage
        polarity= np.sign(self.current_desired_voltage)
        voltage_step = self.slew_rate * self.communication_time * polarity
        set_voltage = self.nanonis_module.Bias.Get()

        while  (current_setpoint - self.get_irec()) > 0:
            set_voltage += voltage_step
            self.nanonis_module.Bias.Set(set_voltage)

        return 0

    
    # helper function to measure the difference between the current voltage and the desired voltage
    def measure_voltage_difference(self, desired_voltage):
        measured_voltage = self.nanonis_module.Bias.Get()
        diff_voltage = measured_voltage - desired_voltage

        return diff_voltage
    
    # function to turn off the z-controller and wait for the height to stabilize
    def turn_off_z_controller_and_wait(self):
        """
        Function to turn off the z-controller and wait for the height to stabilize.
        """
        try:
            self.nanonis_module.ZCtl.OnOffSet(0) # turn off z-controller
        
            while self.nanonis_module.ZCtl.OnOffGet() == 1:
                time.sleep(0.01)

            return 0
        
        except Exception as e:
            print(f"Error while turning off z-controller: {e}. Executing escape routine.")
            self.escape_routine()

    # return to default state
    def return_to_starting_state(self):

        """
        Sets all parameters as before the measurement.
        TODO: check what exactly is needed
        """

        try:
            # TODO: check with Nicolaj for best sequence of operations
    
            # set off delay to initial value
            self.nanonis_module.ZCtl.SwitchOffDelaySet(self.initial_z_controler_switch_off_delay_s)
            self.maneeuver_to_state(self.initial_voltage, self.initial_current_A)

            # restore atom tracking settings
            self.nanonis_module.ATrack.PropsSet(*self.initial_tracking_settings.values())

            return 0

        except Exception as e:
            print(f"Error while returning to starting state: {e}. Executing escape routine.")
            self.escape_routine()      
 
    # if an error occurs, execute this command
    def escape_routine(self):
        """
        Function which is executed if an error occurs. 
        Sets the z-controller value to its default position and activates the controller.
    
        """
        ### REMOVE AFTER TESTING!!! ###
        print("Error occured!")
        exit(1)
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
            self.maneeuver_to_state(self.initial_voltage, self.initial_current_A)

        except Exception as e:

            # debugging ONLY!!!! TODO: Remove this comment
            exit(1)
            print(f"Error during escape routine: {e}. Retrying...")
            self.escape_routine()

        # TODO: what shall happen after recovering?
        self.save_data()
        exit(1)


    # function to go to measurement position and height
    def prepare_measurement(self):
        """
        Function to prepare the measurement.
        """
        try:
            # TODO: verify sequence
        
            if self.atom_tracking_interval < len(self.sweep_frequencies):
                print("Starting atom tracking.")
                self.track_atom()
                
            # ensure the z_off_delay is set to the desired value
            self.nanonis_module.ZCtl.SwitchOffDelaySet(self.height_averaging_time)
            print("Atom tracking finished.")

            # DEBUG ONLY:
            # move to 2.5 V and 10 pA
            voltage = 2.5
            amps = 20e-12
            self.maneeuver_to_state(voltage, amps)

            print("Moved to 2.5 V and 10 pA for testing purposes. Remove this after testing!!!")
            time.sleep(5)
            # move to active state position, if specified
            if self.active_state_voltage is not None and self.active_state_current is not None:
                # update current desired parameters for escape routine
                print(f"Moving to active state given by voltage {self.active_state_voltage} V and current {self.active_state_current} A.")
                self.current_desired_voltage = self.active_state_voltage
                self.current_desired_current = self.active_state_current
                self.maneeuver_to_state(self.active_state_voltage, self.active_state_current)

                print("moved to active state.")
                time.sleep(10)

            # turn off the z-controller to allow for height averaging
            self.turn_off_z_controller_and_wait()

            # set measurement voltage
            self.ramp_bias(self.measurement_voltage)

            return 0
        
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
            
            return 0
        
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
    def get_irec(self, integration_time = None):
        """
        Function to get the current Irec value.

        Returns
        Irec (float): The recorded Irec value in Amperes.
        """
        
        # no integration time -> single shot
        if integration_time is None:
            return self.nanonis_module.Sig.ValGet(signal_index=self.current_index, wait_for_newest_data=True)
        
        readout = self.nanonis_module.Sig.MeasSig(sig_names = ["Current (A)"], averaging_time=integration_time) # returns dictionary with signal names as keys and measured values as values
        
        # if current not in the returned dictionary, raise error
        if "Current (A)" not in readout:
            raise ValueError("Current (A) not found in the measured signals. Check if the channel name is correct and if the signal is properly configured in Nanonis.")
        
        return readout["Current (A)"]

    # record Irec at reference amplitude and frequency
    def record_reference_irec(self):
        """
        Function to record the Irec value at the reference amplitude and frequency. 
        This value isused as a reference for the tuning process.
        """
        try:
            # configure AWG to output the reference signal
            self.awg.configure_continuous_sine_wave(frequency=self.reference_frequency, 
                                                    granularity_frequency=self.granularity_frequency,
                                                    lockin_frequency=self.lockin_frequency, 
                                                    starting_amplitude=self.reference_amplitude)
            # activate the output of the AWG and measure at reference amplitude
            self.awg.start_playing()
            time.sleep(self.awg_settling_time)
            self.reference_i_rec = self.get_irec(integration_time=self.integration_time)
            self.awg.stop_playing()
            
            return 0
            
        except Exception as e:
            print(f"Error while recording reference Irec value: {e}. Executing escape routine.")
            print(f"Error in line {e.__traceback__.tb_lineno}")
            self.escape_routine()
          
    # record Irec vs the reference amplitudes
        """ def record_irec_for_references(self):
        
        """
        #Function to record the Irec value for the reference amplitude, and all sweep amplitudes at the reference frequency.
        """
        try:
            # configure AWG to output the reference signal
            self.awg.configure_continuous_sine_wave(frequency=self.reference_frequency, lockin_frequency=self.lockin_frequency, starting_amplitude=self.reference_amplitude)
            
            # activate the output of the AWG and measure at reference amplitude
            self.awg.start_playing()
            time.sleep(self.awg_settling_time)
            self.reference_i_rec = self.get_irec(integration_time=self.integration_time)
            print("AWG ON")

            for amplitude in self.sweep_amplitudes:

                matched_amplitude = self.awg.update_continuous_sine_wave_amplitude(new_amplitude=amplitude)
                time.sleep(self.awg_settling_time)
                irec = self.get_irec(integration_time=self.integration_time)

                # add data to list
                self.irec_vs_sweep_amplitudes.append((matched_amplitude, irec))

            # stop the AWG output
            self.awg.stop_playing()
            print("AWG OFF")

        except Exception as e:
            print(f"Error while recording Irec for reference amplitudes: {e}. Executing escape routine.")
            self.escape_routine
            """
   
    # function to tune awg amplitude for a specific frequency to match reference irec
    def tune_awg_amplitude_for_frequency(self, frequency, starting_amplitude = 0.1,
                                         tolerance=0.01, max_iterations=2):
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
                                                    granularity_frequency=self.granularity_frequency,
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

            i_rec = self.get_irec(integration_time=self.integration_time)

            while ((i_rec > upper_bound_irec or i_rec < lower_bound_irec) 
                    and iteration < max_iterations):
                
                tuned_amplitude = self.tuning_controller.update(I_ref=self.reference_i_rec, I_meas=i_rec)
                
                self.awg.update_continuous_sine_wave_amplitude(new_amplitude=tuned_amplitude)
                time.sleep(self.awg_settling_time)
                i_rec = self.get_irec(integration_time=self.integration_time)
                iteration += 1
                
                """                
                # TODO: find better tuning strategy, e.g. proportional control based on the difference between recorded Irec and default Irec, instead of just increasing or decreasing by a fixed percentage
                # or linear sweep
                # CAUTION: too much current rips the sample/tip
                # Nicolaj mentioned something as the momentum method
                if self.get_irec(integration_time=self.integration_time) > upper_bound_irec:
                    tuned_amplitude *= 0.9 # decrease amplitude by 10%
                else:
                    tuned_amplitude *= 1.1 # increase amplitude by 10%

                # the awg function will clip to the resolution and return the applied value
                tuned_amplitude = self.awg.update_continuous_sine_wave_amplitude(new_amplitude=tuned_amplitude)
                time.sleep(self.integration_time)
                iteration +=1 
                """

            # turn off the AWG output
            self.awg.stop_playing()

            # log the result
            print(f"Tuned amplitude for frequency {frequency} Hz: {tuned_amplitude} V")
            
            return tuned_amplitude
    
        except Exception as e:
            print(f"Error while tuning AWG amplitude for frequency {frequency} Hz: {e} in line {e.__traceback__.tb_lineno}. Executing escape routine.")
            self.escape_routine()
    
    # function to estimate the starting amplitude for the tuning process
    def estimate_starting_amplitude_for_frequency(self, frequency, mode):
        """
        Function to estimate the starting amplitude for the tuning process based on the recorded Irec values for the reference amplitudes. The function can use different strategies for the estimation, e.g. linear interpolation between the two reference amplitudes that are closest to the desired frequency, or using the transfer function at the default frequency to extrapolate the expected Irec value at the desired frequency and then find the corresponding amplitude.

        Args:
            - frequency (float): The frequency for which to estimate the starting amplitude.
            - mode (str): The strategy to use for the estimation. 
                Options are:
                - "known": use the known transfer function value for the desired frequency
                - "closest": use the recorded Irec values for the reference amplitudes to find the reference frequency that is closest to the desired frequency, and perform a linear interpolation to estimate the Irec value at the desired frequency, then find the corresponding amplitude.
                - "half": assume 0.5 transmission 
                - "fixed": use a fixed value for the starting amplitude, currently 0.2 V
                - idea: 50% transmission -> problem: what is the desired amplitude at the output of the channel?
        """
        if mode not in ["known", 
                        "closest", 
                        "half", 
                        "fixed"]:
            raise ValueError(f"Invalid mode for estimating starting amplitude: {mode}. Valid options are 'known', 'closest', 'half', 'fixed'.")
        
        starting_amplitude = 0
        if mode == "known":
            # could be modified to not only use the same but a close old frequency

            # check if the array exists
            if self.old_compensation_amplitudes is None:
                raise ValueError("No old compensation amplitudes provided for 'known' mode. Please provide a list of previously evaluated frequencies and their corresponding compensation amplitudes.")
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

        if mode == "half":
            starting_amplitude = self.reference_STM_amplitude * 2 # assume 50% transmission
            
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
            values = self.nanonis_module.Sig.MeasSig(sig_names = self.nanonis_channels, averaging_time=self.integration_time) # returns dictionary with signal names as keys and measured values as values

            data_list = [frequency, tuned_amplitude]
            
            for channel in self.nanonis_channels:
                # if channel not in the returned dictionary, raise error
                if channel not in values:
                    raise ValueError(f"{channel} not found in the measured signals. Check if the channel name is correct and if the signal is properly configured in Nanonis.")
                data_list.append(values[channel])

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

    """# function to save the reference Irec values for the reference amplitudes
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
    # CAUTION: Shall not be used, replaced by json_dump of individual dictionaries
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
"""
    # function to save the recorded data
    def save_data(self):
        # get current date and time
        current_time = time.strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"{self.session_path}/{self.filename}_{current_time}.json"

        #print(f"Saving data to {filename}...")
        #print("Data headers:")
        #print(self.recorded_data_headers)   

        #print("Data to be saved:")
        #print(self.recorded_data_values)

        
        # merge all dictionaries into one for easier dumping
        data_to_dump = {
            "type": "dummy_data",
            "version": self.version,
            "header": f"{self.header}",
            "start_time": f"{self.start_time}",
            "end_time": f"{current_time}",
        }

        data_to_dump["nanonis_measurement_settings"] = self.nanonis_settings
        data_to_dump["nanonis_pre_settings"] = self.nanonis_parameters
        data_to_dump["awg_settings"] = self.awg_settings
        data_to_dump["tuning_settings"] = self.tuning_settings
        data_to_dump["reference_i_rec"] = self.reference_i_rec

        # TODO: also save settings
        with open(filename, 'w') as f:

  

            f.write(json.dumps(data_to_dump, indent=4)[:-1]+",\n") # remove the last closing curly brace to add the data

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
    

    # function to read all parameters from an old logging file
    def read_parameters_from_old_measurement(self, filepath):
        """
        Function to read all parameters from an old logging file. This can be used to read the settings and parameters from a previous measurement and use them for the current measurement, e.g. for the tuning process or for the estimation of the starting amplitude for the tuning process.
        Args:
            - filepath (str): The path to the old logging file.
        Returns:
            - parameters (dict): A dictionary containing all the parameters read from the old logging file.
        """
        # TODO: Implement missing parameters (the dictionaries are only for reference as of now)
        with open(filepath, "r") as f:
            data = json.load(f)

        parameters = {
            "nanonis_measurement_settings": data.get("nanonis_measurement_settings", {}),
            "nanonis_pre_settings": data.get("nanonis_pre_settings", {}),
            "awg_settings": data.get("awg_settings", {}),
            "tuning_settings": data.get("tuning_settings", {}),
        }

        return parameters