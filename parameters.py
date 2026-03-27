# define classes to hold parameters for the awg, the nanonis and others

# general parameter class
class Parameter:
    def __init__(self, param, max=None, min=None, name=None, use_in_dict=True):
        self.param = param
        self.max = max
        self.min = min
        self.name = name if name is not None else str(param)
        self.use_in_dict = use_in_dict

    def check_bounds(self, value):
        if self.min is None and self.max is None:
            return True
        
        if value < self.min or value > self.max:
            raise ValueError(f"Value {value} for parameter {self.name} is out of bounds ({self.min}, {self.max})")
        return True
    
# parent-class for a list of parameters
class ParameterList:
    def __init__(self, params):
        self.params = params

    def load_from_dict(self, param_dict):
        for param in self.params:
            if not param.use_in_dict:
                continue
            if param.name in param_dict:
                value = param_dict[param.name]
                param.check_bounds(value)
                param.param = value
            else:
                raise KeyError(f"Parameter {param.name} not found in input dictionary")
            
    def return_dict(self):
        return {param.name: param.param for param in self.params}

# general parameters
class GeneralParameters(ParameterList):
    def __init__(self, 
                 version,                 
                 new_logging_file,
                 header,
                 old_logging_file,
    ):
        params = [
            Parameter(version, name="version"),
            Parameter(new_logging_file, name="new_logging_file", use_in_dict=False),
            Parameter(old_logging_file, name="old_logging_file"),
            Parameter(header, name="header"),
            Parameter(old_logging_file, name="old_logging_file")
        ]
        super().__init__(params)
# AWG parameters
class AWGParameters(ParameterList):
    def __init__(self, reference, 
                 settling_time, 
                 granularity_frequency, 
                 lockin_frequency,
                 ):
        params = [
            Parameter(reference, name="AWG reference", use_in_dict=False),
            Parameter(settling_time, name="settling_time"),
            Parameter(granularity_frequency, name="granularity_frequency"),
            Parameter(lockin_frequency, name="lockin_frequency")
        ]

        super().__init__(params)

# Nanonis measurement parameters
class NanonisMeasurementParameters(ParameterList):
    def __init__(self,
                 nanonis_module,                 
                 #use_active_state,
                 active_state_current,
                 active_state_voltage,
                 integration_time,
                 reference_transmission,
                 reference_STM_amplitude,
                 reference_frequency,
                 sweep_frequencies,
                 atom_tracking_time,
                 atom_tracking_interval,
                 # start_frequency,
                 # stop_frequency,
                                 
                 
                 ):
        params = [
            Parameter(nanonis_module, name="nanonis_module", use_in_dict=False),
            #Parameter(use_active_state, name="use_active_state"),
            Parameter(active_state_current, name="active_state_current"),
            Parameter(active_state_voltage, name="active_state_voltage"),
            Parameter(integration_time, max=100.0, min=1e-3, name="integration_time"),
            Parameter(reference_transmission, max=1.0, min=0.0, name="reference_transmission"),
            Parameter(reference_STM_amplitude, max=1.0, min=0.0, name="reference_STM_amplitude"),
            Parameter(reference_frequency, max=1e9, min=1e-3, name="reference_frequency"),
            Parameter(sweep_frequencies, name="sweep_frequencies"),
            Parameter(atom_tracking_time, max=100.0, min=0.0, name="atom_tracking_time"),
            Parameter(atom_tracking_interval, max=100.0, min=0.0, name="atom_tracking_interval"),
            #Parameter(start_frequency, max=1e9, min=1e-3, name="start_frequency"),
            #Parameter(stop_frequency, max=1e9, min=1e-3, name="stop_frequency")

        ]
        super().__init__(params)

# Nanonis parameters before measurement
class NanonisPreMeasurementParameters(ParameterList):
    def __init__(self, desired_current):
        params = [
            Parameter(desired_current, max=1e-9, min=1e-12, name="desired_current")
        ]
        super().__init__(params)

# Regulator parameters
class AmplitudeTuningParameters(ParameterList):
    def __init__(self, 
                 amplitude_guess_mode,
                 i_rec_tolerance,
                 Kp,
                 Ti,
                 dt, 
                 V_min,
                 V_max,
                 ):
        params = [
        ]
        super().__init__(params)


# class to hold measurement parameters
class MeasurementValues:
    def __init__(self, labels):
        self.labels = labels
        self.values = []

    def add_values(self, new_values):
        if len(new_values) != len(self.labels):
            raise ValueError(f"Number of new values {len(new_values)} does not match number of labels {len(self.labels)}")
        self.values.append(new_values)

    def get_values(self):
        return self.values
    
    def get_labels(self):
        return self.labels