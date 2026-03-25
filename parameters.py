# define classes to hold parameters for the awg, the nanonis and others

# general parameter class
class Parameter:
    def __init__(self, param, max, min, name=None):
        self.param = param
        self.max = max
        self.min = min
        self.name = name if name is not None else str(param)

    def check_bounds(self, value):
        if value < self.min or value > self.max:
            raise ValueError(f"Value {value} for parameter {self.name} is out of bounds ({self.min}, {self.max})")
        return True
    
# parent-class for a list of parameters
class ParameterList:
    def __init__(self, params):
        self.params = params

    def load_from_dict(self, param_dict):
        for param in self.params:
            if param.name in param_dict:
                value = param_dict[param.name]
                param.check_bounds(value)
                param.param = value
            else:
                raise KeyError(f"Parameter {param.name} not found in input dictionary")
            
    def return_dict(self):
        return {param.name: param.param for param in self.params}

# AWG parameters
class AWGParameters(ParameterList):
    def __init__(self, frequency, amplitude, phase):
        params = [
            Parameter(frequency, max=5e9, min=1e3, name="frequency"),
            Parameter(amplitude, max=2.0, min=0.01, name="amplitude"),
            Parameter(phase, max=360.0, min=0.0, name="phase")
        ]
        super().__init__(params)

# Nanonis measurement parameters
class NanonisMeasurementParameters(ParameterList):
    def __init__(self, desired_current, frequency, amplitude, phase):
        params = [
            Parameter(desired_current, max=1e-9, min=1e-12, name="desired_current"),
            Parameter(frequency, max=5e9, min=1e3, name="frequency"),
            Parameter(amplitude, max=2.0, min=0.01, name="amplitude"),
            Parameter(phase, max=360.0, min=0.0, name="phase")
        ]
        super().__init__(params)

# Nanonis parameters before measurement
class NanonisPreMeasurementParameters(ParameterList):
    def __init__(self, desired_current):
        params = [
            Parameter(desired_current, max=1e-9, min=1e-12, name="desired_current")
        ]
        super().__init__(params)

