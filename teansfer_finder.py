# module to run a measurement that allows computing the transfer function for some frequencies
class transfer_finder:
    def __init__(self, safe_voltage, safe_current):
        self.safe_voltage = safe_voltage
        self.safe_current = safe_voltage


    # if an error occurs, execute this command
    def escape_routine(self):
        """
        Function which is executed if an error occurs. 
        Sets the z-controller value to its default position and activates the controller.
    
        """


        print("Error occured!")