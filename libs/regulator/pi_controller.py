class PIController:
    def __init__(self, Kp, Ti, dt, V_min=0.0, V_max=1.0):

        """
        PI-Controller for current regulation in an STM system.
        
        Parameters
            - Kp (float): Proportional gain
            - Ti (float): Integrator time constant in seconds
            - dt (float): Time step for the controller update in seconds
            - V_min (float): Minimum output voltage (default: 0.0 V)
            - V_max (float): Maximum output voltage (default: 1.0 V)
        """

        self.Kp = Kp
        self.Ti = Ti
        self.Ki = Kp / Ti
        self.dt = dt
        self.last_error = 0.0 # use trapezoidal rule for integration (error at t and t-1)

        self.V_min = V_min
        self.V_max = V_max

        self.integral = 0.0
        self.V_out = 0.0

    def reset(self):
        self.integral = 0.0
        self.V_out = 0.0

    def update(self, I_ref, I_meas):
        """
        Update the PI controller.
        
        Parameters
            - I_ref (float): Reference current
            - I_meas (float): Measured current

        Returns
            - float: Output voltage (for the next time step)
        """
        error = I_ref - I_meas

        # integration
        self.integral += self. dt * (error + self.last_error) / 2.0  # trapezoidal rule
        self.last_error = error
        V = self.Kp * error + self.Ki * self.integral

        # Anti-Windup 
        V = max(self.V_min, min(self.V_max, V))

        self.V_out = V
        return V