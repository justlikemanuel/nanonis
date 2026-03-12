class PIStromReglerMinimal:
    def __init__(self, Kp, Ti, dt, V_min=0.0, V_max=5.0):
        self.Kp = Kp
        self.Ti = Ti
        self.Ki = Kp / Ti   # automatisch aus Ti berechnet
        self.dt = dt

        self.V_min = V_min
        self.V_max = V_max

        self.integral = 0.0
        self.V_out = 0.0

    def reset(self):
        self.integral = 0.0
        self.V_out = 0.0

    def update(self, I_ref, I_meas):
        error = I_ref - I_meas

        # Integrator
        self.integral += error * self.dt

        # PI-Regelung
        V = self.Kp * error + self.Ki * self.integral

        # Anti-Windup über Begrenzung
        V = max(self.V_min, min(self.V_max, V))

        self.V_out = V
        return V