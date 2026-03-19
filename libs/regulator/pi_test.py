import numpy as np
import matplotlib.pyplot as plt
from pi_controller import PIStromRegler

# ==========================================
# HF-Strecke (gedämpftes Resonanzsystem + STM)
# ==========================================
class HFStreckeSTM:
    def __init__(self, K=0.8, omega0=40.0, zeta=0.2, dt=0.001):
        self.K = K
        self.omega0 = omega0
        self.zeta = zeta
        self.dt = dt
        self.I = 0.0
        self.I_dot = 0.0

    def update(self, V_in):
        # 2. Ordnung: d²I/dt² + 2ζω₀ dI/dt + ω₀² I = K ω₀² V_in
        I_ddot = self.K * self.omega0**2 * V_in - 2 * self.zeta * self.omega0 * self.I_dot - self.omega0**2 * self.I
        self.I_dot += I_ddot * self.dt
        self.I += self.I_dot * self.dt

        # STM-HF-Komponente
        stm = 0.05 * np.sin(2*np.pi*50*self.dt)  # 50 Hz HF
        return self.I + stm

# ==========================================
# Simulationseinstellungen
# ==========================================
dt = 0.001
sim_time = 0.5
steps = int(sim_time / dt)
I_ref = 1.0
Kp = 2.0

# Verschiedene Integrator-Zeitkonstanten Ti
Ti_values = np.linspace(0.5, 5.0, 5)  # von 0.5s bis 5s

plt.figure(figsize=(10,5))

for Ti in Ti_values:
    controller = PIStromRegler(Kp=Kp, Ti=Ti, dt=dt, V_min=0.0, V_max=5.0)
    plant = HFStreckeSTM(K=0.8, omega0=40.0, zeta=0.2, dt=dt)

    I_values = []

    for step in range(steps):
        I_meas = plant.I
        V_in = controller.update(I_ref, I_meas)
        I_out = plant.update(V_in)
        I_values.append(I_out)

    plt.plot(np.linspace(0, sim_time, steps), I_values, label=f"Ti={Ti}s")

# Groundtruth
plt.axhline(I_ref, linestyle="--", color="gray", label="Groundtruth")

plt.xlabel("Zeit [s]")
plt.ylabel("Strom [A]")
plt.title("HF-Strecke: Einfluss der Integrator-Zeitkonstante Ti")
plt.legend()
plt.grid(True)
plt.show()