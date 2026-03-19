# go one folder up
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from libs.AWG_M8195A_interface.M8195A_transfer import M8195A_transfer
from libs.pyNanonisMeasurements.nanonisTCP.nanonisTCP import nanonisTCP
# create an instance of the M8195A_transfer class
awg = M8195A_transfer()

############## generate allowed frequencies between 1 MHz and 10 MHz with 10 frequencies in total
# "Simple" method:
allowed_freqs, _ = awg.compute_allowed_frequencies_simple(start_freq_Hz=1e6, stop_freq_Hz=10e6, num_freqs=10, base_frequency_Hz=1e5)

#"Fancy" method 
# CAUTION: Suffers from rounding errors
#allowed_freqs, corresponding_num_periods = awg.compute_allowed_frequencies_lcm(start_freq_Hz=1e6, stop_freq_Hz=10e6, num_freqs=10)
print("Allowed frequencies (Hz) for continuous sine wave generation:")

# Dumb method (checking integer wise the number of periods)
corresponding_num_periods = []
for i, freq in enumerate(allowed_freqs):
    num_periods = awg.find_num_periods_for_granularity(freq)
    corresponding_num_periods.append(num_periods)
    timespan = num_periods / freq
    print(f"Frequency: {freq:.2f} Hz, Number of periods: {num_periods}, Timespan: {timespan*1e6:.2f} microseconds")

# test the frequency verification function
results = awg.test_frequencies(allowed_freqs, num_periods=corresponding_num_periods)
num_failed = 0
print("\nFrequency testing results:")
for freq, (meets_requirement, samples_per_period, total_samples, remainder) in results.items():
    if not meets_requirement:
        num_failed += 1
    num_periods = corresponding_num_periods[allowed_freqs.index(freq)]
    print(f"Frequency: {freq:.2f} Hz, num_periods: {num_periods}, samples_per_period: {samples_per_period:.2f}, total_samples: {total_samples}, remainder: {remainder:.2f}, Meets granularity requirement: {meets_requirement}")
print(f"Number of failed frequencies: {num_failed} out of {len(allowed_freqs)}")