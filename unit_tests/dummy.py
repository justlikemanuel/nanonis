dummydict = {
    "nanonis": {
        "i_rec_integration_time_s": 0.1,
        "atom_tracking_settings": {
            "Igain": 5.7e-10,
            "Frequency": 10.0,
        },
    },
    "AWG": {
        "awg_settling_time": 0.1,
        "lockin_frequency": 1000.0,
        "sweep_frequencies": [1000000.0, 2000000.0, 3000000.0, 4000000.0, 5000000.0],
        "reference_frequency": 10000.0,
        "max_reference_amplitude": 0.5,
        "num_amplitudes": 5
    }
}

# save this dict to a json file
import json
with open("dummy_data.json", "w") as f:
    json.dump(dummydict, f, indent=4)