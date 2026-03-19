import json
import random

header = {
    "name": "test",
    "date": "2024-06-30",
}

n_samples = 5

with open("data.json", "w") as f:
    f.write("{\n")

    # Header schreiben (ohne äußere Klammern)
    header_str = json.dumps(header, indent=2)[1:-1]
    f.write(header_str + ",\n")

    # Samples
    f.write('  "samples": [\n')

    for i in range(n_samples):
        sample = [random.randint(0,100) for _ in range(5)]
        line = json.dumps(sample)  # bleibt in EINER Zeile

        if i < n_samples - 1:
            line += ","

        f.write("    " + line + "\n")

    f.write("  ]\n")
    f.write("}\n")