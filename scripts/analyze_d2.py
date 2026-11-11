
import csv
import statistics
import matplotlib.pyplot as plt

input_file = "results/part_d/d2_raw.csv"

groups = {}

with open(input_file, newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:
        t = int(row["threads"])
        groups.setdefault(t, []).append(float(row["seconds"]))

baseline = statistics.median(groups[1])

# Read empirical serial fraction from D1
d1_fractions = []

with open("results/part_d/d1_summary.csv", newline="") as f:
    for row in csv.DictReader(f):
        value = row["karp_flatt"]
        if value:
            d1_fractions.append(float(value))

serial_fraction = statistics.mean(d1_fractions)

summary = []

for t in sorted(groups):

    times = groups[t]
    median_time = statistics.median(times)

    # Weak scaling efficiency
    efficiency = baseline / median_time

    # Scaled speedup
    scaled_speedup = t * efficiency

    # Gustafson's law prediction
    gustafson = t - serial_fraction * (t - 1)

    summary.append({
        "threads": t,
        "median_seconds": median_time,
        "min_seconds": min(times),
        "max_seconds": max(times),
        "spread_seconds": max(times) - min(times),
        "efficiency": efficiency,
        "scaled_speedup": scaled_speedup,
        "gustafson_prediction": gustafson
    })

with open("results/part_d/d2_summary.csv", "w", newline="") as f:

    fields = list(summary[0].keys())
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summary)

threads = [r["threads"] for r in summary]
efficiencies = [r["efficiency"] for r in summary]
measured = [r["scaled_speedup"] for r in summary]
predicted = [r["gustafson_prediction"] for r in summary]

plt.figure()
plt.plot(threads, efficiencies, marker="o")
plt.axhline(1.0, linestyle="--", label="Ideal efficiency")
plt.xlabel("Number of threads")
plt.ylabel("Efficiency")
plt.title("D2 - Weak Scaling Efficiency")
plt.legend()
plt.grid(True)
plt.savefig("results/part_d/d2_efficiency.png", dpi=160)

plt.figure()
plt.plot(threads, measured, marker="o", label="Measured scaled speedup")
plt.plot(threads, predicted, linestyle="--", label="Gustafson prediction")
plt.xlabel("Number of threads")
plt.ylabel("Scaled speedup")
plt.title("D2 - Gustafson's Law")
plt.legend()
plt.grid(True)
plt.savefig("results/part_d/d2_gustafson.png", dpi=160)

print("\nD2 SUMMARY")

for r in summary:
    print(
        f"T={r['threads']:2d} "
        f"Median={r['median_seconds']:.3f}s "
        f"Efficiency={r['efficiency']:.3f} "
        f"Scaled speedup={r['scaled_speedup']:.3f} "
        f"Gustafson={r['gustafson_prediction']:.3f}"
    )

print("\nSaved results/part_d/d2_summary.csv")
print("Saved D2 graphs")

