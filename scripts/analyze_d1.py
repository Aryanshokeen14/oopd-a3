# ARYAN SHOKEEN MT26116
import csv
import statistics
import matplotlib.pyplot as plt

input_file = "results/part_d/d1_raw.csv"
summary_file = "results/part_d/d1_summary.csv"

groups = {}

with open(input_file, newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:
        threads = int(row["threads"])
        seconds = float(row["seconds"])

        groups.setdefault(threads, []).append(seconds)

if 1 not in groups:
    raise ValueError("T=1 baseline is missing")

baseline = statistics.median(groups[1])

summary = []

for threads in sorted(groups):

    times = groups[threads]

    median_time = statistics.median(times)
    min_time = min(times)
    max_time = max(times)
    spread = max_time - min_time

    speedup = baseline / median_time
    efficiency = speedup / threads

    if threads == 1:
        karp_flatt = ""
    else:
        karp_flatt = (
            (1.0 / speedup - 1.0 / threads)
            / (1.0 - 1.0 / threads)
        )

    summary.append({
        "threads": threads,
        "median_seconds": median_time,
        "min_seconds": min_time,
        "max_seconds": max_time,
        "spread_seconds": spread,
        "speedup": speedup,
        "efficiency": efficiency,
        "karp_flatt": karp_flatt
    })

with open(summary_file, "w", newline="") as f:

    fields = [
        "threads",
        "median_seconds",
        "min_seconds",
        "max_seconds",
        "spread_seconds",
        "speedup",
        "efficiency",
        "karp_flatt"
    ]

    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summary)

threads = [r["threads"] for r in summary]
speedups = [r["speedup"] for r in summary]
efficiencies = [r["efficiency"] for r in summary]

kf_values = [
    r["karp_flatt"]
    for r in summary
    if r["karp_flatt"] != ""
]

average_f = statistics.mean(kf_values)

amdahl = [
    1.0 / (average_f + (1.0 - average_f) / p)
    for p in threads
]

ideal = threads

# Speedup graph
plt.figure()
plt.plot(threads, speedups, marker="o", label="Measured speedup")
plt.plot(threads, ideal, linestyle="--", label="Ideal speedup")
plt.plot(threads, amdahl, linestyle=":", label="Amdahl prediction")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.title("D1 - Strong Scaling Speedup")
plt.legend()
plt.grid(True)
plt.savefig("results/part_d/d1_speedup.png", dpi=160)

# Efficiency graph
plt.figure()
plt.plot(threads, efficiencies, marker="o")
plt.xlabel("Number of threads")
plt.ylabel("Parallel efficiency")
plt.title("D1 - Strong Scaling Efficiency")
plt.grid(True)
plt.savefig("results/part_d/d1_efficiency.png", dpi=160)

print("\nD1 SUMMARY")
print("Average Karp-Flatt serial fraction:", average_f)

for row in summary:
    print(
        f"T={row['threads']:2d} "
        f"Median={row['median_seconds']:.3f}s "
        f"Spread={row['spread_seconds']:.3f}s "
        f"Speedup={row['speedup']:.3f} "
        f"Efficiency={row['efficiency']:.3f} "
        f"Karp-Flatt={row['karp_flatt']}"
    )

print("\nSaved:")
print(summary_file)
print("results/part_d/d1_speedup.png")
print("results/part_d/d1_efficiency.png")

