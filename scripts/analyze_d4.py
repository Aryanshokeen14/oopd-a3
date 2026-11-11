
import csv
import statistics
import matplotlib.pyplot as plt

groups = {}

with open("results/part_d/d4_batches.csv") as f:
    for row in csv.DictReader(f):
        key = (int(row["R"]), int(row["T"]))
        groups.setdefault(key, []).append(
            (float(row["wall_seconds"]),
             float(row["runs_per_second"]))
        )

summary = []

for (R, T), values in sorted(groups.items()):

    times = [v[0] for v in values]
    rates = [v[1] for v in values]

    summary.append({
        "R": R,
        "T": T,
        "median_wall_seconds": statistics.median(times),
        "spread_seconds": max(times) - min(times),
        "median_runs_per_second": statistics.median(rates)
    })

with open("results/part_d/d4_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

labels = [f"{r['R']}x{r['T']}" for r in summary]

plt.figure()
plt.bar(labels, [r["median_wall_seconds"] for r in summary])
plt.xlabel("Concurrent runs x threads per run")
plt.ylabel("Median total time (seconds)")
plt.title("D4 - Total Execution Time")
plt.grid(axis="y")
plt.savefig("results/part_d/d4_total_time.png", dpi=160)

plt.figure()
plt.bar(labels, [r["median_runs_per_second"] for r in summary])
plt.xlabel("Concurrent runs x threads per run")
plt.ylabel("Median runs per second")
plt.title("D4 - Throughput")
plt.grid(axis="y")
plt.savefig("results/part_d/d4_throughput.png", dpi=160)

for r in summary:
    print(
        f"R={r['R']} T={r['T']} "
        f"Median={r['median_wall_seconds']:.3f}s "
        f"Spread={r['spread_seconds']:.3f}s "
        f"Throughput={r['median_runs_per_second']:.3f}"
    )

best = max(summary, key=lambda x: x["median_runs_per_second"])

print(f"\nBest configuration: R={best['R']}, T={best['T']}")

