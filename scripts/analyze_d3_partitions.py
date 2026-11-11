
import csv
import statistics
import matplotlib.pyplot as plt

groups = {}

with open("results/part_d/d3_partitions_raw.csv") as f:
    for row in csv.DictReader(f):
        name = row["partition"]
        groups.setdefault(name, []).append(float(row["seconds"]))

summary = []

for name in ["strips", "blocks", "interleaved"]:

    times = groups[name]

    summary.append({
        "partition": name,
        "median_seconds": statistics.median(times),
        "min_seconds": min(times),
        "max_seconds": max(times),
        "spread_seconds": max(times) - min(times)
    })

with open("results/part_d/d3_partitions_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

plt.figure()
plt.bar(
    [r["partition"] for r in summary],
    [r["median_seconds"] for r in summary]
)
plt.xlabel("Partitioning strategy")
plt.ylabel("Median execution time (seconds)")
plt.title("D3.2 - Partitioning Strategy Comparison")
plt.grid(axis="y")
plt.savefig("results/part_d/d3_partitions.png", dpi=160)

for r in summary:
    print(
        f"{r['partition']}: "
        f"Median={r['median_seconds']:.3f}s "
        f"Spread={r['spread_seconds']:.3f}s"
    )

best = min(summary, key=lambda x: x["median_seconds"])
print("\nBest partition:", best["partition"])

