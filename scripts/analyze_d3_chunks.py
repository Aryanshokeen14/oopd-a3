
import csv
import statistics
import matplotlib.pyplot as plt

groups = {}

with open("results/part_d/d3_chunks_raw.csv") as f:
    for row in csv.DictReader(f):
        c = int(row["chunks"])
        groups.setdefault(c, []).append(float(row["seconds"]))

summary = []

for c in sorted(groups):
    times = groups[c]

    summary.append({
        "threads": 8,
        "chunks": c,
        "median_seconds": statistics.median(times),
        "min_seconds": min(times),
        "max_seconds": max(times),
        "spread_seconds": max(times) - min(times)
    })

best = min(summary, key=lambda x: x["median_seconds"])

with open("results/part_d/d3_chunks_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

plt.figure()
plt.plot(
    [r["chunks"] for r in summary],
    [r["median_seconds"] for r in summary],
    marker="o"
)
plt.xlabel("Number of chunks (C)")
plt.ylabel("Median execution time (seconds)")
plt.title("D3.1 - Chunk Size Comparison")
plt.grid(True)
plt.savefig("results/part_d/d3_chunks.png", dpi=160)

for r in summary:
    print(
        f"C={r['chunks']:3d} "
        f"Median={r['median_seconds']:.3f}s "
        f"Spread={r['spread_seconds']:.3f}s"
    )

print("\nBest chunk count:", best["chunks"])
print("Best median time:", best["median_seconds"])

