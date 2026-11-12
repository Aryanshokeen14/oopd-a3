# ARYAN SHOKEEN MT26116
import csv
import statistics
import matplotlib.pyplot as plt

groups = {}

with open("results/part_d/d3_small_grids_raw.csv") as f:
    for row in csv.DictReader(f):
        key = (int(row["N"]), row["mode"])
        groups.setdefault(key, []).append(float(row["seconds"]))

modes = ["sequential", "thread-per-day", "thread-pool"]
sizes = [50, 100, 200, 500, 1000]

summary = []

for n in sizes:
    for mode in modes:

        times = groups[(n, mode)]

        summary.append({
            "N": n,
            "mode": mode,
            "median_seconds": statistics.median(times),
            "min_seconds": min(times),
            "max_seconds": max(times),
            "spread_seconds": max(times) - min(times)
        })

with open("results/part_d/d3_small_grids_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

for mode in modes:
    print(f"\n{mode}")

    for r in summary:
        if r["mode"] == mode:
            print(
                f"N={r['N']:4d} "
                f"Median={r['median_seconds']:.6f}s "
                f"Spread={r['spread_seconds']:.6f}s"
            )

sequential = {
    r["N"]: r["median_seconds"]
    for r in summary if r["mode"] == "sequential"
}

for mode in ["thread-per-day", "thread-pool"]:

    crossover = next(
        (
            n for n in sizes
            if next(
                r["median_seconds"]
                for r in summary
                if r["N"] == n and r["mode"] == mode
            ) < sequential[n]
        ),
        None
    )

    print(f"\nSmallest N where {mode} beats sequential: {crossover}")

plt.figure()

for mode in modes:
    values = [
        next(
            r["median_seconds"]
            for r in summary
            if r["N"] == n and r["mode"] == mode
        )
        for n in sizes
    ]

    plt.plot(sizes, values, marker="o", label=mode)

plt.xlabel("Grid dimension N")
plt.ylabel("Median execution time (seconds)")
plt.title("D3.3 - Small Grid Performance Comparison")
plt.legend()
plt.grid(True)
plt.savefig("results/part_d/d3_small_grids.png", dpi=160)

print("\nSaved D3.3 summary and graph.")

