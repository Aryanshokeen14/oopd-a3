
import csv
import statistics
import matplotlib.pyplot as plt

RAW = "results/part_d/d5_raw.csv"

groups = {}

with open(RAW) as f:
    for row in csv.DictReader(f):

        key = (
            row["mode"],
            int(row["K"]),
            int(row["latency_ms"]),
            int(row["Q"])
        )

        groups.setdefault(key, []).append(row)


def median_value(rows, field):
    values = [float(r[field]) for r in rows if r[field] != ""]
    return statistics.median(values) if values else ""


baseline = statistics.median([
    float(r["wall_seconds"])
    for r in groups[("none", 0, 0, 0)]
])

summary = []

for key, rows in groups.items():

    mode, K, latency, Q = key

    wall = median_value(rows, "wall_seconds")
    loop = median_value(rows, "simulation_seconds")
    memory = median_value(rows, "max_rss_mb")

    hidden_seconds = ""
    hidden_percent = ""

    if mode == "async":

        blocking = groups[("blocking", K, latency, 0)]
        blocking_wall = median_value(blocking, "wall_seconds")

        extra_blocking_time = blocking_wall - baseline

        hidden_seconds = blocking_wall - wall

        if extra_blocking_time > 0:
            hidden_percent = (
                hidden_seconds / extra_blocking_time
            ) * 100

    summary.append({
        "mode": mode,
        "K": K,
        "latency_ms": latency,
        "Q": Q,
        "median_wall_seconds": wall,
        "spread_seconds": (
            max(float(r["wall_seconds"]) for r in rows)
            - min(float(r["wall_seconds"]) for r in rows)
        ),
        "median_simulation_seconds": loop,
        "median_max_rss_mb": memory,
        "hidden_seconds": hidden_seconds,
        "hidden_percent": hidden_percent
    })


fields = list(summary[0].keys())

with open("results/part_d/d5_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summary)


# Total time graph
for K in [10, 20]:

    plt.figure()

    for mode, Q, label in [
        ("blocking", 0, "Blocking"),
        ("async", 1, "Async Q=1"),
        ("async", 4, "Async Q=4")
    ]:

        times = []

        for latency in [0, 10, 50]:

            key = (mode, K, latency, Q)
            times.append(
                statistics.median([
                    float(r["wall_seconds"])
                    for r in groups[key]
                ])
            )

        plt.plot([0, 10, 50], times, marker="o", label=label)

    plt.axhline(
        baseline,
        linestyle="--",
        label="No snapshots"
    )

    plt.xlabel("I/O latency (ms)")
    plt.ylabel("Total execution time (seconds)")
    plt.title(f"D5 - Snapshot Performance K={K}")
    plt.legend()
    plt.grid(True)

    plt.savefig(
        f"results/part_d/d5_time_K{K}.png",
        dpi=160
    )


# Memory graph for async snapshots
plt.figure()

for K in [10, 20]:

    for latency in [10, 50]:

        values = []

        for Q in [1, 4]:
            values.append(
                statistics.median([
                    float(r["max_rss_mb"])
                    for r in groups[("async", K, latency, Q)]
                    if r["max_rss_mb"] != ""
                ])
            )

        plt.plot(
            [1, 4], values,
            marker="o",
            label=f"K={K}, latency={latency}ms"
        )

plt.xlabel("Pending-write limit Q")
plt.ylabel("Maximum resident memory (MB)")
plt.title("D5 - Async Snapshot Memory")
plt.legend()
plt.grid(True)

plt.savefig("results/part_d/d5_memory.png", dpi=160)


print("\nD5 SUMMARY")

for r in summary:
    if r["mode"] == "none":
        continue

    print(
        f"Mode={r['mode']:8s} "
        f"K={r['K']:2d} "
        f"Latency={r['latency_ms']:2d}ms "
        f"Q={r['Q']} "
        f"Time={r['median_wall_seconds']:.3f}s "
        f"Memory={r['median_max_rss_mb']} MB "
        f"Hidden={r['hidden_percent']}"
    )

print("\nNo-snapshot baseline:", round(baseline, 4), "seconds")
print("Saved D5 summary and graphs.")

