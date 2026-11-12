# ARYAN SHOKEEN MT26116
import csv
import os
import re
import subprocess
import time

os.makedirs("results/part_d/d5_snapshots", exist_ok=True)

subprocess.run(["make"], check=True)

N = 1000
DAYS = 200
T = 8
C = 64
REPEATS = 5

K_VALUES = [10, 20]
LATENCIES = [0, 10, 50]
Q_VALUES = [1, 4]

RAW = "results/part_d/d5_raw.csv"

cases = [
    {"mode": "none", "K": 0, "latency": 0, "Q": 0}
]

for K in K_VALUES:
    for latency in LATENCIES:

        cases.append({
            "mode": "blocking",
            "K": K,
            "latency": latency,
            "Q": 0
        })

        for Q in Q_VALUES:
            cases.append({
                "mode": "async",
                "K": K,
                "latency": latency,
                "Q": Q
            })

fields = [
    "mode", "K", "latency_ms", "Q", "repetition",
    "wall_seconds", "simulation_seconds", "max_rss_mb"
]

with open(RAW, "w", newline="") as f:

    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()

    total = len(cases) * REPEATS
    completed = 0

    for case in cases:

        for repetition in range(1, REPEATS + 1):

            command = [
                "./pandemic_sim",
                "-N", str(N),
                "-D", str(DAYS),
                "-s", "42",
                "-T", str(T),
                "-C", str(C),
                "--pool",
                "--partition", "interleaved",
                "--fixed-days",
                "--benchmark",
                "-K", str(case["K"]),
                "--io-latency", str(case["latency"])
            ]

            if case["mode"] != "none":
                command.extend([
                    "--snapshot-prefix",
                    "results/part_d/d5_snapshots/snapshot"
                ])

            if case["mode"] == "async":
                command.extend([
                    "--async-snapshots",
                    "-Q", str(case["Q"])
                ])

            start = time.perf_counter()

            result = subprocess.run(
                ["/usr/bin/time", "-l"] + command,
                capture_output=True,
                text=True,
                check=True
            )

            wall_seconds = time.perf_counter() - start

            match = re.search(
                r"Simulation loop time:\s*([0-9.]+)",
                result.stdout
            )

            memory_match = re.search(
                r"(\d+)\s+maximum resident set size",
                result.stderr
            )

            if not match:
                raise RuntimeError("Simulation timing missing")

            simulation_seconds = float(match.group(1))

            max_rss_mb = (
                int(memory_match.group(1)) / (1024 * 1024)
                if memory_match else ""
            )

            writer.writerow({
                "mode": case["mode"],
                "K": case["K"],
                "latency_ms": case["latency"],
                "Q": case["Q"],
                "repetition": repetition,
                "wall_seconds": wall_seconds,
                "simulation_seconds": simulation_seconds,
                "max_rss_mb": max_rss_mb
            })

            f.flush()

            completed += 1

            print(
                f"{completed}/{total} | "
                f"Mode={case['mode']} "
                f"K={case['K']} "
                f"Latency={case['latency']}ms "
                f"Q={case['Q']} "
                f"Time={wall_seconds:.3f}s "
                f"Memory={max_rss_mb} MB"
            )

print("\nD5 benchmark completed.")
print("Results:", RAW)

