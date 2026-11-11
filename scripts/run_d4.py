
import concurrent.futures
import subprocess
import csv
import time
import re
import os

os.makedirs("results/part_d", exist_ok=True)
subprocess.run(["make"], check=True)

P = 8
N = 1000
DAYS = 200
REPEATS = 5

configs = [(1, 8), (2, 4), (4, 2), (8, 1)]

betas = [0.05, 0.08, 0.10, 0.12]
vaccinations = [0.0, 0.2, 0.4, 0.6]
seeds = [42, 43]

specs = []
run_id = 1

for beta in betas:
    for vaccination in vaccinations:
        for seed in seeds:
            specs.append({
                "run_id": run_id,
                "beta": beta,
                "vaccination": vaccination,
                "seed": seed
            })
            run_id += 1

assert len(specs) == 32

raw_rows = []
batch_rows = []

def execute_run(spec, T):

    command = [
        "./pandemic_sim",
        "-N", str(N),
        "-D", str(DAYS),
        "-s", str(spec["seed"]),
        "-b", str(spec["beta"]),
        "-v", str(spec["vaccination"]),
        "-T", str(T),
        "-C", str(8 * T),
        "--pool",
        "--partition", "interleaved",
        "--benchmark"
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True
    )

    match = re.search(
        r"Simulation loop time:\s*([0-9.]+)\s+seconds",
        result.stdout
    )

    if not match:
        raise RuntimeError("Timing not found")

    return float(match.group(1))


for R, T in configs:

    for repetition in range(1, REPEATS + 1):

        print(f"\nR={R}, T={T}, repetition={repetition}")

        start = time.perf_counter()
        completed = 0

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=R
        ) as executor:

            futures = {
                executor.submit(execute_run, spec, T): spec
                for spec in specs
            }

            for future in concurrent.futures.as_completed(futures):

                spec = futures[future]
                seconds = future.result()
                completed += 1

                raw_rows.append({
                    "R": R,
                    "T": T,
                    "repetition": repetition,
                    "run_id": spec["run_id"],
                    "beta": spec["beta"],
                    "vaccination": spec["vaccination"],
                    "seed": spec["seed"],
                    "simulation_seconds": seconds
                })

                print(f"{completed}/32 completed")

        wall_time = time.perf_counter() - start

        batch_rows.append({
            "R": R,
            "T": T,
            "repetition": repetition,
            "total_runs": 32,
            "wall_seconds": wall_time,
            "runs_per_second": 32 / wall_time
        })

        print(f"Total time: {wall_time:.3f}s")
        print(f"Throughput: {32 / wall_time:.3f} runs/sec")


for filename, rows in [
    ("results/part_d/d4_raw.csv", raw_rows),
    ("results/part_d/d4_batches.csv", batch_rows)
]:

    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=rows[0].keys()
        )
        writer.writeheader()
        writer.writerows(rows)

print("\nD4 FIVE-REPETITION BENCHMARK COMPLETED")


