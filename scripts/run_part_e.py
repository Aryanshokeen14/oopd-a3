
import csv
import subprocess
from pathlib import Path

N = 500
MAX_DAYS = 5000
SEEDS = range(42, 52)

ROOT = Path("results/part_e")
RUNS = ROOT / "runs"
RUNS.mkdir(parents=True, exist_ok=True)

experiments = []

def add_experiment(scenario, setting, beta, vaccination,
                   mask, lockdown, theta_on=0.02,
                   theta_off=0.005):

    for seed in SEEDS:
        experiments.append({
            "scenario": scenario,
            "setting": setting,
            "beta": beta,
            "vaccination": vaccination,
            "mask": mask,
            "lockdown": lockdown,
            "theta_on": theta_on,
            "theta_off": theta_off,
            "seed": seed
        })


# E1: Default parameters
add_experiment(
    "default", "default",
    0.08, 0.0, 1.0, 0.3
)

# E2: Transmission probability
for beta in [0.01, 0.02, 0.03, 0.05, 0.08, 0.10]:
    add_experiment(
        "beta", f"beta_{beta}",
        beta, 0.0, 1.0, 1.0
    )

# E3: Vaccination
for vaccination in [
    0.0, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9
]:
    add_experiment(
        "vaccination", f"v_{vaccination}",
        0.08, vaccination, 1.0, 1.0
    )

# E4: Masks and lockdown
interventions = [
    ("no_intervention", 1.0, 1.0),
    ("masks", 0.5, 1.0),
    ("lockdown", 1.0, 0.3),
    ("masks_lockdown", 0.5, 0.3)
]

for name, mask, lockdown in interventions:
    add_experiment(
        "intervention", name,
        0.08, 0.0, mask, lockdown
    )

# E5: Grid run without interventions for SIR comparison
add_experiment(
    "sir_grid", "no_intervention",
    0.08, 0.0, 1.0, 1.0
)


summary_file = ROOT / "run_summaries.csv"
summary_rows = []

total = len(experiments)

print("Total runs:", total)

for index, spec in enumerate(experiments, 1):

    folder = RUNS / spec["scenario"]
    folder.mkdir(parents=True, exist_ok=True)

    safe_setting = spec["setting"].replace(".", "p")
    output = folder / f"{safe_setting}_seed_{spec['seed']}.csv"

    command = [
        "./pandemic_sim",
        "-N", str(N),
        "-D", str(MAX_DAYS),
        "-s", str(spec["seed"]),
        "-b", str(spec["beta"]),
        "-v", str(spec["vaccination"]),
        "-m", str(spec["mask"]),
        "-l", str(spec["lockdown"]),
        "--theta-on", str(spec["theta_on"]),
        "--theta-off", str(spec["theta_off"]),
        "-o", str(output)
    ]

    subprocess.run(
        command,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True
    )

    with open(output, newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise RuntimeError(f"Empty result: {output}")

    peak = max(rows, key=lambda r: int(r["I"]))
    final = rows[-1]

    total_cells = N * N

    total_infected = (
        total_cells
        - int(final["S"])
        - int(final["V"])
    )

    summary_rows.append({
        "scenario": spec["scenario"],
        "setting": spec["setting"],
        "seed": spec["seed"],
        "beta": spec["beta"],
        "vaccination": spec["vaccination"],
        "mask": spec["mask"],
        "lockdown": spec["lockdown"],
        "peak_infected": int(peak["I"]),
        "peak_infected_fraction": int(peak["I"]) / total_cells,
        "day_of_peak": int(peak["day"]),
        "total_infected": total_infected,
        "total_infected_fraction": total_infected / total_cells,
        "total_deaths": int(final["D"]),
        "days": int(final["day"]),
        "epidemic_ended": int(final["I"]) == 0,
        "output_file": str(output)
    })

    print(f"{index}/{total} completed: {spec['scenario']} / {spec['setting']} / seed {spec['seed']}")


with open(summary_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=summary_rows[0].keys()
    )
    writer.writeheader()
    writer.writerows(summary_rows)

print("\nPART E EXPERIMENTS COMPLETED")
print("Runs:", total)
print("Summary:", summary_file)

not_ended = [
    r for r in summary_rows if not r["epidemic_ended"]
]

print("Runs reaching maximum days:", len(not_ended))

