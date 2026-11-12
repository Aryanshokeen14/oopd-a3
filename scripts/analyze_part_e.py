
import csv
import statistics
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path("results/part_e")

with open(ROOT / "run_summaries.csv") as f:
    summaries = list(csv.DictReader(f))


def mean_std(rows, field):
    values = [float(r[field]) for r in rows]
    mean = statistics.mean(values)
    std = statistics.stdev(values) if len(values) > 1 else 0.0
    return mean, std


def write_group_summary(filename, scenario, setting_field):

    groups = {}

    for row in summaries:
        if row["scenario"] == scenario:
            key = row[setting_field]
            groups.setdefault(key, []).append(row)

    output = []

    for key, rows in groups.items():

        infected, infected_std = mean_std(
            rows, "total_infected_fraction"
        )

        peak, peak_std = mean_std(
            rows, "peak_infected_fraction"
        )

        deaths, deaths_std = mean_std(
            rows, "total_deaths"
        )

        peak_day, peak_day_std = mean_std(
            rows, "day_of_peak"
        )

        output.append({
            setting_field: key,
            "mean_total_infected_fraction": infected,
            "std_total_infected_fraction": infected_std,
            "mean_peak_fraction": peak,
            "std_peak_fraction": peak_std,
            "mean_deaths": deaths,
            "std_deaths": deaths_std,
            "mean_peak_day": peak_day,
            "std_peak_day": peak_day_std
        })

    output.sort(key=lambda r: float(r[setting_field]))

    filename = ROOT / filename

    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=output[0].keys()
        )
        writer.writeheader()
        writer.writerows(output)

    return output


def get_curve(scenario, setting):

    selected = [
        r for r in summaries
        if r["scenario"] == scenario
        and r["setting"] == setting
    ]

    runs = []

    for summary in selected:

        with open(summary["output_file"]) as f:
            runs.append(list(csv.DictReader(f)))

    if not runs:
        raise ValueError(f"No runs for {scenario}/{setting}")

    max_day = max(int(run[-1]["day"]) for run in runs)

    result = []

    for day in range(max_day + 1):

        record = {"day": day}

        for field in ["S", "I", "R", "D", "V"]:

            values = [
                float(run[min(day, len(run) - 1)][field])
                for run in runs
            ]

            record["mean_" + field] = statistics.mean(values)
            record["std_" + field] = statistics.stdev(values)

        result.append(record)

    return result


def save_curve(curve, filename):

    fields = ["day"]

    for field in ["S", "I", "R", "D", "V"]:
        fields.extend(["mean_" + field, "std_" + field])

    with open(ROOT / filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(curve)


# E1: Default epidemic curves
default_curve = get_curve("default", "default")

save_curve(default_curve, "e1_default_timeseries.csv")

days = [r["day"] for r in default_curve]
population = 500 * 500

plt.figure()

for field in ["S", "I", "R", "D"]:

    means = [r["mean_" + field] / population for r in default_curve]
    stds = [r["std_" + field] / population for r in default_curve]

    plt.plot(days, means, label=field)
    plt.fill_between(
        days,
        [max(0, m-s) for m, s in zip(means, stds)],
        [min(1, m+s) for m, s in zip(means, stds)],
        alpha=0.15
    )

plt.xlabel("Day")
plt.ylabel("Population fraction")
plt.title("E1 - Default Epidemic Spread")
plt.legend()
plt.grid(True)
plt.savefig(ROOT / "e1_default_epidemic.png", dpi=160)


# E2: Transmission probability
beta_summary = write_group_summary(
    "e2_beta_summary.csv", "beta", "beta"
)

betas = [float(r["beta"]) for r in beta_summary]
infected = [float(r["mean_total_infected_fraction"]) for r in beta_summary]
peaks = [float(r["mean_peak_fraction"]) for r in beta_summary]

plt.figure()
plt.plot(betas, infected, marker="o")
plt.xlabel("Transmission probability beta")
plt.ylabel("Mean total infected fraction")
plt.title("E2 - Effect of Transmission Probability")
plt.grid(True)
plt.savefig(ROOT / "e2_beta_outbreak.png", dpi=160)

plt.figure()
plt.plot(betas, peaks, marker="o")
plt.xlabel("Transmission probability beta")
plt.ylabel("Mean peak infected fraction")
plt.title("E2 - Peak Infection vs Beta")
plt.grid(True)
plt.savefig(ROOT / "e2_beta_peak.png", dpi=160)


# E3: Vaccination
vacc_summary = write_group_summary(
    "e3_vaccination_summary.csv",
    "vaccination",
    "vaccination"
)

vaccinations = [
    float(r["vaccination"]) for r in vacc_summary
]

infected = [
    float(r["mean_total_infected_fraction"])
    for r in vacc_summary
]

plt.figure()
plt.plot(vaccinations, infected, marker="o")
plt.xlabel("Vaccinated fraction")
plt.ylabel("Mean total infected fraction")
plt.title("E3 - Vaccination and Herd Immunity")
plt.grid(True)
plt.savefig(ROOT / "e3_vaccination.png", dpi=160)


# E4: Masks and lockdown
interventions = [
    "no_intervention",
    "masks",
    "lockdown",
    "masks_lockdown"
]

intervention_summary = []

for name in interventions:

    rows = [
        r for r in summaries
        if r["scenario"] == "intervention"
        and r["setting"] == name
    ]

    peak, peak_std = mean_std(rows, "peak_infected_fraction")
    peak_day, day_std = mean_std(rows, "day_of_peak")
    total, total_std = mean_std(rows, "total_infected_fraction")

    intervention_summary.append({
        "intervention": name,
        "mean_peak_fraction": peak,
        "std_peak_fraction": peak_std,
        "mean_peak_day": peak_day,
        "std_peak_day": day_std,
        "mean_total_infected_fraction": total,
        "std_total_infected_fraction": total_std
    })

with open(ROOT / "e4_interventions_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=intervention_summary[0].keys()
    )
    writer.writeheader()
    writer.writerows(intervention_summary)

plt.figure()

for name in interventions:

    curve = get_curve("intervention", name)

    plt.plot(
        [r["day"] for r in curve],
        [r["mean_I"] / population for r in curve],
        label=name
    )

plt.xlabel("Day")
plt.ylabel("Mean infected fraction")
plt.title("E4 - Effects of Masks and Lockdown")
plt.legend()
plt.grid(True)
plt.savefig(ROOT / "e4_interventions.png", dpi=160)


# E5: Standard SIR model comparison
grid_curve = get_curve("sir_grid", "no_intervention")

# Approximate homogeneous contact rate:
# Eight neighbours, each with transmission probability beta.
beta_sir = 8 * 0.08
gamma = 1.0 / 7.0

S = float(population - 10)
I = 10.0
R = 0.0

dt = 0.01
sir_curve = []

for day in range(len(grid_curve)):

    sir_curve.append({
        "day": day,
        "S": S,
        "I": I,
        "R": R
    })

    for _ in range(100):

        infection = beta_sir * S * I / population
        recovery = gamma * I

        S -= infection * dt
        I += (infection - recovery) * dt
        R += recovery * dt

sir_file = ROOT / "e5_sir_comparison.csv"

with open(sir_file, "w", newline="") as f:

    fields = [
        "day",
        "grid_S_fraction",
        "grid_I_fraction",
        "grid_R_fraction",
        "sir_S_fraction",
        "sir_I_fraction",
        "sir_R_fraction"
    ]

    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()

    for grid, sir in zip(grid_curve, sir_curve):

        writer.writerow({
            "day": grid["day"],
            "grid_S_fraction": grid["mean_S"] / population,
            "grid_I_fraction": grid["mean_I"] / population,
            "grid_R_fraction": grid["mean_R"] / population,
            "sir_S_fraction": sir["S"] / population,
            "sir_I_fraction": sir["I"] / population,
            "sir_R_fraction": sir["R"] / population
        })

plt.figure()

for field in ["S", "I", "R"]:

    plt.plot(
        days if len(default_curve) == len(grid_curve)
        else [r["day"] for r in grid_curve],
        [r["mean_" + field] / population for r in grid_curve],
        label="Grid " + field
    )

    plt.plot(
        [r["day"] for r in sir_curve],
        [r[field] / population for r in sir_curve],
        linestyle="--",
        label="SIR " + field
    )

plt.xlabel("Day")
plt.ylabel("Population fraction")
plt.title("E5 - Grid Simulation vs Standard SIR")
plt.legend()
plt.grid(True)
plt.savefig(ROOT / "e5_sir_comparison.png", dpi=160)


print("\nPART E ANALYSIS COMPLETED")

print("\nE2 - Beta experiment")
for r in beta_summary:
    print(r["beta"], r["mean_total_infected_fraction"])

print("\nE3 - Vaccination experiment")
for r in vacc_summary:
    print(r["vaccination"], r["mean_total_infected_fraction"])

print("\nE4 - Interventions")
for r in intervention_summary:
    print(
        r["intervention"],
        "Peak:", r["mean_peak_fraction"],
        "Peak day:", r["mean_peak_day"]
    )

print("\nE5 SIR approximation:")
print("beta_sir =", beta_sir)
print("gamma =", gamma)

print("\nAll graphs and summary CSV files saved in", ROOT)

