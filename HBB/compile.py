import os
import pandas as pd
import numpy as np
from glob import glob

# Function to calculate bootstrap error
def bootstrap_error(data, n_iterations=1000):
    means = [np.mean(np.random.choice(data, size=len(data), replace=True)) for _ in range(n_iterations)]
    return np.std(means)

# Define expected temperatures to enforce proper order
expected_temps = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95]

# Sort dist.* files in expected temperature order
dist_files = sorted(glob("dist.*"), key=lambda x: expected_temps.index(int(x.split('.')[-1])))

# Initialize containers
step_set = set()
raw_data = {}
avg_row = {"Step": "AVG"}
err_row = {"Step": "ERR"}
summary_data = []

# Process each dist file
for file in dist_files:
    temp = int(file.split('.')[-1])
    with open(file) as f:
        lines = f.readlines()

    # Filter and clean data
    lines = [line.strip() for line in lines if line.strip() and not line.startswith('-')]
    lines = lines[1:]  # Skip header line

    steps, ww, wh_top, wh_bottom = [], [], [], []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            step = int(parts[0])
            step_set.add(step)
            steps.append(step)
            ww.append(float(parts[1]))
            wh_top.append(float(parts[2]))
            wh_bottom.append(float(parts[3]))

    # Store for sideways data
    raw_data[f"{temp}_HB_WW"] = dict(zip(steps, ww))
    raw_data[f"{temp}_HB_WH_top"] = dict(zip(steps, wh_top))
    raw_data[f"{temp}_HB_WH_bottom"] = dict(zip(steps, wh_bottom))
    
    avg_row[f"{temp}_HB_WW"] = np.mean(ww)
    avg_row[f"{temp}_HB_WH_top"] = np.mean(wh_top)
    avg_row[f"{temp}_HB_WH_bottom"] = np.mean(wh_bottom)

    err_row[f"{temp}_HB_WW"] = bootstrap_error(ww)
    err_row[f"{temp}_HB_WH_top"] = bootstrap_error(wh_top)
    err_row[f"{temp}_HB_WH_bottom"] = bootstrap_error(wh_bottom)

    # Store for compact summary
    summary_data.append([
        temp,
        np.mean(ww), bootstrap_error(ww),
        np.mean(wh_top), bootstrap_error(wh_top),
        np.mean(wh_bottom), bootstrap_error(wh_bottom)
    ])

# Construct sideways DataFrame
final_rows = [avg_row, err_row]
for step in sorted(step_set):
    row = {"Step": step}
    for key in raw_data:
        row[key] = raw_data[key].get(step, "")
    final_rows.append(row)

df_wide = pd.DataFrame(final_rows)
df_wide.to_csv("hbond_summary_sideways.csv", index=False)

# Construct compact summary DataFrame
df_summary = pd.DataFrame(
    summary_data,
    columns=[
        "Temp", 
        "OA_HB_WW", "Error_WW", 
        "OA_HB_WH_top", "Error_WH_top", 
        "OA_HB_WH_bottom", "Error_WH_bottom"
    ]
)
df_summary = df_summary.sort_values("Temp")
df_summary.to_csv("hbond_summary_compact.csv", index=False)

print("✅ Output files written:")
print(" - hbond_summary_sideways.csv")
print(" - hbond_summary_compact.csv")
