
"""
1_get_csv_file_from_configs.py

Extracts parameter values from config JSON files and saves them as CSV and Excel files for further analysis.

USAGE EXAMPLES:
---------------
1. Extract from ddemodel/config/parameters_*.json (default):
    python analysis/code/1_get_csv_file_from_configs.py

2. Extract from ddemodel/config/parameters_sleep_*.json:
    python analysis/code/1_get_csv_file_from_configs.py --pattern parameters_sleep_

3. Extract from analysis/data/extracted_configs/parameters_sleep*.json:
    python analysis/code/1_get_csv_file_from_configs.py --config_dir analysis/data/extracted_configs --pattern parameters_sleep

OPTIONS:
--------
--config_dir : Directory containing config files (default: ddemodel/config)
--pattern    : Prefix pattern for config files (e.g., parameters_ or parameters_sleep_)
"""

import json
import os
import pandas as pd
import argparse

# === CONFIGURATION & ARGUMENTS ===

parser = argparse.ArgumentParser(
    description="Extract parameters from config files to CSV/XLSX. "
                "See script docstring for usage examples.")
parser.add_argument('--config_dir', type=str, default="model/config/base",
                    help='Directory with config files (default: ddemodel/model/config)')
parser.add_argument('--pattern', type=str, default="parameters_",
                    help='Prefix pattern for config files (e.g., parameters_ or parameters_sleep_)')
args = parser.parse_args()

CONFIG_DIR = args.config_dir
CONFIG_PATTERN = args.pattern
OUTPUT_DIR = "ddebiftool/helpers/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === MAIN LOGIC ===
data = []

# Only process files directly in CONFIG_DIR, not in subdirectories
for file in os.listdir(CONFIG_DIR):
    if file.startswith(CONFIG_PATTERN) and file.endswith(".json"):
        file_path = os.path.join(CONFIG_DIR, file)
        with open(file_path, "r") as f:
            content = json.load(f)
            participant = content.get("participant")
            parameters = content.get("parameters", {})
            row = {"Participant": participant, **parameters}
            data.append(row)

df = pd.DataFrame(data)

excel_file = os.path.join(OUTPUT_DIR, "parameters.xlsx")
csv_file = os.path.join(OUTPUT_DIR, "parameters.csv")
df.to_excel(excel_file, index=False)
df.to_csv(csv_file, index=False)

print(f"Saved Excel file to: {os.path.abspath(excel_file)}")
print(f"Saved CSV file to:   {os.path.abspath(csv_file)}")

