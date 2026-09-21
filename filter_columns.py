import sys
import pandas as pd
from pathlib import Path
# check a csv

proj_dir = Path.cwd()
input_file = str(sys.argv[1])
input_regex = str(sys.argv[2])

with open(input_file, 'r') as data:
    df = pd.read_csv(data, sep='\t')
    csv_df = df.filter(regex=fr"{input_regex}")
    csv_df.to_csv(f"{proj_dir}/{input_file[:len(input_file)-4]}_filtered.csv", sep="\t")
    