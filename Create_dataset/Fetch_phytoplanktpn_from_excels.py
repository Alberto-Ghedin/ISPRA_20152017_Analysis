import pymysql
from pymysql.constants import CLIENT
import os
import json
import pandas as pd
import string
import numpy as np
import argparse
from datetime import datetime


excel_path = "/mnt/d/PHD/MSFD/Data/Modulo1/"

excel_file_names = filter(lambda x:  "Modulo_1" in x, os.listdir(excel_path))

df = df = pd.DataFrame({"NationalStationID": [], "Year": [], "Month": [], "Day": [], "Time": [
    ], "SampleDepth": [], "Gruppo": [], "Taxon": [], "NuovoTaxon": [], "Autore": [], "Num_cell_l": []})
for excel_file_name in excel_file_names:
    print("_"*10 + f"reading {excel_file_name}" + "_"*10)
    # INSERTING DATA
    df = pd.concat([
        df,
        pd.read_excel(excel_path + excel_file_name, "Fitoplancton").loc[:, ["NationalStationID",
                                                                    "Year",
                                                                    "Month",
                                                                    "Day",
                                                                    "Time",
                                                                    "SampleDepth",
                                                                    "Gruppo",
                                                                    "Taxon",
                                                                    "NuovoTaxon",
                                                                    "Autore",
                                                                    "Num_cell_l"
                                                                    ]]
    ]
    )

print(df.describe())

df.to_csv(excel_path + "/phyto_abund_raw.csv")
