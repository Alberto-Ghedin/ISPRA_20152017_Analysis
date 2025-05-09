import os
import pandas as pd
import string
from datetime import datetime

_HOME_ = os.path.expanduser("~")
excel_path = _HOME_ + "/PHD/MSFD/Data/Modulo1/"

excel_file_names = filter(lambda x:  "Modulo_1" in x, os.listdir(excel_path))

df = df = pd.DataFrame({
    "NationalStationID": [], 
    "Year": [], 
    "Month": [], 
    "Day": [], 
    "Time": [], 
    "SampleID": [], 
    "SampleDepth": [], 
    "FitoZoo": [], 
    "ClasseDim": [], 
    "Abbondanza_CD": []
})
for excel_file_name in excel_file_names:
    print("_"*10 + f"reading {excel_file_name}" + "_"*10)
    # INSERTING DATA
    df = pd.concat([
        df,
        pd.read_excel(excel_path + excel_file_name, "Plancton_CD").loc[:, [
    "NationalStationID", 
    "Year", 
    "Month", 
    "Day", 
    "Time", 
    "SampleID", 
    "SampleDepth", 
    "FitoZoo", 
    "ClasseDim", 
    "Abbondanza_CD"]].assign(file_name = excel_file_name)
    ]
    )

print(df.describe())

df.to_csv(excel_path + "/nano_micro_phyt.csv", index = 0)
