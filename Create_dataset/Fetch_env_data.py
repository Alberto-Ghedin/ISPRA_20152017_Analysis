import os
import pandas as pd
import string
from datetime import datetime
from os import path
import json

with open(path.expanduser("~") + "/sys_specific.json") as file: 
    params = json.load(file)
_HOME_ = params["home"] + "/PHD"

excel_path = _HOME_ + "/MSFD/Data/Modulo1/"

excel_file_names = filter(lambda x:  "Modulo_1" in x, os.listdir(excel_path))

df = df = pd.DataFrame({"NationalStationID": [], "Year": [], "Month": [], "Day": [], "Time": [
    ], "SampleDepth": [], "Determinand_Nutrients": [], "Determinand_Nutrients": [], "Concentration": [], "LOD_LOQ_Flag" : [], "Remarks": [], "file_name" : []})
for excel_file_name in excel_file_names:
    print("_"*10 + f"reading {excel_file_name}" + "_"*10)
    # INSERTING DATA
    temp = pd.read_excel(excel_path + excel_file_name, "Chimico-fisici")
    if "Remarks*" in temp.columns: 
        temp.rename(columns= {"Remarks*" : "Remarks"}, inplace=True)
    df = pd.concat([
        df,
        temp.loc[:, ["NationalStationID",
                                                                    "Year",
                                                                    "Month",
                                                                    "Day",
                                                                    "Time",
                                                                    "SampleDepth",
                                                                    "Determinand_Nutrients",
                                                                    "Concentration", 
                                                                    "LOD_LOQ_Flag",
                                                                    "Remarks"
                                                                    ]].assign(file_name = excel_file_name)
    ]
    )

print(df.describe())

df.to_csv(excel_path + "/chemical_fis_raw.csv", index = 0)
