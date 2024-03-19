import json
import pandas as pd 
import requests

with open("/mnt/d/PhD_coding/SQL/data_MSFD_Marche_MOD1E.json") as file:
    data = json.loads(file.read())
    df = pd.DataFrame.from_dict(data)

cols = ["code", "fullpath"]

#name_path = df[cols]
 #& name_path["code"].str.contains("Veneto")
#Region_df = name_path[name_path["code"].str.contains("Modulo_1") & name_path["code"].str.contains("Veneto")]

inital_path = "http://www.db-strategiamarina.isprambiente.it/api/v1.0.0/download/"    

for code, fullpath in df[cols].values: 
    print(code, fullpath)
    req = requests.get(inital_path + fullpath, allow_redirects=True)
    open("/mnt/d/PHD/MSFD/Data/Modulo1/" + code + f".{fullpath.rsplit('.')[1]}", "wb").write(req.content)
# extension = Veneto.iloc[0]["fullpath"].rsplit(".")[1]

# print("/mnt/d/PHD/MSFD/Data/Modulo1/" + Veneto.iloc[0]["code"] + f".{extension}")
# 