import json
import pandas as pd 
import requests
import argparse

def main(): 

    parser = argparse.ArgumentParser("Script download files from the MSFD web page")
    parser.add_argument("path_datafile", help="Path of the json file containing the js object from which extract links (e.g. /mnt/d/PhD_coding/data_MSFD_Veneto.json)", type=str) 
    parser.add_argument("path_directory", help="Path of the directory where to store the files (e.g. /mnt/d/PHD/MSFD/Data/Modulo1/)", type=str)
    args = parser.parse_args()

    with open(args.path_datafile) as file:
        data = json.loads(file.read())
        df = pd.DataFrame.from_dict(data)

    cols = ["code", "fullpath"]

    inital_path = "http://www.db-strategiamarina.isprambiente.it/api/v1.0.0/download/"    

    for code, fullpath in df[cols].values: 
        print(code, fullpath)
        req = requests.get(inital_path + fullpath, allow_redirects=True)
        open(args.path_directory + code + f".{fullpath.rsplit('.')[1]}", "wb").write(req.content)


if __name__ == "__main__": 

    main()