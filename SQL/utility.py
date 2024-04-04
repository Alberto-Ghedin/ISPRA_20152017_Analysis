import json
from openpyxl import load_workbook
import xlrd
import os

def clean_dict(dic, key_2_delete): 

    for key in dic.keys():
        if isinstance(dic[key],dict): 
            clean_dict(dic[key], key_2_delete)
        elif dic[key] == key_2_delete: 
            del dic[key]
    
    return dic

def remove_additional_white_spaces(string): 
    return " ".join(string.strip().split())