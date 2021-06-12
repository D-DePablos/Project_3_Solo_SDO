"""This script creates the cases for constant velocity and accelerated cases
"""
from datetime import datetime
import pickle

cases = [
    {
        "aiaTime": datetime(2020, 5, 28, 11, 0),
        "matchTime": datetime(2020, 5, 31, 20, 0),
        "soloDurn": 1,
        "caseName": "AIA_28th_T11",
        "MARGINHOURSSOLO": 12,
    },
    {
        "aiaTime": datetime(2020, 5, 28, 9, 0),
        "matchTime": datetime(2020, 5, 31, 18, 0),
        "soloDurn": 2,
        "caseName": "AIA_28th_T09"
    },
    {
        "aiaTime": datetime(2020, 5, 28, 7, 0),
        "matchTime": datetime(2020, 5, 31, 16, 0),
        "soloDurn": 2,
        "caseName": "AIA_28th_T07"
    },
    {
        "aiaTime": datetime(2020, 5, 28, 5, 0),
        "matchTime": datetime(2020, 5, 31, 14, 0),
        "soloDurn": 2,
        "caseName": "AIA_28th_T05"
    },
    {
        "aiaTime": datetime(2020, 5, 28, 3, 0),
        "matchTime": datetime(2020, 5, 31, 12, 0),
        "soloDurn": 2,
        "caseName": "AIA_28th_T03"
    },
    {
        "aiaTime": datetime(2020, 5, 28, 1, 0),
        "matchTime": datetime(2020, 5, 31, 10, 0),
        "soloDurn": 2,
        "caseName": "AIA_28th_T01"
    },
    {
        "aiaTime": datetime(2020, 5, 27, 23, 30),
        "matchTime": datetime(2020, 5, 31, 9, 0),
        "soloDurn": 1,
        "caseName": "AIA_27th_T23"
    },
]

accCases = [
    {
        "aiaTime": datetime(2020, 5, 27, 8, 0),
        "matchTime": datetime(2020, 5, 31, 20, 0),
        "soloDurn": 2,
        "caseName": "ACC_AIA_27th_T08",
        "MARGINHOURSSOLO": 12,
    },
    {
        "aiaTime": datetime(2020, 5, 27, 11, 0),
        "matchTime": datetime(2020, 5, 31, 23, 0),
        "soloDurn": 1,
        "caseName": "ACC_AIA_27th_T11"
    },
]

with open(
        "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/EMDComparison/pickleCases/consCases.pickle",
        "wb") as f:
    pickle.dump(cases, f)

with open(
        "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/EMDComparison/pickleCases/accCases.pickle",
        "wb") as f:
    pickle.dump(accCases, f)