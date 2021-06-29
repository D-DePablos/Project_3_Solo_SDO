"""
Generate the cases for all possible SolO - AIA times (every hour)
"""

from datetime import datetime, timedelta
import pickle

MARGINHOURSSOLO = 25

start = datetime(2020, 5, 30, 12)
end = datetime(2020, 6, 1, 13, 32)

midSOLO = start + timedelta(seconds=(end - start).total_seconds() / 2)


def main():

    AIA_base = datetime(2020, 5, 27, 19, 30)
    AIA_max = datetime(2020, 5, 28, 13, 1)
    AIA_dt = 1  # How many hours to advance

    # Around the base case, MARGINHOURSSOLO is utilised
    SOLO_base = datetime(2020, 5, 31, 20, 0)

    cases = []
    i = 0
    tAIA = AIA_base
    while tAIA <= (AIA_max - timedelta(hours=1)):
        tAIA = AIA_base + timedelta(hours=AIA_dt) * i
        tSOLO = SOLO_base

        cases.append({
            "aiaTime": tAIA,
            "matchTime": tSOLO,
            "soloDurn": 1,
            "caseName": f"AIA_{tAIA.day}_T{tAIA.hour:02d}",
            "MARGINHOURSSOLO": MARGINHOURSSOLO,
        })

        i += 1

    with open(
            "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/EMDComparison/pickleCases/newCases_ALLSOLO.pickle",
            "wb") as f:
        pickle.dump(cases, f)


if __name__ == "__main__":
    main()
