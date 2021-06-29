from datetime import datetime, timedelta
import pickle

# Info about IS data
# startIS = datetime(2018, 10, 31, 8,0)
# endIS = datetime(2018, 11, 2, 8)

MARGIN_IS = 24
MARGIN_HMI = 4


def main():
    HMI_base = datetime(2018, 10, 28, 0, 0)
    HMI_max = datetime(2018, 10, 30, 7, 48)

    IS_base = datetime(2018, 11, 1, 8)

    cases = []
    i = 0
    tHMI = HMI_base

    while tHMI <= (HMI_max - timedelta(hours=MARGIN_HMI)):
        tHMI = HMI_base + timedelta(hours=MARGIN_HMI) * i
        tIS = IS_base

        cases.append({
            "aiaTime": tHMI,
            "matchTime": tIS,
            "soloDurn": 1,
            "caseName": f"HMI_{tHMI.day}_T{tHMI.hour:02d}",
            "MARGINHOURSSOLO": MARGIN_IS
        })

        i += 1

    with open("/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/ISSI/data/HMIcases.pickle", "wb") as f:
        pickle.dump(cases, f)


if __name__ == "__main__":
    main()
