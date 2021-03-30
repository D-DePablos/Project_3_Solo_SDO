from glob import glob
import pandas as pd
import numpy as np
from datetime import datetime


class LcurveManager:
    def __init__(
        self,
        csvPath="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/",
        objCad=60,
        wavelength=193,
    ):
        Lc_csv = glob(f"{csvPath}{wavelength}_complete_lcurves.csv")[0]
        self.df = pd.read_csv(Lc_csv)
        self.df.index = pd.to_datetime(self.df["Time"])
        del self.df["Time"]
        for column in self.df.columns:
            # Clear out columns which are empty
            if self.df[column].mean() == 0:
                del self.df[column]

        self.df = self.df.resample(f"{objCad}s").mean()


if __name__ == "__main__":
    lcDF = LcurveManager()