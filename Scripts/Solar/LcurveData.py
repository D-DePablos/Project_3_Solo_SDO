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
        """Manages getting lightcurves to correct format

        Args:
            csvPath (str, optional): Path to Lightcurve csv. Defaults to "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/".
            objCad (int, optional): objective Cadence in seconds. Defaults to 60.
            wavelength (int, optional): Wavelength to investigate. Defaults to 193.
        """
        Lc_csv = glob(f"{csvPath}complete_lcurves_{wavelength}.csv")[0]
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