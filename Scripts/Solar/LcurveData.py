from glob import glob
import pandas as pd


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

    def psf_test(
        self,
        csvFile="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/PSF_193.csv",
        objCad=60,
    ):
        """
        Docstring
        """
        Lc_csv = pd.read_csv(csvFile)
        self.df = Lc_csv
        self.df.index = pd.to_datetime(self.df["Time"])
        del self.df["Time"]

        for column in self.df.columns:
            # Clear out columns which are empty
            if self.df[column].mean() == 0:
                del self.df[column]

        self.df = self.df.resample(f"{objCad}s").mean()
        return self.df


def plotLcurves(lcurveDic, column="11:13_21:23", vspanStartEnd=(None, None)):
    """
    Plots lightcurves in fancy manner
    """

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(3, 1)

    axs[0].plot(lcurveDic["94"].df[f"{column}"])
    axs[1].plot(lcurveDic["193"].df[f"{column}"])
    axs[2].plot(lcurveDic["211"].df[f"{column}"])

    if vspanStartEnd[0] != None:
        for ax in axs:
            ax.vspan(vspanStartEnd)

    plt.show()


if __name__ == "__main__":
    lcDic = {}
    for WVL in [94, 193, 211]:
        lcDic[f"{WVL}"] = LcurveManager(wavelength=WVL)

    plotLcurves(lcDic)
