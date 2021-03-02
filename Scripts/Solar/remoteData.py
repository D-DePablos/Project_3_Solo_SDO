import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.map import Map
from sunpy.time import TimeRange
from datetime import datetime
import numpy as np
import pandas as pd
from os import makedirs

PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/remoteData/"


class RemoteManager:
    """
    General Functions for AIA and GONG
    This class loads up all of the SDOAIA files and allows for selection of specific ones
    self.filesdf = pd dataframe with time as extracted from AIA filename as well as filename
    Has function that finds closest AIA observation to a given date and throws exception if too far
    """
    def __init__(self,
                 times=("2020/05/25", "2020/05/27"),
                 cadence=1 * u.hour,
                 instrument=None,
                 path=PATH):

        self.attrsTime = a.Time(TimeRange(times))
        self.cadence = a.Sample(cadence)
        self.instrument = a.Instrument(instrument)
        self.path = path
        self.dfiles = None

    @staticmethod
    def findClosestFile(dfiles, dateFind=datetime(2020, 5, 27, 1)):
        idx = dfiles.index[dfiles.index.get_loc(dateFind, method='nearest')]
        print(f"{dfiles.loc[idx].name - dateFind} difference to relevant date")
        return dfiles.loc[idx]

    @staticmethod
    def plot(map, title, **kwargs):
        """
        Pass a map, optionally seeds and field lines
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1, projection=map)
        map.plot(ax)

        if "flines" in kwargs:
            for fline in kwargs["flines"]:
                ax.plot_coord(
                    fline.coords,
                    color='blue',
                    linewidth=1,
                    alpha=0.4,
                )
                ax.plot_coord(fline.solar_footpoint,
                              color="red",
                              marker="x",
                              linewidth=0,
                              markersize=5)
            title = title + "| flines (blue)"

        ax.set_title(title)
        plt.show()


class GONGManager(RemoteManager):
    """
    Simple function that allows to get the GONG synoptic map from Sunpy as requested
    """
    def __init__(self,
                 times=("2020/05/27", "2020/05/27 00:20"),
                 instrument="GONG",
                 gongPath=f"{PATH}GONG/"):
        super().__init__(times=times, instrument=instrument, path=gongPath)
        makedirs(self.path, exist_ok=True)

    def downloadData(self):
        results = Fido.search(self.attrsTime, self.instrument)
        self.file = Fido.fetch(results[0], path=self.path)  # Fetching only one
        self.map = Map(self.file)
        self.map = Map(self.map.data - np.mean(self.map.data), self.map.meta)
        return self.map


class SDOAIAManager(RemoteManager):
    """
    Simple function to get SDOAIA data from sunpy into a map
    TODO : Add AIA plotting
    """
    def __init__(self,
                 times,
                 cadence=1 * u.hour,
                 instrument="AIA",
                 wavelength=193 * u.Angstrom,
                 aiaPath=f"{PATH}AIA/"):
        super().__init__(times=times,
                         cadence=cadence,
                         instrument=instrument,
                         path=aiaPath)
        self.wavelength = a.Wavelength(wavelength)
        self.path = f"{self.path}{int(wavelength.value)}/"
        makedirs(self.path, exist_ok=True)

    def downloadData(self):
        results = Fido.search(self.attrsTime, self.instrument, self.cadence,
                              self.wavelength)
        files = sorted(Fido.fetch(results, path=self.path))
        fileTime = [
            datetime.strptime(file[-39:-20], "%Y_%m_%dt%H_%M_%S")
            for file in files
        ]
        self.dfiles = pd.DataFrame({"file": files}, index=fileTime)

    @staticmethod
    def aiaprep(fitsFile):
        map = Map(fitsFile)
        from aiapy.calibrate import register, update_pointing, normalize_exposure
        m_updated_pointing = update_pointing(map)
        m_registered = register(m_updated_pointing)
        m_normalized = normalize_exposure(m_registered)
        return m_normalized


if __name__ == "__main__":
    timesGONG = ("2020/05/27", "2020/05/27 00:20")
    timesAIA = ("2020/05/27", "2020/05/27 3:20")
    gong = GONGManager(times=timesGONG)
    gong.downloadData()
    sdoaia = SDOAIAManager(timesAIA)
    sdoaia.downloadData()
    sdoaia.findClosestFile(dfiles=sdoaia.dfiles)
