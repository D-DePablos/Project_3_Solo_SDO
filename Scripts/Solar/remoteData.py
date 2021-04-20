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
        print(f"Using {self.path}")
        self.dfiles = None

    def findClosestFile(
            self,
            dfiles,
            dateFind=datetime(2020, 5, 27, 1),
            margin=3,
    ):
        idx = dfiles.index[dfiles.index.get_loc(dateFind, method='nearest')]
        timediff = np.abs(
            (dfiles.loc[idx].name - dateFind).total_seconds() / 3600)
        if timediff >= margin:
            from datetime import timedelta
            timeStrStart = dateFind.strftime("%Y/%m/%d %H:%M")
            timeStrEnd = (dateFind +
                          timedelta(hours=margin)).strftime("%Y/%m/%d %H:%M")
            self.attrsTime = a.Time(TimeRange(timeStrStart, timeStrEnd))
            self.downloadData(force=True)
            idx = dfiles.index[dfiles.index.get_loc(dateFind,
                                                    method='nearest')]
            timediff = np.abs(
                (dfiles.loc[idx].name - dateFind).total_seconds() / 3600)

        return dfiles.loc[idx], timediff

    @staticmethod
    def plot(map, title, margin=3, **kwargs):
        """
        Pass a map, optionally field lines and another map projected on top
        """
        import matplotlib.pyplot as plt
        plt.figure()
        ax = plt.subplot(1, 1, 1, projection=map)
        map.plot(ax)

        if "extraMap" in kwargs:
            overlay = ax.get_coords_overlay(
                kwargs["extraMap"].coordinate_frame)
            overlay.grid()

        if "flines" in kwargs and "flineTimes" in kwargs and "flineOriginTime" in kwargs:
            with open(
                    "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/logsCoords/logs.txt",
                    "a") as f:
                f.write(f"\n### \nAIA MAP {map.date}\n")
                for time, fline, spcTime in zip(kwargs["flineTimes"],
                                                kwargs["flines"],
                                                kwargs["flineOriginTime"]):
                    if margin >= np.abs(
                        (map.date.datetime - time.value).total_seconds() /
                            3600):
                        # For each of the field lines, compare to the time of the original map

                        ax.plot_coord(
                            fline.coords.helioprojective,
                            # color='blue',
                            linewidth=1,
                            alpha=0.4,
                            label=f"{spcTime.value}",
                        )
                        f.write(
                            f"{fline.solar_footpoint.helioprojective.Tx} | {fline.solar_footpoint.helioprojective.Ty} | {fline.solar_footpoint.lon} \n"
                        )
                        color = "white" if np.abs(
                            fline.solar_footpoint.heliographic_stonyhurst.lon.
                            value) < 90 else "red"
                        ax.plot_coord(fline.solar_footpoint.helioprojective,
                                      color=color,
                                      marker="x",
                                      linewidth=0,
                                      markersize=5)
                plt.legend()
                title = title + f"| flines (SPC time) -> Margin {margin} hours vs AIA"
                f.close()

        ax.set_title(title)
        plt.savefig(
            f"{kwargs['savePath']}AIA{int(map.wavelength.value)}{map.date.datetime.strftime('%Y-%m-%d_%H-%M')}.png"
        )
        plt.close()


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
        gongmap = Map(self.file)
        self.map = Map(gongmap.data - np.mean(gongmap.data), gongmap.meta)
        return self.map


class SDOAIAManager(RemoteManager):
    """
    Simple function to get SDOAIA data from sunpy into a map
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

    def downloadData(self, force=False):
        results = None
        if not force:
            try:
                from glob import glob
                files = sorted(glob(f"{self.path}*"))

            except FileNotFoundError:
                results = Fido.search(self.attrsTime, self.instrument,
                                      self.cadence, self.wavelength)
                files = sorted(Fido.fetch(results, path=self.path))

                print(results)

        elif force:

            results = Fido.search(self.attrsTime, self.instrument,
                                  self.cadence, self.wavelength)

            print(results)
        response = input(
            f"Would you like to download the above shown files into {self.path}?"
        )

        if response.lower() == "y":
            print(self.path)
            files = sorted(Fido.fetch(results, path=self.path))

        else:
            raise ValueError(
                f"{response} is not 'y' and files are not downloaded. Stopping"
            )

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
    #     timesGONG = ("2020/05/27", "2020/05/27 00:20")
    #     gong = GONGManager(times=timesGONG)
    #     gong.downloadData()

    # AIA Times
    timesAIA = ("2020/05/27", "2020/05/27 3:20")
    sdoaia = SDOAIAManager(timesAIA)
    sdoaia.downloadData()
    sdoaia.findClosestFile(dfiles=sdoaia.dfiles)
