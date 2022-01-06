from cgitb import small
import astropy.units as u
import sunpy
from sunpy.net import Fido, attrs as a
from sunpy.map import Map
from sunpy.time import TimeRange
from datetime import datetime
import numpy as np
import pandas as pd
from os import makedirs

import Solar.Proj_1_Imports.Maptools as Mtool

# Base path for Project 3
PATH = "/Users/ddp/Documents/PhD/solo_sdo/unsafe/remoteData/"


class RemoteManager:
    """
    General Functions for AIA and GONG
    This class loads up all of the SDOAIA files and allows for selection of specific ones
    self.filesdf = pd dataframe with time as extracted from AIA filename as well as filename
    Has function that finds closest AIA observation to a given date and throws exception if too far
    """

    def __init__(
        self,
        times=("2020/05/25", "2020/05/27"),
        cadence=1 * u.hour,
        instrument=None,
        path=PATH,
    ):

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
        idx = dfiles.index[dfiles.index.get_loc(dateFind, method="nearest")]
        timediff = np.abs((dfiles.loc[idx].name - dateFind).total_seconds() / 3600)
        if timediff >= margin:
            from datetime import timedelta

            timeStrStart = dateFind.strftime("%Y/%m/%d %H:%M")
            timeStrEnd = (dateFind + timedelta(hours=margin)).strftime("%Y/%m/%d %H:%M")
            self.attrsTime = a.Time(TimeRange(timeStrStart, timeStrEnd))
            self.downloadData(force=True)
            idx = dfiles.index[dfiles.index.get_loc(dateFind, method="nearest")]
            timediff = np.abs((dfiles.loc[idx].name - dateFind).total_seconds() / 3600)

        return dfiles.loc[idx], timediff

    @staticmethod
    def plot(map, title, margin=3, **kwargs):
        """
        Pass a map, optionally field lines and another map projected on top
        map: AIA or GONG map
        margin: Number of hours allowed between map and flineTimes


        KWARGS:
            extraMap: Map object to plot on top
            show: bool to show plot
            flines: What flines to plot
            flineTimes: What times are relevant for flines
            flineOriginTime: What time fline reaches spacecraft given speed
        """
        import matplotlib.pyplot as plt

        if "Regions" in kwargs and "refMap" in kwargs:
            refMap = kwargs["refMap"]  # Assuming already is a sunpy Map
            customRegion = kwargs["Regions"]["customRegion"]
            all_hpc = sunpy.map.all_coordinates_from_map(refMap)
            all_hgs = all_hpc.transform_to("heliographic_stonyhurst")

            if refMap.coordinate_frame.obstime != map.coordinate_frame.obstime:
                print("Rotating Sun...")
                large_roi = Mtool.calculateLons(
                    refMap, map, Coords=kwargs["ref_large_roi"]
                )

            else:
                large_roi = kwargs["ref_large_roi"]
            small_rois = Mtool.splitCoords(large_roi, nrowcol=5)

            # To show analysis regions
            customRegion = {
                "lon": (
                    small_rois[customRegion[0]][1]["lon"][0],
                    small_rois[customRegion[1]][1]["lon"][1],
                ),
                "lat": (
                    small_rois[customRegion[2]][1]["lat"][0],
                    small_rois[customRegion[1]][1]["lat"][1],
                ),
            }

            # Remove a sixth of the size
            customSegment = Mtool.getSegment(customRegion, all_hgs_coords=all_hgs)

            for roi_index, coords in small_rois:
                _cregion = {
                    "lon": (coords["lon"][0], coords["lon"][1]),
                    "lat": (coords["lat"][0], coords["lat"][1]),
                }
                _segment = Mtool.getSegment(_cregion, all_hgs_coords=all_hgs)

                borders = extractSegmentBorders(_segment, 500)
                for border in borders:
                    map.data[border] = 0

            # borders = extractSegmentBorders(customSegment)
            # for border in borders:
            #     map.data[border] = 0

        plt.figure(figsize=(9, 9))
        ax = plt.subplot(1, 1, 1, projection=map)
        map.plot(ax)
        if "extraMap" in kwargs:
            overlay = ax.get_coords_overlay(kwargs["extraMap"].coordinate_frame)
            overlay.grid()

        if (
            "flines" in kwargs
            and "flineTimes" in kwargs
            and "flineOriginTime" in kwargs
        ):
            with open(
                "/Users/ddp/Documents/PhD/solo_sdo/unsafe/logsCoords/logs.txt",
                "a",
            ) as f:
                f.write(f"\n### \nAIA MAP {map.date}\n")
                for time, fline, spcTime in zip(
                    kwargs["flineTimes"], kwargs["flines"], kwargs["flineOriginTime"]
                ):
                    if margin >= np.abs(
                        (map.date.datetime - time.value).total_seconds() / 3600
                    ):
                        # For each of the field lines, compare to the time of the original map

                        # Field line
                        ax.plot_coord(
                            fline.coords.helioprojective,
                            color="blue",
                            linewidth=5,
                            alpha=0.4,
                            label=f"{spcTime.value}",
                        )
                        f.write(
                            f"{fline.solar_footpoint.helioprojective.Tx} | {fline.solar_footpoint.helioprojective.Ty} | {fline.solar_footpoint.lon} \n"
                        )

                        # Footpoints  - White if in front of disk
                        color = (
                            "black"
                            if np.abs(
                                fline.solar_footpoint.heliographic_stonyhurst.lon.value
                            )
                            < 90
                            else "red"
                        )
                        ax.plot_coord(
                            fline.solar_footpoint.helioprojective,
                            color=color,
                            marker="x",
                            linewidth=5,
                            markersize=15,
                        )
                plt.legend()
                # title = title + f"| flines (SPC time) -> Margin {margin} hours vs AIA"
                f.close()

        plt.xlim(2200, 3800)
        plt.ylim(2200, 3800)
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(f"{kwargs['savePath']}{kwargs['index']:03d}.png")
        if "show" in kwargs:
            if kwargs["show"]:
                plt.show()
        plt.close()


class GONGManager(RemoteManager):
    """
    Simple function that allows to get the GONG synoptic map from Sunpy as requested
    """

    def __init__(
        self,
        times=("2020/05/27", "2020/05/27 00:20"),
        instrument="GONG",
        gongPath=f"{PATH}GONG/",
    ):
        super().__init__(times=times, instrument=instrument, path=gongPath)
        makedirs(self.path, exist_ok=True)

    def downloadData(self):
        results = Fido.search(self.attrsTime, self.instrument)

        closest = results["gong"][0]
        print(closest)
        self.file = Fido.fetch(closest, path=self.path)  # Fetching only one
        gongmap = Map(self.file)
        self.map = Map(gongmap.data - np.mean(gongmap.data), gongmap.meta)
        return self.map


class SDOAIAManager(RemoteManager):
    """
    Simple function to get SDOAIA data from sunpy into a map
    """

    def __init__(
        self,
        times,
        cadence=1 * u.hour,
        instrument="AIA",
        wavelength=193 * u.Angstrom,
        aiaPath=f"{PATH}AIA/",
    ):
        super().__init__(
            times=times, cadence=cadence, instrument=instrument, path=aiaPath
        )
        self.wavelength = a.Wavelength(wavelength)
        self.path = f"{self.path}{int(wavelength.value)}/"
        makedirs(self.path, exist_ok=True)

    def downloadData(self, force=False):
        results = None
        files = []

        if not force:
            try:
                from glob import glob

                files = sorted(glob(f"{self.path}*"))

            except FileNotFoundError:
                results = Fido.search(
                    self.attrsTime, self.instrument, self.cadence, self.wavelength
                )
                files = sorted(Fido.fetch(results, path=self.path))

                print(results)

        elif force:

            results = Fido.search(
                self.attrsTime, self.instrument, self.cadence, self.wavelength
            )

            files = sorted(Fido.fetch(results, path=self.path))

        assert len(files) > 0, "Fits files for given params. not found"
        fileTime = [
            datetime.strptime(file[-39:-20], "%Y_%m_%dt%H_%M_%S") for file in files
        ]
        self.dfiles = pd.DataFrame({"file": files}, index=fileTime)

    @staticmethod
    def aiaprep(fitsFile):
        map = Map(fitsFile)
        from aiapy.calibrate import register, update_pointing, normalize_exposure

        m_updated_pointing = update_pointing(map)

        # Quick hack to deal with broken data
        try:
            m_registered = register(m_updated_pointing)
        except ValueError as e:
            print(e)
            return m_updated_pointing

        m_normalized = normalize_exposure(m_registered)
        return m_normalized


def getSideBorder(lons, lats, margin: int, side: str):
    """
    Get the border of the side of the map
    lons: longitudes for a portion of the map
    lats: latitudes for left or right
    margin: how many different longitudes to use
    """
    lon_min = np.min(lons)
    lon_max = np.max(lons)

    assert lon_max - lon_min > margin, "Margin too large. Select more longitudes"

    resultLons, resultLats = [], []
    for i in range(margin):
        if side == "left":
            lon = lon_min + i
        elif side == "right":
            lon = lon_max - i

        for lat in lats:
            resultLons.append(lon)
            resultLats.append(lat)

    return resultLons, resultLats


def gettopBottom(
    Lons,
    leftLats,
    rightLats,
    margin: int,
):
    latMinRight, latMaxRight = min(rightLats), max(rightLats)
    latMinLeft, latMaxLeft = min(leftLats), max(leftLats)

    # For each of the longitudes in the map, get approximate lat
    lonsList = list(range(Lons[0], Lons[-1]))
    latsTop = np.linspace(latMaxRight, latMaxLeft, len(lonsList))
    latsBot = np.linspace(latMinRight, latMinLeft, len(lonsList))

    resultLatsTop, resultLatsBot = [], []
    resultLons = []
    for i in range(margin):
        latsTop = latsTop - i
        latsBot = latsBot + i
        for j, lon in enumerate(lonsList):
            resultLons.append(lon)
            resultLatsTop.append(latsTop[j])
            resultLatsBot.append(latsBot[j])

    return resultLons, resultLatsTop, resultLatsBot
    # Need to select top and bottom latitudes, from each border, and then get all lons


def extractSegmentBorders(seg, marginCut=1500):
    # marginCut = len(customSegment[0]) // 6
    # Extract a Margin from the large region
    x_cut = (
        seg[0][: marginCut * 15],
        seg[0][-marginCut * 15 :],
    )
    y_cut = (
        seg[1][:marginCut],
        seg[1][-marginCut:],
    )

    leftLons, leftLats = getSideBorder(x_cut[0], y_cut[0], 5, "left")
    rightLons, rightLats = getSideBorder(x_cut[1], y_cut[1], 5, "right")
    lons, latsTop, latsBot = gettopBottom(seg[0], rightLats, leftLats, 2)

    return (
        [leftLons, leftLats],
        [rightLons, rightLats],
        [np.int0(lons), np.int0(latsBot)],
        [np.int0(lons), np.int0(latsTop)],
    )


if __name__ == "__main__":
    #     timesGONG = ("2020/05/27", "2020/05/27 00:20")
    #     gong = GONGManager(times=timesGONG)
    #     gong.downloadData()

    # AIA Times
    timesAIA = ("2020/05/27", "2020/05/27 3:20")
    sdoaia = SDOAIAManager(timesAIA)
    sdoaia.downloadData()
    sdoaia.findClosestFile(dfiles=sdoaia.dfiles)
