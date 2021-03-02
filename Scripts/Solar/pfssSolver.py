# Astropy for physical constants, unit and coordinate handling
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sunpy.map import Map

import pfsspy
import pfsspy.tracing as tracing
from Solar.remoteData import GONGManager


class PFSSSolver:
    def __init__(self,
                 gongMap,
                 SSradius: float = 2.5,
                 nrho: int = 25,
                 tracer=tracing.FortranTracer()):
        """
        PFSS solution to solar magnetic field is calculated on a 3D grid (phi, s, rho).
        rho = ln(r), and r is the standard spherical radial coordinate.
        """
        self.gongMap = gongMap
        self.SSradius = SSradius
        self.SSradiusRsun = SSradius * const.R_sun
        self.nrho = nrho

        self.input = pfsspy.Input(gongMap, nrho, SSradius)
        self.output = pfsspy.pfss(self.input)
        self.tracer = tracer

    @staticmethod
    def createSeeds(lons, lats, radius, frame):
        return SkyCoord(lons, lats, radius, frame=frame)

    def createSeedsFromList(self, coordList):
        lons, lats, rads, obstimes = [], [], [], []
        for coord in coordList:
            assert (
                type(coord) == SkyCoord
            ), f"Please ensure you are using a list of SkyCoords. Type {type(coord)} not valid"

            lons.append(coord.lon)
            lats.append(coord.lat)
            rads.append(coord.radius)
            obstimes.append(coord.obstime)

        # use the default createSeeds call. Note that it allows for only one frame
        self.seeds = self.createSeeds(lons,
                                      lats,
                                      rads,
                                      frame=self.gongMap.coordinate_frame)

        # Store the time information for each seed on dfseeds
        # Seeds have a time but field lines do not
        self.dfseeds = pd.DataFrame({"seed": self.seeds}, index=obstimes)

    def seedMesh(self,
                 sinLat=(0.25, 0.5),
                 long=(60, 100),
                 rootPoints=5,
                 seedHeight=None):
        s = np.linspace(sinLat[0], sinLat[1], rootPoints)
        phi = np.linspace(long[0], long[1], rootPoints)
        s, phi = np.meshgrid(s, phi)

        # Unpack the lon and lat, add units
        lat = (np.arcsin(s) * u.rad).ravel()
        lon = (phi * u.deg).ravel()
        radius = self.SSradiusRsun if seedHeight == None else seedHeight

        seeds = self.createSeeds(lons=lon,
                                 lats=lat,
                                 radius=radius,
                                 frame=self.gongMap.coordinate_frame)

        self.seeds = seeds
        return seeds

    def traceFlines(self):
        self.flines = self.tracer.trace(self.seeds, self.output)

    def plotMG(self, **kwargs):
        """
        Plot the synoptic magnetogram and seeds or field lines if required
        """
        m = self.input.map
        fig = plt.figure()
        ax = plt.subplot(projection=m)
        m.plot()
        # plt.colorbar()
        title = f" {m.date.value} | Synoptic GONG magnetogram |"

        if "seeds" in kwargs:
            ax.plot_coord(kwargs["seeds"],
                          color='black',
                          marker='o',
                          linewidth=0,
                          markersize=5)
            title = title + "| seeds (black) "

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


if __name__ == "__main__":

    # Test that the
    gong = GONGManager()
    gongMap = gong.downloadData()
    SSradius = 3
    pfss = PFSSSolver(gongMap=gongMap, SSradius=SSradius)
    pfss.seedMesh()
    pfss.traceFlines()
    pfss.plotMG(flines=pfss.flines, seeds=pfss.seeds)
