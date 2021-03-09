import cdflib
from glob import glob
import pandas as pd
import numpy as np
from datetime import datetime
import astropy.units as u
import astropy.constants as const
from sunpy.coordinates import frames

SOLROTRATE = 24.47 / 86400
cdfEpoch = cdflib.cdfepoch()


class SoloManager:
    def __init__(
        self,
        times=(),
        objCad=3600,
        cdfPath="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/"
    ):
        """ 
        Variables within the cdf:
            [
                'Epoch', 
                'Half_interval', 
                'SCET', 
                'Info', 
                'validity', 
                'N', 
                'V_SRF', 
                'V_RTN', 
                'P_SRF', 
                'P_RTN', 
                'TxTyTz_SRF', 
                'TxTyTz_RTN', 
                'T'
            ] 
        """

        # Must download myself!
        for index, file in enumerate(sorted(glob(f"{cdfPath}*.cdf"))):
            cdf = cdflib.CDF(file)
            time = cdfEpoch.to_datetime(cdf["Epoch"])
            _df = pd.DataFrame({}, index=time)

            for i in ("V_RTN", "N", "T", "validity"):
                if i == "V_RTN":  # Allow for multidimensional
                    for n, arg in zip(
                        (0, 1, 2), ("_R", "_T", "_N")):  # R is radial velocity
                        _df[f"V{arg}"] = cdf[i][:, n]

                else:
                    _df[i] = cdf[i]

            # Join the dataframes only after the first instance
            if index == 0:
                _swe_df = _df
            else:
                _swe_df = _swe_df.append(_df)

        # Mask values outside of times to nan
        mask = (_swe_df.index > times[0]) & (_swe_df.index <= times[1])
        _swe_df = _swe_df.loc[mask]

        # NAN and interpolate values where validity < 3
        _swe_df[_swe_df["validity"] < 3] = np.nan
        del _swe_df["validity"]
        _swe_df = _swe_df.resample(f"{objCad}s").mean()
        self.df = _swe_df

    def extractSolOrbit(self):
        """
        Get the SPICE orbit data for SolO and put into coordsCarrington
        """
        import heliopy.data.spice as spicedata
        import heliopy.spice as spice

        spicedata.get_kernel("solo")
        solo = spice.Trajectory("Solar Orbiter")
        times = list(self.df.index.to_pydatetime())
        solo.generate_positions(times, 'Sun', 'IAU_SUN')  # Is in KM
        self.coordsCarrington = solo.coords.transform_to(
            frame=frames.HeliographicCarrington)

    def backmapSun(self, ssRadius=2.5 * const.R_sun.to(u.km), accelerated=1):
        """
        Backmap to Source Surface and save seeds as self.seeds 
        :param ssRadius: Source surface Radius in meters
        :param accelerated: Acceleration factor. Is 4/3 in classic 1 AU
        """
        from datetime import timedelta
        from astropy.coordinates import SkyCoord
        flineCoordsList = []

        for i, coord in enumerate(self.coordsCarrington):
            vSW = self.df["V_R"][i]
            lon, lat, r = coord.lon, coord.lat, coord.radius
            dt = accelerated * ((r - ssRadius) / (vSW * u.km / u.s)).value
            timeSPC = coord.obstime.value
            timeSun = timeSPC - timedelta(seconds=dt)

            rotationDeg = SOLROTRATE * r.value / vSW
            lonSS = lon + rotationDeg * u.deg
            lonSteps = np.arange(lon.value, lonSS.value, step=2) * u.deg
            latSteps = np.repeat(lat.value, len(lonSteps)) * u.deg
            stepRadius = (r.value - ssRadius.value) / len(lonSteps)
            rSteps = np.arange(ssRadius.value, r.value, step=stepRadius)[::-1]
            stepTime = (timeSPC - timeSun).total_seconds() / len(lonSteps)
            flTimes = [(timeSPC - timedelta(seconds=stepTime * lonstep))
                       for lonstep in range(len(lonSteps))]
            rSteps = rSteps[:-1] if (len(rSteps) > len(lonSteps)) else rSteps

            flineCoords = SkyCoord(lon=lonSteps,
                                   lat=latSteps,
                                   radius=rSteps * u.km,
                                   obstime=flTimes,
                                   frame="heliographic_carrington")
            flineCoordsList.append(flineCoords)

        self.flineCoordsList = flineCoordsList
        self.SSFootPoints = [coord[-1] for coord in flineCoordsList]
        self.sourcePointsSPC = self.coordsCarrington.obstime


if __name__ == "__main__":
    start = datetime(2020, 5, 30, 12)
    end = datetime(2020, 6, 1)
    solo = SoloManager(times=(start, end))
    solo.extractSolOrbit()
    solo.backmapSun()
    print(solo.SSFootPoints)
