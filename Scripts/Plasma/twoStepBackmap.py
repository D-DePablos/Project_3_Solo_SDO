# Generate maps of the solar atmosphere overlaid with the field lines - Figure 4.1 in the thesis

BASE_PATH = "/Users/ddp/Documents/PhD/solo_sdo/"
from sys import path

path.append(f"{BASE_PATH}Scripts/")
from os import makedirs
from os.path import exists

# from sunpy.time.timerange import TimeRange
from EMDComparison.LcurveSolOEMD import BASE_PATH
from Plasma.SoloData import SoloManager
from Solar.remoteData import SDOAIAManager, GONGManager
from Solar.pfssSolver import PFSSSolver
import pandas as pd
from astropy import units as u
from Solar.Proj_1_Imports import Maptools as Mtool

from datetime import datetime, timedelta


# Get SolO data
ssRadius = 2.5
# accFactor = 4 / 3
accFactor = 1
MARGINAIA, MARGINFLINE = 0.2, 0.5
FORCE = True  # Whether to force creation of new images
timeRangeFreq = "30min"  # Default 30min
Region = {
    "FOV_lon_lat": [30, 27, 35, 35],
    "extractRegions": [11, 12, 13, 16, 17, 18, 21, 22, 23],
    "customRegion": (11, 13, 21, 23),
}
# The first AIA file is aia_lev1_193a_2020_05_27t20_57_04_85z_image_lev1.fits 27.05.2020 20:57:04
# The last AIA file is aia_lev1_193a_2020_05_28t13_59_04_84z_image_lev1.fits 28.05.2020 13:59:04

# Start and end times for SolO
startSolo = datetime(2020, 5, 30, 12, 0)
endSolo = datetime(2020, 6, 1, 13, 30, 0)
# startSolo = datetime(2020, 7, 7, 12)
# endSolo = datetime(2020, 7, 8)

# Set AIA date, possibly change it within loop?
# dateAIA = gongTimeStart  # - timedelta(days=1)

# Choose HI. LO. Or measured for vSW
solo = SoloManager(times=(startSolo, endSolo), objCad=3600)
SOLOHI, SOLOLO, MEAN, MEASURED = (
    int(solo.df["V_R"].max()),
    int(solo.df["V_R"].min()),
    int(solo.df["V_R"].mean()),
    None,
)
CUSTOMSPEED = MEAN
aia = SDOAIAManager(
    times=("2020/5/27 20:30", "2020/5/28 14:00"), cadence=0.5 * u.hour
)  # Default times to download and prep. AIA files

# Extract orbit and backmap to Sun
solo.extractSolOrbit()
solo.backmapSun(accelerated=accFactor, customSpeed=CUSTOMSPEED)
# print(solo.SSFootPoints)

# Overplot Field lines on SDOAIA map
# used to be - timedelta(days=2, hours =13) - To line up with 23:00
time_1 = startSolo - timedelta(days=3)
time_2 = time_1 + timedelta(hours=60)
aia.downloadData()
timeRange = [d for d in pd.date_range(time_1, time_2, freq=timeRangeFreq)]

# DONE: Add plot of the solo data
if CUSTOMSPEED == None:
    spPath = ""
    savePath = f"{BASE_PATH}unsafe/logsCoords/{accFactor:.2f}/{spPath}"
    solo.plot(savePath=savePath)  # Plot only if not using custom speed

else:
    spPath = f"Speed_{CUSTOMSPEED}kms/"
    savePath = f"{BASE_PATH}unsafe/logsCoords/{accFactor:.2f}/{spPath}"

makedirs(savePath, exist_ok=True)

for index, date_AIA in enumerate(timeRange):
    if exists(f"{savePath}{index:03d}.png") and not FORCE:
        print(
            f"Image {date_AIA} at {index:03d} already exists, set Force to True to recreate"
        )

    else:
        # derive a PFSS for each AIA date
        gongTimeStart = date_AIA
        gongTimeEnd = gongTimeStart + timedelta(hours=13)
        gongTime = (
            f"{gongTimeStart.year}/{gongTimeStart.month}/{gongTimeStart.day} {gongTimeStart.hour}:{gongTimeStart.minute}",
            f"{gongTimeEnd.year}/{gongTimeEnd.month}/{gongTimeEnd.day} {gongTimeEnd.hour}:{gongTimeEnd.minute}",
        )
        gong = GONGManager(times=gongTime)
        gongMap = gong.downloadData()
        pfss = PFSSSolver(gongMap, SSradius=ssRadius)
        pfss.createSeedsFromList(solo.SSFootPoints)
        seedTimes = pfss.traceFlines(seedtimes=pfss.dfseeds.index)

        # pfss.plotMG(seeds=pfss.seeds, flines=pfss.flines)
        print(date_AIA)
        aiaFile, timeDiff = aia.findClosestFile(
            dfiles=aia.dfiles, dateFind=date_AIA, margin=MARGINAIA
        )
        aiaMap = aia.aiaprep(aiaFile.file)
        if index == 0:
            refMap = aiaMap
            refRoi = Mtool.ROI_defCoords(Region["FOV_lon_lat"]).copy()
            # Can try to calculate one time instead of multiple
        aia.plot(
            aiaMap,
            title=f"SDO AIA 193 - {aiaMap.date.datetime.strftime('%m-%d %H:%M')} \n  GONG White Light - {gongMap.date.datetime.strftime('%m-%d %H:%M')}",
            margin=MARGINFLINE,
            flines=pfss.flines,
            flineTimes=seedTimes,
            flineOriginTime=solo.sourcePointsSPC,
            savePath=savePath,
            index=index,
            force=FORCE,
            Regions=Region,
            refMap=refMap,
            ref_large_roi=refRoi,
            show=False,
            # show=True,
        )
