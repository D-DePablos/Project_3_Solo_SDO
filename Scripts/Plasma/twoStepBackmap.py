BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"
from sys import path

path.append(f"{BASE_PATH}Scripts/")
from os import makedirs
from sunpy.time.timerange import TimeRange
from EMDComparison.LcurveSolOEMD import BASE_PATH
from Plasma.SoloData import SoloManager
from Solar.remoteData import SDOAIAManager, GONGManager
from Solar.pfssSolver import PFSSSolver
from astropy import constants as const
from astropy import units as u
import pandas as pd

from datetime import datetime, timedelta
# Get SolO data
ssRadius = 2.5
# accFactor = 4 / 3
accFactor = 1
MARGINAIA, MARGINFLINE = 0.2, 0.5

startSolo = datetime(2020, 5, 30, 12, 0)
endSolo = datetime(2020, 6, 1, 13, 30, 0)
# startSolo = datetime(2020, 7, 7, 12)
# endSolo = datetime(2020, 7, 8)

# Set AIA date, possibly change it within loop?
# dateAIA = gongTimeStart  # - timedelta(days=1)

# Choose HI. LO. Or measured for vSW
solo = SoloManager(times=(startSolo, endSolo), objCad=3600)
SOLOHI, SOLOLO, MEAN, MEASURED = (int(solo.df["V_R"].max()),
                                  int(solo.df["V_R"].min()),
                                  int(solo.df["V_R"].mean()), None)
CUSTOMSPEED = MEAN
aia = SDOAIAManager(
    times=("2020/5/26",
           "2020/5/28 18:00"))  # Default times to download and prep. AIA files

# Extract orbit and backmap to Sun
solo.extractSolOrbit()
solo.backmapSun(accelerated=accFactor, customSpeed=CUSTOMSPEED)
# print(solo.SSFootPoints)

# Overplot Field lines on SDOAIA map
# used to be - timedelta(days=2, hours =13) - To line up with 23:00
time_1 = (startSolo - timedelta(days=3))
time_2 = (time_1 + timedelta(hours=60))
aia.downloadData()
timeRange = [d for d in pd.date_range(time_1, time_2, freq="30min")]

#DONE: Add plot of the solo data
if CUSTOMSPEED == None:
    spPath = ""
    savePath = f"{BASE_PATH}unsafe/logsCoords/{accFactor:.2f}/{spPath}"
    solo.plot(savePath=savePath)  # Plot only if not using custom speed

else:
    spPath = f"Speed_{CUSTOMSPEED}kms/"
    savePath = f"{BASE_PATH}unsafe/logsCoords/{accFactor:.2f}/{spPath}"

makedirs(savePath, exist_ok=True)

for index, date_AIA in enumerate(timeRange):

    # derive a PFSS for each AIA date
    gongTimeStart = date_AIA
    gongTimeEnd = gongTimeStart + timedelta(hours=13)
    gongTime = (
        f"{gongTimeStart.year}/{gongTimeStart.month}/{gongTimeStart.day} {gongTimeStart.hour}:{gongTimeStart.minute}",
        f"{gongTimeEnd.year}/{gongTimeEnd.month}/{gongTimeEnd.day} {gongTimeEnd.hour}:{gongTimeEnd.minute}"
    )
    gong = GONGManager(times=gongTime)
    gongMap = gong.downloadData()
    pfss = PFSSSolver(gongMap, SSradius=ssRadius)
    pfss.createSeedsFromList(solo.SSFootPoints)
    seedTimes = pfss.traceFlines(seedtimes=pfss.dfseeds.index)

    # pfss.plotMG(seeds=pfss.seeds, flines=pfss.flines)
    print(date_AIA)
    aiaFile, timeDiff = aia.findClosestFile(dfiles=aia.dfiles,
                                            dateFind=date_AIA,
                                            margin=MARGINAIA)
    aiaMap = aia.aiaprep(aiaFile.file)
    aia.plot(
        aiaMap,
        title=
        f"SDO AIA 193 - {aiaMap.date.datetime.strftime('%m-%d %H:%M')} \n  GONG White Light - {gongMap.date.datetime.strftime('%m-%d %H:%M')}",
        margin=MARGINFLINE,
        flines=pfss.flines,
        flineTimes=seedTimes,
        flineOriginTime=solo.sourcePointsSPC,
        savePath=savePath,
        index=index,
        # show=True,
    )
