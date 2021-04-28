from os import makedirs
from sunpy.time.timerange import TimeRange
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
MARGINAIA, MARGINFLINE = 1, 1

startSolo = datetime(2020, 5, 30, 12)
endSolo = datetime(2020, 6, 1)
# startSolo = datetime(2020, 7, 7, 12)
# endSolo = datetime(2020, 7, 8)

gongTimeStart = startSolo - timedelta(days=3)
gongTimeEnd = gongTimeStart + timedelta(minutes=20)
gongTime = (
    f"{gongTimeStart.year}/{gongTimeStart.month}/{gongTimeStart.day} {gongTimeStart.hour}:{gongTimeStart.minute}",
    f"{gongTimeEnd.year}/{gongTimeEnd.month}/{gongTimeEnd.day} {gongTimeEnd.hour}:{gongTimeEnd.minute}"
)

# Set AIA date, possibly change it within loop?
# dateAIA = gongTimeStart  # - timedelta(days=1)

# Choose HI. LO. Or measured for vSW
SOLOHI, SOLOLO, MEASURED = (285, 255, None)
CUSTOMSPEED = SOLOHI
solo = SoloManager(times=(startSolo, endSolo), objCad=3600)
gong = GONGManager(times=gongTime)
aia = SDOAIAManager(
    times=("2020/5/26",
           "2020/5/28"))  # Default times to download and prep. AIA files

# Extract orbit and backmap to Sun
solo.extractSolOrbit()
solo.backmapSun(accelerated=accFactor, customSpeed=CUSTOMSPEED)
# print(solo.SSFootPoints)

# Start up PFSS and track field lines
gongMap = gong.downloadData()
pfss = PFSSSolver(gongMap, SSradius=ssRadius)
pfss.createSeedsFromList(solo.SSFootPoints)
seedTimes = pfss.traceFlines(seedtimes=pfss.dfseeds.index)
# pfss.plotMG(seeds=pfss.seeds, flines=pfss.flines)

# Overplot Field lines on SDOAIA map
aia.downloadData()
time_1 = (gongTimeStart - timedelta(hours=24))
time_2 = (gongTimeStart + timedelta(hours=24))
timeRange = [d for d in pd.date_range(time_1, time_2, freq="1h")]

#DONE: Add plot of the solo data
if CUSTOMSPEED == None:
    spPath = ""
    savePath = f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/logsCoords/{accFactor:.2f}/{spPath}"
    solo.plot(savePath=savePath)  # Plot only if not using custom speed

else:
    spPath = f"Speed_{CUSTOMSPEED}kms/"
    savePath = f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/logsCoords/{accFactor:.2f}/{spPath}"

makedirs(savePath, exist_ok=True)

for index, date_AIA in enumerate(timeRange):
    aiaFile, timeDiff = aia.findClosestFile(dfiles=aia.dfiles,
                                            dateFind=date_AIA,
                                            margin=MARGINAIA)
    aiaMap = aia.aiaprep(aiaFile.file)
    aia.plot(
        aiaMap,
        title=f"SDO AIA 93 {aiaMap.date.datetime.strftime('%m-%d %H:%M')}",
        margin=MARGINFLINE,
        flines=pfss.flines,
        flineTimes=seedTimes,
        flineOriginTime=solo.sourcePointsSPC,
        savePath=savePath,
        index=index,
    )
