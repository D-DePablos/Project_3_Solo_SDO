from Plasma.SoloData import SoloManager
from Solar.remoteData import SDOAIAManager, GONGManager
from Solar.pfssSolver import PFSSSolver
from astropy import constants as const
from astropy import units as u

from datetime import datetime, timedelta

# Get SolO data
ssRadius = 2.5
accFactor = 4 / 3

startSolo = datetime(2020, 5, 30, 12)
endSolo = datetime(2020, 6, 1)

gongTimeStart = startSolo - timedelta(days=3)
gongTimeEnd = gongTimeStart + timedelta(minutes=20)
gongTime = (
    f"{gongTimeStart.year}/{gongTimeStart.month}/{gongTimeStart.day} {gongTimeStart.hour}:{gongTimeStart.minute}",
    f"{gongTimeEnd.year}/{gongTimeEnd.month}/{gongTimeEnd.day} {gongTimeEnd.hour}:{gongTimeEnd.minute}"
)

solo = SoloManager(times=(startSolo, endSolo))
gong = GONGManager(times=gongTime)
# aia = SDOAIAManager(times=("2020/5/26", "2020/5/28"))
aia = SDOAIAManager(times=("2020/5/16", "2020/5/18"))

# Extract orbit and backmap to Sun
solo.extractSolOrbit()
solo.backmapSun(accelerated=accFactor)
# print(solo.SSFootPoints)

# Start up PFSS and track field lines
gongMap = gong.downloadData()
pfss = PFSSSolver(gongMap, SSradius=ssRadius)
pfss.createSeedsFromList(solo.SSFootPoints)
pfss.traceFlines()
# pfss.plotMG(seeds=pfss.seeds, flines=pfss.flines)

# Overplot Field lines on SDOAIA map
aia.downloadData()
aiaFile = aia.findClosestFile(dfiles=aia.dfiles,
                              dateFind=datetime(2020, 5, 16, 12))
aiaMap = aia.aiaprep(aiaFile.file)
aia.plot(aiaMap, title=f"SDO AIA 193 {aiaMap.date.value}", flines=pfss.flines)

# TODO Check that there is a coordinate conversion from Heliographic Carrington to SDOAIA coords
