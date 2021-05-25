BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"
UNSAFE_EMD_DATA_PATH = f"{BASE_PATH}unsafe/EMD_Data/"

from sys import path

path.append(f"{BASE_PATH}Scripts/")
"""Main routine to compare remote and in-situ observations"""
from os import makedirs
from signalHelpers import compareTS, new_plot_format
import numpy as np
from datetime import datetime, timedelta
from Solar.LcurveData import LcurveManager
from Plasma.SoloData import SoloManager
from astropy import constants as const
from astropy import units as u

# Set parameters here
objCad = 60  # Objective cadence in seconds for comparisons
WVLLIST = [94, 193, 211]
DELETE = False
SHOWFIG = True
FILTERP = True
PLOT_ALL_TOGETHER = True
accelerated = 1
# accelerated = 4 / 3
PERIODMINMAX = [3, 20]

# Solar Orbiter Data requires start, end
start = datetime(2020, 5, 30, 12)
end = datetime(2020, 6, 2)

# Lcurve regions
lcRegs = ["11", "12", "13", "16", "17", "18", "21", "22", "23", "11:13_21:23"]

# lcRegs = ["11:13_21:23"]


# Import the following functions into the AnySpacecraft_data script
def compareLcurvesToSolO(
    remoteObj,
    insituObj,
    remVars,
    insituVars,
    remoteTimes,
    insituTimes,
    remoteName,
    insituName,
    remoteCad,
    insituCad,
    objDirExt,
    expectedLocationList=False,
    PeriodMinMax=[1, 20],
    filterPeriods=False,
    delete=True,
    showFig=True,
    renormalize=False,
    DETREND_BOX_WIDTH=None,
    soloHI=None,
    soloLO=None,
):
    """
    Feed in the PSP Spacecraft and insituObjectSpc object
    Self is expected to be solar orbiter
    Other is expected to be the lcurves
    """

    assert soloHI != None, "No High speed set"
    assert soloLO != None, "No Low speed set"

    # Set header of directories
    makedirs(UNSAFE_EMD_DATA_PATH, exist_ok=True)

    ### Directory structure
    # Specific folder to have all extracted datasets and plots
    mainDir = f"{UNSAFE_EMD_DATA_PATH}{objDirExt}/"
    makedirs(mainDir, exist_ok=True)

    # Set the Self and Other dataframe to those within the Spacecraft object
    dfRemote = remoteObj.df[remVars]
    dfRemote.columns = [f"{remoteName}_{i}"
                        for i in remVars]  # Rename the columns

    dfInsitu = insituObj.df[insituVars]
    dfInsitu.columns = [f"{insituName}_{i}" for i in insituVars]

    # Cut down the self and other dataseries
    dfRemote = dfRemote[remoteTimes[0]:remoteTimes[1]]
    dfInsitu = dfInsitu[insituTimes[0]:insituTimes[1]]
    cadSelf = remoteCad
    cadOther = insituCad

    compareTS(
        dfRemote,
        dfInsitu,
        cadSelf,
        cadOther,
        labelOther=insituName,
        winDispList=[60],
        corrThrPlotList=np.arange(0.65, 1, 0.05),
        PeriodMinMax=PeriodMinMax,
        filterPeriods=filterPeriods,
        savePath=mainDir,
        useRealTime=True,
        expectedLocationList=expectedLocationList,
        detrend_box_width=DETREND_BOX_WIDTH,
        delete=delete,
        showFig=showFig,
        renormalize=renormalize,
        showSpeed=True,
        LOSPEED=soloHI,
        HISPEED=soloLO,
    )


def extractDiscreteExamples(Caselist, margin, AIAduration=1):
    """Extracts disxrete AIA - insituObject pairs

    Args:
        Caselist ([type]): [description]
        margin ([type]): [description]
        AIAduration (int, optional): [description]. Defaults to 1.
        noColumns (bool, optional): Whether to skip plotting backmapped time columns
    """
    def _constructExpectedLocation(_times, _color="blue", _label="BBMatch"):
        """Construct the Expected Location dic

        Args:
            _times ([type]): Tuple of two datetimes, start and end
            _color (str, optional): [description]. Defaults to "orange".
            label (str, optional): [description]. Defaults to "BBMatch".

        Returns:
            Dictionary with proper formatting
        """

        return {
            "start": _times[0],
            "end": _times[1],
            "color": _color,
            "label": _label
        }

    aiaTimes = []
    matchTimes = []
    insituTimes = []
    caseNames = []
    refLocations = []

    # Open each of the list dictionaries
    for case in Caselist:
        # Get the AIA start and end
        aiaStart = case["aiaTime"]
        aiaEnd = aiaStart + timedelta(hours=AIAduration)
        aiaTimes.append((aiaStart, aiaEnd))

        # Get the match, which is used for reference later
        matchStart = case["matchTime"]

        # insituObjectDurn gives reference to amount of insituObject datapoints that match
        matchEnd = matchStart + timedelta(hours=case["soloDurn"])
        matchAvg = matchStart + (matchEnd - matchStart) / 2
        matchTimes.append((matchStart, matchEnd))

        # Get the Solar Orbiter measurements
        insituObjectStart = matchAvg - timedelta(hours=margin)
        insituObjectEnd = matchAvg + timedelta(hours=margin)
        insituTimes.append((insituObjectStart, insituObjectEnd))

        # Get the specific case Name
        caseNames.append(case["caseName"])

        refLocations.append(
            _constructExpectedLocation(_times=(matchStart, matchEnd)))

    return aiaTimes, insituTimes, caseNames, refLocations


def deriveAndPlotSeparately():
    # The main program uses each of the wavelengths separately
    for wv in WVLLIST:
        WAVELENGTH = wv

        insituObject = SoloManager(
            times=(start, end),
            objCad=objCad,
            cdfPath=
            "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/",
        )
        insituObject.df = insituObject.df.interpolate()  # Fill gaps
        # Velocities are modified with 4/3 factor. Gives slightly better idea
        soloHI, soloLO, MEAN = (int(
            insituObject.df["V_R"].max() /
            accelerated), int(insituObject.df["V_R"].min() / accelerated),
                                int(insituObject.df["V_R"].mean() /
                                    accelerated))

        # Calculate mass flux
        Vx = (insituObject.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
        mp = const.m_p
        N = (insituObject.df["N"].values * (u.cm**(-3))).to(u.m**(-3))
        insituObject.df["Mf"] = (N * mp * Vx).value

        # Variables in situ
        insituObjectVars = ["N", "T", "V_R", "Mf"]
        # insituObjectVars = ["V_R"]

        # Light Curves
        lc = LcurveManager(
            csvPath=
            "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/",
            objCad=objCad,
            wavelength=WAVELENGTH,
        )
        lc.df = lc.df.interpolate()  # Interpolate after forming lc object

        # Open the cases file
        caseName = "accCases" if accelerated == 4 / 3 else "consCases"
        with open(
                f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/EMDComparison/{caseName}.pickle",
                "rb") as f:
            import pickle
            cases = pickle.load(f)

        # We set a margin around original obs.
        aiaTimesList, soloTimesList, caseNamesList, refLocations = extractDiscreteExamples(
            cases,
            margin=12,
        )

        for index, aiaTimes in enumerate(aiaTimesList):
            dirName = f"""{caseNamesList[index]}"""
            compareLcurvesToSolO(
                remoteObj=lc,
                insituObj=insituObject,
                remVars=lcRegs,
                insituVars=insituObjectVars,
                remoteTimes=aiaTimes,
                insituTimes=soloTimesList[index],
                remoteName=f"{WAVELENGTH}",
                insituName="SolO",
                remoteCad=objCad,
                insituCad=objCad,
                objDirExt=dirName,
                filterPeriods=FILTERP,
                PeriodMinMax=PERIODMINMAX,
                delete=DELETE,
                showFig=SHOWFIG,
                expectedLocationList=[refLocations[index]],
                renormalize=False,
                soloHI=soloHI,
                soloLO=soloLO,
            )


def combinedPlot():
    lcDic = {}
    for _wvl in WVLLIST:
        lcDic[f"{_wvl}"] = LcurveManager(objCad=objCad, wavelength=_wvl)
        lcDic[f"{_wvl}"].df = lcDic[f"{_wvl}"].df.interpolate()
        del lcDic[f"{_wvl}"].df["Unnamed: 0"]

    # Solar Orbiter Data
    insituObject = SoloManager(
        times=(start, end),
        objCad=objCad,
        cdfPath=
        "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/",
    )
    insituObject.df = insituObject.df.interpolate()  # Fill gaps
    # Velocities are modified with 4/3 factor. Gives slightly better idea
    soloHI, soloLO = (
        int(insituObject.df["V_R"].max() / accelerated),
        int(insituObject.df["V_R"].min() / accelerated),
    )

    # Calculate mass flux
    Vx = (insituObject.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
    mp = const.m_p
    N = (insituObject.df["N"].values * (u.cm**(-3))).to(u.m**(-3))
    insituObject.df["Mf"] = (N * mp * Vx).value

    # Variables in situ
    insituObjectVars = ["V_R", "Mf", "N", "T"]
    # insituObjectVars = ["V_R"]

    # Open the cases file
    caseName = "accCases" if accelerated == 4 / 3 else "consCases"
    with open(
            f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/EMDComparison/{caseName}.pickle",
            "rb") as f:
        import pickle
        cases = pickle.load(f)

    # We set a margin around original obs.
    aiaTimesList, soloTimesList, caseNamesList, refLocations = extractDiscreteExamples(
        cases,
        margin=12,
    )

    for index, aiaTimes in enumerate(aiaTimesList):
        # Need to cut up dataframes
        isTimes = soloTimesList[index]
        dfInsituCut = insituObject.df[isTimes[0]:isTimes[1]]
        dfInsituCut = dfInsituCut[insituObjectVars]

        lcDicCut = {}
        for _wvl in lcDic:
            lcDicCut[f"{_wvl}"] = lcDic[_wvl].df[aiaTimes[0]:aiaTimes[1]]

        dirExtension = f"{caseNamesList[index]}"
        base_folder = f"{UNSAFE_EMD_DATA_PATH}{dirExtension}/"
        new_plot_format(dfInsitu=dfInsituCut,
                        lcDic=lcDicCut,
                        regions=lcDicCut[list(lcDicCut.keys())[0]].columns,
                        base_folder=base_folder,
                        period=PERIODMINMAX,
                        addResidual=False,
                        SPCKernelName="solo",
                        spcSpeeds=(soloLO, soloHI),
                        showFig=SHOWFIG)


if __name__ == "__main__":
    if not PLOT_ALL_TOGETHER:
        deriveAndPlotSeparately()

    else:
        combinedPlot()
