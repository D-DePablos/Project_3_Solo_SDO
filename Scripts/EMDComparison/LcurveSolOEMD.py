# Set up UNSAFE_EMD_DATA_PATH: global variable
from sys import path

BASE_PATH = "/Users/ddp/Documents/PhD/solo_sdo/"
path.append(f"{BASE_PATH}Scripts/")

from Plasma.SoloData import SoloManager
from collections import namedtuple
from os import makedirs
from EMDComparison.signalHelpers import compareTS, new_plot_format, plot_super_summary
import numpy as np
from datetime import datetime, timedelta
from Solar.LcurveData import LcurveManager
from astropy import constants as const
from astropy import units as u


"""Main routine to compare remote and in-situ observations"""

newCases = True
UNSAFE_EMD_DATA_PATH = f"{BASE_PATH}unsafe/EMD_Data/"
UNSAFE_EMD_DATA_PATH = (
    UNSAFE_EMD_DATA_PATH if not newCases else f"{UNSAFE_EMD_DATA_PATH}newCases_allSOLO/"
)
makedirs(UNSAFE_EMD_DATA_PATH, exist_ok=True)

# Set parameters here
objCad = 60  # Objective cadence in seconds for comparisons
WVLLIST = [94, 193, 211]
PERIODMINMAX = [3, 20]
SHOWFIG = False
FILTERP = True

# Plot all in-situ variables?
PLOT_ALL_TOGETHER = False

# Add residual to non-super summary?
ADDRESIDUAL = False

# Plot summary?
SUPER_SUMMARY_PLOT = True

# Accelerated cases?
accelerated = 1
# accelerated = 4 / 3

# Solar Orbiter Data requires start, end
start = datetime(2020, 5, 30, 12)
end = datetime(2020, 6, 1, 13, 32)

# Lcurve regions
# lcRegs = ["11", "12", "13", "16", "17", "18", "21", "22", "23", "11:13_21:23"]

# Do only 9 for square plots
lcRegs = ["11", "12", "13", "16", "17", "18", "21", "22", "23"]
# lcRegs = ["11"]

# Open the cases file
caseName = "accCases" if accelerated == 4 / 3 else "consCases"
caseName = "newCases_ALLSOLO" if newCases else caseName

with open(
    f"{BASE_PATH}Scripts/EMDComparison/pickleCases/{caseName}.pickle",
    "rb",
) as f:
    import pickle

    cases = pickle.load(f)

MARGINHOURSSOLO = cases[0]["MARGINHOURSSOLO"]


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

    # Directory structure
    # Specific folder to have all extracted datasets and plots
    mainDir = f"{UNSAFE_EMD_DATA_PATH}{objDirExt}/"
    makedirs(mainDir, exist_ok=True)

    # Set the Self and Other dataframe to those within the Spacecraft object
    dfRemote = remoteObj.df[remVars]
    dfRemote.columns = [f"{remoteName}_{i}" for i in remVars]  # Rename the columns

    dfInsitu = insituObj.df[insituVars]
    dfInsitu.columns = [f"{insituName}_{i}" for i in insituVars]

    # Cut down the self and other dataseries
    dfRemote = dfRemote[remoteTimes[0] : remoteTimes[1]]
    dfInsitu = dfInsitu[insituTimes[0] : insituTimes[1]]
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
        showFig=showFig,
        renormalize=renormalize,
        showSpeed=True,
        LOSPEED=soloHI,
        HISPEED=soloLO,
    )


def extractDiscreteExamples(Caselist, margin, AIAduration=1):
    """Extracts discrete AIA - insituObject pairs

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

        return {"start": _times[0], "end": _times[1], "color": _color, "label": _label}

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

        refLocations.append(_constructExpectedLocation(_times=(matchStart, matchEnd)))

    return aiaTimes, insituTimes, caseNames, refLocations


def first_DeriveAndPlotSeparately():
    # The main program uses each of the wavelengths separately
    for wv in WVLLIST:
        WAVELENGTH = wv

        insituObject = SoloManager(
            times=(start, end),
            objCad=objCad,
        )
        insituObject.df = insituObject.df.interpolate()  # Fill gaps
        # Velocities are modified with 4/3 factor. Gives slightly better idea
        soloHI, soloLO, MEAN = (
            int(insituObject.df["V_R"].max() / accelerated),
            int(insituObject.df["V_R"].min() / accelerated),
            int(insituObject.df["V_R"].mean() / accelerated),
        )

        # Calculate mass flux
        Vx = (insituObject.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
        mp = const.m_p
        N = (insituObject.df["N"].values * (u.cm ** (-3))).to(u.m ** (-3))
        insituObject.df["Mf"] = (N * mp * Vx).value

        # Variables in situ
        insituObjectVars = ["N", "T", "V_R", "Mf"]
        # insituObjectVars = ["V_R"]

        # Light Curves
        lc = LcurveManager(
            objCad=objCad,
            wavelength=WAVELENGTH,
        )
        lc.df = lc.df.interpolate()  # Interpolate after forming lc object

        # We set a margin around original obs.
        (
            aiaTimesList,
            soloTimesList,
            caseNamesList,
            refLocations,
        ) = extractDiscreteExamples(
            cases,
            margin=MARGINHOURSSOLO,
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
                showFig=SHOWFIG,
                expectedLocationList=[refLocations[index]],
                renormalize=False,
                soloHI=soloHI,
                soloLO=soloLO,
            )


def combinedPlot(superSummaryPlot=False):
    lcDic = {}
    for _wvl in WVLLIST:
        lcDic[f"{_wvl}"] = LcurveManager(
            objCad=objCad,
            wavelength=_wvl,
            csvPath=f"{BASE_PATH}sharedData/",
        )
        lcDic[f"{_wvl}"].df = lcDic[f"{_wvl}"].df.interpolate()
        try:
            del lcDic[f"{_wvl}"].df["Unnamed: 0"]
        except KeyError:
            pass

    # Solar Orbiter Data
    insituObject = SoloManager(
        times=(start, end),
        objCad=objCad,
    )
    insituObject.df = insituObject.df.interpolate()  # Fill gaps
    # Velocities are modified with 4/3 factor. Gives slightly better idea
    soloHI, soloLO, soloAVG = (
        int(insituObject.df["V_R"].max() / accelerated),
        int(insituObject.df["V_R"].min() / accelerated),
        int(insituObject.df["V_R"].mean() / accelerated),
    )

    # Calculate mass flux
    Vx = (insituObject.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
    mp = const.m_p
    N = (insituObject.df["N"].values * (u.cm ** (-3))).to(u.m ** (-3))
    insituObject.df["Mf"] = (N * mp * Vx).value

    # Variables in situ
    insituObjectVars = ["V_R", "Mf", "N", "T"]
    insituObject.df = insituObject.df[insituObjectVars]

    figName = "accelerated" if accelerated == 4 / 3 else "constant"

    # We set a margin around original obs.
    aiaTimesList, soloTimesList, caseNamesList, refLocations = extractDiscreteExamples(
        cases,
        margin=MARGINHOURSSOLO,
    )

    # When necessary to make summary of all summaries
    if superSummaryPlot:

        # Possibly this is not great
        soloStendTotal = (
            insituObject.df.index[0].to_pydatetime(),
            insituObject.df.index[-1].to_pydatetime(),
        )

        wvlList = WVLLIST

        allCases = []
        Casetuple = namedtuple("Case", ["dirExtension", "isStend_t", "rsStend_t"])
        for index, aiaTimes in enumerate(aiaTimesList):
            _isT = soloTimesList[index]
            dirExtension = f"{caseNamesList[index]}"
            allCases.append(
                Casetuple(dirExtension, (_isT[0], _isT[1]), (aiaTimes[0], aiaTimes[1]))
            )

        # Figure out whether to show yellow bar - DONE

        insituObject.df.columns = ["SolO_" + param for param in insituObject.df.columns]

        for insituParam in insituObject.df.columns:
            plot_super_summary(
                allCasesList=allCases,
                longSpan=soloStendTotal,
                wvlList=wvlList,
                insituParam=insituParam,
                regions=lcRegs,
                unsafeEMDDataPath=UNSAFE_EMD_DATA_PATH,
                period=PERIODMINMAX,
                SPCKernelName="solo",
                speedSuper=soloHI,
                speedSuperLow=soloLO,
                speedAVG=soloAVG,
                showFig=SHOWFIG,
                figName=figName,
            )

    else:
        for index, aiaTimes in enumerate(aiaTimesList):

            if aiaTimes[0] == datetime(2020, 5, 28, 1, 30):
                # Need to cut up dataframes
                isTimes = soloTimesList[index]
                dfInsituCut = insituObject.df[isTimes[0] : isTimes[1]]
                dfInsituCut = dfInsituCut[insituObjectVars]

                lcDicCut = {}
                for _wvl in lcDic:
                    lcDicCut[f"{_wvl}"] = (
                        lcDic[_wvl].df[aiaTimes[0] : aiaTimes[1]].copy()
                    )

                dirExtension = f"{caseNamesList[index]}"
                base_folder = f"{UNSAFE_EMD_DATA_PATH}{dirExtension}/"
                new_plot_format(
                    dfInsitu=dfInsituCut,
                    lcDic=lcDicCut,
                    regions=lcRegs,
                    base_folder=base_folder,
                    period=PERIODMINMAX,
                    addResidual=ADDRESIDUAL,
                    SPCKernelName="solo",
                    spcSpeeds=(soloLO, soloHI),
                    showFig=SHOWFIG,
                )


if __name__ == "__main__":
    if not PLOT_ALL_TOGETHER:
        first_DeriveAndPlotSeparately()

    else:
        combinedPlot(superSummaryPlot=SUPER_SUMMARY_PLOT)
