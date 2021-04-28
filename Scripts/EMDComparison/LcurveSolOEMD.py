BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"

from sys import path

path.append(f"{BASE_PATH}Scripts/")
"""Main routine to compare remote and in-situ observations"""
from os import makedirs
from signalHelpers import compareTS
import numpy as np
from datetime import datetime, timedelta
from Solar.LcurveData import LcurveManager
from Plasma.SoloData import SoloManager


# Import the following functions into the AnySpacecraft_data script
def compareLCurvesToSolo(
    selfObj_Short,
    otherObj_Long,
    selfVars,
    otherVars,
    selfTimes,
    otherTimes,
    selfName,
    otherName,
    selfCad,
    otherCad,
    objDirExt,
    expectedLocationList=False,
    PeriodMinMax=[1, 180],
    filterPeriods=False,
    delete=True,
    showFig=True,
    renormalize=False,
    DETREND_BOX_WIDTH=None,
):
    """
    Feed in the PSP Spacecraft and SolOSpc object
    Self is expected to be solar orbiter
    Other is expected to be the lcurves
    """
    # Set header of directories
    general_directory = f"{BASE_PATH}unsafe/EMD_Data/"
    makedirs(general_directory, exist_ok=True)

    ### Directory structure
    # Specific folder to have all extracted datasets and plots
    mainDir = f"{general_directory}{objDirExt}/"
    makedirs(mainDir, exist_ok=True)

    # Set the Self and Other dataframe to those within the Spacecraft object
    dfSelf = selfObj_Short.df[selfVars]
    dfSelf.columns = [f"{selfName}_{i}"
                      for i in selfVars]  # Rename the columns

    dfOther = otherObj_Long.df[otherVars]
    dfOther.columns = [f"{otherName}_{i}" for i in otherVars]

    # Cut down the self and other dataseries
    dfSelf = dfSelf[selfTimes[0]:selfTimes[1]]
    dfOther = dfOther[otherTimes[0]:otherTimes[1]]
    cadSelf = selfCad
    cadOther = otherCad

    compareTS(
        dfSelf,
        dfOther,
        cadSelf,
        cadOther,
        labelOther=otherName,
        winDispList=[60],
        corrThrPlotList=np.arange(0.65, 1, 0.05),
        PeriodMinMax=PeriodMinMax,  # TODO : Figure out best period
        filterPeriods=filterPeriods,
        savePath=mainDir,
        useRealTime=True,
        expectedLocationList=expectedLocationList,
        detrend_box_width=DETREND_BOX_WIDTH,
        delete=delete,
        showFig=showFig,
        renormalize=renormalize,
    )


def extractDatetimePairs(Case, soloTesting=False):
    """
    For a given dictionary, create a list of datetime pairs for SolO and AIA
    """

    aiaMargin = Case["margin"]
    aiaTimes, soloTimes = [], []

    for margin in range(-aiaMargin, aiaMargin + 1):
        aiaStart = Case["AIATimes"] + timedelta(hours=margin)
        aiaEnd = aiaStart + timedelta(hours=1)
        aiaTimes.append((aiaStart, aiaEnd))

    if soloTesting:
        for margin in range(0, Case["soloDurn"]):
            soloStart = Case["soloTimes"] + timedelta(hours=margin)
            soloEnd = soloStart + timedelta(hours=Case["soloDurn"])
            soloTimes.append((soloStart, soloEnd))

    else:
        soloStart = Case["soloTimes"]
        soloEnd = soloStart + timedelta(hours=Case["soloDurn"])
        soloTimes.append((soloStart, soloEnd))

    return aiaTimes, soloTimes


if __name__ == "__main__":
    # import matplotlib.pyplot as plt
    from astropy import constants as const
    from astropy import units as u

    objCad = 60  # Objective cadence in seconds for comparisons

    DELETE = False
    SHOWFIG = False
    FILTERP = True
    PERIODMINMAX = [3, 20]
    WAVELENGTH = 94
    SOLOHI, SOLOLO, MEAN = 285, 255, 271

    # Solar Orbiter Data requires start, end
    start = datetime(2020, 5, 30, 23)
    end = datetime(2020, 6, 1)
    solo = SoloManager(
        times=(start, end),
        objCad=objCad,
        cdfPath=
        "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/",
    )
    solo.df = solo.df.interpolate()  # Fill gaps
    """
                                V_R        V_T       V_N          N         T
        2020-05-30 12:00:00  274.279297 -20.320356 -5.449193  43.087570  1.601013
        2020-05-30 12:01:00  274.317902 -20.294191 -5.402593  43.328201  1.592807
        2020-05-30 12:02:00  274.356537 -20.268026 -5.355993  43.568836  1.584601
        2020-05-30 12:03:00  274.395142 -20.241861 -5.309393  43.809467  1.576395
        2020-05-30 12:04:00  274.433746 -20.215696 -5.262794  44.050102  1.568189
        ...                         ...        ...       ...        ...       ...
        2020-05-31 22:56:00  273.647736 -45.103485 -4.894425  59.370270  2.708913
        2020-05-31 22:57:00  273.763336 -44.869247 -4.896916  59.102261  2.712214
        2020-05-31 22:58:00  273.878967 -44.635010 -4.899406  58.834251  2.715516
        2020-05-31 22:59:00  273.994568 -44.400772 -4.901897  58.566242  2.718817
        2020-05-31 23:00:00  274.110168 -44.166534 -4.904387  58.298233  2.722118
    """

    # Calculate mass flux
    Vx = (solo.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
    mp = const.m_p
    N = (solo.df["N"].values * (u.cm**(-3))).to(u.m**(-3))
    solo.df["Mf"] = (N * mp * Vx).value

    # Variables in situ
    soloVars = ["N", "T", "V_R", "Mf"]

    # Light Curves
    lc = LcurveManager(
        csvPath="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/",
        objCad=objCad,
        wavelength=WAVELENGTH,
    )
    lc.df = lc.df.interpolate()  # Interpolate after forming lc object
    """
                                 16          17          18          21          22          23
        Time                                                                                       
        2020-05-27 00:00:00  329.801805  607.425191  613.393714  621.675743  720.785361  623.828448
        2020-05-27 00:01:00  326.327427  604.308078  614.007916  618.429190  718.380477  621.880857
        2020-05-27 00:02:00  325.820149  605.108373  614.363865  615.858402  717.268210  625.081203
        2020-05-27 00:03:00  324.287363  605.050359  614.706581  612.511004  717.981663  627.441994
        2020-05-27 00:04:00  324.318297  605.238100  615.577086  611.048792  717.956587  628.376957
        ...                         ...         ...         ...         ...         ...         ...
        2020-05-28 10:41:00  381.538028  618.639577  664.715820  637.554964  661.129843  574.473581
        2020-05-28 10:42:00  382.460766  617.707741  664.790350  637.610576  658.359670  572.943388
        2020-05-28 10:43:00  384.760020  618.726230  667.828106  639.802201  660.028347  571.286007
        2020-05-28 10:44:00  387.059273  619.744718  670.865861  641.993826  661.697023  569.628625
        2020-05-28 10:45:00  386.551898  623.982588  672.670828  641.362444  665.743955  570.076031
    """
    # lcRegs = ["16"]
    lcRegs = [
        "11", "12", "13", "16", "17", "18", "21", "22", "23",
        "Summary_regions_11:13_21:23"
    ]

    # Cadence setup
    selfCad = int(lc.df.index.freq.delta.total_seconds())
    otherCad = int(solo.df.index.freq.delta.total_seconds())

    # Set up either case that we need to splice dframes
    # These are fo
    soloCon = [datetime(2020, 5, 31, 11, 0), datetime(2020, 5, 31, 13, 0)]
    soloAcc = [datetime(2020, 5, 31, 17, 0), datetime(2020, 5, 31, 21, 0)]

    # Dictionaries to extract data from timeseries
    accCase = {
        "AIATimes":
        datetime(2020, 5, 27, 5, 0),
        "soloTimes":
        datetime(2020, 5, 31, 12, 0),
        "margin":
        1,
        "soloDurn":
        12,
        "relevantTimes": [
            {
                "start": soloAcc[0],
                "end": soloAcc[1],
                "color": "blue",
                "label": "Acc.",
            },
        ],
    }
    # soloTimes should be 12:00

    conCase = {
        "AIATimes":
        datetime(2020, 5, 28, 5, 0),
        "soloTimes":
        datetime(2020, 5, 31, 5, 0),
        "margin":
        1,
        "soloDurn":
        8,
        "relevantTimes": [
            {
                "start": datetime(2020, 5, 31, 10),
                "end": datetime(2020, 5, 31, 11),
                "color": "yellow",
                "label": "Con",
            },
        ],
        "label":
        f"{WAVELENGTH}_Con"
    }

    # Relevant times actually change. Could be good to show minimum and maximum SolO speeds?
    multiAIACase = {
        "AIATimes": datetime(2020, 5, 28, 2),
        "soloTimes": datetime(2020, 5, 31, 0, 0),
        "margin": 7,
        "soloDurn": 24,
        "relevantTimes": [],
        "label": f"{WAVELENGTH}_R",
    }

    # Select a case
    selectedCase = multiAIACase
    aiaTimesList, soloTimesList = extractDatetimePairs(selectedCase,
                                                       soloTesting=False)
    # Invert N timeseries
    # solo.df["N"] = solo.df["N"].iloc[::-1].values

    for aiaTimes in aiaTimesList:
        if len(soloTimesList) > 1:
            for soloTimes in soloTimesList:
                dirName = f"""AIA_D{aiaTimes[0].day}_H{aiaTimes[0].hour}:{aiaTimes[1].hour}___________SoloD{soloTimes[0].day}_H{soloTimes[0].hour}:D{soloTimes[0].day}_H{soloTimes[1].hour}"""

                compareLCurvesToSolo(
                    selfObj_Short=lc,
                    otherObj_Long=solo,
                    selfVars=lcRegs,
                    otherVars=soloVars,
                    selfTimes=aiaTimes,
                    otherTimes=soloTimes,
                    selfName=selectedCase["label"],
                    otherName="Solo",
                    selfCad=60,
                    otherCad=60,
                    objDirExt=dirName,
                    filterPeriods=FILTERP,
                    PeriodMinMax=PERIODMINMAX,
                    delete=DELETE,
                    showFig=SHOWFIG,
                    expectedLocationList=selectedCase["relevantTimes"],
                    renormalize=False,
                )
        else:
            soloTimes = soloTimesList[0]
            dirName = f"""AIA_D{aiaTimes[0].day}_H{aiaTimes[0].hour}:{aiaTimes[1].hour}___________SoloD{soloTimes[0].day}_H{soloTimes[0].hour}:D{soloTimes[1].day}_H{soloTimes[1].hour}"""

            compareLCurvesToSolo(
                selfObj_Short=lc,
                otherObj_Long=solo,
                selfVars=lcRegs,
                otherVars=soloVars,
                selfTimes=aiaTimes,
                otherTimes=soloTimes,
                selfName=selectedCase["label"],
                otherName="Solo",
                selfCad=60,
                otherCad=60,
                objDirExt=dirName,
                filterPeriods=FILTERP,
                PeriodMinMax=PERIODMINMAX,
                delete=DELETE,
                showFig=SHOWFIG,
                expectedLocationList=selectedCase["relevantTimes"],
                renormalize=False,
            )
