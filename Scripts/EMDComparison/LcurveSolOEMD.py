BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"

from sys import path
from typing import Dict

path.append(f"{BASE_PATH}Scripts/")
"""Main routine to compare remote and in-situ observations"""
from os import makedirs
from signalHelpers import compareTS
import numpy as np
from datetime import datetime, timedelta
from Solar.LcurveData import LcurveManager
from Plasma.SoloData import SoloManager

# TODO: Add signal comparing structure from ISSI work


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
):
    """
    Feed in the PSP Spacecraft and SolOSpc object
    Self is expected to be solar orbiter
    Other is expected to be the lcurves
    """
    # Set header of directories
    general_directory = f"{BASE_PATH}unsafe/EMD_Data/"
    # TODO: Delete only the specific folder
    # if delete:
    #     from shutil import rmtree
    #     rmtree(general_directory, ignore_errors=True)
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
        corrThrPlotList=np.arange(0.4, 1, 0.05),
        PeriodMinMax=PeriodMinMax,  # TODO : Figure out best period
        filterPeriods=filterPeriods,
        savePath=mainDir,
        useRealTime=True,
        expectedLocationList=expectedLocationList,
        detrend_box_width=None,
        delete=delete,
        showFig=showFig,
    )


def extractDatetimePairs(Case: Dict):
    """
    For a given dictionary, create a list of datetime pairs for SolO and AIA
    """

    aiaMargin = Case["margin"]
    aiaTimes = []

    for margin in range(-aiaMargin, aiaMargin + 1):
        aiaStart = Case["AIATimes"] + timedelta(hours=margin)
        aiaEnd = aiaStart + timedelta(hours=Case["margin"])
        aiaTimes.append((aiaStart, aiaEnd))

    soloStart = Case["SolOTimes"]
    soloEnd = soloStart + timedelta(hours=Case["soloDurn"])
    soloTimes = (soloStart, soloEnd)

    return aiaTimes, soloTimes


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    objCad = 60
    # Solar Orbiter Data requires start, end
    start = datetime(2020, 5, 30, 12)
    end = datetime(2020, 6, 1)
    solo = SoloManager(
        times=(start, end),
        objCad=objCad,
        cdfPath=
        "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/")
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
    soloVars = ["N"]
    # otherVars = ["B_R", "B_T", "B_N", "N", "T"]

    lc = LcurveManager(
        csvPath="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/",
        objCad=objCad,
    )
    # TODO : Check that lightcurves are correctly derived
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
    lcRegs = ["16", "17", "18", "21", "22", "23"]

    # Set up either case that we need to splice dframes
    accCase = {
        "AIATimes": datetime(2020, 5, 27, 5, 0),
        "SolOTimes": datetime(2020, 5, 31, 12, 0),
        "margin": 1,
        "soloDurn": 10,
        "relevantTimes": False,
    }

    conCase = {
        "AIATimes": datetime(2020, 5, 28, 5, 0),
        "SolOTimes": datetime(2020, 5, 31, 5, 0),
        "margin": 1,
        "soloDurn": 6,
        "relevantTimes": False,
    }

    soloCon = [datetime(2020, 5, 31, 11, 0), datetime(2020, 5, 31, 13, 0)]
    soloAcc = [datetime(2020, 5, 31, 17, 0), datetime(2020, 5, 31, 21, 0)]

    # List of relevant times

    longSolOCase = {
        "AIATimes":
        conCase["AIATimes"],
        "SolOTimes":
        datetime(2020, 5, 30, 12, 0),
        "margin":
        1,
        "soloDurn":
        100,
        "relevantTimes": [{
            "start": soloCon[0],
            "end": soloCon[1],
            "color": "yellow",
            "label": "Con",
        }, {
            "start": soloAcc[0],
            "end": soloAcc[1],
            "color": "blue",
            "label": "Acc.",
        }],
    }

    selectedCase = longSolOCase
    aiaTimesList, soloTimes = extractDatetimePairs(selectedCase)
    # TODO : indicate which of the cases is being tested

    selfCad = int(lc.df.index.freq.delta.total_seconds())
    otherCad = int(solo.df.index.freq.delta.total_seconds())
    # Invert N timeseries
    # solo.df["N"] = solo.df["N"].iloc[::-1].values

    for aiaTimes in aiaTimesList:

        dirName = f"AIA_D{aiaTimes[0].day}_H{aiaTimes[0].hour}:{aiaTimes[1].hour}"

        compareLCurvesToSolo(
            selfObj_Short=lc,
            otherObj_Long=solo,
            selfVars=lcRegs,
            otherVars=soloVars,
            selfTimes=aiaTimes,
            otherTimes=soloTimes,
            selfName="Lcurve",
            otherName="Solo",
            selfCad=60,
            otherCad=60,
            objDirExt=dirName,
            filterPeriods=True,
            PeriodMinMax=[20, 80],
            delete=False,
            showFig=False,
            expectedLocationList=selectedCase["relevantTimes"],
        )

        # TODO : Perform analysis (Add relevant coloured columns)
        # TODO  : Add 211 and 94 Angstrom
        # break
