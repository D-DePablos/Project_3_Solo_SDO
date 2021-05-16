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
from astropy import constants as const
from astropy import units as u


def plot_all(
    soloTimes,
    solo,
    soloVars,
    lcTimes,
    lc94,
    lc193,
    lc211,
    regions,
):
    # First get all dataframes with good column names and copped to relTimes
    dfsolo = solo.df[soloVars]
    dfsolo.columns = [f"SolO_{i}" for i in soloVars]  # Rename the columns
    dfsolo = dfsolo[soloTimes[0]:soloTimes[1]]

    # Construct lc dictionary + rename lcurves
    lcDic = {}
    for lcObj, _wvl in zip((lc94, lc193, lc211), ("94", "193", "211")):
        lcObj.df = lcObj.df[regions]
        lcObj.df.columns = [f"{_wvl}_{i}" for i in regions]
        lcObj.df = lcObj.df[lcTimes[0]:lcTimes[1]]
        lcDic[_wvl] = lcObj.df

    # The new plot format overplots coloured columns (depending on wavelength) on top of SolO observations
    # It then shows the 3 lightcurves with extracted IMFs
    # This makes the ~11 panel figure a 7 panel figure!
    """
    Collect remote, collect in situ
    
    Select a region
    
    For all wavelengths, collect all of the correlation matrices
    
    """

    # TODO: Complete this!
    def construct_plot(dfsolo, lcDic, region):
        from glob import glob

        def _find_corr_mat(
                region,
                wvl,
                insituParams,
                _base_folder="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/ISSI/SB_6789/",
                windDisp="60s",
                period="5 - 20"):
            """Finds correlation matrices for all In-situ variables and a given wavelength

            Args:
                region (str): Which remote sensing region to explore
                wvl (int): Which wavelength is relevant
                insituParams (list): In situ parameters to find correlation matrices for
            """

            resultingMatrices = {}
            for isparam in insituParams:
                # TODO: Need to set this to load in the specific Corr_matrix, short signal and time
                _subfolder = f"{_base_folder}{wvl}_{region}/*{isparam}/{windDisp}/{period}/"
                foundMatrix = glob(f"{_subfolder}IMF/Corr_matrix_all.npy")
                short_D = glob(f"{_subfolder}IMF/short*.npy")
                short_T = glob(f"{_subfolder}IMF/time*.npy")
                resultingMatrices["isparam"] = foundMatrix

        # Select the correlation matrix for each solo variable, given region, given lcurve

        for wvl in lcDic:
            # _find_corr_mat()

        pass

    # Each of the times are selected above
    for region in regions:
        construct_plot(
            dfsolo,
            lcDic,
            region,
        )

        pass


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
    SOLOHI=None,
    SOLOLO=None,
    plot_all_together=False,
):
    """
    Feed in the PSP Spacecraft and SolOSpc object
    Self is expected to be solar orbiter
    Other is expected to be the lcurves
    """

    assert SOLOHI != None, "No High speed set"
    assert SOLOLO != None, "No Low speed set"

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
        LOSPEED=SOLOHI,
        HISPEED=SOLOLO,
    )


def extractDiscreteExamples(Caselist, margin, AIAduration=1):
    """Extracts disxrete AIA - SolO pairs

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
    soloTimes = []
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

        # SoloDurn gives reference to amount of SolO datapoints that match
        matchEnd = matchStart + timedelta(hours=case["soloDurn"])
        matchAvg = matchStart + (matchEnd - matchStart) / 2
        matchTimes.append((matchStart, matchEnd))

        # Get the Solar Orbiter measurements
        soloStart = matchAvg - timedelta(hours=margin)
        soloEnd = matchAvg + timedelta(hours=margin)
        soloTimes.append((soloStart, soloEnd))

        # Get the specific case Name
        caseNames.append(case["caseName"])

        refLocations.append(
            _constructExpectedLocation(_times=(matchStart, matchEnd)))

    return aiaTimes, soloTimes, caseNames, refLocations


if __name__ == "__main__":
    objCad = 60  # Objective cadence in seconds for comparisons

    DELETE = False
    SHOWFIG = True
    FILTERP = True
    PLOT_ALL_TOGETHER = True
    accelerated = 1
    # accelerated = 4 / 3
    PERIODMINMAX = [3, 20]

    for wv in [94, 193, 211]:
        WAVELENGTH = wv

        # Solar Orbiter Data requires start, end
        start = datetime(2020, 5, 30, 12)
        end = datetime(2020, 6, 2)
        solo = SoloManager(
            times=(start, end),
            objCad=objCad,
            cdfPath=
            "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/soloData/",
        )
        solo.df = solo.df.interpolate()  # Fill gaps
        # Velocities are modified with 4/3 factor. Gives slightly better idea
        SOLOHI, SOLOLO, MEAN = (int(solo.df["V_R"].max() / accelerated),
                                int(solo.df["V_R"].min() / accelerated),
                                int(solo.df["V_R"].mean() / accelerated))
        """
                                    V_R        V_T       V_N          N         T
            2020-05-30 12:00:00  274.279297 -20.320356 -5.449193  43.087570  1.601013
            ...                         ...        ...       ...        ...       ...
            2020-05-31 23:00:00  274.110168 -44.166534 -4.904387  58.298233  2.722118
        """

        # Calculate mass flux
        Vx = (solo.df["V_R"].values * (u.km / u.s)).to(u.m / u.s)
        mp = const.m_p
        N = (solo.df["N"].values * (u.cm**(-3))).to(u.m**(-3))
        solo.df["Mf"] = (N * mp * Vx).value

        # Variables in situ
        soloVars = ["N", "T", "V_R", "Mf"]
        # soloVars = ["V_R"]

        # Light Curves
        lc = LcurveManager(
            csvPath=
            "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/sharedData/",
            objCad=objCad,
            wavelength=WAVELENGTH,
        )
        lc.df = lc.df.interpolate()  # Interpolate after forming lc object
        """
                                    16          17          18          21          22          23
            Time                                                                                       
            2020-05-27 00:00:00  329.801805  607.425191  613.393714  621.675743  720.785361  623.828448
            ...                         ...         ...         ...         ...         ...         ...
            2020-05-28 10:45:00  386.551898  623.982588  672.670828  641.362444  665.743955  570.076031
        """
        lcRegs = [
            "11", "12", "13", "16", "17", "18", "21", "22", "23", "11:13_21:23"
        ]

        # Cadence setup
        selfCad = int(lc.df.index.freq.delta.total_seconds())
        otherCad = int(solo.df.index.freq.delta.total_seconds())

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

            if not PLOT_ALL_TOGETHER:
                compareLCurvesToSolo(
                    selfObj_Short=lc,
                    otherObj_Long=solo,
                    selfVars=lcRegs,
                    otherVars=soloVars,
                    selfTimes=aiaTimes,
                    otherTimes=soloTimesList[index],
                    selfName=f"{WAVELENGTH}",
                    otherName="Solo",
                    selfCad=60,
                    otherCad=60,
                    objDirExt=dirName,
                    filterPeriods=FILTERP,
                    PeriodMinMax=PERIODMINMAX,
                    delete=DELETE,
                    showFig=SHOWFIG,
                    expectedLocationList=[refLocations[index]],
                    renormalize=False,
                    SOLOHI=SOLOHI,
                    SOLOLO=SOLOLO,
                )
            elif PLOT_ALL_TOGETHER:
                pass
