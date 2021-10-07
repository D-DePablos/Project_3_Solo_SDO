BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"

from sys import path

path.append(f"{BASE_PATH}Scripts/")
path.append(
    f"/home/diegodp/Documents/PhD/Paper_2/InsituEMDCorrelation/Scripts/EMD/")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import makedirs
from EMDComparison import LcurveSolOEMD as lc

from datetime import datetime
from importsProj3.signalAPI import (
    emdAndCompareCases,
    caseCreation,
)
import idlsave
from collections import namedtuple


def main(plotSeparately=True, forceCreateCases=False, multiCPU=4):
    # Set the unsafe, target safe, and dataFolder
    unsafe_dir = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/"
    saveFolder = f"{unsafe_dir}ISSI/New_Method/"
    dataFolder = f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/ISSI/data/"

    # Parameters for showing FIG
    SHOWFIG = False

    # We set a large possible set of periodicities
    PeriodMinMax = [5, 30]
    makedirs(saveFolder, exist_ok=True)

    # IN SITU DATA
    df_is = pd.read_csv(f"{dataFolder}small_ch_in_situ.csv")
    df_is.index = pd.to_datetime(df_is["Time"])
    del df_is["Time"]

    insituParams = ["Vr", "Mf", "Np", "T", "Br"]
    df_is = df_is[insituParams]

    # Set up the dataframes with proper cadence, etc.
    # Attempt to read in dataframes
    try:
        df_171 = pd.read_csv(
            f"{dataFolder}small_ch_171_lc_in.csv", index_col="Time")
        df_193 = pd.read_csv(
            f"{dataFolder}small_ch_193_lc_in.csv", index_col="Time")
        df_flux = pd.read_csv(f"{dataFolder}ch_flux.csv", index_col="Time")
        print("Loaded csv successfully")

        for _df in (df_171, df_193, df_flux):
            _df.index = pd.to_datetime(_df.index)

    # If unable to load CSVs, generate them from base data
    except FileNotFoundError:
        # TODO: Make into function (make ISSI csv)
        # REMOTE DATA
        rs_171 = idlsave.read(
            f"{dataFolder}small_ch_171_lc_in.sav", verbose=False)
        rs_193 = idlsave.read(
            f"{dataFolder}small_ch_193_lc_in.sav", verbose=False)
        ch_flux = idlsave.read(f"{dataFolder}chflux.sav", verbose=False)

        # 171 and 193 observations
        time_array = rs_171.date_obs_171.copy()
        time_array = [t.decode() for t in list(time_array)]

        df_171 = pd.DataFrame(
            {
                "plume": rs_171.lc_171_plume_in,
                "cbpoint": rs_171.lc_171_bp_in,
                "chplume": rs_171.lc_171_ch_plume_in,
                "chole": rs_171.lc_171_ch_in,
                "qsun": rs_171.lc_171_qs_in,
            },
            index=pd.to_datetime(time_array),
        )

        df_193 = pd.DataFrame(
            {
                "plume": rs_193.lc_193_plume_in,
                "cbpoint": rs_193.lc_193_bp_in,
                "chplume": rs_193.lc_193_ch_plume_in,
                "chole": rs_193.lc_193_ch_in,
                "qsun": rs_193.lc_193_qs_in,
            },
            index=pd.to_datetime(time_array),
        )

        # Open and Bright point flux
        flux_time = ch_flux.hmitimes.copy()
        flux_time = [t.decode() for t in list(flux_time)]

        df_flux = pd.DataFrame(
            {
                "ch_open_flux": ch_flux.chofluxes,
                "ch_bpoint_flux": ch_flux.chbpfluxes,
            },
            index=pd.to_datetime(flux_time, format="%Y.%m.%d_%H:%M:%S_TAI"),
        )

        df_171.to_csv(f"{dataFolder}small_ch_171_lc_in.csv",
                      index_label="Time")
        df_193.to_csv(f"{dataFolder}small_ch_193_lc_in.csv",
                      index_label="Time")
        df_flux.to_csv(f"{dataFolder}ch_flux.csv", index_label="Time")

    # # Create test cases
    """
	Generate the cases for all possible SolO - SHORT times (every hour)
	"""
    AIACases = {
        "shortTimes": (datetime(2018, 10, 29, 16), datetime(2018, 10, 30, 23, 50)),
        "longTimes": (datetime(2018, 10, 31, 8), datetime(2018, 11, 2, 8)),
        "shortDuration": 3,
        "caseName": "SDO_AIA",
        "shortDisplacement": 3,
        "MarginHours": 24,
        "savePicklePath": "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/ISSI/cases/AIAcases.pickle",
        "forceCreate": forceCreateCases,
    }

    # Get the cases and put them together with respective AIA observations in Dic
    cases = caseCreation(**AIACases)
    AIACase = namedtuple("AIACase", ["name", "df", "regions", "cases"])
    LongCase = namedtuple("LongCase", ["name", "df"])

    shortDFDic = [
        AIACase(171, df_171.copy(), df_171.columns, cases),
        AIACase(193, df_193.copy(), df_193.columns, cases),
    ]
    longDF = LongCase("PSP", df_is.copy())

    if plotSeparately:

        emdAndCompareCases(
            shortDFDic,
            longDF,
            saveFolder=saveFolder,
            PeriodMinMax=PeriodMinMax,
            showFig=SHOWFIG,
            detrendBoxWidth=200,
            corrThrPlotList=np.arange(0.65, 1, 0.05),
            multiCPU=multiCPU,
        )
    else:
        raise NotImplementedError("Still unable to plot all together")


if __name__ == "__main__":
    main()

    # TODO: Super summary plot - what do they look like for these cases?
