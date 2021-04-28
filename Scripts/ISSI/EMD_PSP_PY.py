BASE_PATH = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/"

from sys import path

path.append(f"{BASE_PATH}Scripts/")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import makedirs

from datetime import datetime, timedelta
import PyEMD
from EMDComparison.signalHelpers import Signal, SignalFunctions, compareTS
from astropy.convolution import Box1DKernel, convolve
from sunpy.time import parse_time
import idlsave

# Set the unsafe, target safe, and dataFolder
unsafe_dir = "/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/"
saveFolder = f"{unsafe_dir}ISSI/"
dataFolder = f"/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/Scripts/ISSI/data/"
makedirs(saveFolder, exist_ok=True)

# IN SITU DATA
df_is = pd.read_csv(f"{dataFolder}small_ch_in_situ.csv")
df_is.index = pd.to_datetime(df_is["Time"])
del df_is["Time"]

# REMOTE DATA
rs_171 = idlsave.read(f'{dataFolder}small_ch_171_lc_in.sav', verbose=False)
rs_193 = idlsave.read(f'{dataFolder}small_ch_193_lc_in.sav', verbose=False)
ch_flux = idlsave.read(f'{dataFolder}chflux.sav', verbose=False)

# 171 and 193 observations
time_array = rs_171.date_obs_171.copy()
time_array = [t.decode() for t in list(time_array)]

df_171 = pd.DataFrame(
    {
        'plume': rs_171.lc_171_plume_in,
        'cbpoint': rs_171.lc_171_bp_in,
        'chplume': rs_171.lc_171_ch_plume_in,
        'chole': rs_171.lc_171_ch_in,
        'qsun': rs_171.lc_171_qs_in,
    },
    index=pd.to_datetime(time_array))

df_193 = pd.DataFrame(
    {
        'plume': rs_193.lc_193_plume_in,
        'cbpoint': rs_193.lc_193_bp_in,
        'chplume': rs_193.lc_193_ch_plume_in,
        'chole': rs_193.lc_193_ch_in,
        'qsun': rs_193.lc_193_qs_in,
    },
    index=pd.to_datetime(time_array))

# Open and Bright point flux
flux_time = ch_flux.hmitimes.copy()
flux_time = [t.decode() for t in list(flux_time)]

df_flux = pd.DataFrame(
    {
        "ch_open_flux": ch_flux.chofluxes,
        "ch_bpoint_flux": ch_flux.chbpfluxes,
    },
    index=pd.to_datetime(flux_time, format="%Y.%m.%d_%H:%M:%S_TAI"))

# # Functions


def plot_variables(df, BOXWIDTH=200):
    plt.figure(figsize=(20, 20))
    for index, parameter in enumerate(list(df)):
        signal = df[parameter].copy()
        csignal = convolve(signal, Box1DKernel(BOXWIDTH), boundary="extend")
        smp_signal = 100 * (signal - csignal) / csignal
        plt.subplot(5, 1, index + 1)
        plt.plot(smp_signal, color='black')
        plt.ylabel(parameter)
    plt.show()
    plt.close()


# # Times of relevance

PeriodMinMax = [5, 60]

timeInsitu = (datetime(2018, 10, 31, 12), datetime(2018, 10, 31, 22))
timeInsitu = (datetime(2018, 10, 31, 12), datetime(2018, 10, 31, 15))

# All of 171 and all of 193 to correlate?
time171 = (df_171.head(1).index.to_pydatetime()[0],
           df_171.tail(1).index.to_pydatetime()[0])
time193 = (df_193.head(1).index.to_pydatetime()[0],
           df_193.tail(1).index.to_pydatetime()[0])

timeFlux = (parse_time("2018-10-29 16:00:09.350").datetime,
            parse_time("2018-10-30 23:49:57.350").datetime)

dfRemotes = {
    "df_171": df_171[time171[0]:time171[1]],
    "df_193": df_171[time193[0]:time193[1]],
    "df_flux": df_flux[timeFlux[0]:timeFlux[1]],
}

# Coloured columns for reference
acc_column = {
    "start": datetime(2018, 10, 30, 2),
    "end": datetime(2018, 10, 30, 19),
    "label": "4/3 Acc. Vsw",
    "color": "yellow",
    "height": 1,
}

con_column = {
    "start": datetime(2018, 10, 30, 11),
    "end": datetime(2018, 10, 31, 2),
    "label": "Con. Vsw",
    "color": "blue",
    "height": 1,
}


)
