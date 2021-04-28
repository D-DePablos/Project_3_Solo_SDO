"""Helper functions"""
from matplotlib import rc
from os import makedirs

# import h5py
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy import signal as scipy_sig
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from glob import glob
import matplotlib
from PyEMD import EMD, Visualisation  # This import works so let's just leave it
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates
from datetime import timedelta
from copy import deepcopy

# Quick fix for cross-package import

emd = EMD()
vis = Visualisation()

font = {"family": "DejaVu Sans", "weight": "normal", "size": 25}
rc("font", **font)


# Local helper functions
def normalize_signal(s: np.ndarray):
    """
    Normalised an input signal
    :param s: Input signal
    """

    _min = np.min(s)
    _max = np.max(s)

    for i, x_i in enumerate(s):
        s[i] = (x_i - _min) / (_max - _min)

    return s


def transformTimeAxistoVelocity(timeAxis, originTime, SPCKernelName=None):
    """
    Gives a corresponding velocity axis for a time axis and originTime

    Args:
        timeAxis: X axis with time values
        originTime: Time that is to be compared (hour of AIA images)
        SPCKernelName: Kernel name for spacecraft
    """
    import heliopy.data.spice as spicedata
    import heliopy.spice as spice

    if SPCKernelName == "solo":
        spicedata.get_kernel("solo")
        sp_traj = spice.Trajectory("Solar Orbiter")

    elif SPCKernelName == "psp":

        spicedata.get_kernel("psp")
        sp_traj = spice.Trajectory("SPP")
        pass

    else:
        raise ValueError(
            f"{SPCKernelName} is not a valid spacecraft kernel, please choose one from ['solo'] "
        )

    times = list(timeAxis)
    sp_traj.generate_positions(times, 'Sun', 'IAU_SUN')  # Is in Km
    R = sp_traj.coords.radius.mean()

    dtAxis = [(t - originTime).total_seconds() for t in timeAxis]
    vSwAxis = R.value / dtAxis  # In km / s

    return vSwAxis


def emd_and_save(s, t, saveFolder, save_name, plot=False):
    """
    Generate ALL EMDs for a relevant timeseries and save them in a numpy file for later use
    Checks if the file exists already to not repeat efforts

    Input params:
    :param s: signal (data)
    :param t: long_signal_time array
    :param saveFolder: Path to save relevant long or short
    :param save_name: Name to save (generally uses long_signal_time references)

    """
    makedirs(saveFolder, exist_ok=True)
    saved_npy = f"{saveFolder}{save_name}.npy"

    try:
        imfs = np.load(saved_npy)
        # print(f"loading imfs {saved_npy}")
        return imfs
    except FileNotFoundError:
        pass

    # Will always use EMD. Will always get imfs and residue separately
    imfs = emd.emd(S=s, T=t)

    if plot:
        plt.figure()
        for imf in imfs:
            plt.plot(imf)
        # plt.show()

    # Find relevant timeID
    time_id = f"time_{t[0]:08d}-{t[-1]:08d}"
    np.save(f"{saveFolder}{time_id}.npy", t)
    np.save(saved_npy, imfs)  # Keeps noise as we can do filtering later

    return imfs


def check_imf_periods(t,
                      imfs,
                      filterPeriods=False,
                      pmin=None,
                      pmax=None,
                      filter_low_high=(0, 0)):
    """
    Check period of several IMFs in minutes. Uses values from SignalFilter class definition
    """
    if filterPeriods:
        if pmin == None or pmax == None:
            raise ValueError(
                f"pmin, pmax value not valid ({pmin}, {pmax}) if IMF number filtering is set to None"
            )

    # Important to drop t to keep same N
    N = t[-1] - t[0]
    Period_valid = np.ndarray((len(imfs), 2))

    if pmin is not None and filter_low_high == (0, 0):
        for i, imf in enumerate(imfs):
            # Maybe not the nicest use of find_extrema
            n_extrema = len(emd.find_extrema(T=t, S=imf)[0])

            if n_extrema != 0 and np.std(imf) > 0.00001:
                P = 2 * N / n_extrema / 60  # Period of each of the IMFs in minutes
                Period_valid[i, 0] = True if pmin < P < pmax else False
                Period_valid[i, 1] = P

            else:
                Period_valid[i, 0] = False
                Period_valid[i, 1] = np.nan

    else:
        if len(imfs) < (filter_low_high[0] + filter_low_high[1]):
            raise ValueError("Filter larger than IMF length")

        Period_valid[:, 0] = False
        # print(f"Invalidating {filter_low_high[0]}:{len(imfs)} - {filter_low_high[1]}")
        Period_valid[filter_low_high[0]:len(imfs) - filter_low_high[1],
                     0] = True
        Period_valid[:, 1] = np.nan

    return Period_valid


def plot_imfs(time_array, vector, imfs, title, savePath, width, show):

    # Find optimal spacing
    spacing = 0  # Start  0 base spacing
    for curr_imf in imfs:
        # Update spacing
        if spacing < abs(curr_imf.max() - curr_imf.mean()):
            spacing = abs(curr_imf.max() - curr_imf.mean())
        if spacing < abs(curr_imf.mean() - curr_imf.min()):
            spacing = abs(curr_imf.mean() - curr_imf.min())

    ax = plt.figure(figsize=(10, 9))
    plt.suptitle(f"{title}")

    for index, curr_imf in enumerate(imfs):
        ax.add_subplot(len(imfs) + 1, 1, len(imfs) - index)
        ax_frame = plt.gca()

        if type(time_array) is np.ndarray:
            plt.plot(time_array, curr_imf, color="black")
        else:
            plt.plot(curr_imf, color="black")
        leg = plt.legend([f"IMF {index}"], loc=1)
        leg.get_frame().set_linewidth(0.0)
        plt.ylim(curr_imf.mean() - spacing - width,
                 curr_imf.mean() + spacing + width)
        ax_frame.axes.get_xaxis().set_visible(False)

    ax.add_subplot(len(imfs), 1, len(imfs))

    if type(time_array) is np.ndarray:
        plt.plot(time_array, vector, color="blue")
        plt.plot(time_array, imfs[-1], color="red")
    else:
        plt.plot(vector, color="blue")
        plt.plot(imfs[-1], color="red")

    plt.legend((["Signal", "Residual"]), frameon=False, loc=1, ncol=2)
    plt.savefig(f"{savePath}.pdf")

    if show:
        plt.show()

    plt.close("all")


def heatmap(
        data,
        row_labels,
        col_labels,
        valid_data=None,
        ax=None,
        cbar_kw={},
        cbarlabel="",
        vmin_vmax=(-1, 1),
        **kwargs,
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    valid_data
        Data which needs to be highlighted
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    vmin_vmax
    |   Minimum and maximum values for image (Hence colourbar)
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    vmin, vmax = vmin_vmax

    if not ax:
        ax = plt.gca()

    # Plot the heatmap

    im = ax.imshow(data, vmin=vmin, vmax=vmax, **kwargs)
    if valid_data is not None:
        ax.imshow(valid_data, vmin=0, vmax=1, alpha=0.8, cmap="Greys")

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(),
             rotation=-30,
             ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)

    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
        im,
        data=None,
        valfmt="{x:.2f}",
        textcolors=("black", "white"),
        threshold=None,
        **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **textkw
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """
    import matplotlib

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j,
                                i,
                                valfmt(data[i, j], None),
                                fontsize="large",
                                **kw)
            texts.append(text)

    return texts


# Classes that we need: a Signal and extended Functions to the signal.
class Signal:
    """
    This class creates a Signal and provides functions to deform and extract it.
    """
    def __init__(
        self,
        saveFolder=None,
        duration=0,
        cadence=0,
        mean=0,
        sig_type=(True, "flat"),
        name="Unnamed Signal",
        custom_data=False,
        saveformat="png",
        **kwargs,
    ):
        """
        This function creates a single line with a given mean_1, with a selected observation cadence
        :param duration: Duration of the signal to be created in seconds
        :param cadence: Observation cadence, within duration
        :param mean: Average value
        :param sig_type: Whether Gaussian or not
        :param real_signal: Whether we are using a real signal or not
        :param custom_data: whether custom data or not!
        """
        assert cadence != 0, "Cadence set at 0"
        self.cadence = cadence
        self.saveFolder = saveFolder
        self.imfs = []
        self.residue = []
        self.name = name
        self.data_smoothed = []
        self.location_signal_peak = []
        self.true_time = None
        self.saveformat = saveformat
        # Minimum and Maximum Periods
        self.pmin = None
        self.pmax = None

        if custom_data is False:
            self.duration = duration
            self.mean = mean
            self.noise = 0
            self.number_of_signals = 0
            self.location_signal_peak = [
            ]  # Set the location signal peak as false

            self.long_signal_time = np.arange(0, duration, step=cadence)
            if sig_type[0] is True:
                # data for observations is a flat line that we add on to
                self.data = np.repeat(float(self.mean), duration / cadence)
            else:
                raise ValueError("different types of signal not implemented")

            for key, value in kwargs.items():
                if key == "gaussian":
                    self.gaussian_curve = value
                if key == "corr_mode":
                    self.correlation_mode = value

        else:
            # Check the type of data that is being fed:
            if type(custom_data) == np.ndarray:
                self.data = custom_data  # If array, can set immediately
                self.true_time = None

            elif type(custom_data) == pd.DataFrame or type(
                    custom_data) == pd.Series:

                if type(custom_data.index) == pd.DatetimeIndex:
                    self.true_time = pd.to_datetime(custom_data.index)

                # When using pandas dataframe
                if type(custom_data) is pd.DataFrame:
                    for column in list(custom_data):
                        self.data = custom_data[column].values
                else:
                    self.data = custom_data.values

                assert (
                    self.true_time is not None
                ), "Was unable to set true long_signal_time using either index or Time Column"

                assert (
                    self.data is not None
                ), "Was unable to set true long_signal_time using either index or Time Column"

            else:
                raise ValueError(
                    f"Did not pass a valid type {type(custom_data)}")

            # Generate constant long_signal_time array for IMFs
            self.duration = self.cadence * len(self.data)
            self.long_signal_time = np.arange(0,
                                              self.duration,
                                              step=self.cadence)

    def __repr__(self):
        """
        Good representation of the object
        """

        return f"Signal({self.saveFolder}, cadence={self.cadence})"

    def __str__(self):
        """
        What to print out
        """

        return f"Signal object with a cadence of {self.cadence}s. Length of {self.duration}"

    def detrend(self, box_width=200):
        """
        Detrend a dataseries and express as percentage difference from new mean
        """

        if box_width is not None:
            from astropy.convolution import convolve, Box1DKernel

            x = self.data
            x_conv = convolve(x, Box1DKernel(box_width), boundary="extend")
            x_pdiff = 100 * (x - x_conv) / x_conv

            self.data = x_pdiff

    def decimate_to(self, other):
        """
        Reshape the array into a given cad. Note that total long_signal_time is preserved.
        :param other: other signal object, with objective cad to decimate down to
        """

        # Cut down into given cad
        cad_factor = int(other.cadence / self.cadence)

        if self.true_time is not None:
            self.true_time = self.true_time[::
                                            cad_factor]  # Lose long_signal_time information

        # Get rid of NA values
        _data = self.data
        nans, x = np.isnan(_data), lambda z: z.nonzero()[0]
        _data[nans] = np.interp(x(nans), x(~nans), _data[~nans])

        _data = _data[::cad_factor]
        self.data = _data

        # Cadence and long_signal_time
        self.cadence = other.cadence
        self.long_signal_time = np.arange(0,
                                          len(self.data) * self.cadence,
                                          step=self.cadence)

        # Update name to reflect decimation
        self.name = f"Decimated {self.name}"
        return self

    def add_noise(self, noise_std=0):
        """
        Add white noise with a given standard deviation - Changes data itself!
        :param noise_std: standard deviation for randomly distributed noise - can be linked to how big original data looks
        """
        self.peak_above = round(self.data.max(), 1) - self.mean

        if noise_std != 0:
            self.noise = np.random.normal(0, noise_std, size=len(self.data))
            self.data = self.data + self.noise
            self.noise_perc = f"{( noise_std/self.peak_above ) * 100:.1f}"

        else:
            self.noise_perc = "0"

    def create_gaussian_add(self, where, duration, std):
        """
        Add a Gaussian into data. Changes self.data Used as "signal"
        :param duration: Duration of the signal (generally relative to the data array)
        :param std: Standard dev. for the created signal
        :param where: First value where gaussian included
        """

        # Mask an array from the initial point for a duration with a gaussian
        mask = ma.masked_inside(self.long_signal_time, where, where + duration)
        mask = ma.getmask(mask)

        if len(self.data) != len(mask):
            mask = mask[:-1]
        extent = np.height_nonzero(mask)

        # Create a gaussian to put inside masked array
        gaussian = scipy_sig.general_gaussian(M=extent, p=1, sig=std)
        gaussian = gaussian - (gaussian[0] - self.data.mean()
                               )  # Only works if flat signal!
        #

        # This is the specific signal that has been placed inside.
        self.gaussian_curve = gaussian
        try:
            self.data[mask] = self.gaussian_curve
        except IndexError:
            self.data[mask] = self.gaussian_curve[-1]

        self.location_signal_peak = where

    def mod_gaussian(self, alter_mode, show=False):
        """
        Modify the gaussian curve in a given way
        :param alter_mode: Mode to modify gaussian in. stretch, multiply, add
        :param show: Whether to show
        """
        # Start by setting it, then modify it in each of the steps.
        print(f"Modifying the signal using alteration {alter_mode}")
        gaussian = self.gaussian_curve

        if not alter_mode:
            return gaussian
        else:
            for key in alter_mode:

                # These would have to be applied multiple times
                if key == "stretch":
                    # Stretch the signal by ... Duplicating or repeating values n times?
                    temp_gaussian = []
                    for i, data in enumerate(gaussian):
                        if data == gaussian[-1] and i > (len(gaussian) / 2):
                            # For the last datapoint, have same twice
                            for n_g in range(alter_mode[key]):
                                temp_gaussian.append(gaussian[i])
                        else:
                            dy = gaussian[i + 1] - gaussian[i]
                            dx = (i + 1) - i
                            grad = dy / dx
                            # Append twice, once at same value, once at half the gradient up/down
                            for n_g in range(alter_mode[key]):
                                # Doing it like this means that you never reach the higher value - desired behaviour!
                                temp_gaussian.append(gaussian[i] + (
                                    grad * n_g / len(range(alter_mode[key]))))

                    gaussian = np.array(temp_gaussian)
                elif key == "height_mod":

                    temp_gaussian = []

                    for i in range(len(gaussian)):

                        # Find gradient at each point by subtracting next to current
                        dy = gaussian[i + 1] - gaussian[i]
                        dx = (i + 1) - i
                        grad = dy / dx
                        n = alter_mode[key]

                        if i == len(gaussian) - 1:
                            temp_gaussian.append(gaussian[i])

                        else:
                            temp_gaussian.append(gaussian[i] + n * grad)

                    gaussian = np.array(temp_gaussian)
                elif key == "multiply":

                    temp_gaussian = gaussian * alter_mode[key]

                    return temp_gaussian
                else:
                    print(f"Mode {key} not supported")

        fig = plt.figure(figsize=(10, 6))
        plt.subplot(1, 2, 1)
        plt.plot(self.gaussian_curve, color="black")
        plt.subplot(1, 2, 2)
        plt.plot(gaussian, color="black")

        if show:
            plt.show()

        return gaussian

    def inc_custom_signal(self, signal, where):
        """
        Include a custom signal from another dataset
        :param signal: What signal to include
        :param where: Where in long_signal_time to include signal
        """
        # In the case where there is a signal peak already
        self.location_signal_peak.append(where)
        if self.number_of_signals == 0:
            mask = ma.masked_inside(self.long_signal_time, where,
                                    where + len(signal) * 12 - 1)
            mask = ma.getmask(mask)

            # Only works if flat signal!
            signal = signal - (signal[0] - self.mean)

            self.gaussian_curve = signal
            self.data[mask] = signal

            self.name = (
                # First long_signal_time change the title
                f"Modified {self.name}")

        else:
            # This simply replaces and that is not correct. Should build up
            mask = ma.masked_inside(self.long_signal_time, where,
                                    where + len(signal) * 12 - 1)
            mask = ma.getmask(mask)

            # Only works if flat signal!
            signal = signal - (signal[0] - self.mean)

            self.gaussian_curve = signal
            self.data[mask] = signal + self.data[mask] - self.mean
            self.name = (
                # First long_signal_time change the title
                f"Duplicated {self.name}")

        self.number_of_signals += 1


class SignalFunctions(Signal):
    """
    Class to give functionality to a given signal object. Allows for things like EMD
    """
    def __init__(
        self,
        signal,
        noise=None,
        filterIMFs=False,
        norm=True,
        PeriodMinMax=[5, 50],
        saveformat=None,
    ):
        self.signalObject = signal
        self.cadence = signal.cadence
        self.base_signal = signal.data
        self.pmin = signal.pmin
        self.pmax = signal.pmax
        self.noise = noise
        self.saveformat = signal.saveformat if saveformat is None else saveformat

        if norm:
            print(f"Normalising {self.signalObject.name}")
            self.s = normalize_signal(signal.data.copy())
        else:
            self.s = signal.data.copy()

        self.t = signal.long_signal_time
        self.true_time = signal.true_time

        if filterIMFs:
            self.filtered = True
            self.pmin = PeriodMinMax[0]
            self.pmax = PeriodMinMax[1]

        else:
            self.filtered = False

        self.name = signal.name
        self.saveFolder = signal.saveFolder

        self.filter_low_high = None
        self.path_to_signal = None
        self.path_to_corr_matrix = None
        self.windowDisp = None

        self.hitrate = 0, 0
        self.table = None
        self.no_displacements = 0

    def __str__(self):
        """
        What to return when asking for string description
        """

        return f"{self.name} Signal. Filtering outside of {self.pmin}:{self.pmax}"

    def __repr__(self):
        """
        Representation of Signal Object
        """

        return f"SignalFunctions({self.name})"

    def plot_norm(self,
                  save_to=None,
                  labels=("Time (s)", "Data (arb.units)"),
                  show=True):

        from matplotlib import rc

        title = f"{self.name} signal sampled at {self.cadence}s cad"
        _ = plt.figure(figsize=(10, 8))
        plt.plot(self.t, self.s, color="black")
        plt.title(title)
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])

        if save_to:
            plt.savefig(
                f"{save_to}{self.name}_{self.cadence}s.{self.saveformat}")

        if show:
            plt.show()

        plt.close("all")

    def generate_windows(
            self,
            other,
            windowDisp,
            plot_long_imfs=False,
            long_window_imf_list=[],
            filterIMFs_on_plot=False,
            useRealTime=False,
            filterPeriods=False,
            filter_low_high=(0, 0),
            renormalize=False,
    ):
        """
        Generate array of relevant windows with two timeseries of different length
        """
        # Note: The difference between self.data and self.s is that S may be normalised
        assert (
            self.t[-1] < other.t[-1]
        ), f"Other Signal Object {other.t[-1]} is shorter than {self.t[-1]}. Please alter order"

        # Once asserted that the first signal is smaller, continue
        short = self
        long = other

        self.filter_low_high = filter_low_high
        self.path_to_signal = f"{short.saveFolder}Split_signal_all_IMFs.npy"
        self.path_to_corr_matrix = f"{short.saveFolder}IMF/Corr_matrix_all.npy"
        self.windowDisp = windowDisp

        print(f"Saving IMF files to {short.saveFolder}")

        short_imfs = emd_and_save(
            s=short.s,
            t=short.t,
            saveFolder=f"{short.saveFolder}IMF/",
            save_name=f"short_{short.t[0]:08d}_{short.t[-1]:08d}",
            plot=False,
        )
        self.imfs = short_imfs

        # Find the valid IMFs on short signal
        if filterPeriods:
            valid_imfs_short = check_imf_periods(
                imfs=short_imfs,
                t=short.t,
                pmin=self.pmin,
                pmax=self.pmax,
                filterPeriods=True,
            )

        else:
            valid_imfs_short = check_imf_periods(
                imfs=short_imfs,
                t=short.t,
                filter_low_high=filter_low_high,
                filterPeriods=False,
            )

        # Setup to perform many EMDs on long dataset
        self.no_displacements = int(
            np.floor(
                (long.t[-1] - short.t[-1]) / self.windowDisp))  # In seconds

        print(f"Found {self.no_displacements} displacements possible.")

        # If the correlation matrix is set, skip
        try:
            np.load(self.path_to_corr_matrix)

            # When it loads correctly, it either continues or returns None
            if plot_long_imfs:
                pass

            # Found files, did not want to plot long imfs
            else:
                return None

        # Continue with the rest of the function, calculating the array
        except FileNotFoundError:
            pass

        # Otherwise continue
        height = 0
        # Assume that two datasets are standardised
        complete_array = np.ndarray((3, short.no_displacements, len(short.s)))

        if not useRealTime:
            # generate 10 by 10 matrix for IMF correlation
            corr_matrix = np.zeros(shape=(12, 12, short.no_displacements, 3))

        else:
            # When needed to save mid_point_time
            corr_matrix = np.zeros(shape=(12, 12, short.no_displacements, 4))

        # Bounds to move window
        left_bound = short.t[0]
        right_bound = short.t[-1]

        # While inside the long timeseries
        while height < short.no_displacements:
            # Only do if necessary!
            if (height in long_window_imf_list
                    and plot_long_imfs) or (plot_long_imfs == False):
                # Find and set relevant window
                i = int(np.where(long.t == left_bound)[0])
                j = int(np.where(long.t == right_bound)[0])

                # Make copies of data instead of using directly
                _data_long = deepcopy(long.s[i:j + 1])
                _data_long = _data_long.reshape(len(_data_long))
                _data_long = normalize_signal(
                    _data_long) if renormalize else _data_long
                _time_long = deepcopy(long.t[i:j + 1])

                # Set values for array
                complete_array[0, height, :] = _time_long
                complete_array[1, height, :] = _data_long
                complete_array[2, height, :] = deepcopy(short.s)

                # Derive EMD and save to relevant folder
                _long_imfs = emd_and_save(
                    s=_data_long,
                    t=_time_long,
                    saveFolder=f"{long.saveFolder}IMF/",
                    save_name=f"long_{_time_long[0]:08d}_{_time_long[-1]:08d}",
                    plot=False,
                )

                # Uses pmin and pmax from short dataseries
                if filterPeriods:
                    _valid_imfs_long = check_imf_periods(
                        imfs=_long_imfs,
                        t=_time_long,
                        pmin=self.pmin,
                        pmax=self.pmax,
                        filterPeriods=filterPeriods,
                    )

                else:
                    _valid_imfs_long = check_imf_periods(
                        imfs=_long_imfs,
                        t=_time_long,
                        filter_low_high=filter_low_high,
                        filterPeriods=filterPeriods,
                    )

                # If required to plot specific IMFs
                if plot_long_imfs and height in long_window_imf_list:
                    __long_signal = Signal(
                        cadence=long.cadence,
                        custom_data=_data_long,
                        name=f"{long.name} - Window #{height}",
                        saveFolder="",
                        norm=False,
                    )

                    __long_sf = SignalFunctions(
                        __long_signal,
                        filterIMFs=filterIMFs_on_plot,
                        PeriodMinMax=(long.pmin, long.pmax),
                        norm=False,
                    )

                    __long_sf.true_time = _time_long
                    __long_sf.pmin = long.pmin
                    __long_sf.pmax = long.pmax

                    __long_sf.plot_emd_inst_freq(
                        ncols=2,
                        savepath=f"{long.saveFolder}IMFplots/",
                        save_name=f"{height:08d}",
                        with_residue=True,
                    )

                    print(
                        f"Saved long IMF {height:08d} figures to {long.saveFolder}IMFplots/ \n"
                    )

                # For all of the short, long IMFs
                for i, row in enumerate(short_imfs):
                    short_valid = valid_imfs_short[i, 0]
                    for j, col in enumerate(_long_imfs):
                        long_valid = _valid_imfs_long[j, 0]
                        if short_valid and long_valid:
                            valid = 1
                        else:
                            valid = 0

                        # Contains pearsonR values, SpearmanR values, and validity
                        # Get the p value for each!
                        corr_matrix[i, j, height, 0] = pearsonr(row, col)[0]
                        corr_matrix[i, j, height, 1] = pearsonr(row, col)[1]
                        corr_matrix[i, j, height, 2] = valid

                if useRealTime:  # We only have the real time in some ocasions
                    mid_point_time = np.floor(
                        (_time_long[-1] + _time_long[0]) / 2)
                    corr_matrix[0, 0, height, 3] = mid_point_time

            # Increase height by one before advancing
            height += 1
            left_bound, right_bound = (
                left_bound + self.windowDisp,
                right_bound + self.windowDisp,
            )

        # TODO: Chck the P value for statistical significance
        np.save(self.path_to_signal, complete_array)  # Signal + IMFs
        np.save(self.path_to_corr_matrix, corr_matrix)  # IMF corrs
        return None

    def plot_emd_inst_freq(
        self,
        ncols,
        savepath,
        save_name,
        imfs=None,
        norm=False,
        with_residue=False,
    ):

        # Begin by generating relevant IMFs
        s = self.s if norm else self.base_signal
        t = self.t
        self.with_residue = with_residue
        if imfs == None:
            imfs = self.generate_imfs_entire_signal(imfs)
            res = self.res
            figsize = (16, 2 * (self.N_imfs + 1))

        fig, axs = plt.subplots(
            self.N_imfs + 1,
            ncols,
            figsize=figsize,
            sharex=True,
            constrained_layout=True,
        )

        # Figure
        # plt.subplot2grid((self.N_imfs + 1, ncols), (0, 0), colspan=1)
        ax0 = axs[0][0]
        ax0.plot(t, s, "r", label="Signal")

        if with_residue:
            ax0.plot(t, res, "k--", label="Residue")
        else:
            ax0.plot(t, imfs[-1], "k--", label="Residue")

        ax0.text(t[0], 0.85 * np.max(s), f"Signal", color="r")
        ax0.text(t[0], 0.65 * np.max(s), f"Residue", color="k")

        # ax0.set_title(f"Normalised {self.name}" + r": $S(t)$")
        # ax0.set_xlabel("Time [s]")
        ax0.set_ylabel("S [a.u.]")

        # Now plot another without the Trend
        ax1 = axs[0][1]
        ax1.plot(t, self.base_signal, "r", label="Signal")
        ax1.set_title("Input Signal")

        ax1.text(t[0], 0.9 * np.max(self.base_signal), f"Signal", color="r")
        # ax1.set_xlabel("Time [s]")
        ax1.set_ylabel("True value [a.u.]")

        # Allows to track which height the IMFs should be at
        plot_index = 1

        for index, (imf, (flag, P), inst_freq) in enumerate(
                zip(self.imfs, self.period_valid, self.inst_freq)):
            if self.filtered and flag or not self.filtered:

                ax = axs[plot_index][0]

                if self.with_residue:
                    ax.plot(t, imf + res, "g")
                    ax.plot(t, res, "k--")
                    ylim_0 = np.mean(res) - 1
                    ylim_1 = np.mean(res) + 1

                else:
                    ax.plot(t, imf, "g")
                    ylim_0, ylim_1 = -1, 1

                if plot_index == self.N_imfs:
                    ax.set_xlabel("Time [s]")

                # ax.set_ylabel("A [a.u.]")
                ax.set_ylim(ylim_0, ylim_1)
                ax.text(
                    t[0],
                    0.7 * ylim_1,
                    f"IMF {index + 1} / {len(self.imfs)}",
                    color="green",
                )
                ax.text(t[0], 0.4 * ylim_1, f"P: {P: .02f} min", color="blue")

                # Instant Freqs
                if ncols > 1:
                    axb = axs[plot_index][1]
                    inst_per = (1 / inst_freq) / 60
                    axb.plot(t, np.log(np.abs(inst_per)), "r")
                    # plt.ylim(0,Phigh/4)
                    axb.set_ylabel("Period [min]")

                    if plot_index == self.N_imfs:
                        axb.set_xlabel("Time [s]")

                plot_index += 1

        makedirs(savepath, exist_ok=True)
        plt.savefig(f"{savepath}{save_name}.{self.saveformat}",
                    bbox_inches="tight",
                    dpi=300)
        # plt.show()
        plt.close()
        print(f"Saved IMFs to {savepath} \n")

    def plot_all_results(
            self,
            other,
            Label_long_ts="No name",
            useRealTime=False,
            expectedLocationList=False,
            savePath=None,
            plot_heatmaps=False,
            minCorrThrPlot=0.7,
            corrThrPlotList=np.arange(0.7, 1.01, 0.1),
            margin_hours=0.5,
            bar_width=1.2,
            filterPeriods=True,
            showFig=False,
            showSpeed=True,  # Whether to show speed instead of time
            SPCKernelName="solo"):
        """
        This function plots the number of IMFs with high correlation for all heights
        Takes signal objects.

        Parameters:
        self: SignalFunctions object
        other: SignalFunctions object

        """
        try:
            corr_matrix = np.load(self.path_to_corr_matrix)

        except FileNotFoundError:
            raise InterruptedError(
                "Please use generate windows before plotting results")

        assert self.pmin == other.pmin and self.pmax == other.pmax, (
            f"Unequal periodicity filtering of"
            f"{self.pmin}:{self.pmax} to {other.pmin}:{other.pmax}")

        assert len(self.s) < len(other.s), (
            "Data of other timeseries is smaller than given timeseries."
            "Please swap their positions")

        makedirs(savePath, exist_ok=True)

        # Init the pearson, spearman and approximate location lists
        # pearsonr_array, spearmanr_array, approximate_locations = [], [], []
        corr_locations = np.ndarray(
            (len(corr_matrix[0, 0, :, 0]), len(corrThrPlotList), 3))

        # Derive the relevant correlation matrix for valid values
        # In the case where there is a signal peak in "Other object"
        if other.signalObject.location_signal_peak != []:
            """
            When there is a location of the signal peak,
            check start and end time, and whether
            valid pearson / spearman correlations were found
            """
            # Arrays that contains positive / negative pairs
            pe_sp_pairs = np.ndarray(
                (len(corr_matrix[0, 0, :, 0]), len(corrThrPlotList), 2))

        for height in range(len(corr_matrix[0, 0, :, 0])):
            # Get all pearson, spearman, and valid values
            pearson = corr_matrix[:, :, height, 0]
            spearman = corr_matrix[:, :, height, 1]
            valid = corr_matrix[:, :, height, 2]

            if useRealTime:
                midpoint = corr_matrix[0, 0, height, 3]
                midpoint_time = other.true_time[0] + timedelta(
                    seconds=midpoint)
                time = midpoint_time

            else:
                # Only when using fake time we have location of the peak
                midpoint = self.t[int(len(self.t) / 2)]
                time = other.t[0] + midpoint + (self.windowDisp * height)

                # NUMPY INT64
                if other.signalObject.location_signal_peak:
                    start_inferred = time - midpoint
                    end_inferred = time + midpoint

                    for loc_peak in other.signalObject.location_signal_peak:
                        if start_inferred <= loc_peak <= end_inferred:
                            location = "inside"
                            break
                    else:
                        location = "outside"

            # Get rid of borders
            pearson[pearson == 0] = np.nan
            spearman[spearman == 0] = np.nan
            mask = ~np.isnan(pearson)

            # Maximum number of IMFs for each
            wind_max = 12 - np.isnan(pearson[:, 0]).sum()
            aia_max = 12 - np.isnan(pearson[0, :]).sum()

            # wind_max_sp = 10 - np.isnan(spearman[:, 0]).sum()
            # aia_max_sp = 10 - np.isnan(spearman[0, :]).sum()
            # Filter to only take values where considered valid due to period filtering
            pvalid = pearson[valid == 1]
            rvalid = spearman[valid == 1]

            # For all relevant correlation thresholds, depending on list
            relevant_pearson_index = np.where(
                corrThrPlotList == minCorrThrPlot)[0][0]

            for index, corr_thr in enumerate(corrThrPlotList):
                _number_high_pe = len(pvalid[np.abs(pvalid) >= corr_thr])
                corr_locations[height, index, 1] = _number_high_pe
                _number_high_sp = len(rvalid[np.abs(rvalid) >= corr_thr])
                corr_locations[height, index, 2] = _number_high_sp

                # TODO: If hit rate important, modify
                # if other.signalObject.location_signal_peak:
                #     if location == "inside":
                #         pe_sp_pairs[height, index, 0] = _number_high_pe
                #         pe_sp_pairs[height, index, 1] = _number_high_sp

                #     else:
                #         pe_sp_pairs[height, index, 0] = -_number_high_pe
                #         pe_sp_pairs[height, index, 1] = -_number_high_sp

            # Only generate heatmap when above threshold
            if plot_heatmaps and corr_locations[height, relevant_pearson_index,
                                                1] >= 1:
                # and 4600 < height < 6000  # Last and conditional for heatmaps
                makedirs(f"{self.saveFolder}corr_matrix/", exist_ok=True)

                # Reduce the arrays to get rid of the noise
                pearson_masked = pearson[mask].reshape(wind_max, aia_max)
                valid_masked = valid[mask].reshape(wind_max, aia_max)

                pvalue_masked = spearman[mask].reshape(wind_max, aia_max)

                # Prepare heatmap
                row_labels, col_labels = [], []

                if not filterPeriods:
                    # In this case, need to cut by specific amount of data
                    # print(f"Filtering with {self.filter_low_high}")
                    low = self.filter_low_high[0]
                    high = -self.filter_low_high[1]

                    # Get rid of edges
                    pearson_hmap = pearson_masked[0:-1, 0:-1]
                    valid_hmap = valid_masked[0:-1, 0:-1]
                    pvalue_hmap = pvalue_masked[0:-1, 0:-1]

                    for row in range(pearson_hmap.shape[0]):
                        row_labels.append(
                            f"AIA IMF {0 + (row + 1)}/ {pearson_hmap.shape[0]}"
                        )

                    for col in range(pearson_hmap.shape[1]):
                        col_labels.append(
                            f"WIND IMF {(col) + 1}/ {pearson_hmap.shape[1]}")
                else:
                    low = 0
                    high = -1
                    # When filtering depending on periods, the number of ignored imfs is different per step. Is this a problem?

                    pearson_hmap = pearson_masked[low:high, low:
                                                  high]  # Eliminate residual
                    valid_hmap = valid_masked[low:high, low:high]
                    pvalue_hmap = pvalue_masked[low:high, low:high]

                    for row in range(pearson_hmap.shape[0]):
                        row_labels.append(
                            f"AIA IMF {low + (row + 1)}/ {pearson_hmap.shape[0]}"
                        )

                    for col in range(pearson_hmap.shape[1]):
                        col_labels.append(
                            f"WIND IMF {(col) + 1}/ {pearson_hmap.shape[1]}")

                # Plot Heatmap
                plt.figure(figsize=(12, 12))
                # Using valid_hmap and  pearson_hmap, attempt to plot just the valid numbers
                im, _ = heatmap(
                    pearson_hmap,
                    row_labels,
                    col_labels,
                    valid_data=valid_hmap if (filterPeriods == True) else None,
                    vmin_vmax=[-1, +1],
                    cmap="RdBu",
                    cbarlabel=f"PearsonR correlation",
                )
                _ = annotate_heatmap(im,
                                     valfmt="{x:.2f}",
                                     textcolors=["black", "black"])

                # Title information is contained in filename instead
                # plt.title(f"Correlation matrix at window #{height} {time}")

                if self.signalObject.location_signal_peak != []:
                    fig_locn = f"{self.saveFolder}corr_matrix/IMF_Heatmap{height:08d}_{time}_peak_at_{self.signalObject.location_signal_peak}.{self.saveformat}"

                else:
                    fig_locn = f"{self.saveFolder}corr_matrix/IMF_Heatmap{height:08d}_{time}.{self.saveformat}"

                plt.tight_layout(pad=0.001)
                plt.savefig(fig_locn, bbox_inches="tight", dpi=300)
                print(f"Saved heatmap to {fig_locn}")
                # plt.show()
                plt.close()

                # Plot PVALUE Heatmap
                plt.figure(figsize=(12, 12))
                # Using valid_hmap and  pearson_hmap, attempt to plot just the valid numbers
                im, _ = heatmap(
                    pvalue_hmap,
                    row_labels,
                    col_labels,
                    valid_data=valid_hmap if (filterPeriods == True) else None,
                    cmap="Blues",
                    vmin_vmax=[0, 1],
                    cbarlabel=f"P-values per IMF",
                )
                _ = annotate_heatmap(im,
                                     valfmt="{x:.2f}",
                                     textcolors=["black", "black"])

                # Title information is contained in filename instead
                # plt.title(f"Correlation matrix at window #{height} {time}")

                if self.signalObject.location_signal_peak != []:
                    fig_locn = f"{self.saveFolder}corr_matrix/IMF_Pvalues{height:08d}_{time}_peak_at_{self.signalObject.location_signal_peak}.{self.saveformat}"

                else:
                    fig_locn = f"{self.saveFolder}corr_matrix/IMF_Pvalues{height:08d}_{time}.{self.saveformat}"

                plt.tight_layout(pad=0.001)
                plt.savefig(fig_locn, bbox_inches="tight", dpi=300)
                print(f"Saved pvalue to {fig_locn}")
                plt.close()

                def create_ts_plot(a=self.s, b=other.s, start=height):
                    c = b[start:start + len(a)]
                    t = other.t[start:start + len(a)] / 60

                    plt.figure(figsize=(16, 12))
                    plt.plot(t, a, color="black", label=r"AIA$_{synth}$")
                    plt.plot(t, c, color="blue", label=r"WIND$_{synth}$")
                    plt.ylabel("Normalised value")
                    plt.xlabel("Minutes since start")
                    plt.legend()
                    plt.title(
                        f"Direct Pearson R correlation: {pearsonr(a, c)[0]:.2f}"
                    )
                    save_to = f"{self.saveFolder}corr_matrix/IMF_Heatmap{height:08d}_{time}_plot.{self.saveformat}"
                    plt.savefig(save_to, bbox_inches="tight", dpi=300)
                    # plt.show()
                    plt.close()

                create_ts_plot()

        # # Save the hitrates
        # if other.signalObject.location_signal_peak:
        #     df_pe = pd.DataFrame({})
        #     df_sp = pd.DataFrame({})

        #     for index, corr_thr in enumerate(corrThrPlotList):
        #         df_pe[f"{corr_thr}"] = pe_sp_pairs[:, index, 0]
        #         df_sp[f"{corr_thr}"] = pe_sp_pairs[:, index, 1]

        #     # Then have the hitrate information available
        #     self.hitrate = (df_pe, df_sp)
        #     # hitrate_tables = self.calculate_hitrate(save=True)

        # Establish relevant signal Objects
        # SHORT
        short_signal = self
        short_values = short_signal.s
        short_time = short_signal.t
        # n_short_imfs = len(short_signal.imfs)

        # LONG
        long_signal = other
        long_values = long_signal.s
        long_true_values = long_signal.base_signal
        # long_true_values

        # window_width = max(short_time) * 12 / 60

        if useRealTime:
            short_duration = short_signal.true_time[
                -1] - short_signal.true_time[0]
            time_axis = long_signal.true_time
        else:
            short_duration = (short_signal.t[-1] - short_signal.t[0]) / 60
            time_axis = long_signal.t / 60

        region_string = self.name
        #######################################
        # Figure
        fig, axs = plt.subplots(2, figsize=(20, 10), sharex=True)

        # First plot
        ax = axs[0]
        ax.plot(time_axis,
                long_true_values,
                color="black",
                label=Label_long_ts,
                alpha=1)
        # TODO: Plot detrended dataset

        if filterPeriods:
            ax.set_title(f"{self.pmin} < {r'$P_{IMF}$'} < {self.pmax} min")

        elif not filterPeriods:
            ax.set_title(f"Noise {self.signalObject.noise_perc}%")

        else:
            raise ValueError("Please determine whether using period or not")

        ax.set_ylabel(f"{Label_long_ts}")
        if useRealTime:
            ax.set_xlim(
                time_axis[0] - timedelta(hours=margin_hours),
                time_axis[-1] + timedelta(hours=margin_hours),
            )
            ax.xaxis.grid(True)

            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
            ax.xaxis.set_tick_params(rotation=45)

        else:
            # Set xticks every hour
            xticks = np.arange(time_axis[0], time_axis[-1] + 1, step=180)
            plt.xticks(xticks)

        if useRealTime:
            # Plot the hits
            try:
                # Saved only sometimes!
                true_time_secs = corr_matrix[0, 0, :, 3]
            except IndexError as ind_err:
                raise IndexError(
                    "Please ensure that the true timeseries is saved!")

            ttindex = (true_time_secs / long_signal.cadence).astype(int)
            true_time_datetime = long_signal.true_time[ttindex]
            time_mid_min = true_time_datetime
            bar_width = short_duration

        else:
            time_mid_min = corr_locations[:, 0, 0]
            bar_width = short_duration

        ###################################################
        # BOTTOM PLOT
        ax2 = axs[1]

        #######################################
        # INSET AXIS
        # displays short signal

        # if Label_long_ts == 'Mf': --> Just /inset and change to if True
        in_ax = inset_axes(
            ax,
            width="15%",  # width = 30% of parent_bbox
            height=1,  # height : 1 inch
            loc="upper left",
        )
        in_ax.plot(short_time, short_values, color="black")
        in_ax.set_title(f"{region_string}")
        plt.xticks([])
        plt.yticks([])

        # Adding specific colormap for pairs
        possibleColors = {"1": "green", "2": "magenta", "3": "yellow"}

        for height in range(len(time_mid_min)):
            # Open up pearson and spearman IMF pair list
            pearson_array = corr_locations[height, :, 1]
            # spearman_array = corr_locations[height, :, 2]

            if useRealTime:
                midpoint = corr_matrix[0, 0, height, 3]
                midpoint_time = other.true_time[0] + timedelta(
                    seconds=midpoint)
                time = midpoint_time
                barchart_time = time.to_pydatetime()

            else:
                midpoint = self.t[int(len(self.t) / 2)]
                time = other.t[0] + midpoint + (self.windowDisp * height)
                barchart_time = time / 60
                # NUMPY INT64

            # Bar charts for each of the heights
            for index, corr_label in enumerate(corrThrPlotList):
                if pearson_array[index] != 0:  # If some pairs are found

                    try:
                        _color = possibleColors[f"{int(pearson_array[index])}"]
                    except KeyError:
                        _color = "red"

                    # Bug with periodicity calculation. WTF
                    _alpha = 0.35
                    ax2.bar(
                        barchart_time,
                        corr_label,
                        width=bar_width,
                        color=_color,
                        edgecolor="white",
                        alpha=_alpha,
                        zorder=1,
                    )

                # if spearman_array[index] != 0:
                #     try:
                #         _color = possibleColors[
                #             f"{int(spearman_array[index])}"]
                #     except KeyError:
                #         _color = "red"
                #     _alpha = 0.35 if spearman_array[index] > 0 else 0
                #     ax2.bar(
                #         barchart_time,
                #         corr_label,
                #         width=bar_width,
                #         color=_color,
                #         edgecolor="white",
                #         alpha=_alpha,
                #         zorder=2,
                #     )

        # Columns on bottom plot
        if useRealTime:
            margin_label = timedelta(hours=1)
            # interval = int((other.t[-1] / 3600) / 10)
            interval = 3

            ax2.xaxis.set_minor_locator(mdates.HourLocator(1))
            ax2.xaxis.set_major_locator(mdates.HourLocator(interval=interval))
            ax2.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
            ax2.xaxis.set_tick_params(rotation=10)
            # Set the x limits
            ax2.set_xlim(
                time_axis[0] - timedelta(hours=margin_hours),
                time_axis[-1] + timedelta(hours=margin_hours),
            )
            # If using some expected location dictionary add here
            if expectedLocationList is not False:
                for expected_location in expectedLocationList:
                    start_expected = expected_location["start"]
                    end_expected = expected_location["end"]
                    color = expected_location["color"]
                    label = expected_location["label"]

                    ax2.axvspan(
                        xmin=start_expected,
                        xmax=end_expected,
                        ymin=0,
                        ymax=1,
                        alpha=0.3,
                        color=color,
                    )

                    ax2.text(
                        x=start_expected + timedelta(minutes=15),
                        y=0.95,
                        s=f"{label}",
                    )

        else:
            ax2.xaxis.set_tick_params(rotation=25)
            margin_label = -10
            interval = int((other.t[-1] / 3600) / 10)

        # Set out the legend and info about pairs
        legend_x = time_axis[-1] + margin_label
        # ax2.text(x=legend_x, y=0.95, s="1 pair", color=possibleColors["1"])
        # ax2.text(x=legend_x, y=0.9, s="2 pairs", color=possibleColors["2"])
        # ax2.text(x=legend_x, y=0.85, s="3 pairs", color=possibleColors["3"])

        # Divide by 2 number of ticks and add 1 as highest Corr
        corr_thr_new = []
        for i, value in enumerate(corrThrPlotList):
            if np.mod(i, 2) == 0:
                corr_thr_new.append(value)
        corr_thr_new.append(1)

        ax2.set_yticks(corr_thr_new)  # Ensure ticks are good
        ax2.set_ylim(corrThrPlotList[0], 1.01)  # Allows for any correlation
        ax2.set_ylabel("Highest corr. found")

        if showSpeed:
            # Extract the Velocities
            vSWAxis = transformTimeAxistoVelocity(
                time_axis,
                originTime=short_signal.true_time.mean().to_pydatetime(),
                SPCKernelName=SPCKernelName)

            axV = ax2.twiny()
            axV.plot(vSWAxis, long_values, alpha=0)
            axV.set_title("Implied Vsw (km/s)")

        ax2.grid(True)

        if useRealTime:
            ax2.set_xlabel("Date (dd HH:MM)")
        else:
            ax2.set_xlabel("Minutes since start")

        # Save, show and close
        save_name = f"{long_signal.name}_{region_string}_IMF_Summary"
        plt.tight_layout()
        plt.savefig(
            f"{self.saveFolder}{save_name}.{self.saveformat}",
            dpi=300,
            bbox_inches="tight",
        )
        if showFig:
            plt.show()
        plt.close()

    def calculate_hitrate(self, save=False):
        """
        Creates dataframe with hitrate and prints it out
        """
        assert (
            self.hitrate
            is not False), "Please calculate hitrates with generate_windows"

        self.table = {}
        df_pe, df_sp = self.hitrate

        for label, df_hrate in zip(("pearson", "spearman"), (df_pe, df_sp)):
            results_df = pd.DataFrame({})

            # For the different correlation thresholds
            for _, corr_thr in enumerate(df_hrate):

                # Reset hit-rate to ensure working well
                _hit_rate = 0

                # Find hitrate
                ncorrect = np.sum(df_hrate[corr_thr][df_hrate[corr_thr].gt(0)])
                # nincorrect = np.sum(df_pe[corr_thr][df_pe[corr_thr].lt(0)])
                ntotal = np.sum(abs(df_hrate[corr_thr]))

                if ntotal != 0:
                    _hit_rate = ncorrect / ntotal * 100

                results_df[f"{corr_thr}"] = [f"{_hit_rate:.2f}"]

            self.table[label] = results_df

        if save:
            # This will be one of the results dataframes
            self.table["pearson"].to_csv(
                f"{self.saveFolder}pearson_hitrate.csv")
            self.table["spearman"].to_csv(
                f"{self.saveFolder}spearman_hitrate.csv")

        return self.table

    def plot_single_hitrate(self, mode="pearson", show=False):
        """
        Plots hitrates
        Needs to replicate  TestIMF functionality
        Input of a dictionary with several
        """
        assert (
            self.table is not None
        ), f"Hit-rate table not created. Use {self.calculate_hitrate.__name__}"

        corr_table = self.table
        n = len(list(corr_table))  # All correlation thresholds
        x = range(n)
        # font = {"family": "DejaVu Sans", "weight": "bold", "size": 18}

        matplotlib.rc("font", **font)

        fig = plt.figure(figsize=(12, 8))
        ax = plt.gca()
        ax.set_prop_cycle("color", "blue")

        for corr in list(corr_table[mode]):
            # This is plot of hit rate per correlation threshold, with the windowing as colour
            plt.plot(
                float(corr),
                float(corr_table[mode][corr][0]) + 0.1,
                label=f"{self.windowDisp}s",
                marker="o",
                markersize=14,
                color="black",
            )

        plt.grid("both")
        ax.set_ylabel("Hit rate (%)")
        ax.set_xlabel("Correlation threshold")
        ax.set_ylim(-5, 105)

        plt.title(f"{self.signalObject.noise_perc}% noise")
        fig.subplots_adjust(bottom=0.25)
        if show:
            plt.show()

        plt.savefig(
            f"{self.saveFolder}{self.windowDisp}"
            f"_{self.signalObject.noise_perc}%_hitrates.{self.saveformat}",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()


def compareTS(
    dfSelf,
    dfOther,
    cadSelf,
    cadOther,
    labelOther,
    winDispList=[60],
    corrThrPlotList=[],
    PeriodMinMax=[1, 180],
    filterPeriods=False,
    savePath=None,
    useRealTime=False,
    expectedLocationList=False,
    detrend_box_width=200,
    delete=False,
    showFig=True,
    renormalize=False,
    showSpeed=False,
    SPCKernelName="solo",
):
    """
    Takes two dataframes sampled at same cadence
    
    dfSelf is a dataframe
    dfOther is another dataFrame
    cadSelf, cadOther for cadence
    labelOther is the label to be shown on second dataset
    
    winDispList is a list of window displacements, in seconds
    corrThrPlotList is 

    
    
    
    # """
    # This function takes two dataframes with variables sampled at the same rate, and performs EMD analysis on them
    # dfSelf is the dataframe of an object (Self)
    # dfOther is the dataframe of a different object (Self)
    # cadSelf, cadOther is cadence
    # labelOther = Which label to show for other dataset
    # """

    assert savePath != None, "Please set savePath to store relevant arrays"
    deleted_already = False
    # For all of the lightcurves
    for varOther in list(dfOther):
        otherPath = f"{savePath}{labelOther}/{varOther}/"
        makedirs(otherPath, exist_ok=True)

        signalOther = Signal(
            cadence=cadOther,
            custom_data=dfOther[varOther],
            name=varOther,
            time=dfOther.index,
            saveFolder=otherPath,
        )

        signalOther.detrend(box_width=detrend_box_width)

        # Add functionality, filter and generate IMFs
        otherSigFunc = SignalFunctions(signalOther,
                                       filterIMFs=True,
                                       PeriodMinMax=PeriodMinMax)

        # For each of the in-situ variables
        for varSelf in dfSelf:
            selfPath = f"{savePath}{varSelf}/{varOther}/"
            makedirs(selfPath, exist_ok=True)

            dataSelf = dfSelf[varSelf]
            signalSelf = Signal(
                cadence=cadSelf,
                custom_data=dataSelf,
                saveFolder=selfPath,
                name=varSelf,
            )
            signalSelf.detrend(box_width=detrend_box_width)
            selfSigFunc = SignalFunctions(
                signalSelf,
                filterIMFs=True,
                PeriodMinMax=PeriodMinMax,
            )

            for windowDisp in winDispList:
                otherWinFolder = f"{otherPath}{windowDisp}s/"
                selfWinFolder = (
                    f"{selfPath}{windowDisp}s/{PeriodMinMax[0]} - {PeriodMinMax[1]}/"
                )
                otherSigFunc.saveFolder = otherWinFolder
                selfSigFunc.saveFolder = selfWinFolder

                if delete and not deleted_already:
                    from shutil import rmtree
                    print("Deleting otherSigFunc and selfSigFunc folders")
                    rmtree(otherSigFunc.saveFolder, ignore_errors=True)
                    rmtree(selfSigFunc.saveFolder, ignore_errors=True)

                selfSigFunc.generate_windows(
                    other=otherSigFunc,
                    windowDisp=windowDisp,
                    useRealTime=useRealTime,
                    filterPeriods=filterPeriods,
                    renormalize=renormalize,
                )
                selfSigFunc.plot_all_results(
                    other=otherSigFunc,
                    Label_long_ts=varOther,
                    plot_heatmaps=False,
                    minCorrThrPlot=corrThrPlotList[0],
                    savePath=f"{otherSigFunc.saveFolder}corr_matrix/",
                    corrThrPlotList=corrThrPlotList,
                    bar_width=None,
                    useRealTime=useRealTime,
                    expectedLocationList=expectedLocationList,
                    showFig=showFig,
                    showSpeed=showSpeed,
                    SPCKernelName=SPCKernelName,
                )


def plot_variables(df, BOXWIDTH=200):
    """Plots detrended variables found within a dataframe"""
    from astropy.convolution import convolve, Box1DKernel

    plt.figure(figsize=(20, 20))
    for index, parameter in enumerate(list(df)):
        signal = df[parameter].copy()
        csignal = convolve(signal, Box1DKernel(BOXWIDTH), boundary="extend")
        smp_signal = 100 * (signal - csignal) / csignal
        plt.subplot(5, 1, index + 1)
        plt.plot(smp_signal, color="black")
        plt.ylabel(parameter)
    plt.show()
    plt.close()
