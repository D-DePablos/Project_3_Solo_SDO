from __future__ import division
import glob
import os
from datetime import timedelta
from multiprocessing import Process

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class AstroMap:
    """
    This class is able to download and use fits files for ligthcurves.
    """

    def __init__(
        self,
        global_save_path,
        timestamp,
        wavelength,
        fits_folder,
        n_row_col,
        n_processes,
        relevantCoords,
        exposure,
        sigma,
        showPreview=False,
    ):
        """
        To initialize an AstroMap object, we require a main path,
        :param global_save_path: path for reference for rest of functions. Saving here and below.
        :param fits_folder: Folder to fits files (if existing)
        :param timestamp: Start and endtime to take into account
        :param wavelength: Relevant wavelength
        :param n_row_col: Number of rows and columns to splice array into

        :param n_processes: Number of processes to use
        """

        self.global_save_path = global_save_path
        self.fits_folder = fits_folder

        self.wvlnth = wavelength
        self.nrowcol = n_row_col
        self.n_cpus = n_processes
        self.timestamp = timestamp

        self.relevantCoords = relevantCoords
        self.exposure = exposure

        self.complete_lcurves = []
        self.sigma = sigma
        self.showPreview = showPreview

        if not self.fits_folder:
            print("No fits folder given")

    @staticmethod
    def translate_number(input_region, n_by_n):
        """
        Take an input region and a given n_by_n to return two dimensional response to where we are on array

        :param input_region: Which number we require
        :param n_by_n: Dimensions of array
        :return: Returns height, length at which box is located
        """

        array = np.arange(0, n_by_n ** 2, 1)
        array = np.split(array, n_by_n)
        array = np.vstack(array)

        return [
            np.where(array == input_region)[0][0],
            np.where(array == input_region)[1][0],
        ]

    def download_all(self, down_dir):
        """
        Simple script that downloads fits into download folder
        :param down_dir: Download directory for the fits
        :return:
        """
        from Scripts.Imports.Data_Fetch.Fido_Download import run_fido

        run_fido(
            start=self.timestamp[0],
            end=self.timestamp[1],
            down_dir=down_dir,
            wavelengths=self.wvlnth,
        )

        self.fits_folder = down_dir

    def create_column_plots(self, figsize=(10, 10)):
        """
        This function saves the generated light curves into a given path_gen/lightcurve
        :param diff: Whether to plot diff or not
        :return: Saves complete images to a given folder
        """
        from matplotlib import gridspec as gs

        def proc_multi(split_indices):
            """
            Function created for multiprocessing of getting images and lcurves together
            :param split_indices: Which indexes to use
            :return: Saves figures into given path
            """
            for i in split_indices:
                curr_time = timearr[i]
                for numpycol in range(lcurve_data.shape[0]):
                    fig = plt.figure(figsize=figsize)
                    for numpyrow in range(lcurve_data.shape[1]):
                        # Save each of the images in their own folder -
                        folder_image = (
                            f"{self.global_save_path}{numpycol:02d}_{numpyrow:02d}"
                        )
                        curr_aia = plt.imread(f"{folder_image}/{i:03d}.png")
                        curr_lcurve = lcurve_data[numpyrow, numpycol, :]

                        plt.subplot(grs[numpyrow, 0])
                        plt.imshow(curr_aia)
                        plt.axis("off")

                        # Plot lightcurve and corresponding column
                        plt.subplot(grs[numpyrow, 1])
                        ax = plt.gca()
                        plt.plot(timearr, curr_lcurve, linewidth=2, color="black")
                        ax.axvline(
                            x=curr_time,
                            color="black",
                            linestyle="--",
                            alpha=0.2,
                            linewidth=2,
                        )

                    plt.tight_layout()

                    combo_path = f"{self.global_save_path}Total_combined/"
                    os.makedirs(combo_path, exist_ok=True)
                    plt.savefig(f"{combo_path}{numpycol:02d}_{i:03d}.png")
                    plt.close(fig)

                print(f"Done index : {i}")

        # Get total images, data and time arrays
        combo_images = sorted(glob.glob(f"{self.global_save_path}Combo/*.png"))

        if combo_images == []:
            raise ValueError(
                f"No images were found in {self.global_save_path}Combo/. "
                f"/n Please use the astro_plot function first"
            )

        lcurve_data = np.load(f"{self.global_save_path}lcurves.npy")

        time_npy = np.load(f"{self.global_save_path}times.npy")
        tdf = pd.DataFrame({"Time": time_npy})
        tdf.Time = pd.to_datetime(tdf.Time)
        timearr = tdf.Time

        # Now need to separate lightcurves as required:
        grs = gs.GridSpec(lcurve_data.shape[0], 2, width_ratios=(1, 4))
        index_split = np.array_split(
            np.array(np.arange(len(combo_images))), self.n_cpus
        )

        proc = []

        # Create all the required processes for the given indices
        for index_list in index_split:
            if index_list.size:
                pr = Process(target=proc_multi, args=[index_list])
                proc.append(pr)
                pr.start()

        for process in proc:
            process.join()

    def load_lcurves(self, path):
        """
        This function opens up the lightcurves that were saved prior
        :param path: Absolute path to csv containing complete_lcurves.csv
        """
        df = pd.read_csv(path)
        self.complete_lcurves = df

        return self.complete_lcurves

    def save_numpy(self, delete=True):
        """
        Save numpy array of entire lightcurve set. Required for EMD analysis
        :param delete:
        :return:
        """
        # Now saving each of the files individually, numbered, and stitching them together!

        all_times = sorted(glob.glob(f"{self.global_save_path}times_*.npy"))
        print(self.global_save_path)
        comb_times = np.empty_like(np.load(all_times[0])[0:1])
        for index, file in enumerate(all_times):
            curr_time = np.load(file)
            comb_times = np.concatenate([comb_times, curr_time])

            if delete:
                os.remove(file)

        comb_times = comb_times[1:]
        np.save(f"{self.global_save_path}times.npy", comb_times)

        all_lcurves = sorted(glob.glob(f"{self.global_save_path}lcurves_*.npy"))
        comb_lcurves = np.empty_like(np.load(all_lcurves[0])[:, :, 0:2])

        for index, file in enumerate(all_lcurves):
            curr_lc = np.load(file)
            comb_lcurves = np.block([comb_lcurves, curr_lc])

            if delete:
                os.remove(file)

        comb_lcurves = comb_lcurves[:, :, 2:]
        np.save(f"{self.global_save_path}lcurves.npy", comb_lcurves)

        all_pix = sorted(glob.glob(f"{self.global_save_path}pix_*.npy"))
        comb_pix = np.empty_like(np.load(all_pix[0])[:, :, 0:2])

        for index, file in enumerate(all_pix):
            curr_pix = np.load(file)
            comb_pix = np.block([comb_pix, curr_pix])

            if delete:
                os.remove(file)

        comb_pix = comb_pix[:, :, 2:]
        np.save(f"{self.global_save_path}pix.npy", comb_pix)

        return [
            f"{self.global_save_path}times.npy",
            f"{self.global_save_path}lcurves.npy",
            f"{self.global_save_path}pix.npy",
        ]

    def extract_info_from_fits(self, extract_regions, files_array=None):
        """
        Extract and save lightcurves from the fits files
        :param extract_regions: Regions which are of specific interest
        :param files_array: Array of fits files (already sorted glob.glob)
        """

        from ...Solar.Plotting.helpers import astro_plot
        import multiprocessing

        if not files_array:
            files_list = sorted(glob.glob(f"{self.fits_folder}*.fits"))

            if not files_list:
                raise ValueError("Your fits array does not exist!")

            else:
                files_array = np.array(files_list)

        # Create arrays which are split, as well as corresponding indices to save images
        split_file_arr = np.array_split(files_array, self.n_cpus)
        split_index_array = np.array_split(
            np.array(np.arange(len(files_array))), self.n_cpus
        )

        # Initialize required multiprocessing methods
        # According to some website, Pool is better when more than 10 cores are available. (Overhead big otherwise)
        procs = []

        print(
            f"The length of our array is {len(files_array)}; Split by {self.n_cpus} CPU(s) "
        )

        if self.n_cpus == 1:
            astro_plot(
                split_file_arr[0],
                split_index_array[0],
                files_array,
                self.wvlnth,
                self.relevantCoords,
                self.nrowcol,
                self.global_save_path,
                self.sigma,
                extract_regions,
                self.showPreview,
            )

            # Debugging on one processor

        # instantiating process with arguments
        for index_list, files_list in zip(split_index_array, split_file_arr):

            print(index_list[0:5])
            print(index_list[-5:-1])

            if files_list.size:
                proc = multiprocessing.Process(
                    target=astro_plot,
                    args=(
                        files_list,
                        index_list,
                        files_array,
                        self.wvlnth,
                        self.relevantCoords,
                        self.nrowcol,
                        self.global_save_path,
                        self.sigma,
                        extract_regions,
                        self.showPreview,
                    ),
                )
                procs.append(proc)
                proc.start()

        # Wait for all processes to complete
        for proc in procs:
            proc.join()

        print("Finished images. Now stitching numpy arrays")
        self.save_numpy(delete=True)

        # Save them in the main folder as complete_lcurves.csv
        self.save_lcurves_to_csv()

    def save_lcurves_to_csv(self, redirect_path=False):
        """
        This function saves all lightcurves individually, labeling them from '0' -> end as a single csv file
        :param redirect_path: Which path to look in
        :return: Returns path to saved csv
        """

        if redirect_path:
            path_main = redirect_path
        else:
            path_main = self.global_save_path

        for file_type in ["pix", "lcurves"]:

            existing_csv = glob.glob(f"{path_main}*{file_type}*.csv")

            if existing_csv == []:
                # Open up the file
                numpy_lcurves = np.load(f"{path_main}{file_type}.npy")

                df = pd.DataFrame({})
                index = 0
                df["Time"] = np.load(f"{path_main}times.npy")

                for row in range(len(numpy_lcurves[:, 0, 0])):
                    for col in range(len(numpy_lcurves[0, :, 0])):
                        curr_curve = numpy_lcurves[row, col, :]
                        df[f"{index}"] = curr_curve
                        index += 1

                df.to_csv(f"{path_main}complete_{file_type}.csv", index=False)

            # When csv exists already
            else:

                print(f"Reading an existing csv: {existing_csv}")
                df = pd.read_csv(existing_csv[0])

        return df

    def create_custom_combo(self, sel_regions, labels):
        """
        For all of the selected regions, make combo plot with corresponding lightcurve
        :param sel_regions: Regions selected to create
        :param labels: Labels to give each of the lightcurves
        :return: Returns path to folder
        """
        from matplotlib import gridspec as gs

        # Count total images from e.g Combo -> This may be incorrect
        combo_images = sorted(glob.glob(f"{self.global_save_path}Combo/*.png"))

        lcurve_data = np.load(
            f"{self.global_save_path}lcurves.npy"
        )  # This may not exist
        pix_data = np.load(f"{self.global_save_path}pix.npy")
        time_npy = np.load(f"{self.global_save_path}times.npy")
        tdf = pd.DataFrame({"Time": time_npy})
        tdf.Time = pd.to_datetime(tdf.Time)
        timearr = tdf.Time

        # Now need to separate lightcurves as required:
        width_ratios = []

        for reg in range(len(sel_regions)):
            width_ratios.append(1)

        width_ratios.append(1.25 * len(sel_regions))

        grs = gs.GridSpec(
            len(sel_regions), len(sel_regions) + 1, width_ratios=tuple(width_ratios)
        )  # May look horrible but fast to check

        # First need to check maxima and minima for plotting
        found_max = 0
        found_min = 500

        found_max_pix = 0
        found_min_pix = 1000

        for curr_sec in sel_regions:
            r, c = AstroMap.translate_number(int(curr_sec), self.nrowcol)
            total_line = lcurve_data[r, c, :]

            if total_line.max() > found_max:
                found_max = total_line.max()

            if total_line.min() < found_min:
                found_min = total_line.min()

        for curr_sec in sel_regions:
            r, c = AstroMap.translate_number(int(curr_sec), self.nrowcol)
            total_line = pix_data[r, c, :]

            if total_line.max() > found_max:
                found_max_pix = total_line.max()

            if total_line.min() < found_min:
                found_min_pix = total_line.min()

        def proc_multi(split_indices=None):
            """
            Function created for multiprocessing of getting images and lcurves together
            :param split_indices: Which indexes to use
            :return: Saves figures into given path
            """
            plt.rcParams.update({"font.size": 20})

            for i in split_indices:
                counter = 0
                curr_time = timearr[i]
                fig = plt.figure(figsize=(19, 10))

                for rel_section, label in zip(sel_regions, labels):
                    numpyrow, numpycol = self.translate_number(
                        int(rel_section), self.nrowcol
                    )
                    curr_aia = plt.imread(combo_images[i])
                    curr_lcurve = lcurve_data[numpyrow, numpycol, :]
                    curr_pix = pix_data[numpyrow, numpycol, :]

                    plt.subplot(grs[0 : len(sel_regions), 0 : len(sel_regions)])
                    plt.imshow(curr_aia)
                    plt.title(f"AIA {self.wvlnth} : {curr_time}")
                    plt.axis("off")

                    # Plot lightcurve and corresponding column
                    plt.subplot(grs[counter : counter + 1, len(sel_regions)])
                    ax = plt.gca()
                    plt.plot(timearr, curr_lcurve, linewidth=2, color="black")
                    ax.set_ylabel("dn/s")
                    plt.ylim(found_min, found_max)
                    plt.title(label)
                    ax.axvspan(
                        ymin=0,
                        ymax=1,
                        xmin=curr_time,
                        xmax=curr_time + timedelta(0, 2),
                        color="black",
                        linestyle="--",
                        alpha=0.2,
                        linewidth=2,
                    )
                    plt.gcf().autofmt_xdate()

                    counter += 1

                plt.tight_layout()
                os.makedirs(f"{self.global_save_path}Selected_Regions/", exist_ok=True)
                plt.savefig(f"{self.global_save_path}Selected_Regions/{i:03d}.png")
                plt.close(fig)

                print(f"Done index : {i}")

        # Multiprocessing requirements: Array of arrays and process list
        vector_index_div = np.array_split(
            np.array(np.arange(len(combo_images))), self.n_cpus
        )
        proc = []

        for list_indices_proc in vector_index_div:
            if list_indices_proc.size:
                pr = Process(target=proc_multi, args=([list_indices_proc]))
                proc.append(pr)
                pr.start()

        for process in proc:
            process.join()

    def create_custom_imageless(self, sel_regions, labels):
        """
        For all of the selected regions, make combo plot with corresponding lightcurve
        :param sel_regions: Regions selected to create
        :param labels: Labels to give each of the lightcurves
        :param nrowcol: Number of rows and columns
        :return: Returns path to folder6
        """
        from matplotlib import gridspec as gs

        # Count total images from e.g Combo)
        lcurve_data = self.complete_lcurves  # Dataframe containing everything
        tdf = pd.DataFrame({"Time": lcurve_data["Time"]})
        tdf.Time = pd.to_datetime(tdf.Time)
        timearr = tdf.Time

        itotal = np.arange(len(timearr))  # This will check each of the timesteps

        # This should create a found maxima, found minima
        found_max = 0
        found_min = 500

        found_max_pix = 0
        found_min_pix = 1000

        for curr_sec in sel_regions:

            total_line = lcurve_data[f"{curr_sec}"]

            if total_line.max() > found_max:
                found_max = total_line.max()

            if total_line.min() < found_min:
                found_min = total_line.min()

        # Then cycle through all of the timesteps
        for i in itotal:
            counter = 0
            # This is true
            curr_time = timearr[i]
            fig = plt.figure(figsize=(19, 10))

            for rel_section, label in zip(sel_regions, labels):
                curr_lcurve = lcurve_data[f"{counter}"]
                # Plot lightcurve and corresponding column
                plt.subplot(self.nrowcol, self.nrowcol, counter + 1)
                ax = plt.gca()
                plt.plot(timearr, curr_lcurve, linewidth=2, color="black")
                ax.set_ylabel("dn/s")
                plt.ylim(found_min, found_max)
                plt.title(label)

                ax.axvspan(
                    ymin=0,
                    ymax=1,
                    xmin=curr_time,
                    xmax=curr_time + timedelta(0, 2),
                    color="black",
                    linestyle="--",
                    alpha=0.2,
                    linewidth=2,
                )
                plt.gcf().autofmt_xdate()

                counter += 1
            plt.tight_layout()
            os.makedirs(f"{self.global_save_path}Selected_Regions/", exist_ok=True)
            plt.savefig(f"{self.global_save_path}Selected_Regions/{i:03d}.png")
            plt.close(fig)

            print(f"Done index : {i}")

    def create_emd_from_lightcurve(self, label, figsize=(20, 28)):

        from matplotlib import gridspec as gs
        from PyEMD import EMD
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from pandas.plotting import register_matplotlib_converters

        register_matplotlib_converters()
        import glob
        from datetime import timedelta
        from matplotlib import rc
        from multiprocessing import Process

        emd = EMD()

        font = {"family": "DejaVu", "weight": "normal", "size": 16}

        rc("font", **font)

        numpyrow, numpycol = self.translate_number(int(label), self.nrowcol)
        picture_folder_id = f"{numpycol:02d}_{numpyrow:02d}"

        path_to_lcurves = self.global_save_path
        path_to_picture = f"{self.global_save_path}{picture_folder_id}/"

        dataset = pd.read_csv(f"{path_to_lcurves}complete_lcurves.csv")
        time = pd.to_datetime(dataset["Time"])
        region = dataset[f"{label}"]
        imfs = emd(region.values)

        image_list = sorted(glob.glob(f"{path_to_picture}*.png"))
        vector_index_div = np.array_split(
            np.array(np.arange(len(image_list))), self.n_cpus
        )
        proc = []

        colors = ["black", "blue", "red"]
        fig = plt.figure(figsize=(10, 20))

        time_curr = time[0]
        lineset = {}

        min_height = 0
        max_height = 0

        for index, i in enumerate(imfs):

            if index < len(imfs) - 1:
                if i.min() < min_height:
                    min_height = i.min()
                if i.max() > max_height:
                    max_height = i.max()

        # Draw all IMFs
        for index, i in enumerate(imfs):
            if index != len(imfs) - 1:
                fig.add_subplot(len(imfs), 1, len(imfs) - 1 - index)
                ax = plt.gca()
                plt.plot(time, i, color=colors[0], linewidth=4, label=f"IMF {index}")
                ybot, ytop = ax.get_ylim()
                delta_y = ytop - ybot
                move_by = (max_height - min_height) / 2
                ax.set_ylim(ybot - move_by, ytop + move_by)

                ax.yaxis.set_label_position("right")
                ax.set_ylim(min_height, max_height)
                ax.set(ylabel=f"IMF {index + 1}", xlabel="", xticks=[])
                lineset[f"{index}"] = ax.axvline(
                    x=time_curr, color="black", alpha=0.4, linewidth=5
                )

        # Draw original signal + last IMF (residual)
        fig.add_subplot(len(imfs), 1, len(imfs))
        ax1 = plt.gca()
        plt.plot(time, region, color=colors[1], linewidth=5, label="Signal")
        plt.plot(
            time,
            imfs[len(imfs) - 1],
            color=colors[2],
            linestyle="dashed",
            label="Residual",
            linewidth=4,
        )
        lineset["signal"] = ax1.axvline(
            x=time_curr, color="black", alpha=0.4, linewidth=5
        )
        ax1.set(ylabel=f"Original Signal")

        import matplotlib.dates as mdates
        import matplotlib.units as munits
        import datetime

        converter = mdates.ConciseDateConverter()
        munits.registry[np.datetime64] = converter
        munits.registry[datetime.date] = converter
        munits.registry[datetime.datetime] = converter
        plt.gcf().autofmt_xdate()
        plt.legend()
        ax1.yaxis.set_label_position("right")

        save_emd = f"{path_to_picture}/emd_results/"
        os.makedirs(save_emd, exist_ok=True)

        def multi_proc(index_sublist, save_path=save_emd):
            for curr_obs in index_sublist:
                # Draw all of the lines#
                for k in list(lineset):
                    lineset[k].set_xdata(time[curr_obs])
                fig.canvas.draw()
                plt.savefig(f"{save_path}{curr_obs:03d}")

        for list_indices_proc in vector_index_div:
            if list_indices_proc.size:
                # noinspection PyTypeChecker
                pr = Process(target=multi_proc, args=([list_indices_proc]))
                proc.append(pr)
                pr.start()

        for process in proc:
            process.join()

        # Then combine all of these

        from Scripts.Imports.Data_analysis.Tools import imgrecomb_any

        imgrecomb_any(
            path_1=f"{path_to_picture}",
            path_2=f"{save_emd}",
            save_path=f"{path_to_picture}with_EMD/",
            figsize=figsize,
        )


def load_map(index, aia_folder, prep="yes"):
    """
    :param index: Index from 0 to length of folder of file to process
    :param aia_folder: Folder where AIA files are located
    :param prep: Defaults to 'yes', allows to use sunpy prep
    :return: Returns map structure at the index
    """
    import glob
    import sunpy.map
    from sunpy.instr.aia import aiaprep as pr_aia

    aia_file = (sorted(glob.glob(f"{aia_folder}/*.fits")))[index]
    aiapr = sunpy.map.Map(aia_file)

    if prep == "no":
        aiapr = pr_aia(aiapr)

    return aiapr


def aia_locator_gen(prepped_map):
    """
    This function takes a prepped AIA map and outputs the locator element for use with SSW in IDL
    :param prepped_map: aia map object
    :return: Returns a string containing necessary information to find the map
    """
    aia_obs = prepped_map.coordinate_frame.obstime.datetime

    yr = aia_obs.year
    mth = aia_obs.month
    day = aia_obs.day

    hour = aia_obs.hour
    mn = aia_obs.minute

    if day < 10:
        day = f"0{day}"
    else:
        day = day

    if hour < 10:
        hour = f"0{hour}"
    else:
        hour = hour

    if mn < 10:
        mn = f"0{mn}"
    else:
        mn = mn

    locator = f"{yr}-{mth}-{day}_{hour}:{mn}"

    return locator


def create_point_df(id, selected_files, n_vec_soln):
    """
    :param id: The locator of the text file
    :param selected_files: This function takes the well defined files as an input
    :param n_vec_soln: The number of vectors in the current case
    :return: Returns a dataframe with the correct values
    """

    import pandas as pd
    import numpy as np

    # Now that CSV selected_files are correctly fetched:
    elat = pd.read_csv(selected_files[0], header=None)
    elong = pd.read_csv(selected_files[1], header=None)
    r_ends = pd.read_csv(selected_files[2], header=None)
    start = pd.read_csv(selected_files[3], header=None)
    nstep_file = pd.read_csv(selected_files[4], header=None)

    timelist, nstep_list, st_r, st_lat, st_long, r_endlist, lat_ends, long_ends = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )

    for vec_LINE in range(n_vec_soln):
        # We now have selected the length of points! - Only useful if want to get height profile

        starray = start[vec_LINE : vec_LINE + 1]
        radius = starray.iloc[0, 0]
        start_lat = starray.iloc[0, 1] * (180 / np.pi)  # This looks promising
        start_long = starray.iloc[0, 2] * 180 / np.pi

        # Whole vectors!
        r_row = r_ends[vec_LINE : vec_LINE + 1]
        lat_row = elat[vec_LINE : vec_LINE + 1]
        long_row = elong[vec_LINE : vec_LINE + 1]
        nstep_row = nstep_file[vec_LINE : vec_LINE + 1]

        # Select the first point of each row
        r_point = r_row.iloc[0, 0]
        lat_endpoint = lat_row.iloc[0, 0]
        long_endpoint = long_row.iloc[0, 0]
        curr_n_step = nstep_row.iloc[0, 0]

        timelist.append(id)
        nstep_list.append(curr_n_step)
        st_r.append(radius)
        st_lat.append(start_lat)
        st_long.append(start_long)
        r_endlist.append(r_point)
        lat_ends.append(lat_endpoint)
        long_ends.append(long_endpoint)

    time = np.asarray(timelist)
    nstep_list = np.asarray(nstep_list)
    start_r = np.asarray(st_r)
    start_lat = np.asarray(st_lat)
    start_long = np.asarray(st_long)
    end_r = np.asarray(r_endlist)
    end_lat = np.asarray(lat_ends)
    end_long = np.asarray(long_ends)

    end_lat = -end_lat
    df = pd.DataFrame(
        {
            "obs_time": time,
            "n_steps": nstep_list,
            "start_r": start_r,
            "start_lat": start_lat,
            "start_long": start_long,
            "end_r": end_r,
            "end_lat": end_lat,
            "end_long": end_long,
        }
    )

    return df


def ROI_def(xylnht):
    """
    :param xylnht: X, Y, length, height in arcsec
    :return: Array containing [x0, xf, y0, yf]
    """
    if xylnht[2] < 0:
        raise ValueError("Please input positive length and/or height")

    import astropy.units as u

    x = xylnht[0] * u.arcsec
    y = xylnht[1] * u.arcsec
    lth = xylnht[2] * u.arcsec
    hth = xylnht[3] * u.arcsec
    x0 = x - lth
    xf = x + lth
    y0 = y - hth
    yf = y + hth

    r_array = [x0, xf, y0, yf]

    return r_array


def ROI_defCoords(ROI_Lon_Lat_Wd_Ht):
    """
    Define a ROI on Map Coordinates
    """
    import astropy.units as u

    lon0, lat0, width, height = ROI_Lon_Lat_Wd_Ht
    Coords = {
        "lon": (lon0 * u.deg, (lon0 + width) * u.deg),
        "lat": (lat0 * u.deg, (lat0 + height) * u.deg),
    }

    return Coords


def split_roi(base_ROI, nrowcol, pix=False, header=None):
    """
    This function takes a region of interest and splits it in an equal amount of rows and columns
    :param base_ROI: Region of interest which is being divided
    :param nrowcol: How many chunks
    :param pix: Whether or not to return in pixels
    :param header: If using pixels, need a header
    :return: returns an array of positions
    """
    total_rois = []
    x0, xf, y0, yf = base_ROI[0], base_ROI[1], base_ROI[2], base_ROI[3]

    # Calculate the differences
    dx = (xf - x0) / nrowcol
    dy = (yf - y0) / nrowcol

    # Starts from the bottom. Goes right, then up 1 row
    for row in range(nrowcol):
        # Is there an extra row here?
        y_top = [yf - dy * row]
        y_bot = [yf - dy * (row + 1)]  # 0 to n-1 bottom to top

        for col in range(nrowcol):
            x_left = [x0 + dx * col]  # 0 to n-1 left to right
            x_right = [x0 + dx * (col + 1)]  # 0 to n-1 left to right

            if pix:
                x_left_pix, y_bot_pix = arc2px(x_left[0].value, y_bot[0].value, header)
                x_right_pix, y_top_pix = arc2px(
                    x_right[0].value, y_top[0].value, header
                )

                total_rois.append([x_left_pix, x_right_pix, y_bot_pix, y_top_pix])

            else:
                total_rois.append([x_left, x_right, y_bot, y_top])

    return total_rois


def splitCoords(base_ROI, nrowcol):
    """
    Alternative to splitROI that splits Coords in LON, LAT
    """
    import astropy.units as u

    total_rois = []
    regionN = 0

    lon0, lonF, lat0, latF = (
        base_ROI["lon"][0],
        base_ROI["lon"][1],
        base_ROI["lat"][1],
        base_ROI["lat"][0],
    )

    dlon = (lonF - lon0) / nrowcol
    dlat = (latF - lat0) / nrowcol

    for row in range(nrowcol):
        lat_top = [lat0 + dlat * row]
        lat_bot = [lat0 + dlat * (row + 1)]

        for col in range(nrowcol):
            lon_left = [lon0 + dlon * col]
            lon_right = [lon0 + dlon * (col + 1)]

            total_rois.append(
                [
                    regionN,
                    {
                        "lon": (lon_left * u.deg, lon_right * u.deg),
                        "lat": (lat_bot * u.deg, lat_top * u.deg),
                    },
                ]
            )

            regionN += 1  # Numbering of regions from top left > right > down

    return total_rois


def track_ROI(init_submap, fullmap):
    """
    :param init_submap: Subregion of map to crop
    :param fullmap: Map that is to be cropped
    :return: New coordinates for new map
    """
    from sunpy.physics.differential_rotation import solar_rotate_coordinate

    # Subtract Time of current map to previous map
    rotated_coord = solar_rotate_coordinate(init_submap.center, fullmap.date)
    dlon = rotated_coord.Tx - init_submap.center.Tx
    dlat = rotated_coord.Ty - init_submap.center.Ty

    # Get new coordinates for box
    bl_new_lon = init_submap.bottom_left_coord.Tx + dlon
    bl_new_lat = init_submap.bottom_left_coord.Ty + dlat
    tr_new_lon = init_submap.top_right_coord.Tx + dlon
    tr_new_lat = init_submap.top_right_coord.Ty + dlat

    # Return the new box coordinates
    new_coords = [
        bl_new_lon,
        tr_new_lon,
        bl_new_lat,
        tr_new_lat,
    ]  # These are on arcsec quantities!
    return new_coords


def getSegment(CoordsHGS, all_hgs_coords):
    """
    Pass a single Coords object, return values for it
    """
    latmin, latmax = CoordsHGS["lat"]
    lonmin, lonmax = CoordsHGS["lon"]
    segment = np.logical_and(
        np.logical_and(all_hgs_coords.lon > lonmin, all_hgs_coords.lon < lonmax),
        np.logical_and(all_hgs_coords.lat > latmin, all_hgs_coords.lat < latmax),
    ).nonzero()
    return segment


def calculateLons(initial_submap, curr_submap, Coords):
    """
    Calculates where the coordinates should be given a time difference between current and initial
    Coords is the Coord = {"lon":(lon0,lonf), "lat":(lat0,latf)} object
    """

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import RotatedSunFrame

    # Subtract Time of current map to previous map

    dinit = initial_submap.date
    dcurr = curr_submap.date
    dt = dcurr - dinit

    startlon0, startlonF = Coords["lon"][0], Coords["lon"][1]
    dlon = startlonF - startlon0
    startlat0, startlatF = Coords["lat"][0], Coords["lat"][1]

    start_coord = SkyCoord(
        startlon0,
        Coords["lat"][0],
        frame="heliographic_stonyhurst",
        obstime=dinit,
    )
    rotated_coord = RotatedSunFrame(
        base=start_coord, duration=dt.to(u.second)
    ).transform_to(start_coord.frame)

    rotlon0 = rotated_coord.lon
    rotlonF = rotated_coord.lon + dlon

    rotCoords = {
        "lon": (rotlon0, rotlonF),
        "lat": (startlat0, startlatF),
    }  # Latitude is constant

    return rotCoords


def astro_track_ROI(ref_map, newmap, observer, large_roi):
    """
    This is a knock-off version to calculate movement for a given amount of time. Use SUNPY?
    :param ref_map: Subregion of map to crop
    :param newmap: Map that is to be cropped
    :param large_roi: Large roi_ref
    :return: New coordinates for new map
    """
    import datetime
    from sunpy.physics.differential_rotation import solar_rotate_coordinate
    from astropy.coordinates import SkyCoord

    t_ref = datetime.datetime.strptime(
        ref_map.fits_header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"
    )
    t_new = datetime.datetime.strptime(
        newmap.fits_header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f"
    )

    x_0 = large_roi[0]
    x_f = large_roi[1]
    y_0 = large_roi[2]
    y_f = large_roi[3]

    dx = x_f - x_0
    dy = y_f - y_0

    center_x, center_y = dx / 2 + x_0, dy / 2 + y_0
    start_coord = SkyCoord(
        center_x, center_y, frame="helioprojective", obstime=t_ref, observer=observer
    )
    rotated_center = solar_rotate_coordinate(start_coord, time=t_new)

    new_coords = [
        rotated_center.Tx - dx / 2,
        rotated_center.Tx + dx / 2,
        rotated_center.Ty - dy / 2,
        rotated_center.Ty + dy / 2,
    ]
    return new_coords


def arc2px(x_arc, y_arc, header):
    """Converts x and y arcsec coords into px."""
    import numpy as np

    try:
        xpx = [0] * len(x_arc)
        ypx = [0] * len(y_arc)
    except TypeError:
        xpx = [0]
        ypx = [0]
        x_arc = [x_arc]
        y_arc = [y_arc]

    for i in range(0, len(xpx)):
        xpx[i] = int(
            np.round(
                ((x_arc[i] - header["crval1"]) / header["cdelt1"]) + header["crpix1"]
            )
        )

    for j in range(0, len(ypx)):
        ypx[j] = int(
            np.round(
                ((y_arc[j] - header["crval2"]) / header["cdelt2"]) + header["crpix2"]
            )
        )

    return xpx, ypx


def track_arcsec(refmap, fullmap):
    """
    By Alexander James (UCL)
    Given a cropped reference map, will return the im cropped in a way that
    tracks solar rotation.
    Example:
        refmap = sunpy.map.Map(sfiles[0])
        bottomleft_ref = [-500*u.arcsec, -600*u.arcsec]
        topright_ref   = [0*u.arcsec, -200*u.arcsec]
        refmap = alexpy.submap(refmap, bottomleft_ref, topright_ref)
        fullmap = (sunpy.map.Map(sfile))
        fullmap = track_arcsec(refmap, fullmap)
    """
    from sunpy.physics.differential_rotation import solar_rotate_coordinate
    from astropy.coordinates import SkyCoord

    # Subtract Time of current map to previous map
    rotated_coord = solar_rotate_coordinate(refmap.center, fullmap.date)
    dlon = rotated_coord.Tx - refmap.center.Tx
    dlat = rotated_coord.Ty - refmap.center.Ty
    # Get new coordinated for box
    bl_new_lon = refmap.bottom_left_coord.Tx + dlon
    bl_new_lat = refmap.bottom_left_coord.Ty + dlat
    tr_new_lon = refmap.top_right_coord.Tx + dlon
    tr_new_lat = refmap.top_right_coord.Ty + dlat
    # Make new coordinates into new coord frame
    bl_new = SkyCoord(bl_new_lon, bl_new_lat, frame=fullmap.coordinate_frame)
    tr_new = SkyCoord(tr_new_lon, tr_new_lat, frame=fullmap.coordinate_frame)
    # Crop map to new location
    fullmap = fullmap.submap(bl_new, tr_new)
    # Give new, fixed map
    return fullmap
