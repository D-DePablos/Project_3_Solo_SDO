import os
import sys
from copy import deepcopy

import sunpy.map
from aiapy.calibrate import register, update_pointing, normalize_exposure
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from pandas.plotting import register_matplotlib_converters

# import matplotlib.colors as colors
from astropy import units as u
from astropy.coordinates import SkyCoord

register_matplotlib_converters()
# This path indicates to main python console where to get Scripts!
# sys.path.append('/home/diego/Documents/PhD')
# Assuming cwd is where ~helpers.py is located
module_path = os.path.abspath(os.path.join(f"{os.getcwd()}/../../../"))

for module in [module_path]:
    if module not in sys.path:
        sys.path.append(module)

# sys.path.append("/home/diegodp/Documents/PhD/Paper_1/Proj_1/")
sys.path.append(f"{os.getcwd()}")
from Proj_1_Imports import MapTools as Mtool


def translate_number(input_region, n_by_n):
    """
    Take an input region and a given n_by_n to return two dimensional response to where we are on array

    :param input_region: Which number we require
    :param n_by_n: Dimensions of array
    :return: Returns height, length at which box is located
    """

    array = np.arange(0, n_by_n**2, 1)
    array = np.split(array, n_by_n)
    array = np.vstack(array)

    return [
        np.where(array == input_region)[0][0],
        np.where(array == input_region)[1][0],
    ]


def astro_plot(
        file_list,
        rel_indexes,
        all_images,
        wvln,
        ROI_LATLON,
        n_row_col,
        save_path,
        threshold=(False, 1),
        regions_highlight=(),
        showPreview=False,
):
    """
    This function is analogous to sunpy plotting, but uses astropy maps instead. Optimised for     parallel processing
    :param file_list: List of files from glob
    :param rel_indexes: Which section this process takes care of
    :param all_images: List containing all files
    :param wvln: Wavelength (unitless)
    :param main_roi: Region of interest to be separated: [x0, y0, length, height]
    :param n_row_col: Amount of rows and columns to separate in
    :param vmin_vmax: min and max exposure for the wavelength
    :param save_path: Path to save figures and lightcurves to
    :param threshold: Whether to set a threshold or not. If true, then input?
    :param regions_highlight: Which regions to highlight
    """
    import datetime

    # Need to open the first fits file!
    # Create a n x n x n_obs array
    os.makedirs(save_path, exist_ok=True)

    # Need to define ROI and exposure limits; Create reference submap to then with respect to it
    large_roi = Mtool.ROI_defCoords(ROI_LATLON)
    large_roi_ref = large_roi.copy()
    smap_ref = sunpy.map.Map(all_images[0])  # Fetch the first map
    all_hpc = sunpy.map.all_coordinates_from_map(smap_ref)
    all_hgs = all_hpc.transform_to("heliographic_stonyhurst")

    # Set out empty time and lcurve ndarrays
    sub_time_array = np.array(np.arange(len(all_images)), dtype="str")
    sub_lcurve = np.zeros([n_row_col, n_row_col, len(all_images)])
    sub_pixel_array = np.zeros([n_row_col, n_row_col, len(all_images)])

    # For each of the images, we need to advance the regions
    for relevantIndex, file in zip(rel_indexes, file_list):
        m = sunpy.map.Map(file)
        m_updated_pointing = update_pointing(m)
        m_registered = register(m_updated_pointing)
        aia_curr = normalize_exposure(m_registered)
        aia_curr_base_copy = deepcopy(aia_curr)
        header = aia_curr.fits_header
        sub_time_array[relevantIndex] = datetime.datetime.strptime(
            header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f")

        # If not the first Map, we must track the Region of Interest
        if relevantIndex > 0:
            large_roi = Mtool.calculateLons(smap_ref,
                                            aia_curr,
                                            Coords=large_roi_ref)
            pass

        small_rois = Mtool.splitCoords(large_roi, n_row_col)
        # largeSegment = Mtool.getSegment(large_roi, aia_curr, all_hgs_coords=all_hgs)

        # These are exclusively the images. The order will hold after vector works!
        for regionIndex, segmentCoord in small_rois:
            if regionIndex in regions_highlight:
                smallCurrSegment = Mtool.getSegment(segmentCoord,
                                                    aia_curr,
                                                    all_hgs_coords=all_hgs)
                # We translate the number into new form
                row, col = translate_number(regionIndex, n_row_col)

                # Rotate vector such that it matches the image when 'lower' is used
                smallRegData = aia_curr.data[smallCurrSegment]
                # subreg_init = smallRegData.copy()

                # Calculate some threshold and use it to only display much higher than given
                if threshold[0]:
                    minus_plus_mean = np.array([
                        smallRegData.mean() -
                        threshold[1] * smallRegData.std(),
                        smallRegData.mean() +
                        threshold[1] * smallRegData.std(),
                    ])

                    above_two_sd = smallRegData > minus_plus_mean[1]
                    below_two_sd = smallRegData < minus_plus_mean[1]

                    smallRegData[below_two_sd] = 0
                    ref_mean = smallRegData[above_two_sd].mean()
                    smallRegData[below_two_sd] = np.nan
                    n_pixels = np.count_nonzero(~np.isnan(smallRegData))
                    sub_pixel_array[row, col, relevantIndex] = n_pixels

                elif not threshold[0]:
                    n_pixels = np.count_nonzero(~np.isnan(smallRegData))
                    sub_pixel_array[row, col, relevantIndex] = n_pixels
                    ref_mean = smallRegData.mean()

                sub_lcurve[row, col, relevantIndex] = ref_mean
                aia_curr.data[smallCurrSegment] = 0

                # TODO : Make a figure using a copy of the base AIA plot, with the area that is selected blacked out
                #        Save as row_col
                # curr_subRegion_savePath = f"{save_path}/{col:02d}_{row:02d}/"

        # Showing a preview to check whether working
        # aia_curr.data[segment] = 0
        bottom_left = SkyCoord(250 * u.arcsec,
                               350 * u.arcsec,
                               frame=aia_curr.coordinate_frame)
        submap = aia_curr.submap(bottom_left=bottom_left,
                                 width=500 * u.arcsec,
                                 height=500 * u.arcsec)
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=submap)
        submap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
        plt.savefig(f"{save_path}{relevantIndex:04d}.png")
        # plt.show()
        plt.close()

        # smallRegData = subreg_init
        # TODO : Set up a second map that does not black out
        submap_base_copy = aia_curr_base_copy.submap(bottom_left=bottom_left,
                                                     width=500 * u.arcsec,
                                                     height=500 * u.arcsec)
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=submap_base_copy)
        submap_base_copy.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
        # plt.axis("off")
        # plt.show()
        filledFolder = f"{save_path}Filled/"
        os.makedirs(filledFolder, exist_ok=True)
        plt.savefig(f"{filledFolder}{relevantIndex:04d}.png")
        plt.close()

    # Now select the relevant chunk of the vector and cut it to that size. Then we can append it onto the larger one?
    sub_time_array = sub_time_array[rel_indexes]
    sub_lcurve = sub_lcurve[:, :, rel_indexes]
    sub_pixel_array = sub_pixel_array[:, :, rel_indexes]

    np.save(f"{save_path}lcurves_{rel_indexes[-1]:04d}.npy", sub_lcurve)
    np.save(f"{save_path}times_{rel_indexes[-1]:04d}.npy", sub_time_array)
    np.save(f"{save_path}pix_{rel_indexes[-1]:04d}.npy", sub_pixel_array)


class NamedVector:
    def __init__(self, ndvector, label, cadence, relevant_imfs):
        """
        This class is used to store vectors which contain relevant information about all the data, and the relevant imfs
        These vectors can be IMFs, etc.
        :param ndvector: N-dimensional numphy array to store within this class
        :param label: Label for plotting
        """

        self.ndvector = ndvector.copy()
        self.label = label
        self.selected_imfs = relevant_imfs
        sum_built = np.zeros(shape=self.ndvector[0].shape)
        self.sum_imfs = sum_built
        self.time = np.arange(0, len(sum_built) * cadence, cadence)
        self.cadence = cadence
        self.flag_all = False

        if self.selected_imfs == -2:
            self.flag_all = True

    def __plot__(self, show=True):
        plt.clf()
        plt.figure()
        ax = plt.gca()

        ax.plot(self.time, self.sum_imfs)
        ax.set_title(f"{self.label} : Sum of IMFs {self.selected_imfs}")
        plt.show()

    def extract_relevant_function(self):

        for relevantIndex, timeseries in enumerate(self.ndvector):

            if self.flag_all:
                if relevantIndex == 0:
                    self.selected_imfs = []

                self.selected_imfs.append(f"{relevantIndex}")
                self.sum_imfs = self.sum_imfs + timeseries

            elif relevantIndex + 1 in self.selected_imfs:
                self.sum_imfs = self.sum_imfs + timeseries
