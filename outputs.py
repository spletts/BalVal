"""
Name and write to output files.
E.g. log files, 
region files for objects,
plot filename.

Comments are ABOVE the code they refer to.
"""

from collections import OrderedDict
import csv
import numpy as np
import os
import sys

from set_constants import *


PLOT_TYPES = OrderedDict({'cornerhist2d':CORNER_HIST_2D, 'scatter':SCATTER, 'hist2d':HIST_2D, 'completeness':PLOT_COMPLETENESS, 'cbar_cm_t':CM_T_CBAR, 'cbar_cm_t_err':CM_T_ERR_CBAR, 'hexbin':HEXBIN})
if PLOT_FLUX is False:
    # Find which plot type is `True`. Note if `PLOT_FLUX` these are all `False` #
    IDX = np.where(PLOT_TYPES.values())[0][0]
    DISPLAY = PLOT_TYPES.keys()[IDX]
if PLOT_FLUX:
    DISPLAY = 'histogram'



def get_directory(balrog_run, low_level_dir, match_type, output_directory, realization, tile):
    """Get path to directory where various output files are to be saved.
    For directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    low_level_dir (str OR list of str)
        The directory structure following `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __dir (str)
        Complete path to directory.
    """

    if isinstance(low_level_dir, str):
        __dir = os.path.join(output_directory, 'outputs', balrog_run, match_type, tile, realization, low_level_dir)

    if isinstance(low_level_dir, list):
        if low_level_dir[1] == 'magnitude' and NORMALIZE_MAG:
            __dir = os.path.join(output_directory, 'outputs', balrog_run, match_type, tile, realization, low_level_dir[0], low_level_dir[1], 'normalized')
        else:
            __dir = os.path.join(output_directory, 'outputs', balrog_run, match_type, tile, realization, low_level_dir[0], low_level_dir[1])


    # FOF analysis #
    if RUN_TYPE is not None:
        __dir = os.path.join(output_directory, 'outputs', balrog_run, match_type, tile, realization, low_level_dir, 'fof_analysis')

    # Check for directory existence #
    if os.path.isdir(__dir) is False:
        if NO_DIR_MAKE is False:
            sys.exit('Directory {} does not exist. Change directory structure in `outputs.get_directory()`.\n'.format(__dir))
        if NO_DIR_MAKE:
            if VERBOSE_ING: print 'Making directory {}...'.format(__dir)
            os.makedirs(__dir)

    return __dir




def get_log_filenames(balrog_run, match_type, output_directory, realization, tile):
    """Create filenames for various log files.
     For directory structure and log filenames, see https://github.com/spletts/BalVal/blob/master/README.md#full-output-directory-structure 

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_flag_log (str)
        Complete filename for the log file that records flags if `LOG_FLAGS=True`. 
        If `LOG_FLAGS=False` this file is not created; 'dummy.log' is created.
        Headers: TILE, REALIZATION, FLAG1_HEADER, FLAG2_HEADER,  FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2, RUN_TYPE.

    __fn_mag_err_log (str)
        Complete filename for the log file that records details of the magnitude error calculation.
        Headers: TILE, REALIZATION, BAND, NUM_OBJS_IN_BIN, MAG_{MEAS/TRUE}_BIN_LHS, MAG_{MEAS/TRUE}_BIN_RHS, MEDIAN_HAX_MAG_{MEAS/TRUE}, MEDIAN_MAG_ERROR

    __fn_color_log (str)
        Complete filename for the log file that details the color calculation. 
        Headers: TILE, REALIZATION, COLOR, MAG_{MEAS/TRUE}_BIN, OBJS_IN_MAG_BIN, TOTAL_OBJECTS_NO_FLAGS, RUN_TYPE.

    __fn_mag_diff_outliers_log (str)
        Complete filename for the log file that records the magnitude differences greater than a threshold.
        Threshold is given by `__cutoff_mag_diff_for_log` in `bin_and_cut_measured_magnitude_error()`. 
        Headers: TILE, REALIZATION, BAND, MAG1, MAG2, MAG_DIFF

    __fn_mag_completeness_log (str)
        Complete filename for the log file that details the magnitude completeness calculation. 
        Headers: TILE, REALIZATION, BAND, INJ_PERCENT, TRUTH_MAG_BIN_LHS, TRUTH_MAG_BIN_RHS, MATCH_CAT_OBJS_IN_BIN, TRUTH_CAT_OBJS_IN_BIN

    __fn_gauss_aper_log (str)
        Complete filename for the log file that records details of the Gaussian aperture method for measuring flux. 
        If `PLOT_FLUX=False` and `GAUSS_APER=False` this file is not created; 'dummy.log' is created.
        Headers: TILE, REALIZATION, BAND, FLUX_MEAS, FLUX_TRUE, GAUSS_APER_FLUX_MEAS, GAUSS_APER_FLUX_TRUE, FLUX_ERR_MEAS, GAUSS_APER_FLUX_DIFF/FLUX_ERR_MEAS


    __fn_flux_sigma_clip_log (str)
        Complete filename for the log file that records the details of the `N`-sigma clip of the normalized (to flux error) flux difference.
        Headers: TILE, REALIZATION, BAND, GAUSSIAN_APER_APPLIED, NUM_OBJS_FLAGS_RM, N-SIGMA_CLIP, NUM_OBJS_CLIPPED
        If `PLOT_FLUX=False` and `SIGMA_CLIP_NORM_FLUX_DIFF=False` this file not created; 'dummy.log' is created. 

    __fn_percent_recovered_log (str)
        Complete filename for the log file that records the percent of Balrog-injected objects recovered in the appropriate measured catalog.
        The percent recovered is computed in two ways: (1) without removing flagged* objects and (2) removing flagged* objects.
        For a description of flags, see https://github.com/spletts/BalVal/blob/master/README.md#flagged-objects
        Headers: TILE, REALIZATION, BAND, TOTAL_OBJS_IN_MATCH1AND2, TOTAL_FLAGGED_OBJS_IN_MATCH1AND2, TOTAL_FLAGGED_OBJS_IN_TRUTH, %_RECOVERED_FLAGS_IN, %_RECOVERED_FLAGS_RM, RUN_TYPE

    __fn_num_objs_in_1sig_mag_log (str)
        Complete filename for the log file that records the number of objects that lie within one-sigma of the measured magnitude error.
        If `PLOT_MAG=False`, this file is not created; 'dummy.log' is created.
        Headers: TILE, REALIZATION, BAND, TOTAL_OBJS_IN_MATCH1AND2, TOTAL_FLAGGED_OBJS_IN_MATCH1AND2, TOTAL_OBJS_IN_1SIGMA_MAG_MEAS, NORMALIZED, RUN_TYPE

    __fn_match_log (str)
        Complete filename for the log file that records then number of objects in various matched catalogs (set by STILTS parameter `join`)
        Headers: TILE, REALIZATION, IN1, IN2, TOTAL_MATCHES_1AND2, TOTAL_MATCHES_1NOT2, TOTAL_MATCHES_2NOT1
    """

    logFileDirectory = get_directory(tile=tile, realization=realization, low_level_dir='log_files', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    # Repeated #
    repeated_str = '{}_{}_{}'.format(tile, realization, match_type)

    # FOF analysis specifies 'ok' or 'rerun' #
    if RUN_TYPE is not None:
        repeated_str = '{}_{}_{}_{}'.format(tile, realization, match_type, RUN_TYPE)


    __fn_flag_log = os.path.join(logFileDirectory, '{}_flags.log'.format(repeated_str))

    __fn_mag_err_log = os.path.join(logFileDirectory, '{}_mag_error_computation.log'.format(repeated_str))

    if MATCH_CAT1 in ('gal_truth', 'star_truth') or MATCH_CAT2 in ('gal_truth', 'star_truth'):
        __fn_percent_recovered_log = os.path.join(logFileDirectory, '{}_percent_recovered.log'.format(repeated_str)) 
    if MATCH_CAT1 not in ('gal_truth', 'star_truth') and MATCH_CAT2 not in ('gal_truth', 'star_truth'):
        __fn_percent_recovered_log = os.path.join(logFileDirectory, 'dummy.log')

    if PLOT_COLOR:
        __fn_color_log = os.path.join(logFileDirectory, '{}_color_from_mag.log'.format(repeated_str))
    if PLOT_MAG or PLOT_FLUX:
        __fn_color_log = os.path.join(logFileDirectory, 'dummy.log')

    __fn_mag_diff_outliers_log = os.path.join(logFileDirectory, '{}_outlier_mag_diffs.log'.format(repeated_str))

    if PLOT_MAG:
        __fn_mag_completeness_log = os.path.join(logFileDirectory, '{}_mag_completeness.log'.format(repeated_str))
        __fn_num_objs_in_1sig_mag_log = os.path.join(logFileDirectory, '{}_num_objs_in_1sig_mag.log'.format(repeated_str))
    if PLOT_COLOR or PLOT_FLUX:
        __fn_mag_completeness_log = os.path.join(logFileDirectory, 'dummy.log')
        __fn_num_objs_in_1sig_mag_log = os.path.join(logFileDirectory, 'dummy.log')

    # For Gaussian aperture #
    if PLOT_FLUX:
        __fn_gauss_aper_log = os.path.join(logFileDirectory, '{}_flux_from_gaussian_aperture.log'.format(repeated_str))
        __fn_flux_sigma_clip_log = os.path.join(logFileDirectory, '{}_sigma_clip_flux.log'.format(repeated_str))
    if PLOT_MAG or PLOT_COLOR:
        __fn_mag_completeness_log = os.path.join(logFileDirectory, 'dummy.log')
        __fn_gauss_aper_log = os.path.join(logFileDirectory, 'dummy.log')
        __fn_flux_sigma_clip_log = os.path.join(logFileDirectory, 'dummy.log')

    __fn_match_log = os.path.join(logFileDirectory, '{}_matched_catalogs.log'.format(repeated_str))

    print 'Saving log file for flags as:\n----->{}'.format(__fn_flag_log)
    print 'Saving log file for magnitude error and bins used as:\n----->{}'.format(__fn_mag_err_log)
    if PLOT_COLOR: print 'Saving log file for color plot as:\n----->{}'.format(__fn_color_log)
    print 'Saving log file for outlier MagnitudeDifference as:\n----->{}'.format(__fn_mag_diff_outliers_log)
    print 'Saving log file for completeness plot as:\n----->{}'.format(__fn_mag_completeness_log)
    if PLOT_FLUX: print 'Saving log file for sigma clipping as:\n----->{}'.format(__fn_flux_sigma_clip_log)
    if PLOT_FLUX and PLOT_GAUSS_APER_FLUX: print 'Saving Gaussian aperture calculation details as:\n----->{}'.format(__fn_gauss_aper_log)


    return __fn_flag_log, __fn_mag_err_log, __fn_color_log, __fn_mag_diff_outliers_log, __fn_mag_completeness_log, __fn_gauss_aper_log, __fn_flux_sigma_clip_log, __fn_percent_recovered_log, __fn_num_objs_in_1sig_mag_log, __fn_match_log 




def get_color_plot_filename(balrog_run, match_type, output_directory, realization, tile):
    """Create the filename for the color plot. 
    If `SAVE_PLOT=True` the filename is used in `matplotlib.pyplot.savefig()`.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_plot (str)
        Complete filename for the color plot.
        If `SAVE_PLOT=True` this is used via `matplotlib.pyplot.savefig(__fn_plot)`
    """

    plotDirectory = get_directory(tile=tile, realization=realization, low_level_dir=['plots', 'color'], output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    # Find which plot type is `True` #
    __idx = np.where(PLOT_TYPES.values())[0][0] #FIXME use const at top
    __plot_type = DISPLAY 


    if COLOR_YLOW is None and COLOR_YHIGH is None:
        ylim = ''
    if COLOR_YLOW is not None and COLOR_YHIGH is not None:
        ylim = '{}y{}'.format(COLOR_YLOW, COLOR_YHIGH)

    if RUN_TYPE is None:
        __fn_end = '{}_{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, 'placehold', __plot_type, ylim)
    if RUN_TYPE is not None:
        __fn_end = '{}_{}_{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, RUN_TYPE, 'placehold', __plot_type, ylim)


    ### Get complete filename (including path) ###
    __fn_plot = os.path.join(plotDirectory, __fn_end)

    return __fn_plot




def get_flux_plot_filename(balrog_run, match_type, output_directory, realization, tile):
    """Create the filename for the flux plot. 
    If `SAVE_PLOT=True` the filename is used in `matplotlib.pyplot.savefig()`.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_plot (str)
        Complete filename for the flux plot.
        If `SAVE_PLOT=True` this is used via `matplotlib.pyplot.savefig(__fn_plot)`
    """

    plotDirectory = get_directory(tile=tile, realization=realization, low_level_dir=['plots', 'flux'], output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    if PLOT_GAUSS_APER_FLUX and SIGMA_CLIP_NORM_FLUX_DIFF:
        plot_type = '{}_sigma_clip_gauss_aper_norm_flux_diff_histogram'.format(int(N))
    if PLOT_GAUSS_APER_FLUX and RAW_NORM_FLUX_DIFF:
        plot_type = 'gauss_aper_norm_flux_diff_histogram'
    if PLOT_CM_FLUX and SIGMA_CLIP_NORM_FLUX_DIFF:
        plot_type = '{}_sigma_clip_norm_cm_flux_diff_histogram'.format(int(N))
    if PLOT_CM_FLUX and RAW_NORM_FLUX_DIFF:
        plot_type = 'norm_cm_flux_diff_histogram'


    if FLUX_XLOW is None and FLUX_XHIGH is None:
        xlim = ''
    if FLUX_XLOW is not None and FLUX_XHIGH is not None:
        xlim = '{}y{}'.format(FLUX_XLOW, FLUX_XHIGH)


    if RUN_TYPE is None:
        __fn_end = '{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, plot_type, xlim)
    if RUN_TYPE is not None:
        __fn_end = '{}_{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, RUN_TYPE, plot_type, xlim)


    ### Get complete filename (including path) ###
    __fn_plot = os.path.join(plotDirectory, __fn_end)

    return __fn_plot



def get_magnitude_plot_filename(balrog_run, match_type, output_directory, realization, tile):
    """Create the filename for the magnitude plot. 
    If `SAVE_PLOT=True` the filename is used in `matplotlib.pyplot.savefig()`.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_plot (str)
        Complete filename for the magnitude plot.
        If `SAVE_PLOT=True` this is used via `matplotlib.pyplot.savefig(__fn_plot)`
    """

    plotDirectory = get_directory(tile=tile, realization=realization, low_level_dir=['plots', 'magnitude'], output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    # Find which plot type is `True` #
    __idx = np.where(PLOT_TYPES.values())[0][0]
    __plot_type = DISPLAY 

    if MAG_YLOW is None and MAG_YHIGH is None:
        # Default scale for the vertical axis (vax) #
        ylim = ''
    if MAG_YLOW is not None and MAG_YHIGH is not None:
        ylim = '{}y{}'.format(MAG_YLOW, MAG_YHIGH)


    # `placehold` will be replaced with {band(s)/color} when `savefig()` is called within various functions #
    if RUN_TYPE is None:
        __fn_end = '{}_{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, 'placehold', __plot_type, ylim)
    if RUN_TYPE is not None:
        __fn_end = '{}_{}_{}_{}_{}_{}_{}.png'.format(tile, realization, match_type, RUN_TYPE, 'placehold', __plot_type, ylim)

    ### Get complete filename (including path) ###
    __fn_plot = os.path.join(plotDirectory, __fn_end)

    if NORMALIZE_MAG and PLOT_MAG:
        __fn_plot = os.path.join(plotDirectory, 'norm_{}'.format(__fn_plot))


    return __fn_plot





def get_plot_filename(balrog_run, match_type, output_directory, realization, tile):
    """Create file name of plot.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_plot (str)
        Complete filename for the plot.
        If `SAVE_PLOT=True` this is used via `matplotlib.pyplot.savefig(__fn_plot)`
    """

    if PLOT_COLOR:
        __fn_plot = get_color_plot_filename(balrog_run=balrog_run, match_type=match_type, output_directory=output_directory, realization=realization, tile=tile)


    if PLOT_FLUX:
        __fn_plot = get_flux_plot_filename(balrog_run=balrog_run, match_type=match_type, output_directory=output_directory, realization=realization, tile=tile)


    if PLOT_MAG:
        __fn_plot = get_magnitude_plot_filename(balrog_run=balrog_run, match_type=match_type, output_directory=output_directory, realization=realization, tile=tile)


    '''
    # By construction, completeness plots contain both 10% and 20% injected catalogs. Remove this distinction #
        if PLOT_COMPLETENESS:
        match_type = match_type.replace('10%_', ''); match_type = match_type.replace('20%_', '')
        if PLOT_COMPLETENESS is False: 
        match_type = match_type
    '''

    return __fn_plot




def get_region_filenames(balrog_run, match_type, output_directory, realization, tile):
    """Generate names for region files of different `join` types (`join` is a STILTS parameter). 

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_reg_1and2 (str)
        Complete filename for region file containing regions of objects in both `MATCH_CAT1` and `MATCH_CAT2` (`join=1and2` in `stilts_matcher`).

    __fn_reg_1not2 (str)
        Complete filename for region file containing regions of objects uniquely in `MATCH_CAT1` and not in `MATCH_CAT2` (`join=1not2` in `stilts_matcher`). 

    __fn_reg_2not1 (str)
        Complete filename for region file containing regions of objects uniquely in `MATCH_CAT2` and not in `MATCH_CAT1` (`join=2not1` in `stilts_matcher`). 
    """

    regionFileDirectory = get_directory(tile=tile, realization=realization, low_level_dir='region_files', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    # Repeated #
    repeated_str = '{}_{}_{}'.format(tile, realization, match_type)

    if RUN_TYPE is not None:
        repeated_str = '{}_{}_{}'.format(tile, realization, match_type, RUN_TYPE)

    __fn_reg_1and2 = os.path.join(regionFileDirectory, '{}_match1and2.reg'.format(repeated_str))
    __fn_reg_1not2 = os.path.join(regionFileDirectory, '{}_match1not2.reg'.format(repeated_str))
    __fn_reg_2not1 = os.path.join(regionFileDirectory, '{}_match2not1.reg'.format(repeated_str))


    if MAKE_REGION_FILES:
        print 'Saving region files as...'
        print '----->', __fn_reg_1and2
        print '----->', __fn_reg_1not2
        print '----->', __fn_reg_2not1

    return __fn_reg_1and2, __fn_reg_1not2, __fn_reg_2not1 




def get_matched_catalog_filenames(balrog_run, match_type, output_directory, realization, tile):
    """Create filenames for the matched catalogs that will be created in `stilts_matcher`.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_match_1and2 (str)
        Complete filename for region file containing regions of objects in both `MATCH_CAT1` and `MATCH_CAT2` (`join=1and2` in `stilts_matcher`).

    __fn_match_1not2 (str)
        Complete filename for region file containing regions of objects uniquely in `MATCH_CAT1` and not in `MATCH_CAT2` (`join=1not2` in `stilts_matcher`). 

    __fn_match_2not1 (str)
        Complete filename for region file containing regions of objects uniquely in `MATCH_CAT2` and not in `MATCH_CAT1` (`join=2not1` in `stilts_matcher`). 
    """

    matchedCatalogDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

    repeated_str = '{}_{}_{}'.format(tile, realization, match_type)

    __fn_match_1and2 = os.path.join(matchedCatalogDirectory, '{}_match1and2.csv'.format(repeated_str))
    __fn_match_1not2 = os.path.join(matchedCatalogDirectory, '{}_match1not2.csv'.format(repeated_str))
    __fn_match_2not1 = os.path.join(matchedCatalogDirectory, '{}_match2not1.csv'.format(repeated_str))

    return __fn_match_1and2, __fn_match_1not2, __fn_match_2not1
    



def get_ngmix_compatible_catalog_filename(balrog_run, fn_unmatched_cat, match_type, output_directory, realization, tile):
    """Create filename for the catalog that is reformatted to be compatible with `ngmix`.
    An `ngmix` compatible catalog is required to use https://github.com/esheldon/ngmixer/blob/master/ngmixer/gaussap.py

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    fn_unmatched_cat (str)
        Complete filename for an unmatched catalog.
        This refers to either the catalog described by `MATCH_CAT1, INJ1, INJ1_PERCENT` or `MATCH_CAT2, INJ2, INJ2_PERCENT`;
        this does not refer to catalogs created by `stilts_matcher`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __full_fn_ngmix_cat (str)
        Complete filename for the catalog to be used in `ngmixer.gaussap.get_gauss_aper_flux_cat()`
    """

    matchedCatalogDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type) 

    # '/data/des41.b/data/Balrog/COSMOS_HEX/DES0220-0250/real_0_DES0220-0250_mof.fits' --> 'real_0_DES0220-0250_mof.fits' #
    if fn_unmatched_cat[-1] != '/':
        fn_unmatched_cat = '{}/'.format(fn_unmatched_cat)
    __fn_end = fn_unmatched_cat[fn_unmatched_cat[:-1].rfind('/')+1:-1]

    __fn_cat = 'gauss_ap_ngmixer_{}'.format(__fn_end)

    __full_fn_ngmix_cat = os.path.join(matchedCatalogDirectory, __fn_cat)

    return __full_fn_ngmix_cat




def get_reformatted_coadd_catalog_filename(balrog_run, match_type, output_directory, realization, tile):
    """Create filename for the reformatted coadd catalog.
    g-, r-, i-, and z-band coadd catalogs exist. Catalogs are reformatted so that g-, r-, i-, and z-band relevant observables are in one catalog.
    Several observables are added to the i-band catalog so that matching is performed on the i-band RA and Dec.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    match_type (str)
        Describes the two catalogs matched (via `join=1and2`) in `stilts_matcher`.
        Example: '10%_inj_mof_cat_10%_inj_truth_cat'.
        Note that `match_type` reflects the order in which catalogs were matched. 
        In the example above, the injected MOF catalog was used as the STILTS parameter `in1`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    Returns
    -------
    __fn_coadd_for_matcher (str)
        Complete filename for the reformatted coadd catalog.
    """

    catalogFileDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type) 

    __fn_coadd_for_matcher = os.path.join(catalogFileDirectory, '{}_i_cat_reformatted.fits'.format(tile))

    return __fn_coadd_for_matcher 




def write_log_file_headers(fn_mag_err_log, fn_flag_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log, fn_flux_sigma_clip_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log):
    """Write headers to log files.

    Parameters
    ----------
    fn_flag_log (str)
        Complete filename for the log file that records flags if `LOG_FLAGS=True`. 
        If `LOG_FLAGS=False` this file is not created; 'dummy.log' is created.

    fn_mag_err_log (str)
        Complete filename for the log file that records details of the magnitude error calculation.

    fn_color_log (str)
        Complete filename for the log file that details the color calculation. 

    fn_mag_diff_outliers_log (str)
        Complete filename for the log file that records the magnitude differences greater than a threshold.
        Threshold is given by `__cutoff_mag_diff_for_log` in `bin_and_cut_measured_magnitude_error()`. 

    fn_mag_completeness_log (str)
        Complete filename for the log file that details the magnitude completeness calculation. 

    fn_flux_sigma_clip_log (str)
        Complete filename for the log file that records the details of the `N`-sigma clip of the normalized (to flux error) flux difference.
        If `PLOT_FLUX=False` and `SIGMA_CLIP_NORM_FLUX_DIFF=False` this file not created; 'dummy.log' is created. 

    fn_percent_recovered_log (str)
        Complete filename for the log file that records the percent of Balrog-injected objects recovered in the appropriate measured catalog.
        The percent recovered is computed in two ways: (1) without removing flagged* objects and (2) removing flagged* objects.
        For a description of flags, see https://github.com/spletts/BalVal/blob/master/README.md#flagged-objects

    fn_num_objs_in_1sig_mag_log (str)
        Complete filename for the log file that records the number of objects that lie within one-sigma of the measured magnitude error.
        If `PLOT_MAG=False`, this file is not created; 'dummy.log' is created.

    Returns
    -------
    Identical to Parameters.
    """

    ### Percent recovered log ###
    if MATCH_CAT1 in ('gal_truth', 'star_truth') or MATCH_CAT2 in ('gal_truth', 'star_truth'):
        with open(fn_percent_recovered_log, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['TILE', 'REALIZATION', 'BAND', 'TOTAL_OBJS_IN_MATCH1AND2', 'TOTAL_FLAGGED_OBJS_IN_MATCH1AND2', 'TOTAL_FLAGGED_OBJS_IN_TRUTH', '%_RECOVERED_FLAGS_IN', '%_RECOVERED_FLAGS_RM', 'RUN_TYPE'])
        csvfile.close()


    ### Number of objects in 1sigma_mag_meas ###
    if PLOT_MAG:
        with open(fn_num_objs_in_1sig_mag_log, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['TILE', 'REALIZATION', 'BAND', 'TOTAL_OBJS_IN_MATCH1AND2', 'TOTAL_FLAGGED_OBJS_IN_MATCH1AND2', 'TOTAL_OBJS_IN_1SIGMA_MAG_MEAS', 'NORMALIZED', 'RUN_TYPE'])
        csvfile.close()

    ### Flag log? ###


    ### Log for magnitude error calculation ###
    with open(fn_mag_err_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        if ('truth' in MATCH_CAT1 and SWAP_HAX is False) or ('truth' in MATCH_CAT2 and SWAP_HAX):
            # mag_true #
            writer.writerow(['TILE', 'REALIZATION', 'BAND', 'NUM_OBJS_IN_BIN', 'MAG_TRUE_BIN_LHS', 'MAG_TRUE_BIN_RHS', 'MEDIAN_HAX_MAG_TRUE', 'MEDIAN_MAG_ERROR'])
        else:
            # mag_meas #
            writer.writerow(['TILE', 'REALIZATION', 'BAND', 'NUM_OBJS_IN_BIN', 'MAG_MEAS_BIN_LHS', 'MAG_MEAS_BIN_RHS', 'MEDIAN_HAX_MAG_MEAS', 'MEDIAN_MAG_ERROR'])  
    csvfile.close()


    ### Log file for color calculation ###
    with open(fn_color_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        if 'truth' in MATCH_CAT1:
            writer.writerow(['TILE', 'REALIZATION', 'COLOR', 'MAG_TRUE_BIN', 'OBJS_IN_MAG_BIN', 'TOTAL_OBJECTS_NO_FLAGS', 'RUN_TYPE'])
        if 'truth' not in MATCH_CAT1:
            writer.writerow(['TILE', 'REALIZATION', 'COLOR', 'MAG_MEAS_BIN', 'OBJS_IN_MAG_BIN', 'TOTAL_OBJECTS_NO_FLAGS', 'RUN_TYPE'])
    csvfile.close()


    ### Log file for magnitude outliers ###
    with open(fn_mag_diff_outliers_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['TILE', 'REALIZATION', 'BAND', 'MAG1', 'MAG2', 'MAG_DIFF'])
    csvfile.close()


    ### Log file for magnitude completeness calculation ###
    with open(fn_mag_completeness_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['TILE', 'REALIZATION', 'BAND', 'INJ_PERCENT', 'TRUTH_MAG_BIN_LHS', 'TRUTH_MAG_BIN_RHS', 'MATCH_CAT_OBJS_IN_BIN', 'TRUTH_CAT_OBJS_IN_BIN'])
    csvfile.close()


    ### Log file for sigma clipping done in `()` ###
    with open(fn_flux_sigma_clip_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['TILE', 'REALIZATION', 'BAND', 'GAUSSIAN_APER_APPLIED', 'NUM_OBJS_FLAGS_RM', 'N-SIGMA_CLIP', 'NUM_OBJS_CLIPPED']) 
    csvfile.close()


    return fn_flag_log, fn_mag_err_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log, fn_flux_sigma_clip_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log
