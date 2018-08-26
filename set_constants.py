"""
Set the catalogs to be matched.
Set plot attributes.
Set the color coding of plot(s).
See 'Table of Constants' in README.md (https://github.com/spletts/BalVal/blob/master/README.md)
Catch error from incompatible parameters set at command line.

Comments are ABOVE the code they refer to.
"""

import numpy as np

# From BalVal #
import calculate_injection_percent


##### For catalogs #####
MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'coadd'
INJ1, INJ2 = True, False

STACK_REALIZATIONS = False
STACK_TILES = False

# FOF analysis #
RUN_TYPE = None


##### For plots #####

### Observable ###
PLOT_COLOR = False
PLOT_FLUX = False
PLOT_MAG = True

### Display ###
# Color plots #
COLOR_YLOW, COLOR_YHIGH = None, None
# For both color plots and magnitude plots #
CORNER_HIST_2D = False
SCATTER = False
# Magnitude plots #
HIST_2D = False
HEXBIN = True
PLOT_COMPLETENESS = False
NORMALIZE_MAG = False
PLOT_68P = True
PLOT_34P_SPLIT = True
PLOT_MAG_ERR = True
CENTER_MAG_ERR_ABT_ZERO = False
CM_T_ERR_CBAR = False
CM_T_CBAR = False
MAG_YLOW, MAG_YHIGH = None, None

# Flux plots #
PLOT_GAUSS_APER_FLUX = False
PLOT_CM_FLUX = True
PLOT_GAUSSIAN_FIT = False
PLOT_PEAKS = True
TRIM_NORM_FLUX_DIFF = False
NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY = True
RAW_NORM_FLUX_DIFF = False
SIGMA_CLIP_NORM_FLUX_DIFF = True
N = 3.0
NUM_ITERS = None
FLUX_XLOW, FLUX_XHIGH = -8, 8


### For plots, in general ###
OVERWRITE_AXLABELS = False

SAVE_PLOT = False
SHOW_PLOT = True

SWAP_HAX = False
SWAP_ORDER_OF_SUBTRACTION = False

SUBPLOT = True

##### Region files #####
MAKE_REGION_FILES = True


### Printouts ###
NOTICE = False
VERBOSE_ING = True
VERBOSE_ED = False

NO_DIR_MAKE = True


### Defaults overwritten only if necessary ###

# `Y3_FIT` only used if `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog #
Y3_MODEL, Y3_FIT = 'CM', None
if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
    if 'mof' in (MATCH_CAT1, MATCH_CAT2) or 'deep_sn_mof' in (MATCH_CAT1, MATCH_CAT2):
        Y3_FIT = 'MOF'
    if 'sof' in (MATCH_CAT1, MATCH_CAT2) or 'deep_sn_sof' in (MATCH_CAT1, MATCH_CAT2) :
        Y3_FIT = 'SOF'

# `FOF_FIT` only used if `RUN_TYPE` is not `None` #
if RUN_TYPE is None:
    FOF_FIT = None 

if RUN_TYPE is not None:
    if VERBOSE_ING: print 'Doing FOF analysis. Please set `FOF_FIT`...'
    # Overwrite MATCH_CATs #
    #NOTE set me 
    FOF_FIT = 'mof'
    MATCH_CAT1, MATCH_CAT2 = FOF_FIT, FOF_FIT

# Note: not really being used #
SAVE_BAD_IDX = False




### Percentiles for `corner.hist2d` contours ###
# `levels` parameter in `corner.hist2d` interpreted as percentiles #
CORNER_HIST_2D_LVLS = 1.0 - np.exp(-0.5 * np.array([1.5]) ** 2)
# LVLS = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5, 2.2]) ** 2)
CORNER_HIST_2D_LVLS_CLRS = ['red']
# Point of confusion: `colors` parameter passed to contour_kwargs is applied in reverse order, so reverse order for labels #
CORNER_HIST_2D_LVLS_CLRS_LABEL = CORNER_HIST_2D_LVLS_CLRS[::-1]


### Color coding ###
CMAPS = {'g':'Greens', 'r':'Purples', 'i':'Greys', 'z':'Blues'}
PT_COLORS = {'g':'green', 'r':'purple', 'i':'darkgrey', 'z':'navy'}
# Flux plots #
GAP_FLUX_HIST = {'g':'blue', 'r':'blue', 'i':'blue', 'z':'blue'}
FIT_TO_GAP_FLUX = {'g':'cyan', 'r':'cyan', 'i':'cyan', 'z':'cyan'}
FLUX_HIST = {'g':'green', 'r':'green', 'i':'green', 'z':'green'}
FIT_TO_FLUX = {'g':'lime', 'r':'lime', 'i':'lime', 'z':'lime'}


### For completeness plot ###
# Keep bin boundaries uniform across multiple tiles, important when `STACK_TILES=True` #
COMPLETENESS_PLOT_MAG_BINS = np.arange(15, 30, 0.5)


### User should not need to change the following ###
# Bands #
ALL_BANDS = ['g', 'r', 'i', 'z']
# For computing color #
BANDS_FOR_COLOR = {'g':'r', 'r':'i', 'i':'z'}
# For writing to log files #
WRITE_COLORS = {'g':'g-r', 'r':'r-i', 'i':'i-z'}
BAND_INDEX = {'g':0, 'r':1, 'i':2, 'z':3}




def catch_error(cmd_line_realizations, cmd_line_tiles):
    """Find errors created by setting incompatible parameters and command line arguments.

    Parameters
    ----------
    cmd_line_realizations (list of str)
        Realizations provided at the command line via `$python plotter.py ....`.
    cmd_line_tiles (list of str)
        Tiles provided at the command line via `$python plotter.py ....`

    Returns
    -------
    __err_msg (str or None)
        Error message. Is `None` if there is no error.
    """

    __err_msg = None

    if cmd_line_realizations[0] == 'None' and INJ1 and INJ2: __err_msg = 'Realization=None entered at command line must be used with `INJ1=True` or `INJ2=True`.'
    if INJ1 is False and INJ2 is False and cmd_line_realizations[0] != 'None': __err_msg = 'If `INJ1=False` and `INJ2=False` must have realizations=None at command line.'

    if STACK_REALIZATIONS and len(cmd_line_realizations) == 1: __err_msg = 'STACK_REALIZATIONS must be used with multiple realizations entered at command line'

    if STACK_TILES and ('.dat' not in cmd_line_tiles[0] and len(cmd_line_tiles) == 1): __err_msg = 'STACK_TILES must be used with multiple tiles entered at command line or .dat file provided at command line.'

    if MAG_YLOW is not None and MAG_YHIGH is not None and MAG_YHIGH == MAG_YLOW: __err_msg = '`MAG_YLOW` =/= `MAG_YHIGH`. These two cannot be equal.'

    if NORMALIZE_MAG and PLOT_MAG_ERR is False: __err_msg = 'If NORMALIZE_MAG is True so must be PLOT_MAG_ERR.'

    if NORMALIZE_MAG and PLOT_COLOR: __err_msg = 'Script not equipped to normalize color plots'

    if PLOT_FLUX is False and PLOT_COLOR is False and PLOT_MAG is False: __err_msg = 'Must plot one observable'

    if PLOT_FLUX and PLOT_GAUSS_APER_FLUX:
        if 'star_truth' in (MATCH_CAT1, MATCH_CAT2): 
            __err_msg = 'Star truth catalogs do not have the necessary  headers to perform Gaussian aperture flux measurement'
        if  'coadd' in (MATCH_CAT1, MATCH_CAT2):
            __err_msg = 'Coadd catalogs do not have the necessary headers to perform Gaussian aperture flux measurement'

    # Colorbar errors #
    cbar_counter = 0
    if HEXBIN: cbar_counter += 1
    if CM_T_ERR_CBAR: cbar_counter += 1
    if HIST_2D: cbar_counter += 1
    if CORNER_HIST_2D: cbar_counter += 1
    if CM_T_CBAR: cbar_counter += 1 
    if cbar_counter > 1: __err_msg = 'Only one colorbar can be used. Edit `HEXBIN`, `CM_T_ERR_CBAR`, `HIST_2D`, `CORNER_HIST_2D`, `CM_T_CBAR`, and `CM_T_CBAR`.' 

    if ('truth' in MATCH_CAT1 and INJ1 is False) or ('truth' in MATCH_CAT2 and INJ2 is False): __err_msg = 'Truth catalogs are always injected. Check `INJ1` and `INJ2`.'

    if PLOT_COMPLETENESS and 'truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2: __err_msg = 'Completeness plots must be made with a truth catalog. Check `MATCH_CAT1` and `MATCH_CAT2`.'

    if PLOT_COMPLETENESS and PLOT_MAG is False: __err_msg = 'Script not equipped to make flux completeness plots nor color completeness plots.'

    # Display errors #
    display_mag_counter = 0
    if HIST_2D: display_mag_counter += 1
    if CORNER_HIST_2D: display_mag_counter += 1 
    if HEXBIN: display_mag_counter += 1
    if SCATTER: display_mag_counter += 1
    if PLOT_MAG and display_mag_counter == 0 and PLOT_COMPLETENESS is False: __err_msg = 'Pick a display for magnitude plot via `HIST_2D`, `CORNER_HIST_2D`, `HEXBIN`, or `SCATTER`.'
    if display_mag_counter > 1: __err_msg = 'Pick only one display for magnitude plot.'

    if PLOT_COLOR and SCATTER is False and CORNER_HIST_2D is False: __err_msg = 'Pick a display for color plot via `SCATTER` or `CORNER_HIST_2D`.'
    if PLOT_COLOR and SCATTER and CORNER_HIST_2D:  __err_msg = 'Pick only one display for color plot via `SCATTER` or `CORNER_HIST_2D`'

    # Just a warning #
    if PLOT_FLUX and display_mag_counter != 0: print 'Flux plots are histograms. Ignoring `HIST_2D`, `CORNER_HIST_2D`, `HEXBIN`, and `SCATTER`...' 

    if 'y3_gold' in MATCH_CAT1 and INJ1: __err_msg = 'Y3 Gold catalogs are not Balrog-injected. Edit `INJ1`.'
    if 'y3_gold' in MATCH_CAT2 and INJ2: __err_msg = 'Y3 Gold catalogs are not Balrog-injected. Edit `INJ2`.'

    return __err_msg
