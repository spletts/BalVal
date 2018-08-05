"""
Set constants that determine the catalogs.
Set constants that describe plot attributes.
Set the color coding of plot.
See 'Table of Constants' in README.md (https://github.com/spletts/BalVal/blob/master/README.md)
Catches error.

Comments are ABOVE the code they refer to.
"""

import numpy as np
import sys

# From BalVal #
import calculate_injection_percent

ALL_BANDS = ['g', 'r', 'i', 'z']

### For catalogs ###
MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'mof'
INJ1, INJ2 = True, True 
NOTICE = True

### For plots ###
# Observable #
PLOT_COLOR = False
PLOT_FLUX = True
PLOT_MAG = False

SAVE_PLOT = False
SHOW_PLOT = True

# Display #
HIST_2D = False
HEXBIN = False
# Can be used for color plot and magnitude plot #
CORNER_HIST_2D = False 
SCATTER = False
COLOR_YLOW, COLOR_YHIGH = None, None

# For flux plots #
PLOT_GAUSS_APER_FLUX = True
PLOT_CM_FLUX = True
PLOT_GAUSSIAN_FIT = True
PLOT_PEAKS = False
TRIM_NORM_FLUX_DIFF = False
NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY = True
RAW_NORM_FLUX_DIFF = False
SIGMA_CLIP_NORM_FLUX_DIFF = True
N = 3.0
NUM_ITERS = None #20
FLUX_XLOW, FLUX_XHIGH = -8, 8

# For magnitude plots #
PLOT_COMPLETENESS = False
NORMALIZE = False
PLOT_68P = True
PLOT_34P_SPLIT = True
PLOT_MAG_ERR = True 
CENTER_ERR_ABT_ZERO = False
CM_T_ERR_CBAR = False
CM_T_CBAR = False
MAG_YLOW, MAG_YHIGH = -1, 1


STACK_REALIZATIONS = False
STACK_TILES = False

OVERWRITE_AXLABELS = False

VERBOSE_ING = True
VERBOSE_ED = True

SWAP_HAX = False
SWAP_ORDER_OF_SUBTRACTION = False


# Defaults overwritten if necessary #
Y3_MODEL, Y3_FIT = 'CM', None


### Handle nonsensical combinations ###

# Only used if MATCH_CAT1 or MATCH_CAT2 is 'y3_gold'. Use MOF, SOF observables? # 
if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
	if 'mof' in (MATCH_CAT1, MATCH_CAT2) or 'deep_sn_mof' in (MATCH_CAT1, MATCH_CAT2):
		Y3_FIT = 'MOF'
	if 'sof' in (MATCH_CAT1, MATCH_CAT2) or 'deep_sn_sof' in (MATCH_CAT1, MATCH_CAT2) :
		Y3_FIT='SOF'


SUBPLOT = True
NO_DIR_MAKE = True
MAKE_REG = False


### For FOF analysis. To ignore FOF analysis set `RUN_TYPE=None` ###
RUN_TYPE = None
if RUN_TYPE is None:
	MOF = None 

if RUN_TYPE is not None:
	if VERBOSE_ING: print 'Doing FOF analysis ... \n '
        # Overwrite MATCH_CATs #
	#NOTE set me 
        FOF_FIT = 'mof'
	MATCH_CAT1, MATCH_CAT2 = FOF_FIT, FOF_FIT


# Only refers to printouts within get_floats_from_string(), get_matrix_diagonal_element(), and bin_and_cut_measured_magnitude_error() # 
VERBOSE_MINOR = False

# Not currently in use or under constructrion #
PLOT_FLAGGED_OBJS = True




### Percentiles for `corner.hist2d` contours ###
# Correct 1sigma for a 2D histogram is given by http://corner.readthedocs.io/en/latest/pages/sigmas.html #
# `LVLS` passed to `corner.hist2d` interpreted as percentiles #
LVLS = 1.0 - np.exp(-0.5 * np.array([1.5]) ** 2)
# LVLS = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5, 2.2]) ** 2)
CLRS = ['red']
# CLRS = ['magenta', 'red', 'blue', 'yellow']
# Question: `colors` passed to contour_kwargs are passed/applied in reversed order so reverse order for labels #
CLRS_LABEL = CLRS[::-1]


### Dictionaries for color coding ###
CMAPS = {'g':'Greens', 'r':'Purples', 'i':'Greys', 'z':'Blues'}
PT_COLORS = {'g':'green', 'r':'purple', 'i':'darkgrey', 'z':'navy'}
GAUSS_FIT_COLORS = {'g':'lime', 'r':'magenta', 'i':'dimgrey', 'z':'cyan'}
# 'Opposite' of `PT_COLORS`
GAUSS_FIT_TO_GAP_COLORS = {'g':'crimson', 'r':'red', 'i':'tomato', 'z':'orangered'}
GAUSSIAN_FIT_TO_GAUSSIAN_APERTURE_MEASUREMENT = {'g':'red', 'r':'yellow', 'i':'blue', 'z':'orange'}

FLUX_HIST = {'g':'green', 'r':'green', 'i':'green', 'z':'green'}
FIT_TO_FLUX = {'g':'lime', 'r':'lime', 'i':'lime', 'z':'lime'}
GAP_FLUX_HIST = {'g':'blue', 'r':'blue', 'i':'blue', 'z':'blue'}
FIT_TO_GAP_FLUX = {'g':'cyan', 'r':'cyan', 'i':'cyan', 'z':'cyan'}

# For computing color #
BANDS_FOR_COLOR = {'g':'r', 'r':'i', 'i':'z'}

# For writing to log files #
WRITE_COLORS = {'g':'g-r', 'r':'r-i', 'i':'i-z'}
BAND_INDEX = {'g':0, 'r':1, 'i':2, 'z':3}

### For completeness plot ###
# Keep bin boundaries uniform across multiple tiles, important when `STACK_TILES=True` #
COMPLETENESS_MAG_BINS = np.arange(15, 30, 1)
COMPLETENESS_PLOT_MAG_BINS = []
for i in np.arange(0, len(COMPLETENESS_MAG_BINS)-1):
        COMPLETENESS_PLOT_MAG_BINS.append(np.median([COMPLETENESS_MAG_BINS[i], COMPLETENESS_MAG_BINS[i+1]]))




def catch_error(cmd_line_realizations, cmd_line_tiles):
	"""Find errors created by setting parameters and command line arguments in an incompatible way."""

	__err_msg = None

	if cmd_line_realizations[0] == 'None' and INJ1 and INJ2: __err_msg = 'Realization of None at cmd line must be used with INJ1=True ot INJ2=True'

	if STACK_REALIZATIONS and len(cmd_line_realizations) == 1: __err_msg = 'STACK_REALIZATIONS must be used with multiple realizations'
	#if STACK_REALIZATIONS and realizations[0] != 'all': __err_msg = 'STACK_REALIZATIONS is True must be used with realization = all'
	if STACK_TILES and ('.dat' not in cmd_line_tiles[0] and len(cmd_line_tiles) == 1): __err_msg = 'STACK_TILES must be used with multiple tiles'
	if MAG_YLOW is not None and MAG_YHIGH is not None and MAG_YHIGH == MAG_YLOW: __err_msg = 'MAG_YLOW and MAG_YHIGH cannot be equal'
	if NORMALIZE and PLOT_MAG_ERR is False: __err_msg = 'If NORMALIZE is True so must be PLOT_MAG_ERR'

	if NORMALIZE and PLOT_COLOR: __err_msg = 'Script not equipped to normalize color plots'

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
	#TODO add NORMALIZE
	if cbar_counter > 1: __err_msg = 'Only one colorbar can be used. Edit HEXBIN, CM_T_ERR_CBAR, HIST_2D, CM_T_CBAR'

	if INJ1 is False and INJ2 is False and cmd_line_realizations[0] != 'None': __err_msg = 'If INJ1 and INJ2 are False realizations must be None at cmd line'

	if ('truth' in MATCH_CAT1 and INJ1 is False) or ('truth' in MATCH_CAT2 and INJ2 is False): __err_msg = 'Truth catalogs are injected'

	if PLOT_COMPLETENESS and 'truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2: __err_msg = 'Completeness plots must be made with a truth catalog'

	if PLOT_FLUX and PLOT_COMPLETENESS: __err_msg = 'Script not equipped to make flux completeness plots'

	plt_type_counter = 0
	if PLOT_COMPLETENESS: plt_type_counter += 1
	if PLOT_FLUX: plt_type_counter += 1
	if plt_type_counter > 1: __err_msg = 'Pick one plot type...' #TODO

	display_type_counter = 0
	if HIST_2D: display_type_counter += 1
	if CORNER_HIST_2D: display_type_counter += 1 
	if HEXBIN: display_type_counter += 1
	if SCATTER: display_type_counter += 1
	if PLOT_MAG and display_type_counter == 0 and PLOT_COMPLETENESS is False: __err_msg = 'Pick at least one display for magnitude plot'

	if display_type_counter > 1: __err_msg = 'Pick only one display for plot'

	if PLOT_COLOR and SCATTER is False and CORNER_HIST_2D is False: __err_msg = 'Color plots currently only made for CORNER_HIST_2D' #TODO

	if PLOT_FLUX and CORNER_HIST_2D: __err_msg = 'Flux plots cannot be corner plots'
	#TODO not normalize and plot_color simult.

	if 'y3_gold' in MATCH_CAT1 and INJ1: __err_msg = 'Y3 Gold catalogs are not Balrog-injected'
	if 'y3_gold' in MATCH_CAT2 and INJ2: __err_msg = 'Y3 Gold catalogs are not Balrog-injected'

	return __err_msg
