"""
Set constants that describe the catalogs.
Set constants that determine the plot.
Set color coding of plot.
See 'Table of Constants' in README.md (https://github.com/spletts/BalVal/blob/master/README.md)
"""

import numpy as np

# From BalVal #
import calculate_injection_percent

ALL_BANDS = ['g', 'r', 'i', 'z']

### For catalogs ###
MATCH_CAT1, MATCH_CAT2 = 'deep_sn_sof', 'y3_gold_2_2'
INJ1, INJ2 = False, False 


### For plots ###
# Observable #
PLOT_COLOR = False
PLOT_FLUX = False
PLOT_MAG = True

SAVE_PLOT = False
SHOW_PLOT = True

# Display magnitude #
PLOT_COMPLETENESS = False
HIST_2D = False
HEXBIN = False
# Can be used for color plot and magnitude plot #
CORNER_HIST_2D = True 
SCATTER = False

# For flux plots #
PLOT_GAUSS_APER_FLUX = False
PLOT_CM_FLUX = True
TRIM_NORM_FLUX_DIFF = False
NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY = True
RAW_NORM_FLUX_DIFF = False
SIGMA_CLIP_NORM_FLUX_DIFF = True
N = 3.0

# For magnitude plots #
NORMALIZE = False
PLOT_68P = True
PLOT_34P_SPLIT = True
PLOT_MAG_ERR = True 
CENTER_ERR_ABT_ZERO = False
CM_T_ERR_CBAR = False
CM_T_CBAR = False

MAG_YLOW, MAG_YHIGH = -1, 1 
COLOR_YLOW, COLOR_YHIGH = None, None
FLUX_XLOW, FLUX_XHIGH = -7, 7

STACK_REALIZATIONS = False
STACK_TILES = False

OVERRIDE_AXLABELS = True

VERBOSE_ING = True
VERBOSE_ED = True

# FOF analysis #
RUN_TYPE = None

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
        if MOF:
                MATCH_CAT1, MATCH_CAT2 = 'mof', 'mof'
        if MOF is False:
                MATCH_CAT1, MATCH_CAT2 = 'sof', 'sof'


# Only refers to printouts within get_floats_from_string(), get_matrix_diagonal_element(), and bin_and_cut_measured_magnitude_error() # 
VERBOSE_MINOR = False

# Not currently in use or under constructrion #
LOG_FLAGS = False
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


# Used if LOG_FLAGS is True #
FLAG_IDX = []
