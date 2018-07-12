"""Creates various plots for Balrog validation testing.


To run: $python ms_plotter.py base_path_to_catalogs output_directory realization tile
Example: $python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /Users/mspletts/BalVal/ all DES2247-4414

Relies on ms_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` with the correct path to ms_matcher.
FOF analysis relies on ms_fof_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher` with the correct path to ms_fof_matcher.

Plot attributes are specified with constants (many of them booleans) at the beginning of this script.
See 'Table of Constants' in README.md for descriptions of these constants. 
Docstrings describe function parameters and function returns, but a GoogleSheet is also availabe: https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing

# Comments are ABOVE the code they correspond to (with the exception of FIXMEs and TODOs) #
In general, function input parameters are `named_with_underscores`. 
Variables defined within functions (and that should only be accessed by said function) are `__named_with_leading_underscores`.
Variable `namesWithoutSpaces` are defined by another function via `nameWithoutSpace = function(param1=x)` 

Note that `None` is typically used as a flag value throughout this script.

Megan Splettstoesser mspletts@fnal.gov"""

#TODO Update docstrings

# `astropy` is needed only if analyzing a coadd catalog or making a completeness plot #
from astropy.io import fits
from astropy.table import Table, vstack, Column
from scipy import stats
# $setup ngmix #
from ngmix import gmix
# For Python environment with `corner` run: $source activate des18a #
import corner
import csv
import fileinput
import math
import matplotlib
import matplotlib.pyplot as plt
import ngmix
import numpy as np
import os
import pandas as pd
import subprocess
import sys




### Command line args ###
# Catch error from inadequate number of command line args #
if len(sys.argv) != 5:
        sys.exit("Args: basepath (location of Balrog catalogs), output directory, realizations (can be a list of form: real1,real2,...), tiles (can be a filename.dat, or a list of form: tile1,tile2,...) \n")
# Convert `realizations` and `tiles` into lists (may be one-element list). Note 'None' entered at command line is interpreted as a str #
BASEPATH, OUTDIR, cmd_line_realizations, cmd_line_tiles = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4].split(',')


### For directory structure ###
# '/data/des71.a/data/kuropat/des2247-4414_sof/' --> 'des2247-4414_sof' #
BALROG_RUN = BASEPATH[BASEPATH[:-1].rfind('/')+1:-1]
# Rewrite #
if BALROG_RUN == 'Balrog': BALROG_RUN = 'TAMU_Balrog'
if BALROG_RUN == 'DES0102-4914': BALROG_RUN = 'blank_test_DES0102-4914'

### Constants needed to loop over bands, realizations, and tiles ###
# Bands #
ALL_BANDS = [ 'g', 'r', 'i', 'z' ]

# Realization(s) #
ALL_REALIZATIONS = cmd_line_realizations

# Tile(s) #
if '.dat' not in cmd_line_tiles[0]: ALL_TILES = cmd_line_tiles
# Accept .dat file of tile names at command line #
if '.dat' in cmd_line_tiles[0]:
	ALL_TILES = []
	for line in fileinput.input(cmd_line_tiles):
		# Get rid of newline character \n #
		ALL_TILES.append(line[:-1])





################################################################### Specify plot and catalog attributes ###################################################################

# Please see 'Table of Parameters' in README.md #

MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'mof'
INJ1, INJ2 = True, True
INJ1_PERCENT, INJ2_PERCENT = 20, 20

# Observable #
PLOT_MAG = False
PLOT_COLOR = True
PLOT_FLUX = False
GAUSS_APER = True
TRIM_FLUX = True

SAVE_PLOT = False
SHOW_PLOT = True

# Display observable #
PLOT_COMPLETENESS = False
HIST_2D = True
CORNER_HIST_2D = False
HEXBIN = False 
SCATTER = False
CM_T_ERR_COLORBAR = False
CM_T_COLORBAR = False

# For magnitude plots #
NORMALIZE = False
PLOT_68P = True
PLOT_34P_SPLIT = True
PLOT_MAG_ERR = True 
CENTER_ERR_ABT_ZERO = False

PLOT_DIFF_ON_VAX = True

MAG_YLOW, MAG_YHIGH = None, None 

STACK_REALIZATIONS = False
STACK_TILES = False

# FOF analysis #
RUN_TYPE = None

# Swap horizontal axis? Default is magnitude1. Matching script ms_matcher determines which catalog is 1 and which is 2. Generally SWAP_HAX does not need to be changed unless the truth catalog values are not on the horizontal axis. #
SWAP_HAX = False



### Handle nonsensical combinations ###
# Always non-injections #
if 'y3_gold' in MATCH_CAT1:
	#INJ1_10PERCENT, INJ1_20PERCENT = False, False
	INJ1_PERCENT = None
if 'y3_gold' in MATCH_CAT2:
	INJ2_PERCENT = None
if cmd_line_realizations[0] == 'None':
	INJ1, INJ2 = False, False

# !!!!! Only used if MATCH_CAT1 or MATCH_CAT2 is 'y3_gold'. If False, SOF catalogs exists in subdirectories of BASEPATH #
Y3_MOF = None
if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
	if 'mof' in (MATCH_CAT1, MATCH_CAT2):
		Y3_MOF = True
	if 'sof' in (MATCH_CAT1, MATCH_CAT2):
		Y3_MOF = False
 


# !!!!! Make 2x2 subplots of each griz band? Or make individual plots? #
SUBPLOT = True

# !!!!! If directories do no exist, make them (True) or force sys.exit() to edit dirs within script (False)? # 
NO_DIR_MAKE = True


### For FOF analysis. To ignore FOF analysis set `RUN_TYPE=None` ###
# !!!!! Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog-injection. 'rerun' refers to groups changed after Balrog-injection. #
RUN_TYPE = None

# !!!!! Only used if RUN_TYPE is not None #
MOF = False

if RUN_TYPE is not None:
	print 'Doing FOF analysis ... \n '
        # Overwrite MATCH_CATs #
        if MOF:
                MATCH_CAT1, MATCH_CAT2 = 'mof', 'mof'
        if MOF is False:
                MATCH_CAT1, MATCH_CAT2 = 'sof', 'sof'


### !!!!! Make region files? #
MAKE_REG = False 

### Miscellaneous ###
# Print progress? #
PRINTOUTS = True
# Only refers to printouts within get_floats_from_string(), get_matrix_diagonal_element(), and bin_and_cut_measured_magnitude_error() # 
PRINTOUTS_MINOR = False

# Not currently in use or under constructrion #
LOG_FLAGS = False
PLOT_FLAGGED_OBJS = True
# Use quality cuts introduced by Eric Huff? Link: https://github.com/sweverett/Balrog-GalSim/blob/master/plots/balrog_catalog_tests.py. Can only be performed if catalog has all the necessary headers: cm_s2n_r, cm_T, cm_T_err, and psfrec_T. #




def catch_error():
	"""Find errors created by setting parameters and command line arguments in an incompatible way."""

	__err_msg = None

	if STACK_REALIZATIONS and len(cmd_line_realizations) == 1: __err_msg = 'STACK_REALIZATIONS must be used with multiple realizations'
	#if STACK_REALIZATIONS and realizations[0] != 'all': __err_msg = 'STACK_REALIZATIONS is True must be used with realization = all'
	if STACK_TILES and ('.dat' not in cmd_line_tiles[0] and len(cmd_line_tiles) == 1): __err_msg = 'STACK_TILES must be used with multiple tiles'
	if MAG_YLOW is not None and MAG_YHIGH is not None and MAG_YHIGH == MAG_YLOW: __err_msg = 'MAG_YLOW and MAG_YHIGH cannot be equal'
	if NORMALIZE and PLOT_MAG_ERR is False: __err_msg = 'If NORMALIZE is True so must be PLOT_MAG_ERR'

	if NORMALIZE and PLOT_COLOR: __err_msg = 'Script not equipped to normalize color plots'

	if PLOT_FLUX is False and PLOT_COLOR is False and PLOT_MAG is False: __err_msg = 'Must plot one observable'

	# Colorbar errors #
	cbar_counter = 0
	if HEXBIN: cbar_counter += 1
	if CM_T_ERR_COLORBAR: cbar_counter += 1
	if HIST_2D: cbar_counter += 1
	if CORNER_HIST_2D: cbar_counter += 1
	if CM_T_COLORBAR: cbar_counter += 1 
	#TODO add NORMALIZE
	if cbar_counter > 1: __err_msg = 'Only one colorbar can be used. Edit HEXBIN, CM_T_ERR_COLORBAR, HIST_2D, CM_T_COLORBAR'

	if INJ1 is False and INJ2 is False and cmd_line_realizations[0] != 'None': __err_msg = 'If INJ1 and INJ2 are False realizations must be None at cmd line'

	if ('truth' in MATCH_CAT1 and INJ1 is False) and ('truth' in MATCH_CAT2 and INJ2 is False): __err_msg = 'Truth catalogs are injected'

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
	#TODO not normalize and plot_color simult.
	return __err_msg

if catch_error() is not None:
        sys.exit('Error: ' + catch_error())


lines = [
['\ncatalog1 (name, injected?, injection percent) :\t ' + str(MATCH_CAT1) + ' ' + str(INJ1) + ' ' + str(INJ1_PERCENT)],
['catalog2 (name, injected?, injection percent) :\t ' + str(MATCH_CAT2) + ' ' + str(INJ2) + ' ' + str(INJ2_PERCENT)],
['Observable (magnitude, color, flux) :\t ' + str(PLOT_MAG) + ', ' + str(PLOT_COLOR) + ', ' + str(PLOT_FLUX)],
['Plot type (hist2d, corner.hist2d, completeness, scatter, hexbin, ...) :', str(HIST_2D) + ', ' + str(CORNER_HIST_2D) + ' , ' + str(PLOT_COMPLETENESS) + ', '+ str(SCATTER) + ', ' + str(HEXBIN)],
['Stack matched catalog? (stack realizations, stack tiles) :\t ' + str(STACK_REALIZATIONS) + ', ' + str(STACK_TILES)],
['Show plot? :\t ' + str(SHOW_PLOT)],
['Save plot? :\t ' + str(SAVE_PLOT)],
['Limits for the vertical axis of plot? :\t ' + str(MAG_YLOW) + ', ' + str(MAG_YHIGH)],
['Plot 1sigma_meas curve? :\t ' + str(PLOT_MAG_ERR)],
['If plotting 1sigma_meas curve, center about zero? (else centered about medians) :\t ' + str(CENTER_ERR_ABT_ZERO)],
['Normalize plot to error? :\t ' + str(NORMALIZE)],
['FOF analysis? :\t ' + str(RUN_TYPE)]
]

for line in lines:
	print('{:>12}'.format(*line))

NOTICE = raw_input('\n Check the above before running! \n --> Press enter to proceed, control+c to stop...\n')
#Make region files?
#If NORMALIZE: plot 68th percentile? plot 34th percentile split?




################################################################### Store catalog information ###################################################################

class CoaddCat():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{band}_cat.fits. There is a separate catalog for each band."""

        # Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj (bool)

		inj_percent (int)
			If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog. 

		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# For plot title #
		if inj:
			if inj_percent == 10:
				self.title_piece = '10% Inj Coadd Cat'
			if inj_percent == 20:
				self.title_piece = '20% Inj Coadd Cat'
		if inj is False:
			self.title_piece = 'Base Coadd Cat'

		self.axlabel = 'meas'
		# Magnitude, is one number #
		self.mag_hdr = 'MAG_AUTO' + str(suf)
		self.mag_axlabel = 'MAG_AUTO_meas'
		# For error calculation #
		self.mag_err_hdr = 'MAGERR_AUTO' + str(suf)
		self.cm_flux_hdr = None
		self.cm_flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = 'FLAGS' + str(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'ALPHAWIN_J2000' + str(suf)
		self.dec_hdr = 'DELTAWIN_J2000' + str(suf)
		self.angle = 'THETA_J2000' + str(suf)
		# Units: pixels #
		self.a_hdr = 'A_IMAGE' + str(suf)
		self.b_hdr = 'B_IMAGE' + str(suf)
		# For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = None




class GalTruthCat():
	"""Declare headers for galaxy truth catalogs. Note that as of May 2018 galaxy truth catalogs are created using MOF (thus have the same headers).""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj (bool)

		inj_percent (int)

		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		 # For plot title #
		# `inj` forced True for truth catalogs #
		if inj:
			if inj_percent == 10:
				self.title_piece = '10% Inj Gal Truth Cat'
			if inj_percent == 20:
				self.title_piece = '20% Inj Gal Truth Cat'
		self.axlabel = 'true'
		# Headers are the same as MOFCat class. Reproduced below for clarity in ms_plotter.py #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
	    	self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_true'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size. cm_T units: arcseconds squared. #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_true'
                # Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
		self.angle = None
                # For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf)
                self.psfrec_t_hdr = 'psfrec_T' + str(suf)



	
class SOFCat():
        """Declare headers for SOF catalogs.""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
		inj (bool)

                inj_percent (int)

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

		if inj:
			if inj_percent == 10:
				self.title_piece = '10% Inj SOF Cat'
			if inj_percent == 20:
				self.title_piece = '20% Inj SOF Cat'
		if inj is False:
                        self.title_piece = 'Base SOF Cat'

		self.axlabel = 'meas'
                # Headers are the same as MOFCat class with the exception of cm_mof_flags_hdr. Reproduced below for clarity #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
                self.cm_mof_flags_hdr = None
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
		self.angle = None
                # For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf)
                self.psfrec_t_hdr = 'psfrec_T' + str(suf)




class MOFCat():
	"""Declare headers for MOF catalogs."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_percent (int)

		inj (bool)

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

		# For plot title #
		if inj:
			if inj_percent == 10:
				self.title_piece = '10% Inj MOF Cat'
			if inj_percent == 20:
				self.title_piece = '20% Inj MOF Cat'
		if inj is False:
			self.title_piece = 'Base MOF Cat'
		self.axlabel = 'meas'
		# Magnitude, is string of form (mag_g, mag_r, mag_i, mag_z)  #
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None 
		self.cm_flux_hdr = 'cm_flux' + str(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
		# Size #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
		# For region file #
		self.ra_hdr = 'ra' + str(suf)
		self.dec_hdr = 'dec' + str(suf)
		self.a_hdr, self.b_hdr = None, None
		self.angle = None
		# For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf) 
                self.psfrec_t_hdr = 'psfrec_T' + str(suf) 




class StarTruthCat(): #are there sep headers for MOFStarTruthCat and SOFStarTruthCat?
        """Declare headers for star truth catalogs."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
                """Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_percent (int)

		inj (bool)

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """	

		# `inj` forced True for truth catalogs #	
		if inj:
			if inj_percent == 10:
				self.title_piece = '10% Inj Star Truth Cat'
			if inj_percent == 20:
				self.title_piece = '20% Inj Star Truth Cat'
		self.axlabel = 'true'
		# Magnitude #
		self.mag_hdr = 'g_Corr' + str(suf) 
		self.mag_axlabel = 'mag_true' 
		# For error calculation #
                # Is of form 'PSF_MAG_ERR_{band}' + str(suf) #
		self.mag_err_hdr = 'PSF_MAG_ERR' + str(suf)
                self.cm_flux_hdr = None
                self.cm_flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = None
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'RA_new' + str(suf)
		self.dec_hdr = 'DEC_new' + str(suf)
		self.a_hdr, self.b_hdr = None, None # Use cm_T
		self.angle = None
		# For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = None
		self.psfrec_t_hdr = None




class Y3Gold():
	"""Declare headers for Y3 Gold catalogs (https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Full_bins_of_Y3_GOLD_2_0_Columns).""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
                """Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_percent (int)

                inj (bool)

                suf (str)
                        Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

		if inj is False: sys.exit('error in class Y3Gold')

		if BALROG_RUN in ('grid_bkg', 'grid_bkg_noise', 'seed_test_base', 'seed_test_same_seeds' 'seed_test_same_seeds1', 'seed_test_new_seeds', 'des2247-4414_sof', 'grid_test_noise_sof') or Y3_MOF is False: 
			pref1 = 'SOF_'

		if Y3_MOF is True:
			pref1 = 'MOF_'

		if 'star' in MATCH_CAT1 or 'star' in MATCH_CAT2:
			pref2 = 'PSF_'
		if 'star' not in MATCH_CAT1 and 'star' not in MATCH_CAT2:
			pref2 = 'CM_'

		print 'Using ', pref1, '&', pref2, ' for Y3 Gold catalog ... \n'
	 
		self.title_piece = 'Y3 Gold Cat' 
                self.axlabel = 'meas'

		# Note: band dependent headers -- MAG, MAG_ERR, FLUX, FLUX_COV, PSF_FLAGS #

                # Magnitude # 
                self.mag_hdr = pref1 + pref2 + 'MAG' + str(suf) 
                self.mag_axlabel = 'mag_meas'
                # For error calculation #
                self.mag_err_hdr = pref1 + pref2 + 'MAG_ERR' + str(suf) 
                self.cm_flux_hdr = pref1 + pref2 + 'FLUX' + str(suf) 
                self.cm_flux_cov_hdr = pref1 + pref2 + 'FLUX_COV' + str(suf) 
                # Size #
                self.cm_t_hdr = pref1 + pref2 + 'T' + str(suf) 
                self.cm_t_err_hdr = pref1 + pref2 + 'T_ERR' + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
                # Flags #
                self.flags_hdr = 'FLAGS_GOLD' + str(suf)
                self.obj_flags_hdr = pref1 + 'OBJ_FLAGS' + str(suf)
                self.psf_flags_hdr = pref1 + 'PSF_FLAGS_filt' + str(suf) 
                self.cm_flags_hdr = pref1 + pref2 + 'FLAGS' + str(suf) 
                self.cm_max_flags_hdr = None
		# MOF_CM_MOF_FLAGS has no SOF equivalent #
		# Note: is duplicate #
                self.cm_mof_flags_hdr = pref1 + pref2 + 'FLAGS' + str(suf)
                self.cm_flags_r_hdr = pref1 + pref2 + 'FLAGS_R' + str(suf)
                # For region file #
		# Note: there is also an ALPHAWIN_J2000 and DELTAWIN_J2000 #
		self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr = 'A_IMAGE'+ str(suf)
		self.b_hdr = 'B_IMAGE' + str(suf) 
                self.angle = 'THETA_J2000' + str(suf) 
                # For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = pref1 + 'PSFREC_T' + str(suf) 




def get_class(cat_type, inj, inj_percent, suf):
        """Get the appropriate class for the catalog type.

        Parameters
	----------
	cat_type (str)
		Catalog type. This is set by `MATCH_CAT1` or `MATCH_CAT2`. 

	inj_percent (int)

	inj (bool)

	suf (str)
		Refers to the order in which `cat_type` was matched (via join=1and2) in ms_matcher (order set by `in1` and `in2` in STILTS script). Allowed values: '_1' '_2'.

        Returns
	-------
	__class (class)
		Points to the appropriate class and class constants.
        """

        if cat_type == 'gal_truth':
                __cat_class = GalTruthCat(inj_percent=inj_percent, inj=inj, suf=suf)

        if cat_type == 'mof':
                __cat_class = MOFCat(inj_percent=inj_percent, inj=inj, suf=suf)

        if cat_type == 'star_truth':
                __cat_class = StarTruthCat(inj_percent=inj_percent, inj=inj, suf=suf)

        if cat_type == 'sof':
                __cat_class = SOFCat(inj_percent=inj_percent, inj=inj, suf=suf)

        if cat_type == 'coadd':
                __cat_class = CoaddCat(inj_percent=inj_percent, inj=inj, suf=suf)

	if 'y3_gold' in cat_type:
		__cat_class = Y3Gold(inj_percent=inj_percent, inj=inj, suf=suf)


        return __cat_class 




def get_match_type(title_piece1, title_piece2):
        """Transform two strings of form '10% Inj MOF Cat' and '10% Inj Truth Cat' to '10%_inj_mof_cat_10%_inj_truth_cat'.

        Parameters
	----------
	title_piece1, title_piece2 (str)
		Set by `self.title_piece` in appropriate class. Class set by `MATCH_CAT*` `INJ*PERCENT`. Example: 10% Inj MOF Cat.

        Returns
	-------
	match_type (str)
		Reflects the order in which catalogs were matched (via join=1and2). Example: 10%_inj_mof_cat_10%_inj_truth_cat, where the injected MOF catalog was STILTS parameter `in1`.
        """

	__match_type = '_'.join([title_piece1.lower().replace(' ', '_'), title_piece2.lower().replace(' ', '_')])
	return __match_type




def get_log_filenames(tile, realization):
        """Generate names for log files. Relies on directory structure: /`OUTDIR`/log_files/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/log_files/.

        Parameters
	----------
	tile (str)

	realization (str)

        Returns
	-------
	__fn_flag_log (str)
		Complete filename for flag log file.

	__fn_mag_err_log (str)
		Complete filename for error calculation log file.

	__fn_main_log (str)
		Complete filename for log file containing number of objects matched, number of objects flagged, number of objects within 1sigma_mag, etc.

	__fn_color_log (str)

	__fn_mag_diff_outliers_log (str)

	__fn_mag_completeness_log (str)

	__fn_gauss_aper_log (str)

        """

        # !!!!! User may wish to edit directory structure #
	log_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, 'log_files') 
        if RUN_TYPE is not None:
		log_dir = os.path.join(log_dir, 'fof_analysis')

	### Check for directory existence ###
        if os.path.isdir(log_dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory ' + str(log_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_log_filenames() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', log_dir, '...\n'
                        os.makedirs(log_dir)

	# Repeated #
	fn_rep = '_'.join([tile, realization, MATCH_TYPE])
	__fn_flag_log = os.path.join(log_dir, 'flag_log_'+fn_rep+'.log')
	__fn_mag_err_log = os.path.join(log_dir, 'mag_err_calc_'+fn_rep+'.log')
	__fn_main_log = os.path.join(log_dir, 'main_'+fn_rep+'.log')
	__fn_color_log = os.path.join(log_dir, 'color_plot_'+fn_rep+'.log')
	__fn_mag_diff_outliers_log = os.path.join(log_dir, 'outlier_magnitude_differences_'+fn_rep+'.log')

	# By construction, completeness plots use both 10% and 20% injected catalogs. Remove the % injection given by `MATCH_TYPE`. #
	fn_complete_rep = fn_rep.replace('10%_', '')
	fn_complete_rep = fn_complete_rep.replace('20%_', '')
	if PLOT_COLOR:
		__fn_mag_completeness_log = os.path.join(log_dir, 'dummy_color_completeness_'+fn_complete_rep+'.log')
	if PLOT_MAG:
		__fn_mag_completeness_log = os.path.join(log_dir, 'mag_completeness_'+fn_complete_rep+'.log')
	if PLOT_FLUX:
		__fn_mag_completeness_log = os.path.join(log_dir, 'dummy_flux_completeness_'+fn_complete_rep+'.log')

        if RUN_TYPE is not None:
		__fn_flag_log = '_'.join([__fn_flag_log[:-4], RUN_TYPE, __fn_flag_log[-4:]])
		__fn_mag_err_log = '_'.join([__fn_mag_err_log[:-4], RUN_TYPE, __fn_mag_err_log[-4:]])
		__fn_main_log = '_'.join([__fn_main_log[:-4], RUN_TYPE, __fn_main_log[-4:]])
		__fn_color_log = '_'.join([__fn_color_log[:-4], RUN_TYPE, __fn_color_log[-4:]])
		__fn_mag_diff_outliers_log = '_'.join([__fn_mag_diff_outliers_log[:-4], RUN_TYPE, __fn_mag_diff_outliers_log[-4:]])
		__fn_mag_completeness_log = '_'.join([__fn_mag_completeness_log[:-4], RUN_TYPE, __fn_mag_completeness_log[-4:]])
		#TODO 

	# For Gaussian aperture #
	__fn_gauss_aper_log = os.path.join(log_dir, 'gauss_aper_'+fn_rep+'.csv')

        print '-----> Saving log file for flags as: ', __fn_flag_log, '\n'
        print '-----> Saving log file for magnitude and error bins as: ', __fn_mag_err_log, '\n'
        print '-----> Saving log file for number of flags, number of objects within 1sigma_mag as: ', __fn_main_log, '\n'
	print '-----> Saving log file for color plot as: ', __fn_color_log, '\n' #FIXME does this need to be created each time? Use PLOT_COLOR
	print '-----> Saving log file for outlier DeltaMagnitude as: ', __fn_mag_diff_outliers_log, '\n'
	print '-----> Saving log file for completeness plot as: ', __fn_mag_completeness_log, '\n'
	if GAUSS_APER: print '-----> Saving Gaussian aperture calculation details as: ', __fn_gauss_aper_log, '\n'
	
	return __fn_flag_log, __fn_mag_err_log, __fn_main_log, __fn_color_log, __fn_mag_diff_outliers_log, __fn_mag_completeness_log, __fn_gauss_aper_log 




def get_region_filenames(tile, realization):
	"""Generate names for region files of different join types (join types specified in ms_matcher or ms_fof_matcher). Relies on directory structure `/OUTDIR/outputs/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/region_files/`

        Parameters
	----------
	tile (str)

	realization (str)

        Returns
	-------
	__fn_reg_1and2 (str)
		Complete filename for region file containing regions with join=1and2.

	__fn_reg_1not2 (str)
		Complete filename for region file containing regions with join=1not2.

	__fn_reg_2not1 (str)
		Complete filename for region file containing regions with join=2not1.
        """

        # !!!!! User may wish to edit directory structure #
	reg_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, 'region_files')
        if RUN_TYPE is not None:
		reg_dir = os.path.join(reg_dir, 'fof_analysis')

	### Check for directory existence ###
        if os.path.isdir(reg_dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory ' + str(reg_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_region_filenames() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', reg_dir, '...\n'
                        os.makedirs(reg_dir)

	# Repeated #
        fn_rep = '_'.join([tile, realization, MATCH_TYPE])
	__fn_reg_1and2 = os.path.join(reg_dir, fn_rep+'_match1and2.reg')
	__fn_reg_1not2 = os.path.join(reg_dir, fn_rep+'_match1not2.reg')
	__fn_reg_2not1 = os.path.join(reg_dir, fn_rep+'_match2not1.reg')


	if RUN_TYPE is not None:
		__fn_reg_1and2 = __fn_reg_1and2[:-15] + '_' + str(RUN_TYPE) + __fn_reg_1and2[-15:]
		__fn_reg_1not2= __fn_reg_1not2[:-15] + '_' + str(RUN_TYPE) + __fn_reg_1not2[-15:]
		__fn_reg_2not1 = __fn_reg_2not1[:-15] + '_' + str(RUN_TYPE) + __fn_reg_2not1[-15:]


	return __fn_reg_1and2, __fn_reg_1not2, __fn_reg_2not1 







################################################################### Declare necessary constants ###################################################################

### For data analysis ###
# CLASS1 refers to in1 in ms_matcher. in1 appends _1 to all the headers, hence suf=1. ms_fof_matcher is done such that injected catalogs have no suffix #
if RUN_TYPE is not None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=INJ1_PERCENT, suf='')
if RUN_TYPE is None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=INJ1_PERCENT, suf='_1')
CLASS2 = get_class(cat_type=MATCH_CAT2, inj=INJ2, inj_percent=INJ2_PERCENT, suf='_2')


# Get arguments to pass to ms_matcher. Need to transform header of form 'ra_1' to 'ra', hence [:-2] #
RA_HDR1, RA_HDR2 = CLASS1.ra_hdr[:-2], CLASS2.ra_hdr[:-2]
DEC_HDR1, DEC_HDR2 = CLASS1.dec_hdr[:-2], CLASS2.dec_hdr[:-2]
# For plot labels. `AXLABEL` is either 'true' or 'meas' #
AXLABEL1, AXLABEL2 = CLASS1.axlabel, CLASS2.axlabel

# Magnitudes #
M_HDR1, M_HDR2 = CLASS1.mag_hdr, CLASS2.mag_hdr
M_AXLABEL1, M_AXLABEL2 = CLASS1.mag_axlabel, CLASS2.mag_axlabel

# For magnitude error calculation #
M_ERR_HDR1, M_ERR_HDR2 = CLASS1.mag_err_hdr, CLASS2.mag_err_hdr
CM_FLUX_HDR1, CM_FLUX_HDR2 = CLASS1.cm_flux_hdr, CLASS2.cm_flux_hdr
CM_FLUX_COV_HDR1, CM_FLUX_COV_HDR2 = CLASS1.cm_flux_cov_hdr, CLASS2.cm_flux_cov_hdr

# For signal to noise calculation #
CM_T_HDR1, CM_T_HDR2 = CLASS1.cm_t_hdr, CLASS2.cm_t_hdr
CM_T_ERR_HDR1, CM_T_ERR_HDR2 = CLASS1.cm_t_err_hdr, CLASS2.cm_t_err_hdr
CM_T_S2N_AXLABEL1, CM_T_S2N_AXLABEL2 = CLASS1.cm_t_s2n_axlabel, CLASS2.cm_t_s2n_axlabel

# Flags #
FLAGS_HDR1, FLAGS_HDR2 = CLASS1.flags_hdr, CLASS2.flags_hdr
OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2 = CLASS1.obj_flags_hdr, CLASS2.obj_flags_hdr
# psf_flags is a string of form '(0,0,0,0)'; must pass through get_floats_from_string() if used. #
PSF_FLAGS_HDR1, PSF_FLAGS_HDR2 = CLASS1.psf_flags_hdr, CLASS2.psf_flags_hdr
CM_FLAGS_HDR1, CM_FLAGS_HDR2 = CLASS1.cm_flags_hdr, CLASS2.cm_flags_hdr
CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2 = CLASS1.cm_max_flags_hdr, CLASS2.cm_max_flags_hdr
CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 = CLASS1.cm_flags_r_hdr, CLASS2.cm_flags_r_hdr
CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2 = CLASS1.cm_mof_flags_hdr, CLASS2.cm_mof_flags_hdr

# For quality cuts introduced by Eric Huff #
CM_S2N_R_HDR1, CM_S2N_R_HDR2 = CLASS1.cm_s2n_r_hdr, CLASS2.cm_s2n_r_hdr
PSFREC_T_HDR1, PSFREC_T_HDR2 = CLASS1.psfrec_t_hdr, CLASS2.psfrec_t_hdr

# For region file #
MAJOR_AX_HDR1, MAJOR_AX_HDR2 = CLASS1.a_hdr, CLASS2.a_hdr 
MINOR_AX_HDR1, MINOR_AX_HDR2 = CLASS1.b_hdr, CLASS2.b_hdr
ANGLE1, ANGLE2 = CLASS1.angle, CLASS2.angle

FLAG_HDR_LIST = [ FLAGS_HDR1, FLAGS_HDR2, CM_FLAGS_HDR1, CM_FLAGS_HDR2, CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2, OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2, PSF_FLAGS_HDR1, PSF_FLAGS_HDR2, CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2, CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 ]


### For `corner.hist2d` contours ###
# Percentile levels for contours. Default is 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2) via https://github.com/dfm/corner.py/blob/master/corner/corner.py #
# Correct 1sigma for a 2D histogram is given by http://corner.readthedocs.io/en/latest/pages/sigmas.html #
##FIXME LVLS = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5, 2.2]) ** 2)

# `LVLS` passed to `corner.hist2d` interpreted as percentiles #
LVLS = 1.0 - np.exp(-0.5 * np.array([1.5]) ** 2)
#CLRS = ['magenta', 'red', 'blue', 'yellow']
CLRS = ['red']
# FIXME `colors` passed to contour_kwargs are passed/applied in reversed order so reverse order for labels #
CLRS_LABEL = CLRS[::-1]


### Dictionaries for color coding ###
CMAPS = {'g':'Greens', 'r':'Purples', 'i':'Greys', 'z':'Blues'}
PT_COLORS = {'g':'green', 'r':'purple', 'i':'darkgrey', 'z':'navy'}
HIST_COLORS = {'g':'lime', 'r':'magenta', 'i':'dimgrey', 'z':'cyan'}

# For computing color #
BANDS_FOR_COLOR = {'g':'r', 'r':'i', 'i':'z'}
#TODO replace bandFollow

# For writing to log files #
WRITE_COLORS = {'g':'g-r', 'r':'r-i', 'i':'i-z'}
BAND_INDEX = {'g':0, 'r':1, 'i':2, 'z':3}

### For completeness plot ###
#FIXME to keep bin size uniform across multiple tiles getting rid of: __bins = np.arange(int(np.min(truth_mag)), int(np.max(truth_mag)), __step)
COMPLETENESS_MAG_BINS = np.arange(15, 30, 1)
COMPLETENESS_PLOT_MAG_BINS = []
for i in np.arange(0, len(COMPLETENESS_MAG_BINS)-1):
	COMPLETENESS_PLOT_MAG_BINS.append(np.median([COMPLETENESS_MAG_BINS[i], COMPLETENESS_MAG_BINS[i+1]]))


# Used if LOG_FLAGS is True #
FLAG_IDX = []


### For plot names, plot titles, log file names ###
TITLE_PIECE1, TITLE_PIECE2 = CLASS1.title_piece, CLASS2.title_piece
MATCH_TYPE = get_match_type(title_piece1=TITLE_PIECE1, title_piece2=TITLE_PIECE2)

if STACK_TILES or STACK_REALIZATIONS:
	log_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE)
	if os.path.isdir(log_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(log_dir) + ' does not exist. \n Change directory structure.')
		if NO_DIR_MAKE:
			print 'Making directory ', log_dir, '...\n'
			os.makedirs(log_dir)

	FN_STACK_MAIN_LOG = os.path.join(log_dir, 'log_'+MATCH_TYPE+'.csv')
	with open(FN_STACK_MAIN_LOG, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		# Write headers #
		writer.writerow(['TILE', 'REALIZATION', 'BAND', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJECTS_IN_MATCH', 'TOTAL_FLAGGED_OBJECTS_IN_TRUTH', 'TOTAL_1SIGMA_MAG_MEAS', 'PERCENT_RECOVERED_FLAGS_INCLUDED', 'PERCENT_RECOVERED_FLAGS_REMOVED'])




def write_log_file_headers(fn_main_log, fn_mag_err_log, fn_flag_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log):
	"""Write headers to log files.

        Parameters
        ----------

	Returns
	-------
	"""

	### Main log ###
	with open(fn_main_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		# Write headers #
		writer.writerow(['TILE', 'REALIZATION', 'BAND', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJS_IN_MATCH', 'TOTAL_FLAGGED_OBJS_IN_TRUTH', 'TOTAL_1SIGMA_MAG_MEAS', '%_RECOVERED_FLAGS_IN', '%_RECOVERED_FLAGS_RM', 'RUN_TYPE'])


	### Log for magnitude error calculation ###
	with open(fn_mag_err_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		if 'truth' in MATCH_CAT1 and SWAP_HAX is False:
			# mag_true #
			writer.writerow(['TILE', 'REALIZATION', 'BAND', 'NUM_OBJS_IN_BIN', 'MAG_TRUE_BIN_LHS', 'MAG_TRUE_BIN_RHS', 'MEDIAN_HAX_MAG_TRUE', 'MEDIAN_MAG_ERROR'])
		if 'truth' in MATCH_CAT1 and 'truth' not in MATCH_CAT2 and SWAP_HAX:
			# mag_meas #
			writer.writerow(['TILE', 'REALIZATION', 'BAND', 'NUM_OBJS_IN_BIN', 'MAG_MEAS_BIN_LHS', 'MAG_MEAS_BIN_RHS', 'MEDIAN_HAX_MAG_MEAS', 'MEDIAN_MAG_ERROR'])	


	### Log file for color calculation ###
	with open(fn_color_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		if 'truth' in MATCH_CAT1:
			writer.writerow(['TILE', 'REALIZATION', 'COLOR', 'MAG_TRUE_BIN', 'OBJS_IN_MAG_BIN', 'TOTAL_OBJECTS_NO_FLAGS', 'RUN_TYPE'])
		if 'truth' not in MATCH_CAT1:
			writer.writerow(['TILE', 'REALIZATION', 'COLOR', 'MAG_MEAS_BIN', 'OBJS_IN_MAG_BIN', 'TOTAL_OBJECTS_NO_FLAGS', 'RUN_TYPE'])


	### Log file for magnitude outliers ###
	with open(fn_mag_diff_outliers_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow(['TILE', 'REALIZATION', 'BAND', 'MAG1', 'MAG2', 'MAG_DIFF'])


	### Log file for magnitude completeness calculation ###
	with open(fn_mag_completeness_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow(['TILE', 'REALIZATION', 'BAND', 'INJ_PERCENT', 'TRUTH_MAG_BIN_LHS', 'TRUTH_MAG_BIN_RHS', 'MATCH_CAT_OBJS_IN_BIN', 'TRUTH_CAT_OBJS_IN_BIN'])


	### Log file for flags... ###
	if LOG_FLAGS:
		with open (fn_flag_log, 'wb') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow(['TILE', 'REALIZATION', 'FLAG1_HEADER', 'FLAG2_HEADER', 'FLAG1_VALUE', 'FLAG2_VALUE', 'MAG1', 'MAG2', 'RUN_TYPE'])

	return fn_main_log, fn_mag_err_log, fn_flag_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log




################################################################### Analysis ###################################################################
def get_floats_from_string(df, band, four_elmt_arrs_hdr):
	"""Transform a list of strings of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]' to a list of floats of form '[1,1,...]' (if band="g"), '[2,2,...]' ("r"), '[3,3,...]' ("i"), or '[4,4,...]' ("z"). This is necessary for CSVs created from ms_matcher or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings. 

	Parameters
	----------
        df (pandas DataFrame)
        
	band (str)
		Allowed values: 'g' 'r' 'i' 'z'.
        
	hdr (str)
		Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '(1,2,3,4)'.
        
	Returns
	-------
        __elmts (list of floats)
		Collection of the numbers corresponding to a particular index in a list of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]. 
	"""

	strings = df[four_elmt_arrs_hdr]; __elmts = []

	# Each element (elmt) is of form '(1, 2, 3, 4)' #
	for elmt in strings:

		if band == 'g':
			i = 1
			idx1 = elmt.find('(') + i
			idx2 = elmt.find(',')

		if band == 'r':
			i = 2
			idx1 = elmt.replace(',', ';', 0).find(',') + i
			idx2 = elmt.replace(',', ';', 1).find(',')

		if band == 'i':
			i = 2
			idx1 = elmt.replace(',', ';', 1,).find(',') + i
			idx2 = elmt.replace(',', ';', 2).find(',')

		if band == 'z':
			i = 2
			idx1 = elmt.replace(',', ';', 2).find(',') + i
			idx2 = elmt.find(')')

		__elmts.append(float(elmt[idx1:idx2]))

	if PRINTOUTS_MINOR:
		print 'Got ', four_elmt_arrs_hdr, ' for band ', band, '...'
		print ' Check: ', strings[0], ' & ', __elmts[0], '\n'


	return __elmts




def get_matrix_diagonal_element(df, band, sq_matrices_hdr):
	"""Transforms a list of 4x4 matrices where each element is a string of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' into a list of either the 11 (if band is "g"), 22 ("r"), 33 ("i"), or 44 ("z") matrix elements. This is necessary for CSVs created from ms_matcher or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings.

	Parameters
	----------
	df (pandas DataFrame)

	band (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	sq_matrices_hdr (str)
		Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))'.

        Returns
	-------
	__matrix_diagonals (list of floats)
		Collection of the numbers corresponding to a particular diagonal element in a list of 4-by-4 matrices.
	"""

	matrices = df[sq_matrices_hdr]; __matrix_diagonals = []

	# Each element in `matrices` is a matrix of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' #
	for matrix in matrices:

		if band == 'g':
 			i, j = 2, 0
			idx1 = 0
			idx2 = matrix.find(',')

		if band == 'r':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 4).find(',')
			idx2 = matrix.replace(',', ';', 5).find(',')

		if band == 'i':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 9).find(',')
			idx2 = matrix.replace(',', ';', 10).find(',')

		if band == 'z':
			i, j = 2, -1
			idx1 = matrix.replace(',', ';', 14).find(',')
			idx2 = matrix.replace(',', ';', 15).find(',')

		__matrix_diagonals.append(float(matrix[idx1+i:idx2+j]))

	if PRINTOUTS_MINOR:
		print 'Got ', sq_matrices_hdr, ' for band ', band
		print ' Check: ', matrices[0], ' & ', __matrix_diagonals[0], '\n'


	return __matrix_diagonals




def get_good_indices_using_primary_flags(df_1and2, full_mag1, full_mag2, cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2, band):
	"""Get indices of objects without flags* where flags* used are indicated in README.md. Store flagged indices if PLOT_FLAGGED_OBJS is True.

	Parameters
	----------
	df_1and2 (pandas DataFrame)

	band (str)
		Used if analysing Y3 Gold catalog. Can be `None`.

	full_mag1, full_mag2 (list of floats)
		Uncleaned lists containing magnitudes. 

	Returns
	-------
	__idx_good (list of ints)
	
	__idx_bad (list of ints)
		Is empty if PLOT_FLAGGED_OBJS is False.
	"""

	if cm_flag_hdr2 is None and cm_flag_hdr1 is None and flag_hdr1 is None and flag_hdr2 is None:
                sys.exit('No headers to clean flags with...')

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if cm_flag_hdr2 is None:
                cm_flag_hdr2 = cm_flag_hdr1

        if cm_flag_hdr1 is None:
                cm_flag_hdr1 = cm_flag_hdr2

        if flag_hdr1 is None:
                flag_hdr1 = flag_hdr2

        if flag_hdr2 is None:
                flag_hdr2 = flag_hdr1

        ### Get flags ###
        flag1, flag2 = df_1and2[flag_hdr1], df_1and2[flag_hdr2]
        cm_flag1, cm_flag2 = df_1and2[cm_flag_hdr1], df_1and2[cm_flag_hdr2]


	# Make arrays to take absolute value in next step #
	full_mag1, full_mag2 = np.array(full_mag1), np.array(full_mag2)

	# Get rid of these objects; 37.5 corresponds to a negative flux #
	if 'y3_gold' not in MATCH_CAT1 and 'y3_gold' not in MATCH_CAT2:
		__idx_good = np.where( (abs(full_mag1) != 9999.0) & (abs(full_mag1) != 99.0) & (abs(full_mag1) != 37.5) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & (abs(full_mag2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) )[0]

	# Additional flags for Y3 Gold MOF catalog. # 
	if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
		# Get flag header #
		if 'y3_gold' in MATCH_CAT1:
			suf = '_1'
		if 'y3_gold' in MATCH_CAT2:
			suf = '_2'
		if Y3_MOF:
			hdr = 'MOF_CM_MOF_FLAGS' + suf
		# Duplicate of cm_flag_hdr because no analogous SOF flag for 'MOF_CM_MOF_FLAGS' #
		if Y3_MOF is False:
			hdr = 'SOF_CM_FLAGS' + suf

		__idx_good = np.where( (df_1and2['SEXTRACTOR_FLAGS_'+band.upper()+suf] == 0) & (df_1and2['IMAFLAGS_ISO_'+band.upper()+suf] == 0) & (df_1and2[hdr] == 0) & (abs(full_mag1) != 9999.0) & (abs(full_mag1) != 99.0) & (abs(full_mag1) != 37.5) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & (abs(full_mag2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) )[0]

		
	if PLOT_FLAGGED_OBJS:
		__idx_bad = np.where( (abs(full_mag1) != 9999.0) & (abs(full_mag1) != 99.0) & (abs(full_mag1) != 37.5) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & (abs(full_mag2) != 9999.0) & (abs(full_mag2) != 99.0) & ((flag2 != 0) | (flag1 != 0) | (cm_flag1 != 0) | (cm_flag2 != 0)) )[0]

	if PLOT_FLAGGED_OBJS is False:
		__idx_bad = None

		
        if PRINTOUTS:
		if 'y3_gold' not in MATCH_CAT1 and 'y3_gold' not in MATCH_CAT2:
			print 'Eliminated ', len(full_mag1) - len(__idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ' ... \n'
		# For Y3 #
		if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
			print 'Eliminated ', len(full_mag1) - len(__idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ',', hdr, ', ', 'SEXTRACTOR_FLAGS_'+band.upper(), ', ',  'IMAFLAGS_ISO_'+band.upper(),  ' ... \n'

	return __idx_good, __idx_bad	




def log_flags(df_1and2, flag_hdr1, flag_hdr2, band, full_mag1, full_mag2, realization, tile, fn_flag_log):
	"""Examine a particular flag and write to log file. Can also be used to check all flags in a list of flags.

	Parameters
	----------
	df_1and2 (pandas DataFrame)

	band (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	full_mag1, full_mag2 (numpy.ndarray if directly from `df[hdr]` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.

	realization (str)
		Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.

	tile (str) 

        Returns
	-------
	__idx_good (list of ints)
		Indices of objects with flags values of zero.

	__idx_bad (list of ints)
		Indices of objects with nonzero flag values.
	"""

	__idx_good, __idx_bad = [], []; counter_idx_bad = 0

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if flag_hdr1 is None:
		flag_hdr1 = flag_hdr2

	if flag_hdr2 is None:
		flag_hdr2 = flag_hdr1


	### psf_flags are strings of form '(0,0,0,0)' ###
	if flag_hdr1 is not None and flag_hdr2 is not None and 'psf' not in flag_hdr1 and 'psf' not in flag_hdr2:
		flag1 = df_1and2[flag_hdr1]
		flag2 = df_1and2[flag_hdr2]

	if 'psf' in flag_hdr1 and 'psf' in flag_hdr2:
		flag1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flag_hdr1, band=band)
		flag2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flag_hdr2, band=band)


	### Check for flags ###
	for i in np.arange(0, len(full_mag1)):

		if abs(full_mag1[i]) != 9999.0 and abs(full_mag2[i]) != 9999.0 and full_mag1[i] != 37.5 and full_mag2[i] != 37.5 and full_mag1[i] != 99.0 and full_mag2[i] != 99:
			
			if flag1[i] == 0 and flag2[i] == 0:
				__idx_good.append(i)
	
			if flag1[i] != 0 or flag2[i] != 0:
				__idx_bad.append(i)
                                counter_idx_bad += 1

				# Append to csv #
				with open(fn_flag_log, 'a') as csvfile:
					writer = csv.writer(csvfile, delimiter=',')
					# TILE, REALIZATION, BAND, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2, RUN_TYPE #
					writer.writerow([str(tile), str(realization), str(band), str(flag_hdr1), str(flag_hdr2), str(flag1[i]), str(flag2[i]), str(full_mag1[i]), str(full_mag2[i]), str(RUN_TYPE)])


	if PRINTOUTS:
		print 'For tile: ', str(tile), ' and band: ', str(band), ', checked flags: ', flag_hdr1, ' & ', flag_hdr2, '...'

	### Check if flags were found ###
	if counter_idx_bad > 0 and PRINTOUTS:
		print ' Number of flags for magnitudes values 9999, 99, 37.5 and flags ', str(flag_hdr1), ' and ', str(flag_hdr2), ' : ', counter_idx_bad, '\n'


	return __idx_good, __idx_bad 




def calculate_measured_magnitude_error_from_flux(flux_cov_hdr, df_1and2, band, flux_hdr, idx_good, match_cat):
	"""Calculate the magnitude error via 1.08 * (flux_cov[i][i])^0.5 / flux[i]. Ignore flagged objects in error calculation. This error is a FRACTIONAL error.

	Parameters
	----------
	cov_hdr (str) 
		Header for flux covariance matrix in the matched catalog.

	df_1and2 (pandas DataFrame)

	band (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	flux_hdr (str)
		Headers refer to column names in the matched catalog.

	idx_good (list of ints)
		Indices with flag values equal to zero.

	Returns
	-------
	__mag_err_from_flux (list of floats)
		The magnitude error corresponding to each object.
	"""

	# Uncleaned lists for flux and flux covariance #
	full_flux = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr, band=band)
	full_flux_cov = get_matrix_diagonal_element(df=df_1and2, sq_matrices_hdr=flux_cov_hdr, band=band)

	__mag_err_from_flux, flux, fluxcov = [], [], []; counter_neg = 0

	# 'Clean' indices #
	for idx in idx_good:
		flux.append(full_flux[idx])
		fluxcov.append(full_flux_cov[idx])

	# Calculations #
	for i in np.arange(0, len(flux)):

		# Throw out negative fluxcov (error calculation involves taking the square root of a negative) #
		if fluxcov[i] < 0:
			__mag_err_from_flux.append(0)
			counter_neg += 1

		if fluxcov[i] == 0:
			print 'cm_flux_cov is 0'

		if fluxcov[i] > 0:
			# Pogsons number = 1.08 #
			__mag_err_from_flux.append(1.08*fluxcov[i]**0.5/flux[i])

	if PRINTOUTS:
		print 'Calculated the magnitude error for band: ', band
		print ' Number of negative cm_flux_cov: ', counter_neg, ' / ', len(flux), '\n'
	#TODO do not compute __mag_err_from_flux as above if truth cat
	if 'truth' in match_cat: __mag_err_from_flux = None

	return __mag_err_from_flux 




#TODO mag_bins, rename func get_68percentile_of_magnitude_difference....
def get_68percentile_of_normalized_magnitude_difference(binned_norm_mag_diff, mag_bins_for_mag_err, binned_hax_mag):
	"""Calculate the point on the normalized vertical axis corresponding to the 68th percentile of each bin used in the error calculation. This percentile can be calculated in several ways: centered about zero, centered about the median of each bin, calculated using the absolute value of the data, calculated by examining the 34th percentile of the positive and negative portions of the data separately. 

	Parameters
	----------
	binned_norm_mag_diff (list of list of floats)
		Normalized delta magnitudes. Bin structure preserved.


	binned_hax_mag (list of list of floats)
		Magnitudes on the horizontal axis. Bin structure preserved.

	Returns
	-------
	__norm_mag_diff_68p (list of floats)
		Point on the vertical axis (vax) corresponding to 68 percentile. Each element in the list corresponds to a different bin.

	bins (list of floats)
		Bins used in error calculation.

	__neg_norm_mag_diff_34p (list of floats)

	__pos_norm_mag_diff_34p (list of floats)

	"""

	__norm_mag_diff_68p, __neg_norm_mag_diff_34p, __pos_norm_mag_diff_34p = [], [], []

	# Loop through bins (b) #
	for b in np.arange(0, len(binned_norm_mag_diff)):

		if binned_norm_mag_diff[b] is None:
			__norm_mag_diff_68p.append(None)
			__neg_norm_mag_diff_34p.append(None)
			__pos_norm_mag_diff_34p.append(None)

		if binned_norm_mag_diff[b] is not None:

			### Find 68th percentile about zero ### 	
			# Values in current bin (icb) #
			vax_mag_bins_icb = binned_norm_mag_diff[b]
			# Take absolute value of each point in bin #
			abs_vax_mag_bins_icb = [abs(elmt) for elmt in vax_mag_bins_icb]
			# Percentile sorts the data #
			__norm_mag_diff_68p.append(np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'))	

			if PRINTOUTS_MINOR:
				# Check the percentile because interpolation='lower' was used #
				num = 0
				for j in np.arange(0, len(binned_norm_mag_diff[b])):
					if abs(binned_norm_mag_diff[b][j]) <= np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'):
						num += 1
				print 'Number of objects within 68 percentile via np.percentile(interpolation=lower): ', float(num)/len(binned_norm_mag_diff[b]), '...\n'


			### Find 34th percentile of positive and negative values separately ###
			neg_vax, pos_vax = [], []
			counter_neg, counter_pos = 0, 0
			for j in np.arange(0, len(vax_mag_bins_icb)):
				if vax_mag_bins_icb[j] < 0:
					neg_vax.append(vax_mag_bins_icb[j])
					counter_neg += 1
				if vax_mag_bins_icb[j] > 0:
					pos_vax.append(vax_mag_bins_icb[j])
					counter_pos += 1

			# Check if lists are populated #
			if counter_neg > 0: 
				__neg_norm_mag_diff_34p.append(np.percentile(neg_vax, 34, interpolation='lower'))
			if counter_pos > 0:
				__pos_norm_mag_diff_34p.append(np.percentile(pos_vax, 34, interpolation='lower'))
			if counter_neg  == 0:
				__neg_norm_mag_diff_34p.append(None)
			if counter_pos  == 0:
				__pos_norm_mag_diff_34p.append(None)


			plot_hist = False
			# Plot histogram to see distrubtion of data (data is not normally distributed) #
                        if plot_hist:
                                plt.figure()
                                norm_dm = [abs(elmt) for elmt in binned_norm_mag_diff[b]]
                                plt.hist(norm_dm)
                                plt.title('Bin LHS: ' + str(mag_bins_for_mag_err[b]))
                                plt.xlabel(r'$\Delta M$')
                                plt.axvline(x=0, color='black', linestyle=':', linewidth=0.5)
                                plt.show()


	return __norm_mag_diff_68p, mag_bins_for_mag_err, __neg_norm_mag_diff_34p, __pos_norm_mag_diff_34p




def bin_and_cut_measured_magnitude_error(clean_mag1, clean_mag2, mag_err1, mag_err2, band, tile, realization, fn_mag_err_log, fn_mag_diff_outliers_log):
        """Clean error. Bin error according to horizontal axis of plot. Remove error values corresponding to objects with |DeltaMagnitude|>3. Do not consider error corresponding to empty bins nor bins with a small number of objects.

        Parameters
	----------
	clean_mag1, clean_mag2 (list of floats)
		Objects with flag values of zero and/or quality cuts performed.
                
	mag_err1, mag_err2 (list of floats)
		1 and 2 refer to the matched catalogs. 
        
	Returns
	-------
	binned_hax_mag_median (list of floats)
		List of medians of the horizontal axis magnitude in each bin.
              
	binned_vax_mag_median (list of floats)
		List of medians of the vertical axis magnitude in each bin. Vertical axis is computed via clean_mag1 - clean_mag2.

	binned_err_median (list of floats)
		Median of the error in each bin.

	bins (list of floats)
		Bins used. Binned according to horizontal axis.
        """

	### !!!!! Comment this block out if errors are to be computed using both catalogs regardless of origin (measured catalog or truth catalog) ###
	if 'meas' in AXLABEL1 and 'meas' not in AXLABEL2 and PRINTOUTS:
			print 'Using measured catalog (', MATCH_CAT1, ') for error calculation ... '

	if 'meas' in AXLABEL2 and 'meas' not in AXLABEL1 and PRINTOUTS:
			print 'Using measured catalog (', MATCH_CAT2, ') for error calculation ... '

	if 'meas' in AXLABEL1 and 'meas' in AXLABEL2 and PRINTOUTS:
		print 'Using measured catalogs (', MATCH_CAT1, '&', MATCH_CAT2, ') for error calculation ... '

	if 'true' in AXLABEL1 and 'true' in AXLABEL2:
		sys.exit('Errors are to be computed using the measured catalog(s), not the truth catalog(s).')


	### Define bins ###
        __step = 0.5
        # Find the absolute min and max of the magnitudes in the matched catalog #
        __lim_low1, __lim_low2 = min(clean_mag1), min(clean_mag2)
        __lim_high1, __lim_high2 = max(clean_mag1), max(clean_mag2)
        __lim_low, __lim_high = min([__lim_low1, __lim_low2]), max([__lim_high1, __lim_high2])

	# Define bin limits by ints #
	__lim_low, __lim_high = int(__lim_low), int(__lim_high)
	# Introduce magnitude cutoff to tame errors #
	if 'gal_truth' in (MATCH_CAT1, MATCH_CAT2):
		__lim_high = 26
	if 'star_truth' in (MATCH_CAT1, MATCH_CAT2):
		__lim_high = 24

        if PRINTOUTS:
                print ' Forcing magnitudes to be binned with max ', __lim_high, '...'

	# Include endpoint #
        __mag_bins_for_mag_err = np.arange(__lim_low, __lim_high+(__step*0.1), __step)


	# Stores median of values in each bin #
	__hax_mag_bin_medians, __mag_diff_bin_medians, __mag_err_bin_medians = [], [], []
	# List of lists. Stores all values in each bin #
	__binned_hax_mag, __binned_mag_diff, __binned_mag_err = [], [], []
	# Counter for empty bins and bins with a small number of objects #
        __counter_empty_bin = 0


        # Bin magnitude errors according to the magnitude on the horizontal axis #
        if SWAP_HAX:
                __hax_mag = clean_mag2
        if SWAP_HAX is False:
                __hax_mag = clean_mag1

	# Magnitude on the vertical axis (vax) #
	# __vax_mag --> __mag_diff
	__vax_mag = np.array(clean_mag1) - np.array(clean_mag2)

	__clean_mag_diff, __clean_hax_mag = [], []


	__cutoff_mag_diff = 3
	__cutoff_mag_diff_for_log = 1
	__counter_mag_diff_geq1 = 0

	### Populate each bin ###
	for b in __mag_bins_for_mag_err: 

		__hax_mag_in_bin, __vax_mag_in_bin, __err_mag_in_bin = [], [], []
		__counter = 0

                for i in np.arange(0, len(clean_mag1)):

			if abs(__vax_mag[i]) > __cutoff_mag_diff_for_log:
				__counter_mag_diff_geq1 += 1
				# Append to csv #
				with open(fn_mag_diff_outliers_log, 'a') as csvfile:
					writer = csv.writer(csvfile, delimiter=',')
					# TILE, REALIZATION, BAND, MAG1, MAG2, MAG_DIFF #
					writer.writerow([str(tile), str(realization), str(band), str(clean_mag1[i]), str(clean_mag2[i]), str(__vax_mag[i])])

                        # Do not calculate errors using outlier magnitudes (given by `__cutoff_mag_diff`). Bin magnitude errors according to the magnitude on the horizontal axis of the plot #
			#TODO is is necessary to do this comparison to b and b+__step?
                        if __hax_mag[i] >= b and __hax_mag[i] < b+__step and abs(__vax_mag[i]) < __cutoff_mag_diff: 
				###__err_mag_in_bin.append((mag_err1[i]**2 + mag_err2[i]**2)**0.5)
				if mag_err1 is None: __err_mag_in_bin.append(abs(mag_err2[i]))
				if mag_err2 is None: __err_mag_in_bin.append(abs(mag_err1[i])) 
				if mag_err1 is not None and mag_err2 is not None: __err_mag_in_bin.append((mag_err1[i]**2 + mag_err2[i]**2)**0.5)
				__vax_mag_in_bin.append(__vax_mag[i])
				__hax_mag_in_bin.append(__hax_mag[i])
				__clean_mag_diff.append(__vax_mag[i])
				__clean_hax_mag.append(__hax_mag[i])
                                __counter += 1

		# Written in log file, hence 'minor' #
                if PRINTOUTS_MINOR:
                        print ' For magnitude, number of objects in bin ', round(j, 2), '-', round(j+__step, 2), ': ', __counter, '...'


		### Write to log file ###
		if __counter == 0: write_median, write_err = None, None
		if __counter > 0: write_median, write_err = np.median(__hax_mag_in_bin), np.median(__err_mag_in_bin)
		# Append to csv #
		with open(fn_mag_err_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			# TILE, REALIZATION, BAND, NUM_OBJS_IN_BIN, BIN_LHS, BIN_RHS, MEDIAN_HAXIS_MAG, MEDIAN_ERROR #
			writer.writerow([str(tile), str(realization), str(band), str(__counter), str(round(b, 2)), str(round(b+__step, 2)), str(write_median), str(write_err)])


                ### Tame error calculation and normalization by adding `None` to empty bins and bins with a small number of points ###
		# Define minimum bin population # 
		if STACK_REALIZATIONS: __min_bin_pop = 30
		if STACK_REALIZATIONS is False: __min_bin_pop = 10 

		if __counter <= __min_bin_pop:
                        __counter_empty_bin += 1
                        __mag_err_bin_medians.append(None)
                        __hax_mag_bin_medians.append(None)
                        __mag_diff_bin_medians.append(None)
			# Add to list of lists to keep bin structure #
			__binned_mag_err.append(None)
                        __binned_hax_mag.append(None)
                        __binned_mag_diff.append(None)		

                if __counter > __min_bin_pop:
                        __mag_err_bin_medians.append(np.median(__err_mag_in_bin))
                        __hax_mag_bin_medians.append(np.median(__hax_mag_in_bin))
                        __mag_diff_bin_medians.append(np.median(__vax_mag_in_bin))
			# Add to list of lists to keep bin structure #	
			__binned_mag_err.append(__err_mag_in_bin)
                        __binned_hax_mag.append(__hax_mag_in_bin)
                        __binned_mag_diff.append(__vax_mag_in_bin)


	if PRINTOUTS:
                if SWAP_HAX:
                        print ' Binned clean_mag2 (from', MATCH_CAT2, ') with step size: ', __step, ', and minimum: ', __lim_low, ', and maximum: ', __lim_high, '...'
		if SWAP_HAX is False:
                        print ' Binned clean_mag1 (from', MATCH_CAT1, ') with step size: ', __step, ', and minimum: ', __lim_low, ', and maximum: ', __lim_high, '...'
		print ' Calculated magnitude error using objects where |DeltaM| <', __cutoff_mag_diff, ' ... '
		print ' Excluded ', __counter_empty_bin, ' bins with less than ', __min_bin_pop, ' objects ... \n'
		print ' Number of objects with |DeltaM| > ', __cutoff_mag_diff_for_log, ' : ', __counter_mag_diff_geq1, ' / ', str(len(clean_mag1)), ' ...\n'

        return __hax_mag_bin_medians, __mag_diff_bin_medians, __mag_err_bin_medians, __mag_bins_for_mag_err, __binned_hax_mag, __binned_mag_diff, __clean_hax_mag, __clean_mag_diff




#TODO rename func... normalize to err? must this be a magnitude difference?
# must be a magnitude difference bc mag diff from bin_and_cut...()
def normalize_magnitude_plot_maintain_bin_structure(clean_mag1, clean_mag2, mag_err1, mag_err2, band, tile, realization, fn_mag_err_log, fn_mag_diff_outliers_log):
	"""Normalize the vertical axis using 1sigma_mag and preserve the bin structure of the horizontal axis. 

	Parameters
	----------
	clean_mag1, clean_mag2 (list of floats)

	mag_err1, mag_err2 (list of floats)

	Returns
	-------
	norm_dm_bins (list of list of floats)

	bins (list of floats)
		Bins used in error calculation. 
	"""

	# List of lists. Stores all values in each bin #
	__binned_norm_mag_diff, __binned_hax_mag = [], []
	# Stores the median of each bin #
	__norm_mag_diff_bin_medians = []

	#__mag_err_bin_medians
	haxBinMedian, vaxBinMedian, magErrBinMedians, magBinsForMagErr, haxBins, vaxBins, hax, vax = bin_and_cut_measured_magnitude_error(clean_mag1=clean_mag1, clean_mag2=clean_mag2, mag_err1=mag_err1, mag_err2=mag_err2, band=band, tile=tile, realization=realization, fn_mag_err_log=fn_mag_err_log, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)
 

	# Loop through bins (b) #
	for b in np.arange(0, len(vaxBins)):

		# Normalized Delta-Magnitudes (dm) in current bin (icb) #
		__norm_vax_in_bin, __hax_in_bin = [], []

		# 0 is a placeholder for empty bins and bins with few objects #
		if magErrBinMedians[b] is None:
			__binned_norm_mag_diff.append(None)
			__binned_hax_mag.append(None)
			__norm_mag_diff_bin_medians.append(None)

		#if vax_mag_icb != 0:
		if magErrBinMedians[b] is not None:
			for i in np.arange(0, len(vaxBins[b])):
				__norm_vax_in_bin.append(vaxBins[b][i]/magErrBinMedians[b])
				__hax_in_bin.append(haxBins[b][i])

			# List of lists to keep bin structure #
			__binned_norm_mag_diff.append(__norm_vax_in_bin)
			__binned_hax_mag.append(__hax_in_bin)
			__norm_mag_diff_bin_medians.append(np.median(__norm_vax_in_bin))

	return __binned_norm_mag_diff, magBinsForMagErr, __binned_hax_mag, magErrBinMedians, __norm_mag_diff_bin_medians




def normalize_magnitude_difference_plot(binned_norm_mag_diff, mag_bins_for_mag_err, binned_hax_mag):
	"""Normalize plot to 1sigma_mag curve using tame magnitude errors

	Parameters
	----------
	norm_dm_bins (list of list of floats)
		Normalized delta magnitudes in each bin. 

	mag_bins_for_mag_err (list of floats)
		Magnitude bins used in error calculation.

	hax_mag_bins (list of list of floats)
		Magnitudes on the horizontal axis. Bin structure preserved.

        Returns
	-------
	norm_dm (list of floats)
		Delta-Magnitude normalized by error. Delta-Magnitude computed via magnitude1 - magnitude2. 

        hax_mag (list of floats)
		Magnitude to be plotted on the horizontal axis.
	"""
	print 'len(mag_bins_for_mag_err)', len(mag_bins_for_mag_err)
	### Remove `None` so that lists can be flattened. `None` is a placeholder for missing lists due to empty or small bin. ###
	__idx_no_none = []
        for i in np.arange(0, len(binned_norm_mag_diff)):
                if binned_norm_mag_diff[i] is not None:
                        __idx_no_none.append(i)

	__mag_bins_for_err_no_none, __binned_norm_mag_diff, __binned_hax_mag = np.array(mag_bins_for_mag_err)[__idx_no_none], np.array(binned_norm_mag_diff)[__idx_no_none], np.array(binned_hax_mag)[__idx_no_none] 

	### Flatten lists ###
	__hax_mag = [item for sublist in __binned_hax_mag for item in sublist]
	__norm_mag_diff = [item for sublist in __binned_norm_mag_diff for item in sublist]

	# Get __idx_relevant #	
	__idx_relevant = np.where((__hax_mag >= min(mag_bins_for_mag_err)) & (__hax_mag < max(mag_bins_for_mag_err)))[0]
	__hax_mag, __norm_mag_diff = np.array(__hax_mag)[__idx_relevant], np.array(__norm_mag_diff)[__idx_relevant]

	return __norm_mag_diff, __hax_mag, __mag_bins_for_err_no_none 




#TODO bins--> mag_bins, are any of these lists of lists?
def one_sigma_magnitude_counter(mag_diff, clean_mag1, mag_bins_for_mag_err, hax_mag, binned_mag_err, vax_mag_bin_medians):
	"""Find the number of objects within 1sigma_mag. This function is called if `NORMALIZE` is False.

	Parameters
	----------
	mag_diff (list of floats)
		NON-normalized Delta-Magnitude.

	Returns
	-------
	__num_1sigma_mag (int)
		Number of objects within 1sigma_mag curve.
        """

	if len(mag_bins_for_mag_err) != len(binned_mag_err):
                sys.exit('len(mag_bins_for_mag_err) not equal to len(error)') #TODO

	tot = len(mag_diff)
	__num_1sigma_mag = 0; counter_objs_considered = 0

	
	# Cutoffs were introduced in error calculation. Consider only points not cutoff #
	# Get rid of `None` placeholders #
	idx_good = []
	for i in np.arange(0, len(binned_mag_err)):
		if binned_mag_err[i] is not None:
			idx_good.append(i)
	mag_bins_for_mag_err, binned_mag_err, vax_mag_bin_medians= np.array(mag_bins_for_mag_err)[idx_good], np.array(binned_mag_err)[idx_good], np.array(vax_mag_bin_medians)[idx_good]

	# Examine objects within the relevant bin bounds #
	__idx_relevant = np.where((hax_mag >= min(mag_bins_for_mag_err)) & (hax_mag < max(mag_bins_for_mag_err)))[0]
	hax_mag, mag_diff = np.array(hax_mag)[__idx_relevant], np.array(mag_diff)[__idx_relevant]


	if CENTER_ERR_ABT_ZERO:
		for b in np.arange(0, len(mag_bins_for_mag_err)-1):
			if binned_mag_err[b] is not None: 
				for i in np.arange(0, len(hax_mag)):
					if hax_mag[i] >= mag_bins_for_mag_err[b] and hax_mag[i] < mag_bins_for_mag_err[b+1]:
						counter_objs_considered += 1
						if abs(mag_diff[i]) < binned_mag_err[b]:
							__num_1sigma_mag += 1


	# Center normalization about the median (of vertical axis) of each bin #
	if CENTER_ERR_ABT_ZERO is False:
		print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(mag_bins_for_mag_err)-1):
                        if binned_mag_err[b] is not None:
				print ' Considering objects with magnitudes (on the horizontal axis) in [', mag_bins_for_mag_err[b], ',', mag_bins_for_mag_err[b+1], ')'
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= mag_bins_for_mag_err[b] and hax_mag[i] < mag_bins_for_mag_err[b+1]:
						counter_objs_considered += 1
						if mag_diff[i] < binned_mag_err[b] + vax_mag_bin_medians[b] and mag_diff[i] >= -1.0*binned_mag_err[b] + vax_mag_bin_medians[b]: 
							__num_1sigma_mag += 1

	print ' NOT Normalized '
	print ' Total objects: ', tot
	print ' Number of objects considered:', counter_objs_considered
	print ' Number of objects within 1sigma_mag: ', __num_1sigma_mag 
	print ' Fraction within 1sigma_mag: ', float(num_1sigma_mag)/counter_objs_considered, '\n'

	return __num_1sigma_mag 




#TODO better name for this function, bins-->mag_bins, which of these are lists of lists?
def one_sigma_magnitude_counter_for_normalized_magnitude_difference(norm_mag_diff, clean_mag1, mag_bins_for_mag_err, hax_mag, mag_err_bin_medians, norm_vax_mag_bin_medians):
	"""Find the number of objects within 1sigma_mag using the normalized (to error) DeltaMagnitude. This function is only called if `NORMALIZE` is True.

	Parameters
	----------
	norm_mag_diff (list of floats)
		Normalized DeltaMagnitude. #FIXME what about colors? 

	clean_mag1 (list of floats)
		List of magnitudes with objects with flags* removed.

	mag_err_bin_medians (list) confirmed Jul 3 FIXME

	mag_bins_for_mag_err (list of floats)
		Magnitude bins

        Returns
	-------
        __num_1sigma_mag_norm (int)
		Number of objects within 1sigma_mag curve.
	"""

	tot = len(norm_mag_diff)
	__num_1sigma_mag_norm, counter_objs_considered = 0, 0


        # Cutoffs were introduced in error calculation. Consider only points not cutoff #
        # Get rid of `None` placeholders #
        idx_good = []
        for i in np.arange(0, len(mag_err_bin_medians)):
                if mag_err_bin_medians[i] is not None:
                        idx_good.append(i)
        # Binned values #
        mag_err_bin_medians, mag_bins_for_mag_err = np.array(mag_err_bin_medians)[idx_good], np.array(mag_bins_for_mag_err)[idx_good]
	norm_vax_median = np.array(norm_vax_mag_bin_medians)[idx_good]

        # Examine objects within the relevant bin bounds #
        __idx_relevant = np.where((hax_mag >= min(mag_bins_for_mag_err)) & (hax_mag < max(mag_bins_for_mag_err)))[0]
        hax_mag, norm_mag_diff = np.array(hax_mag)[__idx_relevant], np.array(norm_mag_diff)[__idx_relevant]


	if CENTER_ERR_ABT_ZERO:
                for k in norm_mag_diff:
			counter_objs_considered += 1
                        if abs(k) < 1.0:
                                __num_1sigma_mag_norm += 1
		#counter_objs_considered = len(norm_mag_diff)


	if CENTER_ERR_ABT_ZERO is False:
		print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(mag_bins_for_mag_err)-1):
			if mag_err_bin_medians[b] is not None:
				print ' Considering objects with magnitudes (on the horizontal axis) in [', mag_bins_for_mag_err[b], ',', mag_bins_for_mag_err[b+1], ')'
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= mag_bins_for_mag_err[b] and hax_mag[i] < mag_bins_for_mag_err[b+1]:
                                                counter_objs_considered += 1
                                                if norm_mag_diff[i] < 1 + norm_vax_median[b] and norm_mag_diff[i] >= -1.0 + norm_vax_median[b]: 
							__num_1sigma_mag_norm += 1

	print ' Normalized '
	print ' Total objects: ', tot
        print ' Number of objects considered:', counter_objs_considered
        print ' Number of objects within 1sigma_mag: ', __num_1sigma_mag_norm 
        print ' Fraction within 1sigma_mag: ', float(__num_1sigma_mag_norm)/counter_objs_considered, '\n'

	return __num_1sigma_mag_norm 




#TODO bins-->mag_bins, error-->mag_err?
def main_logger(vax_mag_bin_medians, tile, band, realization, clean_mag1, full_mag1, mag_bins_for_mag_err, hax_mag, fn_main_log, mag_err_bin_medians, vax_mag):
	"""Write to various (but not all) log files. Records the number of objects matched, the number of objects flagged, the number of objects within 1sigma_mag, etc.

	Parameters
	----------
	vax_mag (list of floats)
		Is normalized to 1sigma_mag if `NORMALIZE` is True.

	vax_mag_bin_medians
		Is normalized to 1sigma_mag if `NORMALIZE` is True.

	fn_main_log (file descriptor)

	clean_mag1 (list of floats)
		Magnitude of catalog1 in matched (join=1and2) catalog with ~flagged~ objects removed.

	full_mag1 (list of floats)
		Magnitude read directly from pandas Dataframe of matched (join=1and2) catalog.

	bins

	hax_mag

	error #TODO is this binned? List of lists?

        band (str) 

        realization (str) 

	tile

        Returns
	-------
        __percent_in_1sigma_mag (float)
		Percent of objects within 1sigma_mag.

	percentRecoveredFlagsIncluded (float)
		Percent of Balrog-injected objects recovered, calculated including objects with flags* (see README.md for the definition of flags*). Is `None` if neither `MATCH_CAT1` nor `MATCH_CAT2` is a truth catalog (`gal_truth` or `star_truth`).

	percentRecoveredFlagsRemoved (float)
		Percent of Balrog-injected objects recovered, calculated without including objects with flags* in the matched (via join=1and2) catalog nor the truth catalog. (see README.md for the definition of flags*). Is `None` if neither `MATCH_CAT1` nor `MATCH_CAT2` is a truth catalog (`gal_truth` or `star_truth`).
	"""

	if NORMALIZE:
		num_1sig = one_sigma_magnitude_counter_for_normalized_magnitude_difference(norm_mag_diff=vax_mag, clean_mag1=clean_mag1, mag_bins_for_mag_err=mag_bins_for_mag_err, hax_mag=hax_mag, mag_err_bin_medians=mag_err_bin_medians, norm_vax_mag_bin_medians=vax_mag_bin_medians)
	if NORMALIZE is False:
		num_1sig = one_sigma_magnitude_counter(mag_diff=vax_mag, clean_mag1=clean_mag1, mag_bins_for_mag_err=mag_bins_for_mag_err, hax_mag=hax_mag, binned_mag_err=mag_err_bin_medians, vax_mag_bin_medians=vax_mag_bin_medians) #was nsmedian 

	num_flags = len(full_mag1)-len(clean_mag1)


	__percent_in_1sigma_mag = 100.0*num_1sig/len(clean_mag1)


	# Write to log #
	if 'truth' in MATCH_CAT1:
		# Note that len(full_mag1) = len(full_mag2) & len(clean_mag1) = len(clean_mag2) so it does not matter which is passed to get_percent_recovered() #
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved, numTruthFlags = get_percent_recovered(full_data=full_mag1, clean_data=clean_mag1, inj_percent=INJ1_PERCENT, band=band, tile=tile, realization=realization)	
	if 'truth' in MATCH_CAT2:
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved, numTruthFlags = get_percent_recovered(full_data=full_mag1, clean_data=clean_mag1, inj_percent=INJ2_PERCENT, band=band, tile=tile, realization=realization)

	if 'truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2:
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = None, None


	### Main log ###
	# Append to csv #
        with open(fn_main_log, 'a') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
		# 'TILE', 'REALIZATION', 'BAND', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJECTS_IN_MATCH', 'TOTAL_FLAGGED_OBJECTS_IN_TRUTH', 'TOTAL_1SIGMA_MAG_MEAS', 'PERCENT_RECOVERED_FLAGS_INCLUDED', 'PERCENT_RECOVERED_FLAGS_REMOVED', 'RUN_TYPE' #
		writer.writerow([str(tile), str(realization), str(band), str(len(full_mag1)), str(num_flags), str(numTruthFlags), str(num_1sig), str(percentRecoveredFlagsIncluded), str(percentRecoveredFlagsRemoved), str(RUN_TYPE)])


	### Main log for stacked catalog ###
	if STACK_TILES or STACK_REALIZATIONS:
		# Append to csv #	
		with open(FN_STACK_MAIN_LOG, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			# 'TILE', 'REALIZATION', 'BAND', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJECTS_IN_MATCH', 'TOTAL_FLAGGED_OBJECTS_IN_TRUTH', 'TOTAL_1SIGMA_MAG_MEAS', 'PERCENT_RECOVERED_FLAGS_INCLUDED', 'PERCENT_RECOVERED_FLAGS_REMOVED' #
			writer.writerow([str(tile), str(realization), str(band), str(len(full_mag1)), str(num_flags), str(numTruthFlags), str(num_1sig), str(percentRecoveredFlagsIncluded), str(percentRecoveredFlagsRemoved)])


        return __percent_in_1sigma_mag, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved 




def get_colorbar_for_magnitude_plot_properties(df_1and2, cm_t_hdr, cm_t_err_hdr, idx_good, clean_mag1, clean_mag2, meas_or_true_cat, inj_percent, inj):
	"""Get data that will be used for the colorbar of the plot. This function will return `None` if no colorbar is to be added to the plot.

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.
	
	cm_t_hdr, cm_t_err_hdr (str)
		Matched (join=1and2) catalog headers for cm_T (size squared) and cm_T_err. Can be `None`.

	idx_good (list of ints)
		Indices of objects without flags.		
		
	clean_mag1, clean_mag2 (list of floats)
		Magnitudes ..FIXME rep

	meas_or_true_cat (str)
		Allowed values: 'true' 'meas'	

	Returns
	-------
	cbarData (list of floats)
		Data used to produce colorbar. Can be `None` if not colorbar is to be plotted.

	__cbar_label (str)
		Label for colorbar. Includes LaTeX \bf{} formatting. Can be `None`.
	"""

	if 'true' in meas_or_true_cat:
                sys.exit('ERROR. Colorbars should describe measured catalog values, not truth catalog values.')


	### Colorbar label ###
	# Prefix to labels #
        if inj:
		pref = str(inj_percent) + '%_inj_' #FIXME do this often
	if inj is False:
		if 'y3_gold' not in match_cat:
                        pref = 'base_'
		if 'y3_gold' in match_cat:
                        pref = 'Y3_'

        if CM_T_COLORBAR: __cbar_label = pref + cm_t_hdr[:-2] + '_' + str(meas_or_true_cat)
        if CM_T_ERR_COLORBAR: __cbar_label = pref + cm_t_err_hdr + '_' + str(meas_or_true_cat)


	### Colorbar value ###
        if CM_T_ERR_COLORBAR:
                # For measured catalog, cuts performed on truth catalogs #
                cbarData = get_good_data(df_1and2=df_1and2, hdr=cm_t_err_hdr, idx_good=idx_good, str_of_arr=False, band=None)

        if CM_T_COLORBAR:
                cbarData = get_good_data(df_1and2=df_1and2, hdr=cm_t_hdr, idx_good=idx_good, str_of_arr=False, band=None)


	return cbarData, __cbar_label




def get_magnitude_error(mag_err_hdr, flux_hdr, flux_cov_hdr, df_1and2, band, idx_good, match_cat):
	"""Get magnitude error from the measured catalog. The magnitude error in truth catalogs are not considered.

	Parameters
	----------
	mag_err_hdr (str)
		Matched catalog header for the magnitude error. Can be `None` if error is to be calculated using flux and flux covariance matrix.

	flux_hdr (str)
		Matched catalog header refering to the flux. FIXME is flux_g, flux_r ?? Can be `None` (??). Used if mag_err_hdr is `None`.

	cov_hdr (str)
		Matched catalog header refering to the flux covariance matrix. Can be `None` (??). Used if mag_err_hdr is `None`.

	df_1and2 (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.

	band (str)

	idx_good (list of ints)
		Indices of objects without flags.

	match_cat (str)
		Catalog containing the data. set by `MATCH_CAT1` or `MATCH_CAT2`.

	Returns
	-------
	error (list of floats)
		Error in magnitude. Will be `None` if `PLOT_MAG_ERR` is False.
	"""

        if PLOT_MAG_ERR and 'truth' not in match_cat:
                if mag_err_hdr is None: 
                        __mag_error = calculate_measured_magnitude_error_from_flux(df_1and2=df_1and2, flux_hdr=flux_hdr, flux_cov_hdr=flux_cov_hdr, band=band, idx_good=idx_good, match_cat=match_cat)
                if mag_err_hdr is not None:
			if match_cat == 'coadd':
				__mag_error = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_err_hdr, band=band)
			if match_cat == 'star_truth' or 'y3_gold' in match_cat:
                                __mag_error = df_1and2[str(mag_err_hdr[:-2]) + '_' + band.upper() + str(mag_err_hdr[-2:])]
			# Pass good indices #
			__mag_error = np.array(__mag_error)[idx_good]

        if PLOT_MAG_ERR is False or 'truth' in match_cat:
		__mag_error = None

	return __mag_error 




#TODO this function is not finished
def get_color_plot_error(mag_err_hdr, flux_hdr, flux_cov_hdr, df, band, idx_good, match_cat):
	"""Compute the color error. See @get_magnitude_plot_error docstring """


	magErrBand1 = get_magnitude_error(mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_cov_hdr=flux_cov_hdr, df_1and2=df, band=band, idx_good=idx_good, match_cat=match_cat)

	magErrBand2 = get_magnitude_error(mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_cov_hdr=flux_cov_hdr, df_1and2=df, band=BANDS_FOR_COLOR[band], idx_good=idx_good, match_cat=match_cat)

	'''
	if mag_err_hdr is None:
		magErrorFollow = calculate_measured_magnitude_error_from_flux((df=df, flux_hdr=flux_hdr, flux_cov_hdr=flux_cov_hdr, band=BANDS_FOR_COLOR[band], idx_good=idx_good, match_cat=match_cat)
	if mag_err_hdr is None:
		if match_cat == 'coadd':
			magErrorFollow = get_floats_from_string(df=df, four_elmt_arrs_hdr=mag_err_hdr, band=BANDS_FOR_COLOR[band])
			if match_cat == 'star_truth' or 'y3_gold' in match_cat:
				magErrorFollow = df[str(mag_err_hdr[:-2]) + '_' + BANDS_FOR_COLOR[band].upper() + str(mag_err_hdr[-2:])]

	__color_error = (np.array(magErr)**2 + np.array(magErrorFollow)**2)**0.5
	'''

	__color_error = (np.array(magErrBand1)**2 + np.array(magErrBand2)**2)**0.5

	return __color_error




def get_good_data(df_1and2, hdr, idx_good, str_of_arr, band):
	"""Get the data corresponding to indices of objects without flags* or indices of objects that satisfy quality cuts (if `EH_CUTS=True`).

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.
	
	hdr (str)
		Matched (via join=1and2) catalog header for the desired column of the DataFrame. 

	idx_good (list of ints)
		Indices of objects without flags* or indices of objects that satistfy quality cuts (if `EH_CUTS=True`).
 
	str_of_arr (bool)
		Set to `True` if the desired column of the DataFrame is a string of form '(data_g, data_r, data_i, data_z)'.

	band (str) 
		Can be `None`. Used if `magnitude=True`.

	Returns
	-------
	clean_data (list of floats)
		Contains column of the DataFrame with objects with flags* removed (or objects that do not pass quality cuts removed if `EH_CUTS=True` .
	"""

	if str_of_arr:
		fullData = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=hdr, band=band)
	if str_of_arr is False:
                fullData = df_1and2[hdr]

	__clean_data = np.array(fullData)[idx_good]

	return __clean_data




#TODO clean_mag1_band1 --> clean_mag1_band1, clean_mag2_band1 "band2
def get_color_from_binned_magnitude(df_1and2, mag_hdr1, mag_hdr2, clean_mag1_band1, clean_mag2_band1, band, idx_good):
	"""Get color with magnitudes binned via [21-22), [22-23), [23-24), [24-25). TODO mag_true 

        Parameters
        ----------
        df_1and2 (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.

        mag_hdr1, mag_hdr2 (str) #FIXME --> mag_hdr
		Matched (via join=1and2) catalog headers for magnitude for `MATCH_CAT1`, `MATCH_CAT2`.

        clean_magnitude_a (list of floats)
                List of magnitudes with objects with flags* removed. #TODO quality cuts may be used 

        band (str)

        idx_good (list of ints)
                Indices of objects without flags. #FIXME point to flags used somewhere in README.md.

        Returns
        -------
	"""

	# From `MATCH_CAT1` #
	__clean_mag1_band2 = np.array(get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=BANDS_FOR_COLOR[band]))[idx_good]	
	# From `MATCH_CAT2` #
	__clean_mag2_band2 = np.array(get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=BANDS_FOR_COLOR[band]))[idx_good] 

	__color1 = np.array(clean_mag1_band1) - np.array(__clean_mag1_band2)
	__color2 = np.array(clean_mag2_band1) - np.array(__clean_mag2_band2)

	# Indices in each bin #
	__idx1, __idx2, __idx3, __idx4 = [], [], [], []

	__mag_bins_for_color = [20, 21, 22, 23, 24]

	# Bin magnitude according to `clean_mag1_band1` TODO #
	for i in np.arange(0, len(clean_mag1_band1)):
		if clean_mag1_band1[i] >= __mag_bins_for_color[0] and clean_mag1_band1[i] < __mag_bins_for_color[1]:
			__idx1.append(i) 
		if clean_mag1_band1[i] >= __mag_bins_for_color[1] and clean_mag1_band1[i] < __mag_bins_for_color[2]:
			__idx2.append(i)
		if clean_mag1_band1[i] >= __mag_bins_for_color[2] and clean_mag1_band1[i] < __mag_bins_for_color[3]:
			__idx3.append(i)
		if clean_mag1_band1[i] >= __mag_bins_for_color[3] and clean_mag1_band1[i] < __mag_bins_for_color[4]:
			__idx4.append(i)

	__binned_color1, __binned_color2 = [], []
	# Append. Want list of lists #
	__binned_color1.append(__color1[__idx1]); __binned_color1.append(__color1[__idx2]); __binned_color1.append(__color1[__idx3]); __binned_color1.append(__color1[__idx4])
	__binned_color2.append(__color2[__idx1]); __binned_color2.append(__color2[__idx2]); __binned_color2.append(__color2[__idx3]); __binned_color2.append(__color2[__idx4])

	return __binned_color1, __binned_color2, __mag_bins_for_color 




def get_magnitude_axlabel(inj, mag_hdr, meas_or_true_cat, match_cat, band, inj_percent):
	"""Get labels for the horizontal axis. 
	
	Parameters
	----------
	inj_percent (int)
		
	inj (bool)
	
	mag_hdr, cm_t_hdr, cm_t_err (str) 
		Headers in the matched catalog that refer to the magnitude, cm_T (size), and cm_T error.
	
	match_cat (str)
		Catalog containing the data the hax_label describes. Set by `MATCH_CAT1` or `MATCH_CAT2`. 

	band (str)

	Returns
	-------
		hax_label (str) -- Label for the horizontal axis. Includes LaTeX \bf{} formatting. 
	"""

	### Prefix to labels ###
	if inj and inj_percent == 10:
		__pref = '10%_inj_'
	if inj and inj_percent == 20:
		__pref = '20%_inj_'

	if inj is False:
		if 'y3_gold' not in match_cat:
			__pref = 'base_'
		if 'y3_gold' in match_cat:
			__pref = 'Y3_'


	### Magnitude label ###
	# Transform 'mag_2' to 'mag_true' or 'mag_meas' #
	__mag_axlabel = str(mag_hdr[:-2]) + '_' + str(meas_or_true_cat)
	# 'mag_true' --> 'mag_{band}_true' with {band} bolded #
	__mag_axlabel = __mag_axlabel[:-4] + '$\\bf{' + str(band) + '}$_' + __mag_axlabel[-4:]	
	__mag_axlabel = __pref + __mag_axlabel 
        # Coadd catalogs. Combined to get '(m_g, m_r, m_i, m_z)' then matched. #
	if match_cat == 'coadd':
		__mag_axlabel = 'MAG_AUTO_' + '$\\bf{' + str(band) + '}$_' + str(meas_or_true_cat) 

	return __mag_axlabel



 
def get_color_axlabel(inj, meas_or_true_cat, match_cat, band, inj_percent):
	"""Get axes labels for color plots.

	Parameters
	----------
	meas_or_true_cat (str)
		Allowed values: 'meas' 'true'. Set by `AXLABEL1` or `AXLABEL2`.	
	
	Returns
	-------
	__color_axlabel (str)
		Contains LaTeX formatting. Ex: 'inj_(g-r)_true'.
	"""

	if inj:
		if inj_percent == 10:
			__pref = '10%_inj_'
		if inj_percent == 20:
			__pref = '20%_inj_'
        if inj is False:
		if 'y3_gold' not in match_cat:
                        __pref = 'base_'
		if 'y3_gold' in match_cat:
                        __pref = 'Y3_'

	if band == 'g': __color_axlabel = '$\\bf{(g-r)}$_' + str(meas_or_true_cat)
        if band == 'r': __color_axlabel = '$\\bf{(r-i)}$_' + str(meas_or_true_cat)
        if band == 'i': __color_axlabel = '$\\bf{(i-z)}$_' + str(meas_or_true_cat)
        __color_axlabel = __pref + __color_axlabel	

	return __color_axlabel 




def stacked_magnitude_completeness_subplotter(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, fn_flag_log, fn_mag_completeness_log):
	"""TODO

	Parameters
	----------

	Returns
	-------
	"""

	if STACK_TILES: print 'Stacking completeness for multiple tiles ...'
	if STACK_REALIZATIONS: print 'TODO' #print 'Stacking completeness for multiple realizations ...'

	### Initialize arrays with dimensions: (number_of_tiles, number_of_bands, number_of_magnitude_bins) ###
	# 1 --> 10% #
	__completeness1 = np.empty((len(ALL_TILES), len(ALL_BANDS), len(COMPLETENESS_PLOT_MAG_BINS)))
	# 2 --> 20% # __completeness2 = np.empty((len(ALL_TILES), len(ALL_BANDS), len(COMPLETENESS_PLOT_MAG_BINS)))


	# Fill first axis: tile #
	for i in np.arange(0, len(ALL_TILES)):
		t = ALL_TILES[i]


		# COSMOS `BALROG_RUN`s so far (Jun 21 2018) have only 20% injections #
                if 'COSMOS' not in BALROG_RUN:

			### Get DataFrame for tile with 10% injection via `inj_percent=10` ###
                        __df1 = get_dataframe(realization=realization, tile=t, inj1=True, inj2=True, inj1_percent=10, inj2_percent=10, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

			# Fill second axis: band #
			for j in np.arange(0, len(ALL_BANDS)):
				b = ALL_BANDS[j]

                                magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel= get_magnitude_plot_variables(band=b, df_1and2=__df1, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=t, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

                                completenessFraction1 = get_magnitude_completeness(inj_percent=10, idx_good=idxGood, df_1and2=__df1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=t, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]
			
				# Fill third axis: completeness (list) #
				__completeness1[i][j] = completenessFraction1


		### Get DataFrame for tile with 20% injection ###
		print 'Getting DataFrame for tile: ', t, ' ... '
		__df2 = get_dataframe(realization=realization, tile=t, inj1=True, inj2=True, inj1_percent=20, inj2_percent=20, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

		# Fill second axis: band #
		for j in np.arange(0, len(ALL_BANDS)):
			b = ALL_BANDS[j]

                        magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel= get_magnitude_plot_variables(band=b, df_1and2=__df2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=t, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

                        completenessFraction2 = get_magnitude_completeness(inj_percent=20, idx_good=idxGood, df_1and2=__df2, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=t, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]

			# Fill third axis: completeness (list) #
			__completeness2[i][j] = completenessFraction2


	# One completeness plot per band #
	for k in np.arange(0, len(ALL_BANDS)):

                plt.figure(figsize=(12, 10))
		
		# COSMOS runs have only 20% injections #
                if 'COSMOS' not in BALROG_RUN:
			# Make subplots #
                        plt.subplot(1, 2, 1)
                        plt.plot(COMPLETENESS_PLOT_MAG_BINS,np.nanmean(__completeness1[:,k,:], axis=0) , color='blue')
                        plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
                        plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
                        plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
                        plt.title('10% Injection') #TODO automate this?
                        plt.ylabel('Magnitude Completeness')
                        plt.xlabel('10%_inj_cm_mag_$\\bf{'+str(ALL_BANDS[k])+'}$_true') #TODO automate this?
                        plt.grid(linestyle='dotted')

                        plt.subplot(1, 2, 2)

                # Get completeness mean (amongst all tiles) for each magnitude bin --> `axis=0` #
                # Note: plot() ignores nans so there is no need to remove them from data in sep step #
                plt.plot(COMPLETENESS_PLOT_MAG_BINS, np.nanmean(__completeness2[:,k,:], axis=0), color='green')
		plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
                plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
                plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
                plt.title('20% Injection') #TODO automate this?
                plt.ylabel('Magnitude Completeness')
                plt.xlabel('20%_inj_cm_mag_$\\bf{'+str(ALL_BANDS[k])+'}$_true') #TODO automate this?
                plt.grid(linestyle='dotted')

                if '10% ' in plot_title:
                        plot_title = plot_title.replace('10% ', '')
                if '20% ' in plot_title:
                        plot_title = plot_title.replace('20% ', '')
                plt.suptitle(plot_title, fontweight='bold') #FIXME check if there is a sigma cutoff

		if SHOW_PLOT: plt.show()

	return 0




def magnitude_completeness_subplotter(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, tile, fn_flag_log, fn_mag_completeness_log):
	"""Creates two completeness plots with two panels of completeness for 10% injections (first subplot panel) and 20% injections (second subplot panel).
	This function is only to be called when comparing Balrog-injected catalogs and truth catalogs.

	Parameters
	----------
	
	Returns
	-------
	"""

	__completeness_griz1, __completeness_griz2 = [], []
	__bin_median_griz1, __bin_median_griz2 = [], []


	# COSMOS `BALROG_RUN`s so far (Jun 21 2018) have only 20% injections #
	if 'COSMOS' not in BALROG_RUN: 

		# Get completeness for 10% inj via `inj_percent=10` #
		#FIXME instances where data_frame depends on band:
		__df1 = get_dataframe(realization=realization, tile=tile, inj1=True, inj2=True, inj1_percent=10, inj2_percent=10, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

		for b in ALL_BANDS:
			magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel = get_magnitude_plot_variables(band=b, df_1and2=__df1, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

			completenessFraction1 = get_magnitude_completeness(inj_percent=10, idx_good=idxGood, df_1and2=__df1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]

			__completeness_griz1.append(completenessFraction1)


	# Get completeness for 20% injection via `inj_percent=20` #
	__df2 = get_dataframe(realization=realization, tile=tile, inj1=True, inj2=True, inj1_percent=20, inj2_percent=20, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

	for b in ALL_BANDS:
		magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel = get_magnitude_plot_variables(band=b, df_1and2=__df2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot)

		completenessFraction2 = get_magnitude_completeness(inj_percent=20, df_1and2=__df2, idx_good=idxGood, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]

		print 'Getting completeness for band: ', b, '\n'

		__completeness_griz2.append(completenessFraction2)


	print 'Plotting magnitude completeness...\n'


	for i in np.arange(0, len(ALL_BANDS)):

		plt.figure(figsize=(12, 10))

		# COSMOS runs have only 20% injections #
		if 'COSMOS' not in BALROG_RUN:
			plt.subplot(1, 2, 1)
			plt.plot(COMPLETENESS_PLOT_MAG_BINS, __completeness_griz1[i], color='blue')
			plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
			plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
			plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
			plt.title('10% Injection') #TODO automate this?
			plt.ylabel('Magnitude Completeness')
			plt.xlabel('10%_inj_cm_mag_$\\bf{'+str(ALL_BANDS[i])+'}$_true') #TODO automate this?
			plt.grid(linestyle='dotted')

			plt.subplot(1, 2, 2)

		plt.plot(COMPLETENESS_PLOT_MAG_BINS, __completeness_griz2[i], color='green')
		plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
		plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
		plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
		plt.title('20% Injection') #TODO automate this?
		plt.ylabel('Magnitude Completeness')
		plt.xlabel('20%_inj_cm_mag_$\\bf{'+str(ALL_BANDS[i])+'}$_true') #TODO automate this?
		plt.grid(linestyle='dotted')

		if '10% ' in plot_title:
			plot_title = plot_title.replace('10% ', '')
		if '20% ' in plot_title:
			plot_title = plot_title.replace('20% ', '')
		plt.suptitle(plot_title, fontweight='bold') #FIXME check if there is a sigma cutoff

		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', ALL_BANDS[i]); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

		if SHOW_PLOT: plt.show()


	return 0




#def get_color_completeness(truth_mag, match_mag, band, realization, tile, mag_hdr1, mag_hdr2):
	'''
	Parameters
	----------

	Returns
	-------

	__mag_completeness, __plot_bins = [], []

	__step = 0.25 #FIXME

	# Get truth_color and truth_match_color #
	# Truth catalog #
	if band == 'g': idx1, idx2 = 0, 1
	if band == 'r': idx1, idx2 = 1, 2
	if band == 'i': idx1, idx2 = 2, 3


	### Read truth catalog which has magnitudes in form (m_g, m_r, m_i, m_z) ###
        if 'truth' in MATCH_CAT1:
                fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, realization=realization, tile=tile, band=band, inj=INJ1)
                mag_hdr = mag_hdr1

        if 'truth' in MATCH_CAT2:
                fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT, realization=realization, tile=tile, band=band, inj=INJ2)
                mag_hdr = mag_hdr2

        hdu = fits.open(fn_truth_cat)
        data = hdu[1].data
        truth_mag_griz = data[mag_hdr[:-2]]
        truth_color = []
        for mag_griz in truth_mag_griz:
                truth_color.append(mag_griz[idx1] - mag_griz[idx2])
        truth_color = np.array(truth_color)


	### Get color from `match_mag` and ensure that `match_mag` refers to the truth catalog value #
	sys.exit('how to bin colors for completeness plot?')

	# Match. Note that this introduces a bin cutoff of m=24 that truth_color does not #
	colorBins1, colorBins2, magBins = get_color_from_binned_magnitude(df, hdr1, hdr2, clean_mag1_band1, clean_mag2_band1, band=band, idx_good)
	# Flatten lists #
	colorBins1 = [item for sublist in colorBins1 for item in sublist]; colorBins2 = [item for sublist in colorBins2 for item in sublist]



	__bins = np.arange(int(np.min(truth_color)), int(np.max(truth_color)), __step)


        for b in np.arange(0, len(__bins)-1):

                # Count completeness in each bin #
                __counter_truth, __counter_match = 0, 0

                for i in np.arange(0, len(truth_color)):
                        if truth_color[i] >= __bins[b] and truth_color[i] < __bins[b+1]:
                                __counter_truth += 1

                for j in np.arange(0, len(__truth_mag_in_match1and2)):
                        if truth_match_color[j] >= __bins[b] and truth_match_color[j] < __bins[b+1]:
                                __counter_match += 1

                if __counter_truth != 0:
                        __mag_completeness.append((1.0*__counter_match)/__counter_truth)
                        __plot_bins.append(np.median([__bins[b], __bins[b+1]]))
	return __mag_completeness, __plot_bins




def color_completeness_subplotter(band, df, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title, realization, tile, fd_flag_log, fd_mag_completeness_log, fn_plot):


	### Force order ###


        err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, haxLabel2, vaxLabel = get_magnitude_plot_variables(band=band, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag_log=fd_flag_log, plot_title=plot_title, fn_plot=fn_plot)

	# Photometry cuts applied to `truthMag`, `matchMag` #
        truthMag1, matchMag1 = get_magnitude_completeness(df=df, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=band, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, error1=err1, error2=err2, fd_mag_completeness_log=fd_mag_completeness_log)[2:]


	colorCompleteness1, bins1 = get_color_completeness(truth_mag=truthMag1, match_mag=matchMag1, band=band, realization=realization, tile=tile, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2)



        err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, haxLabel2, vaxLabel = get_magnitude_plot_variables(band=band, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag_log=fd_flag_log, plot_title=plot_title, fn_plot=fn_plot)

        truthMag2, matchMag2 = get_magnitude_completeness(df=df, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=band, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, error1=err1, error2=err2, fd_mag_completeness_log=fd_mag_completeness_log)[2:]


	writeColor = WRITE_COLORS[band]
	if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', writeColor); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot) 

	return 0
'''



#FIXME df_1and2 not used
def get_magnitude_completeness(idx_good, df_1and2, clean_mag1, clean_mag2, full_mag1, full_mag2, tile, realization, band, mag_hdr1, mag_hdr2, mag_err1, mag_err2, fn_mag_completeness_log, inj_percent): #FIXME inj?
	"""TODO after function is done
	Introduce cutoff for bin size.
	"""

	sys.exit('df_1and2 param is not being used....is right df still used?')

	GRID_PLOT = False
	if GRID_PLOT:
		print 'TODO plot_populated_tile_grid()'

	mag_bin_l = 21

	if inj_percent == 10: 
		print 'Calculating magnitude completeness of 10% injected catalogs ... '
		__write_inj_percent = '10%'
        if inj_percent == 20: 
		print 'Calculating magnitude completeness of 20% injected catalogs ... '
		__write_inj_percent = '20%'

	__mag_completeness, __plot_bins = [], []


	### Read truth catalog which has magnitudes in form (m_g, m_r, m_i, m_z) ###
	if 'truth' in MATCH_CAT1:
		fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT1, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ1)	
		mag_hdr = mag_hdr1
		# Use error in measured catalogs only #
		#TODO use None?
		mag_err1 = np.zeros(len(mag_err2))
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
		flag_hdr = FLAGS_HDR1[:-2]
		cm_flag_hdr = CM_FLAGS_HDR1[:-2] 

	if 'truth' in MATCH_CAT2:
		fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT2, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ2)
		mag_hdr = mag_hdr2
		ra_hdr, dec_hdr = RA_HDR2, DEC_HDR2
		# Use error in measured catalogs only #
		mag_err2 = np.zeros(len(mag_err1))
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
                flag_hdr = FLAGS_HDR2[:-2]
                cm_flag_hdr = CM_FLAGS_HDR2[:-2]

	hdu = fits.open(fn_truth_cat)
	data = hdu[1].data
	truth_mag_griz = data[mag_hdr[:-2]]
	truth_mag = []

	for mag_griz in truth_mag_griz:
		__truth_mag.append(mag_griz[BAND_INDEX[band]])
	__truth_mag = np.array(__truth_mag)

	flag = data[flag_hdr]
	cm_flag = data[cm_flag_hdr]

	# Primary flag cuts to truth catalog #
	__idx_good = np.where( (abs(__truth_mag) != 9999.0) & (abs(__truth_mag) != 99.0) & (abs(__truth_mag) != 37.5) & (flag == 0) & (cm_flag == 0) )[0]
	
	__truth_mag = __truth_mag[__idx_good]


	### Photometry cuts on matched catalog (|DeltaM| < 3sigma) where sigma refers to measured catalog only ###
        match_mag1, match_mag2 = [], []

	__cutoff_mag_diff = 2
	print 'Using only points with |DeltaMagnitude| <', __cutoff_mag_diff, ' ... '
	
        for k in np.arange(0, len(clean_mag1)):
                #if abs(clean_mag1[k] - clean_mag2[k]) < 3.0 * (mag_err1[k]**2 + mag_err2[k]**2)**0.5:
		if abs(clean_mag1[k] - clean_mag2[k]) < __cutoff_mag_diff:
                        match_mag1.append(clean_mag1[k])
                        match_mag2.append(clean_mag2[k])

	if 'truth' in MATCH_CAT1:
                __truth_mag_in_match1and2 = match_mag1
		other_match_mag = match_mag2
	if 'truth' in MATCH_CAT2:
                __truth_mag_in_match1and2 = match_mag2
		other_match_mag = match_mag1

	__bin_size_cutoff = 10

	__mag_completeness = []
	numer = 1.0*np.histogram(__truth_mag_in_match1and2, COMPLETENESS_MAG_BINS)[0]
	denom = 1.0*np.histogram(__truth_mag, COMPLETENESS_MAG_BINS)[0]
	# Bin size cutoff: replace small bins with `nan` #
	numer[numer < __bin_size_cutoff] = np.nan
	denom[denom < __bin_size_cutoff] = np.nan

	__mag_completeness = numer/denom

	print ' Including bins with more than ', __bin_size_cutoff, ' objects ... \n'

	# TILE, REALIZATION, BAND, INJ_PERCENT, TRUTH_MAG_BIN_L, TRUTH_MAG_BIN_R, MATCH_CAT_OBJS_IN_BIN, TRUTH_CAT_OBJS_IN_BIN #
	for i in np.arange(0, len(COMPLETENESS_MAG_BINS)-1):
		if numer[i] != 0 and denom[i] != 0: #TODO better way to eliminate the large range of COMPLETENESS_MAG_BINS
			# Append to csv #
			with open(fn_mag_completeness_log, 'a') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow([str(tile), str(realization), str(band), str(inj_percent), str(COMPLETENESS_MAG_BINS[i]), str(COMPLETENESS_MAG_BINS[i+1]), str(numer[i]), str(denom[i])])


	return __mag_completeness, __truth_mag, __truth_mag_in_match1and2 




def gaussian(x, mu, sig):
	"""Normalized Gaussian function. Normalization as defined here means that the area under the curve is 1.

	Parameters
	----------
	x (array)
		Range over which to compute the normalized Gaussian function.
	
	mu (float)
		Mean of the distribution.

	sig (float)
		Standard deviation of the distrubution.

	Returns
	-------
	__gaussian (array)
		Contains points that lie on the specified Gaussian distribution. 	
	"""

	__gaussian = 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2.0)
	return __gaussian 




def measure_flux_using_gaussian_aperture(fn_gauss_aper_log, tile, realization, cm_tdbyte_meas, cm_tdbyte_true, idx_good, band, cm_g_meas, cm_g_1_true, cm_g_2_true, cm_flux_meas, cm_flux_true, fracdev_meas, fracdev_true, box_size_meas, box_size_true, cm_t_meas, cm_t_true, flux_cov_meas, t_thresh=0.001, psf_t=0.05, aper_t=2.25, vb=False, N='all', pixscale=0.263):
	"""Measure flux using a Gaussian aperture. Reference: Spencer Everett.

	Parameters
	----------
	idx_good (list of ints)
		Contains indices without flags*.

	fn_gauss_aper_log (str)
		Complete file name to log the results of this function.

	cm_tdbyte_meas, cm_tdbyte_true (array-like)
		'cm_TdByTe' for the measured and truth catalogs. 'cm_TdByTe' is the ratio of sizes Tdev/Texp.

	cm_t_meas, cm_t_true (array-like)
		'cm_T' for the measured and truth catalog. 'cm_T' is the size squared of the object. Indices with flags* have NOT been removed from this array.

	cm_fracdev_meas, cm_fracdev_true (array-like)
		Composite model (CM) fraction De Vaucouleurs to exponential. TODO verify this. flags* have NOT been removed from this array.

	flux_cov_meas (array-like)
		Flux covariance diagonals corresponding to `band` for the measured catalog. flags* have NOT been removed from this array. 

	box_size_meas, box_size_true (array-like)
		Size of the postage stamp for the measured and truth catalog. flags* have NOT been removed from this array. 

	cm_flux_meas, cm_flux_true (array-like)
		Flux for the measured and truth catalogs. flags* have NOT been removed from this array. 

	cm_g_meas (array-like with nested arrays)
		Shape components for the measured catalog. Each element is array-like and of form `(cm_g_1, cm_g_2)`. flags* have NOT been removed from this array.

	cm_g_1_true, cm_g_2_true (array-like)
		Shape components 1 and 2 for the truth catalog. flags* have NOT been removed from this array.

	band (str)

	tile (str)

        realization (str)	

	Returns
	-------
	__clean_gauss_aper_norm_flux_diff (list of floats)
		Contains DeltaFlux/SigmaFlux (where DeltaFlux is computed via FluxMeasured-FluxTrue) for objects without flags* and without aper_flags. aper_flags correspond to any or all of the following: 'cm_TdByTe' is 0 for the truth catalog or the measured catalog, `gmix.GMixCM()` fails, or `gm_true.convolve()` fails. 
	"""
	print type(cm_g_meas[0]), cm_g_meas[0]; sys.exit('cm_g_meas[0] chk')
	
	print 'Measuring flux using Gaussian aperture for band:', band, '...'

	all_gauss_aper_flux_chi, index, flag, all_gauss_aper_flux_true, all_gauss_aper_flux_meas, all_gauss_aper_flux_diff, flux_per_err = [], [], [], [], [], [], []

	# Question: use cm_t from true? #
	if N == 'all':
		idx_good_t = np.where(cm_t_true > t_thresh)[0]
	else:
		idx_good_t = random.sample(np.where(cm_t_true > t_thresh)[0], N)

	# Combine indices without flags* to indices which fit cm_T criteria #
	idx = []
	if len(idx_good) > len(idx_good_t): idx_compare1 = idx_good; idx_compare2 = idx_good_t
	if len(idx_good_t) > len(idx_good): idx_compare1 = idx_good_t; idx_compare2 = idx_good
	for c in idx_compare1:
		if c in idx_compare2:
			idx.append(c)

	print '{} objects passed cm_T clipping AND flag cuts ... '.format(len(idx))

	count = 0	
	count_none1, count_none2, count_none3 = 0, 0, 0

	for i in idx: 
		count += 1

		# Ratio of sizes Tdev/Texp #
		if cm_tdbyte_true[i] == 0 or cm_tdbyte_meas[i] == 0:
			count_none1 += 1
			index.append(i)
			flag.append(1)
			all_gauss_aper_flux_true.append(None)
			all_gauss_aper_flux_meas.append(None)
			all_gauss_aper_flux_diff.append(None)
			flux_per_err.append(None)
			all_gauss_aper_flux_chi.append(None)
			continue

		### Parameters for `gmix`: cen1, cen2, g1, g2, T, flux ###
		params_true = [0, 0, cm_g_1_true[i], cm_g_2_true[i], cm_t_true[i], cm_flux_true[i]]

		# Question: no band dep?
		# '(cm_g_1, cm_g_2)' --> cm_g_1 #
		cm_g_1_meas = float(cm_g_meas[i][1:cm_g_meas[i].index(',')])
		# '(cm_g_1, cm_g_2)' --> cm_g_2 #
		cm_g_2_meas = float(cm_g_meas[i][cm_g_meas[i].index(',')+2:-1])

		params_meas = [0, 0, cm_g_1_meas, cm_g_2_meas, cm_t_meas[i], cm_flux_meas[i]] 

		#print 'Performing `GMixCM` ...'
		gm_true = gmix.GMixCM(fracdev_true[i], cm_tdbyte_true[i], params_true)
		gm_meas = gmix.GMixCM(fracdev_meas[i], cm_tdbyte_meas[i], params_meas)

		# PSF = single Gaussian #
		params_psf = [0.0, 0.0, 0.0, 0.0, psf_t, 1.0]
		gm_psf = ngmix.GMixModel(params_psf, 'gauss')

		# Jacobians #
		dims_true = [box_size_true[i], box_size_true[i]]
		dims_meas = [box_size_meas[i], box_size_meas[i]]

		# Check that dimensions of measure and true are equal # 
		if dims_true != dims_meas:
			if vb: print 'Dimensions are not consistent. Forcing equality using max()'
			dims_true, dims_meas = max(dims_true, dims_meas), max(dims_true, dims_meas)
			if vb: print 'dims_true: {}\ndims_meas: {}'.format(dims_true, dims_meas)

		rc_true, cc_true = (dims_true[0]-1.0)/2.0, (dims_true[1]-1.0)/2.0
		jacob_true = ngmix.DiagonalJacobian(scale=pixscale, row=rc_true, col=cc_true)
		rc_meas, cc_meas = (dims_meas[0]-1.0)/2.0, (dims_meas[1]-1.0)/2.0
		jacob_meas = ngmix.DiagonalJacobian(scale=pixscale, row=rc_meas, col=cc_meas)

		# Convolution with PSF #
		try: 
			cgm_true = gm_true.convolve(gm_psf)
			cgm_meas = gm_meas.convolve(gm_psf)
		except:
			count_none3 += 1
			index.append(i)
                        flag.append(1)
                        all_gauss_aper_flux_true.append(None)
                        all_gauss_aper_flux_meas.append(None)
                        all_gauss_aper_flux_diff.append(None)
                        flux_per_err.append(None)
			all_gauss_aper_flux_chi.append(None)
                        continue			

		# Render simulated images #
		im_true = gm_true.make_image(dims_true, jacobian=jacob_true)
		im_meas = gm_meas.make_image(dims_meas, jacobian=jacob_meas)

		# Gaussian aperture #
		params_aper = [0.0, 0.0, 0.0, 0.0, aper_t, 1.0]
		aperture = ngmix.GMixModel(params_aper, 'gauss')

		im_weight_true = aperture.make_image(dims_true, jacobian=jacob_true)
		im_weight_true *= 1.0 / im_weight_true.max()
		im_weight_meas = aperture.make_image(dims_meas, jacobian=jacob_meas)
		im_weight_meas *= 1.0 / im_weight_meas.max()

		# Measure flux #
		gauss_aper_flux_true = (im_true * im_weight_true).sum()
		gauss_aper_flux_meas = (im_meas * im_weight_meas).sum()
		# meas minus true #
		gauss_aper_flux_diff = gauss_aper_flux_meas - gauss_aper_flux_true
		per_err = 100.0 * (gauss_aper_flux_diff / gauss_aper_flux_true)
		__flux_err = flux_cov_meas[i]**0.5

		if vb:
			print 'flux_true: {}\nflux_meas: {}\ndiff: {}'.format(gauss_aper_flux_true, gauss_aper_flux_meas, gauss_aper_diff)
			print 'Percent error: {:.2}%'.format(per_err)

		index.append(i)
		flag.append(0)
		all_gauss_aper_flux_true.append(gauss_aper_flux_true)
		all_gauss_aper_flux_meas.append(gauss_aper_flux_meas)
		all_gauss_aper_flux_diff.append(gauss_aper_flux_diff)
		flux_per_err.append(per_err)
		all_gauss_aper_flux_chi.append(gauss_aper_flux_diff/__flux_err)

		# Append to csv #
		with open(fn_gauss_aper_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			# TILE, REALIZATION, BAND, FLUX_MEAS, FLUX_TRUE, GAUSS_APER_FLUX_MEAS, GAUSS_APER_FLUX_TRUE, FLUX_ERR_MEAS, GAUSS_APER_FLUX_DIFF/FLUX_ERR_MEAS #
			writer.writerow([str(tile), str(realization), str(band), str(cm_flux_meas[i]), str(cm_flux_true[i]), str(gauss_aper_flux_meas), str(gauss_aper_flux_true), str(__flux_err), str(gauss_aper_flux_diff/__flux_err)])


	print 'Number of flags from cm_TdByTe==0:', count_none1
	print 'Number of flags from `gm_{true/meas}.convolve(gm_psf)`, likely with Error: gauss2d det too low:', count_none3

	# Get rid of flagged points #
	__clean_gauss_aper_norm_flux_diff = []
	__clean_gauss_aper_flux_diff = []
	for j in all_gauss_aper_flux_chi:
		if j is not None:
			__clean_gauss_aper_norm_flux_diff.append(j)
	for k in all_gauss_aper_flux_diff:
		if k is not None:
			__clean_gauss_aper_flux_diff.append(k)


	return __clean_gauss_aper_norm_flux_diff, __clean_gauss_aper_flux_diff 




def normalized_flux_histogram_subplotter(fn_gauss_aper_log, fn_plot, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title, tile, realization):
	"""Create 2x2 plot grid of normalized (area under the curve equals 1) 1D DeltaFlux/SigmaFlux histograms. One subplot is created for each griz band. Each plot includes a standard Gaussian (mean=0, standard_deviation=1) and the best-fit Gaussian to DeltaFlux/SigmaFlux. 

	Parameters
	----------
	fn_gauss_aper_log (str)
		Complete file name to log the results of `measure_flux_using_gaussian_aperture()`

	fn_plot (str)
		Complete file name for the plot. Used if `SAVE_PLOT=True`.

	
	
	Returns
	-------
	0
	"""


	### Write headers to csv outside of loop over bands ###
        with open(fn_gauss_aper_log, 'wb') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                # Write headers #
                writer.writerow(['TILE', 'REALIZATION', 'BAND', 'FLUX_MEAS', 'FLUX_TRUE', 'GAUSS_APER_FLUX_MEAS', 'GAUSS_APER_FLUX_TRUE', 'FLUX_ERR_MEAS', 'GAUSS_APER_FLUX_DIFF/FLUX_ERR_MEAS'])


	### Create 4-by-4 subplot ###
        counter_subplot = 1
        # Figure size units: inches #
        plt.figure(figsize=(12, 10))


	### Read columns from DataFrames ###
	if GAUSS_APER:
		print 'Using Gaussian aperture to measure flux...'
		# Columns that are not dependent on band # 
		if 'truth' in MATCH_CAT1: __flux_cov_hdr = CM_FLUX_COV_HDR2; suf_true = '_1'; suf_meas = '_2'
                if 'truth' in MATCH_CAT2: __flux_cov_hdr = CM_FLUX_COV_HDR1; suf_true = '_2'; suf_meas = '_1'
                __cm_tdbyte_meas, __cm_tdbyte_true = df_1and2['cm_TdByTe'+suf_meas], df_1and2['cm_TdByTe'+suf_true]
                # `cm_g_meas` is of form '(cm_g_1, cm_g_2)' #
                __cm_g_meas, __cm_g_1_true, __cm_g_2_true = df_1and2['cm_g_2a'], df_1and2['cm_g_1'+suf_true], df_1and2['cm_g_2'+suf_true]
                __fracdev_meas, __fracdev_true = df_1and2['cm_fracdev'+suf_meas], df_1and2['cm_fracdev'+suf_true]
                __box_size_meas, __box_size_true = df_1and2['box_size'+suf_meas], df_1and2['box_size'+suf_true]
                __cm_t_meas, __cm_t_true = df_1and2['cm_T'+suf_meas], df_1and2['cm_T'+suf_true]
	if GAUSS_APER is False:
		print 'Not using Gaussian aperture to measure flux...'
		# Values will not be used #
		__cm_tdbyte_meas, __cm_tdbyte_true, __cm_g_meas, __cm_g_1_true, __cm_g_2_true, __fracdev_meas, __fracdev_true, __box_size_meas, __box_size_true, __cm_t_meas, __cm_t_true = None, None, None, None, None, None, None, None, None, None

        ### Create one subplot for each griz band ###
        for b in ALL_BANDS:
		plt.subplot(2, 2, counter_subplot)

		# Read columns from DataFrame that are dependent on band #
		if GAUSS_APER:
			__cm_flux_meas, __cm_flux_true = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr='cm_flux'+suf_meas, band=b), get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr='cm_flux'+suf_true, band=b)
			__flux_cov_meas = get_matrix_diagonal_element(df=df_1and2, band=b, sq_matrices_hdr=__flux_cov_hdr)
		if GAUSS_APER is False:
			__cm_flux_meas, __cm_flux_true, __flux_cov_meas = None, None, None

		# Create one subplot #
		normalized_flux_histogram_plotter(tile=tile, realization=realization, band=b, df_1and2=df_1and2, cm_tdbyte_meas=__cm_tdbyte_meas, cm_tdbyte_true=__cm_tdbyte_true, flux_hdr1=flux_hdr1, flux_hdr2=flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=plot_title, cm_g_meas=__cm_g_meas, cm_g_1_true=__cm_g_1_true, cm_g_2_true=__cm_g_2_true, cm_flux_meas=__cm_flux_meas, cm_flux_true=__cm_flux_true, fracdev_meas=__fracdev_meas, fracdev_true=__fracdev_true, box_size_meas=__box_size_meas, box_size_true=__box_size_true, cm_t_meas=__cm_t_meas, cm_t_true=__cm_t_true, flux_cov_meas=__flux_cov_meas, fn_gauss_aper_log=fn_gauss_aper_log)
		counter_subplot += 1

	plt.suptitle(plot_title, fontweight='bold')
	plt.subplots_adjust(hspace=0.3)
	#plt.tight_layout(pad=1, h_pad=1, w_pad=1)

	if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

	if SHOW_PLOT: plt.show()


	return 0




def normalized_flux_histogram_plotter(fn_gauss_aper_log, tile, realization, band, df_1and2, cm_tdbyte_meas, cm_tdbyte_true, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title, cm_g_meas, cm_g_1_true, cm_g_2_true, cm_flux_meas, cm_flux_true, fracdev_meas, fracdev_true, box_size_meas, box_size_true, cm_t_meas, cm_t_true, flux_cov_meas):
	"""Create a normalized (area under curve equals 1) 1D histogram of DeltaFlux/SigmaFlux for a particular `band`. Plot includes a standard Gaussian (mean=0, standard_deviation=1) and the best-fit Gaussian to DeltaFlux/SigmaFlux.

	Parameters
	----------
	TODO repeats

	Returns
	-------
	0
	"""

	### Plot DeltaFlux/SigmaFlux as normalized histogram ###
	plotData, idxGood, fluxDiff = get_flux_plot_variables(band=band, df_1and2=df_1and2, flux_hdr1=flux_hdr1, flux_hdr2=flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)
	__flux_hist_bins = 10**2


	if GAUSS_APER:
		#plotData --> more descriptive because this is not arbitrary. normFluxDiff, fluxDiff = ...
		plotDataGauss, fluxDiffGaussAper = measure_flux_using_gaussian_aperture(fn_gauss_aper_log=fn_gauss_aper_log, tile=tile, realization=realization, idx_good=idxGood, band=band, cm_g_meas=cm_g_meas, cm_g_1_true=cm_g_1_true, cm_g_2_true=cm_g_2_true, cm_flux_meas=cm_flux_meas, cm_flux_true=cm_flux_true, fracdev_meas=fracdev_meas, fracdev_true=fracdev_true, box_size_meas=box_size_meas, box_size_true=box_size_true, cm_t_meas=cm_t_meas, cm_t_true=cm_t_true, flux_cov_meas=flux_cov_meas, cm_tdbyte_meas=cm_tdbyte_meas, cm_tdbyte_true=cm_tdbyte_true)
		gauss_aper_label = 'Gauss aper'

		if TRIM_FLUX:
			__val_2pg, __val_98pg = np.percentile(plotDataGauss, [2, 98], interpolation='nearest')
			plotDataGauss = np.sort(plotDataGauss)
			idx1, idx2 = plotDataGauss.tolist().index(__val_2pg), plotDataGauss.tolist().index(__val_98pg)
			plotDataGauss = plotDataGauss[idx1:idx2]
			gauss_aper_label = 'trim Gauss aper'	


	'''
	### For DeltaFlux plots (test) TODO remove ###
	p1, p2 = np.percentile(fluxDiff, [2, 98], interpolation='nearest')
	fluxDiff = np.sort(fluxDiff)
	idx1, idx2 = fluxDiff.tolist().index(p1), fluxDiff.tolist().index(p2)
	fluxDiff = fluxDiff[idx1:idx2]

	p3, p4 = np.percentile(fluxDiffGaussAper, [2, 98], interpolation='nearest')
	fluxDiffGaussAper = np.sort(fluxDiffGaussAper)
	idx3, idx4 = fluxDiffGaussAper.tolist().index(p3), fluxDiffGaussAper.tolist().index(p4) 
	fluxDiffGaussAper = fluxDiffGaussAper[idx3:idx4]

	plt.hist(fluxDiff, label='DeltaFlux', histtype='step', bins=10**2)
	plt.hist(fluxDiffGaussAper, label='DeltaFluxGauss', histtype='step', bins=10**2)
	plt.legend()

	plt.show()
	sys.exit()
	'''

	if TRIM_FLUX:
		# Returns value (not index) of 2nd and 98th percentiles #
		__val_2p, __val_98p = np.percentile(plotData, [2, 98], interpolation='nearest')
		plotData = np.sort(plotData)
		idx1, idx2 = plotData.tolist().index(__val_2p), plotData.tolist().index(__val_98p) 

		# Check for non-unique idx #
		c1, c2 = 0, 0
		for p in plotData:
			if p == __val_2p:	
				c1 += 1
			if p == __val_98p:
				c2 += 1
		if c1 > 1 or c2 > 1:
			print 'ERROR: non-unique index.'

		plotData = plotData[idx1:idx2]
		data_label = 'trim data'

	if TRIM_FLUX is False: data_label = 'data'

	
	### Plot standard normalized Gaussian and vertical line through mean to guide eye ###
        gx = np.linspace(np.min(plotData), np.max(plotData), 1000) 
        gy = gaussian(x=gx, mu=0.0, sig=1.0)
        plt.plot(gx, gy, linestyle='--', color='black', label=r'$\mu=0 \,, \, \sigma=1$', linewidth=0.7) #maroon
        plt.axvline(x=0, linestyle='--', color='black', linewidth=0.7)

	### Fit Gaussian to DeltaFlux/SigmaFlux http://danielhnyk.cz/fitting-distribution-histogram-using-python/ ###
        #FIXME is this fit normalized?
        mean, stddev = stats.norm.fit(plotData)
        fit_data = stats.norm.pdf(gx, mean, stddev)
        plt.plot(gx, fit_data, color='orangered', label=r'fit: $\mu=$'+str(round(mean, 2))+', $\sigma=$'+str(round(stddev, 2))) #FIXME cyan fit?

	if GAUSS_APER:
		gax = np.linspace(np.min(plotDataGauss), np.max(plotDataGauss), 1000) 
		ga_mean, ga_stddev = stats.norm.fit(plotDataGauss)
		ga_fit_data = stats.norm.pdf(gax, ga_mean, ga_stddev)
		plt.plot(gax, ga_fit_data, color='maroon', label=r'ga_fit: $\mu=$'+str(round(ga_mean, 2))+', $\sigma=$'+str(round(ga_stddev, 2)))

	### Plot histogram of data. Normalize by counts via `density=True` (area under histogram equals 1) ###
	plt.hist(plotData, __flux_hist_bins, density=True, histtype='step', color=PT_COLORS[band], label=data_label, rwidth=1)
	if GAUSS_APER: plt.hist(plotDataGauss, __flux_hist_bins, density=True, histtype='step', color='red', label=gauss_aper_label, rwidth=1)


	### Get median of data. Get x-value corresponding to peak of data distribution ###
	__bin_size, __bin_edges = np.histogram(plotData, __flux_hist_bins, density=True)
	# Find bin with max counts #
	__max_val = max(__bin_size)
	for i in np.arange(0, len(__bin_size)):
		if __bin_size[i] == __max_val:
			__ymax_idx = i
	print 'band:', band
	print 'Corresponding hax value of max vax value:', __bin_edges[__ymax_idx]
	print 'Data median:', np.median(plotData)
	print ' \n'
	data_median = np.median(plotData)
	peak = __bin_edges[__ymax_idx]


	'''	
	plt.plot(peak, __max_val, marker='v', color=c1, label='peak: '+str(round(peak, 2)))
	plt.plot(data_median, __max_val, marker='^', color=c1, label='peak: '+str(round(peak, 2)))
	'''
	#TODO make labels clear that this is NOT for gauss aper
	#plt.axvline(x=data_median, color=HIST_COLORS[band], linestyle=':', label='data median: '+str(round(data_median, 2)))
	#plt.axvline(x=peak, color=HIST_COLORS[band], linestyle='-.', label='data peak: '+str(round(peak, 2)))


	### Labels ###
	haxLabel = get_flux_histogram_haxlabel(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, band=band)
	plt.xlabel(haxLabel)
	plt.ylabel('Count')
	plt.legend(fontsize=10).draggable()

	return 0




def get_flux_histogram_haxlabel(mag_hdr1, mag_hdr2, band):
	"""Get the horizontal axis label for the normalized flux histogram.

	Parameters
	----------
	mag_hdr1, mag_hdr2 (str)
		Magnitude header for `MATCH_CAT1`, `MATCH_CAT2`. Headers refer to matched catalogs.

	band (str)
		Refers to one of the bands: g, r, i, z.

	Returns
	-------
	 __flux_haxlabel (str)
		Label for the horizontal axis of plots.
	"""

	if 'cm_mag' in mag_hdr1 and 'cm_mag' in mag_hdr2:
		 __flux_haxlabel = 'cm_flux_$\\bf{'+str(band)+'}$ meas$-$true) / $\sigma_{flux\_meas}$'

	#TODO will we deal with psf mag ever?

	if INJ1 and INJ2:
		if INJ1_PERCENT == 10 and INJ2_PERCENT == 10:
			 __flux_haxlabel = '(10%_inj_'+ __flux_haxlabel
		if INJ1_PERCENT == 20 and INJ2_PERCENT == 20:
			 __flux_haxlabel = '(20%_inj_'+ __flux_haxlabel
	

	if INJ1 is False and INJ2 is False:
		 __flux_haxlabel = '(base_'+ __flux_haxlabel

	return  __flux_haxlabel




def color_subplotter(band, df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization, tile, fn_flag_log, plot_title, fn_plot, fn_color_log):
	"""Plot 2x2 grid of plots. Each plot is a `corner.hist2d` of the color corresponding to different magnitude bins (see `get_color_from_binned_magnitude()` for the magnitude bins used). Get plot labels for the horizontal and vertical axes.

	Parameters
        ----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.

	mag_hdr1, mag_hdr2 (str)

	mag_err_hdr1, mag_err_hdr2 (str)

	fn_plot (str)
		Complete file name for the plot. Used if `SAVE_PLOT=True`.

	plot_title (str)
		Suptitle for the plot.

	bands (str)

	realization (str)

	tile (str)

	Returns
	-------
	0
	"""	

	# For log file #
	writeColor = WRITE_COLORS[band]

	axLabel1 = get_color_axlabel(inj_percent=INJ1_PERCENT, inj=INJ1, meas_or_true_cat=AXLABEL1, match_cat=MATCH_CAT1, band=band)
	axLabel2 = get_color_axlabel(inj_percent=INJ2_PERCENT, inj=INJ2, meas_or_true_cat=AXLABEL2, match_cat=MATCH_CAT2, band=band)

	vaxDeltaColorLabel = get_color_difference_axlabel(color_haxlabel1=axLabel1, color_haxlabel2=axLabel2, band=band)

	### Get plot data ###
	errMag1, errMag2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel = get_magnitude_plot_variables(band=band, df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot)

	#TODO get_magnitude_error is called twice
	#TODO get error or simply use corner.hist2d
	'''
	if 'truth' not in MATCH_CAT1:
		err1 = get_color_plot_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, flux_cov_hdr=CM_FLUX_COV_HDR1, df=df, band=band, idx_good=idxGood, match_cat=MATCH_CAT1)
	if 'truth' not in MATCH_CAT2:
		err2 = get_color_plot_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, flux_cov_hdr=CM_FLUX_COV_HDR2, df=df, band=band, idx_good=idxGood, match_cat=MATCH_CAT2)
	'''

	if band != 'z':
		cleanColor1, cleanColor2, magBinsForColor = get_color_from_binned_magnitude(df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, band=band, idx_good=idxGood, clean_mag1_band1=cleanMag1, clean_mag2_band1=cleanMag2)


	if PLOT_DIFF_ON_VAX is False:
		__vax_color = cleanColor2

	### Plot ###
        if SWAP_HAX is False:
                __hax_label = axLabel1 
		__vax_label = axLabel2
		__subtitle_pref = magAxLabel1
        if SWAP_HAX:
                __hax_label = axLabel2 
		__vax_label = axLabel1
		__subtitle_pref = magAxLabel2

	if PLOT_DIFF_ON_VAX:
                __vax_color = np.array(cleanColor1) - np.array(cleanColor2)
                __vax_label = vaxDeltaColorLabel 


	__lw = 1.1

        ### Create 4-by-4 subplot. Figure size units: inches ###
	if PLOT_DIFF_ON_VAX is False:
		plt.figure(figsize=(8, 12))

	if PLOT_DIFF_ON_VAX:
		plt.figure(figsize=(12, 10))


        for i in np.arange(0, len(cleanColor1)):

		# Write to log file with headers 'TILE REALIZATION COLOR MAG_BIN OBJS_IN_MAG_BIN TOTAL_OBJS_FLAGS_INC RUN_TYPE #
		__write_bin = '['+str(magBinsForColor[i])+','+str(magBinsForColor[i+1])+')'
		# Append to csv #
		with open(fn_color_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([str(tile), str(realization), str(writeColor), str(__write_bin), str(len(cleanColor1[i])), str(len(cleanMag1)), str(RUN_TYPE)])

		print 'Plotting ', len(cleanColor1[i]), ' objects ... \n'
                plt.subplot(2, 2, i+1)

		if CORNER_HIST_2D:
			# Create bins of 1/4 magnitude for hax # 
			__bin_x = np.linspace(np.min(cleanColor1[i]), np.max(cleanColor1[i]), math.ceil(4.0*(np.max(cleanColor1[i]) - np.min(cleanColor1[i]))))
			# Create bins of 1/20 magnitude for vax #
			__ylow = np.min(np.array(__vax_color[i])) # -1*np.max(abs(np.array(__vax_color[i]))) ?
			__yhigh = np.max(np.array(__vax_color[i])) # np.max(abs(np.array(__vax_color[i]))) ?
			__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))

			# Force symmetric vertical axis #
			__sym_ylim = np.mean([abs(__ylow), __yhigh])
			if abs(abs(__ylow) - __yhigh) > 1:
				__sym_ylim = np.min([abs(__ylow), __yhigh])

			# Plot density bins #
                        #corner.hist2d(cleanColor1[i], __vax_color[i], bins=np.array([__bin_x, __bin_y]), no_fill_contours=False, color=get_color(band=band)[0], levels=LVLS, contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw})

			# Plot points and not density bins #
			corner.hist2d(cleanColor1[i], __vax_color[i], plot_density=False, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, color=PT_COLORS[band], levels=LVLS, contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw}, data_kwargs={'alpha':0.35, 'ms':1.75})

			if MAG_YLOW is None and MAG_YHIGH is None:
				# Force symmetric vertical axis #
				plt.ylim([-1*__sym_ylim, __sym_ylim])
			if MAG_YLOW is not None and MAG_YHIGH is not None:
				plt.ylim([MAG_YLOW, MAG_YHIGH])

			if PLOT_DIFF_ON_VAX is False:
				plt.axis('scaled') #plt.gca().set_aspect('equal')
				# Plot x=y line to guide eye #
				lims = [np.min([plt.xlim(), plt.ylim()]), np.max([plt.xlim(), plt.ylim()])] 
				plt.plot(lims, lims, color='k', linestyle=':', linewidth=0.7)
			if PLOT_DIFF_ON_VAX:
				plt.axhline(y=0, color='black', linestyle=':', linewidth=0.7)
				plt.subplots_adjust(hspace=0.6)

			# Work-around for contour label #
			for j in np.arange(0, len(LVLS)):
				plt.plot([0.5, 0.5], [0.5, 0.5], color=CLRS_LABEL[j], label='$P_{'+str(round(LVLS[j], 2))[2:]+'}$', linewidth=__lw)
			plt.legend().draggable()

                plt.ylabel(__vax_label)
                plt.xlabel(__hax_label)
                plt.title(__subtitle_pref+' bins: ' + __write_bin) 

        plt.subplots_adjust(hspace=0.4)
        plt.suptitle(plot_title, fontweight='bold')

        if SAVE_PLOT:
		fn_plot = fn_plot.replace('placehold', writeColor) 
                print '-----> Saving plot as: ', fn_plot 
                plt.savefig(fn_plot)

        if SHOW_PLOT: plt.show()

	return 0 




def get_flux_plot_variables(band, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2): 
	"""Get variables needed for flux histogram.
	
	Parameters
        ----------
	df_1and2 (pandas DataFrame)
		Contains DeltaFlux/SigmaFlux where DeltaFlux is FluxMeasured-FluxTrue

	flux_hdr1, flux_hdr2 (str)
		Headers for flux. Headers refer to the matched (via join=1and2).

	mag_hdr1, mag_hdr2 (str)
		Headers for magnitude for `MATCH_CAT1`, `MATCH_CAT2`. Headers refer to the matched (via join=1and2). 

	mag_err_hdr1, mag_err_hdr2 (str)
		Headers for the magnitude error for `MATCH_CAT1`, `MATCH_CAT2`. Headers refer to the matched (via join=1and2).

	band

        Returns
        -------
	__norm_flux_diff (list of floats)
		Contains DeltaFlux/SigmaFlux with flags* removed. DeltaFlux is given by FluxMeasured-FluxTrue. SigmaFlux is given by 'cm_flux_cov_{band}_{band}'^0.5

	idxGood (list of ints)

	cleanFlux2-cleanFlux1 (list of floats)
		Contains DeltaFlux where DeltaFlux is FluxMeasured-FluxTrue.
	"""

	fullMag1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=band)
        fullMag2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=band)

        idxGood = get_good_indices_using_primary_flags(df_1and2=df_1and2, full_mag1=fullMag1, full_mag2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, band=band)[0]


	### Get flux ###
	fullFlux1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr1, band=band)
        fullFlux2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr2, band=band)
        cleanFlux1 = np.array(fullFlux1)[idxGood]; cleanFlux2 = np.array(fullFlux2)[idxGood]

        ### Flux error from measured catalog only ###
        if 'truth' in MATCH_CAT1:
		# Use `MATCH_CAT2` error #
		__flux_cov_hdr = CM_FLUX_COV_HDR2
		__flux_s2n_hdr = 'cm_flux_s2n_2'
		__flux = cleanFlux2

        if 'truth' in MATCH_CAT2:
		__flux_cov_hdr = CM_FLUX_COV_HDR1
		__flux_s2n_hdr = 'cm_flux_s2n_1'
		__flux = cleanFlux1

	# Flux covariance matrix #
	__flux_cov = np.array(get_matrix_diagonal_element(df=df_1and2, band=band, sq_matrices_hdr=__flux_cov_hdr))[idxGood]

	# Method 1: flux_error = sqrt(flux_cov_diagonal) #
	__flux_err = __flux_cov**0.5

	# Method 2: flux_err = flux/flux_s2n. Gives same results as Method 1 #
	#__flux_s2n = np.array(get_floats_from_string(df=df, hdr=__flux_s2n_hdr, band=band))[idxGood]
	#__flux_err = __flux/__flux_s2n 
	

	# Subtract: {mof/sof/coadd} - truth #
	if 'truth' in MATCH_CAT1:
                __norm_flux_diff = (cleanFlux2-cleanFlux1)/__flux_err
		__flux_diff = cleanFlux2-cleanFlux1
        if 'truth' in MATCH_CAT2:
                __norm_flux_diff = (cleanFlux1-cleanFlux2)/__flux_err
		__flux_diff = cleanFlux1-cleanFlux2

	return __norm_flux_diff, idxGood, cleanFlux2-cleanFlux1




def get_magnitude_plot_variables(band, df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization, tile, fn_flag_log, plot_title, fn_plot, mag_axlabel1, mag_axlabel2):
	"""Get variables needed for magnitude plots.

	Parameters
	----------
	
	Returns
	-------
	"""

	### Axes labels ###
        magAxLabel1 = get_magnitude_axlabel(inj_percent=INJ1_PERCENT, inj=INJ1, mag_hdr=mag_hdr1, meas_or_true_cat=AXLABEL1, match_cat=MATCH_CAT1, band=band)
        magAxLabel2= get_magnitude_axlabel(inj_percent=INJ2_PERCENT, inj=INJ2, mag_hdr=mag_hdr2, meas_or_true_cat=AXLABEL2, match_cat=MATCH_CAT2, band=band)
        magVaxLabel = get_magnitude_difference_axlabel(mag_haxlabel1=magAxLabel1, mag_haxlabel2=magAxLabel2, band=band)


        # Magnitudes with no flags* removed # 
        fullMag1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=band)
        fullMag2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=band)

        ### Remove objects with flags* or perform quality cuts ###
	idxGood = get_good_indices_using_primary_flags(df_1and2=df_1and2, full_mag1=fullMag1, full_mag2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, band=band)[0]

        # Magnitudes with flags* removed #
        cleanMag1 = get_good_data(df_1and2=df_1and2, hdr=mag_hdr1, idx_good=idxGood, str_of_arr=True, band=band)
        cleanMag2 = get_good_data(df_1and2=df_1and2, hdr=mag_hdr2, idx_good=idxGood, str_of_arr=True, band=band)


        ### Calculate errors. get_magnitude_error() will return array of zeros for truth catalogs. ###
        magErr1 = get_magnitude_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, flux_cov_hdr=CM_FLUX_COV_HDR1, df_1and2=df_1and2, band=band, idx_good=idxGood, match_cat=MATCH_CAT1)
        magErr2 = get_magnitude_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, flux_cov_hdr=CM_FLUX_COV_HDR2, df_1and2=df_1and2, band=band, idx_good=idxGood, match_cat=MATCH_CAT2)


        ### Write flags to log file ###
        if LOG_FLAGS:
                for i in np.arange(0, len(FLAG_HDR_LIST), 2):
                        # Bad index #
                        temp_idx = log_flags(df_1and2=df_1and2, band=band, realization=realization, flag_hdr1=FLAG_HDR_LIST[i], flag_hdr2=FLAG_HDR_LIST[i+1], full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, fn_flag_log=fn_flag_log)[1]
                #FLAG_IDX.append(temp_idx)
                FLAG_IDX.extend(temp_idx)

	#TODO err--> magErr
	return magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel 




def get_color_difference_axlabel(color_haxlabel1, color_haxlabel2, band):
	"""Shorten axis label for ease of readability. Example: '10%_inj_(g-r)_true - 10%_inj_(r-i)_meas' --> '10%_inj_cm_mag_g true-meas'

	Parameters
        ----------
        color_haxlabel1, color_haxlabel2 (str)
		Labels from the horizontal axis. Contain LaTeX \bf{} formatting.

        band (str)

        Returns
        -------
        __short_color_diff_axlabel (str)
		Shortened label for vertical axis. Contains LaTeX \bf{} formatting.
        """

	__short_color_diff_axlabel = ''
	__short_color_diff_axlabel = ''

	if '10%_inj' in color_haxlabel1 and '10%_inj' in color_haxlabel2: __short_color_diff_axlabel += '10%_inj_'
        if '20%_inj' in color_haxlabel1 and '20%_inj' in color_haxlabel2: __short_color_diff_axlabel += '20%_inj_'
        if '10%_inj' in color_haxlabel1 and '20%_inj' in color_haxlabel2: __short_color_diff_axlabel += '(10%$-$20%)_inj_'
        if '20%_inj' in color_haxlabel1 and '10%_inj' in color_haxlabel2: __short_color_diff_axlabel += '(20%$-$10%)_inj_'
        if 'base' in color_haxlabel1 and 'base' in color_haxlabel2: __short_color_diff_axlabel += 'base_'

	if band == 'g': __short_color_diff_axlabel += '$\\bf{(g-r)}$_'
        if band == 'r': __short_color_diff_axlabel += '$\\bf{(r-i)}$_' 
        if band == 'i': __short_color_diff_axlabel += '$\\bf{(i-z)}$_'

        if 'true' in color_haxlabel1 and 'true' in color_haxlabel2: __short_color_diff_axlabel += 'true'
        if 'meas' in color_haxlabel1 and 'meas' in color_haxlabel2: __short_color_diff_axlabel += 'meas'

        ### Label prefix ###
        if 'inj' not in __short_color_diff_axlabel and 'base' not in __short_color_diff_axlabel:
                if 'inj' in color_haxlabel1: pref1 = 'inj'
                if 'inj' in color_haxlabel2: pref2 = 'inj'

                if 'base' in color_haxlabel1: pref1 = 'base'
                if 'base' in color_haxlabel2: pref2 = 'base'

                if 'Y3' in color_haxlabel1: pref1 = 'Y3'
                if 'Y3' in color_haxlabel2: pref2 = 'Y3'

        if 'pref1' not in locals(): pref = None
        if 'pref1' in locals(): pref = pref1 + '$-$' + pref2

        ### Label suffix ###
        if 'meas' not in __short_color_diff_axlabel and 'true' not in __short_color_diff_axlabel: suf = color_haxlabel1[-4:] + '$-$' + color_haxlabel2[-4:]
        if 'meas' in __short_color_diff_axlabel or 'true' in __short_color_diff_axlabel: suf = None

        if pref is not None: __short_color_diff_axlabel= pref + '  ' + __short_color_diff_axlabel 
        if suf is not None: __short_color_diff_axlabel = __short_color_diff_axlabel[:-1] + '  ' + suf


	return __short_color_diff_axlabel 




def get_magnitude_difference_axlabel(mag_haxlabel1, mag_haxlabel2, band):
	"""Shorten vertical axis label for ease of readability. Example: '10%_inj_cm_mag_g_true - 10%_inj_cm_mag_g_meas' --> '10%_inj_cm_mag_g true-meas'

	Parameters
	----------
	mag_haxlabel1, mag_haxlabel2 (str)
		Labels from the horizontal axis. Contain LaTeX \bf{} formatting.

	band (str)

	Returns
	-------
	__short_mag_diff_axlabel (str)
		Shortened label for vertical axis. Contains LaTeX \bf{} formatting.
	"""

	__short_mag_diff_axlabel = ''
	__short_mag_diff_axlabel = ''

	### Look for similarities in vertical axis labels and build a shared label ###
	# Inj #
        if '10%_inj' in mag_haxlabel1 and '10%_inj' in mag_haxlabel2: __short_mag_diff_axlabel += '10%_inj_'
	if '20%_inj' in mag_haxlabel1 and '20%_inj' in mag_haxlabel2: __short_mag_diff_axlabel += '20%_inj_'
	if '10%_inj' in mag_haxlabel1 and '20%_inj' in mag_haxlabel2: __short_mag_diff_axlabel += '(10%$-$20%)_inj_'
	if '20%_inj' in mag_haxlabel1 and '10%_inj' in mag_haxlabel2: __short_mag_diff_axlabel += '(20%$-$10%)_inj_'
	if 'base' in mag_haxlabel1 and 'base' in mag_haxlabel2: __short_mag_diff_axlabel += 'base_'

	# Magnitude type #
        if 'cm_mag' in mag_haxlabel1 and 'cm_mag' in mag_haxlabel2: __short_mag_diff_axlabel += 'cm_mag_$\\bf{' + band + '}$_'

	# Catalog type #
        if 'true' in mag_haxlabel1 and 'true' in mag_haxlabel2: __short_mag_diff_axlabel += 'true'
        if 'meas' in mag_haxlabel1 and 'meas' in mag_haxlabel2: __short_mag_diff_axlabel += 'meas'

	### Label prefix ###
        if 'inj' not in __short_mag_diff_axlabel:
		if '10%_inj' in mag_haxlabel1: pref1 = '10%_inj'
		if '20%_inj' in mag_haxlabel1: pref1 = '20%_inj'
		if '10%_inj' in mag_haxlabel2: pref2 = '10%_inj'
		if '20%_inj' in mag_haxlabel2: pref2 = '20%_inj'
		
		if 'base' in mag_haxlabel1: pref1 = 'base'
		if 'base' in mag_haxlabel2: pref2 = 'base'

		if 'Y3' in mag_haxlabel1: pref1 = 'Y3'
		if 'Y3' in mag_haxlabel2: pref2 = 'Y3'

	if 'pref1' not in locals(): pref = None
	if 'pref1' in locals(): pref = pref1 + '$-$' + pref2

	### Label suffix ###
        if 'meas' not in __short_mag_diff_axlabel and 'true' not in __short_mag_diff_axlabel: suf = mag_haxlabel1[-4:] + '$-$' + mag_haxlabel2[-4:]
        if 'meas' in __short_mag_diff_axlabel or 'true' in __short_mag_diff_axlabel: suf = None

        if pref is not None: __short_mag_diff_axlabel = pref + '  ' + __short_mag_diff_axlabel
        if suf is not None: __short_mag_diff_axlabel = __short_mag_diff_axlabel[:-1] + '  ' + suf

	return __short_mag_diff_axlabel




def normalized_magnitude_difference_plotter(mag_hdr1, mag_hdr2, cbar_data, mag_err1, mag_err2, band, clean_mag1, full_mag1, mag_axlabel1, clean_mag2, mag_axlabel2, plot_title, realization, tile, cbar_label, fn_plot, fn_mag_err_log, vax_label, fn_main_log, fn_mag_diff_outliers_log):
	"""Produce a mangitude plot normalized to 1sigma_mag. If specified with `PLOT_68P` and `PLOT_34P_SPLIT` the 68th percentiles of each magnitude bin (see `normalize_magnitude_plot_maintain_bin_structure()` for bins) are plotted.

	Parameters
        ----------

	band

	tile

	realization (str)

	Returns
	-------
	percentRecoveredFlagsIncluded (float)
		Percent of Balrog-injected objects recovered, calculated including objects with flags* (see README.md for a definition of flags*). Is `None` if neither `MATCH_CAT1` nor `MATCH_CAT2` is a truth catalog (`gal_truth` or `star_truth`). 
	"""

	# Args needed to call normalize_magnitude_difference_plot() #
	normVaxBins, initialBins, haxBins, magErrBinMedians, vaxBinMedian = normalize_magnitude_plot_maintain_bin_structure(clean_mag1=clean_mag1, clean_mag2=clean_mag2, mag_err1=mag_err1, mag_err2=mag_err2, band=band, tile=tile, realization=realization, fn_mag_err_log=fn_mag_err_log, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)


	# Percentiles #

	if PLOT_MAG_ERR and CORNER_HIST_2D is False:

		### Plot 1sigma_mag curve ###
		if CENTER_ERR_ABT_ZERO:
			plt.axhline(y=1.0, color='red', linestyle='--', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
			plt.axhline(y=-1.0, color='red', linestyle='--', linewidth=0.7)

		if CENTER_ERR_ABT_ZERO is False:

			counter_legend0 = 0; color0 = 'red'; lw = 0.7

			for b in np.arange(0, len(initialBins)-1):

				if vaxBinMedian[b] is not None and vaxBinMedian[b+1] is not None:
					# Horizontal bar bounds #
					x_hbound = np.array([initialBins[b], initialBins[b+1]])
					y_hbound = np.array([vaxBinMedian[b], vaxBinMedian[b]]) + 1
					# Vertical bar bounds #
					x_vbound1, x_vbound2 = np.array([initialBins[b], initialBins[b]]), np.array([initialBins[b+1], initialBins[b+1]])
					y_vbound = np.array([-1*vaxBinMedian[b]-1, vaxBinMedian[b]+1])

					# Plot legend once #
					if counter_legend0 == 0:
						plt.plot(x_hbound, y_hbound, color=color0, label=r'$1 \sigma_{mag\_meas}$ center: $\tilde{y}$', linewidth=lw)
						counter_legend0 += 1

					if counter_legend0 == 1:
						# Horizontal bar #
						plt.plot(x_hbound, y_hbound, color=color0, linewidth=lw, linestyle='--')
						plt.plot(x_hbound, -1.0*y_hbound, color=color0, linewidth=lw, linestyle='--')
						# Vertical bar #
						plt.plot(x_vbound1, y_vbound, color=color0, linewidth=lw, linestyle=':')
						plt.plot(x_vbound2, y_vbound, color=color0, linewidth=lw, linestyle=':')



		### Percentiles ###
		if PLOT_68P or PLOT_34P_SPLIT:
			vax_68percentile_bins, percentileBins, neg_vax_34percentile, pos_vax_34percentile = get_68percentile_of_normalized_magnitude_difference(binned_norm_mag_diff=normVaxBins, mag_bins_for_mag_err=initialBins, binned_hax_mag=haxBins)

		### Plot +/-34 percentile of each bin ###

		# Line width for top and sides of bins #
		lwt = 1.1; lws = 0.7

		if PLOT_34P_SPLIT:
			counter_legend1 = 0; color1 = 'cyan'

			for b in np.arange(0, len(neg_vax_34percentile)-1):
				# Horizontal bar bounds #
				x_hbound = np.array([percentileBins[b], percentileBins[b+1]])
				x_vbound1 = np.array([percentileBins[b], percentileBins[b]])
				x_vbound2 = np.array([percentileBins[b+1], percentileBins[b+1]])

				if neg_vax_34percentile[b] is not None:
					# Horizontal bar bounds #
					neg_y_hbound = np.array([neg_vax_34percentile[b], neg_vax_34percentile[b]])
					# Vertical bar bounds #
					y_vbound = np.array([neg_vax_34percentile[b], 0])
					# Plot #
					# Plot legend once #
					if counter_legend1 == 0:
						plt.plot(x_hbound, neg_y_hbound, color=color1, linewidth=lwt, label='$\pm P_{34}$')
						counter_legend1 = 1
					if counter_legend1 == 1:
						plt.plot(x_hbound, neg_y_hbound, color=color1)
						plt.plot(x_vbound1, y_vbound, color=color1, linewidth=lws, linestyle=':')
						plt.plot(x_vbound2, y_vbound, color=color1, linewidth=lws, linestyle=':')


				if pos_vax_34percentile[b] is not None:
					# Horizontal bar bounds #
					pos_y_hbound = np.array([pos_vax_34percentile[b], pos_vax_34percentile[b]])
					# Vertical bar bounds #
					y_vbound = np.array([0, pos_vax_34percentile[b]])
					# Plot #
					plt.plot(x_hbound, pos_y_hbound, color=color1, linewidth=lwt)
					plt.plot(x_vbound1, y_vbound, color=color1, linewidth=lws, linestyle=':')
					plt.plot(x_vbound2, y_vbound, color=color1, linewidth=lws, linestyle=':')


		### Plot 68 percentile of each bin ###
		if PLOT_68P:

			counter_legend2 = 0; color2 = 'fuchsia'

			for b in np.arange(0, len(vax_68percentile_bins)-1):

				if vax_68percentile_bins[b] is not None:

					# Horizontal bar bounds #
					x_hbound = np.array([percentileBins[b], percentileBins[b+1]])
					y_hbound = np.array([vax_68percentile_bins[b], vax_68percentile_bins[b]])
					# Vertical bar bounds #
					x_vbound1, x_vbound2 = np.array([percentileBins[b], percentileBins[b]]), np.array([percentileBins[b+1], percentileBins[b+1]])
					y_vbound = np.array([-1*vax_68percentile_bins[b], vax_68percentile_bins[b]])

					# Plot legend once #
					if counter_legend2 == 0:
						plt.plot(x_hbound, y_hbound, color=color2, label='$P_{68}$', linewidth=lwt)
						counter_legend2 += 1

					if counter_legend2 == 1:
						# Horizontal bar #
						plt.plot(x_hbound, y_hbound, color=color2, linewidth=lwt)
						plt.plot(x_hbound, -1.0*y_hbound, color=color2, linewidth=lwt)
						# Vertical bar #
						plt.plot(x_vbound1, y_vbound, color=color2, linewidth=lws, linestyle=':')
						plt.plot(x_vbound2, y_vbound, color=color2, linewidth=lws, linestyle=':')


	# Logger #
	plotDeltaMag, plotHaxMag, cleanBins = normalize_magnitude_difference_plot(binned_norm_mag_diff=normVaxBins, mag_bins_for_mag_err=initialBins, binned_hax_mag=haxBins)

	percent_1sig, percentRecoveredFlagsIncluded, percentRecoveredFlagsNotIncluded = main_logger(vax_mag=plotDeltaMag, band=band, clean_mag1=clean_mag1, full_mag1=full_mag1, realization=realization, tile=tile, mag_bins_for_mag_err=initialBins, hax_mag=plotHaxMag, fn_main_log=fn_main_log, mag_err_bin_medians=magErrBinMedians, vax_mag_bin_medians=vaxBinMedian)

	# One colorbar at a time. This error is caught at beginning of script #
        if SCATTER:
                plt.scatter(plotHaxMag, plotDeltaMag, color=PT_COLORS[band], alpha=0.25, s=0.25)


        if CM_T_ERR_COLORBAR or CM_T_COLORBAR:
                plt.scatter(plotHaxMag, plotDeltaMag, c=cbar_data, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
                plt.colorbar(label=cbar_label)


        if HEXBIN:
		grid = (100, 1000)
		print ' Normalized hexbin has a large number of grid cells. Will take a moment to plot ... \n'
                plt.hexbin(plotHaxMag, plotDeltaMag, gridsize=grid, cmap=CMAPS[band], bins='log')
                plt.colorbar(label='log(counts)')



	if HIST_2D:
                # 1/10 the bin size of that used in error calculation #
                bin_x = np.arange(min(plotHaxMag), max(plotHaxMag), 0.5/10)
                if MAG_YLOW is not None and MAG_YHIGH is not None:
                        # Somewhat using reported 1% error in magnitude #
                        bin_y = np.arange(MAG_YLOW, MAG_YHIGH, (MAG_YHIGH-MAG_YLOW)*0.01)
                if MAG_YLOW is None and MAG_YHIGH is None:
                        bin_y = np.arange(min(plotDeltaMag), max(plotDeltaMag), (max(plotDeltaMag)-min(plotDeltaMag))*0.01)
                plt.hist2d(plotHaxMag, plotDeltaMag, bins=[bin_x, bin_y], cmap=CMAPS[band], norm=matplotlib.colors.LogNorm())
                plt.colorbar()


        if CORNER_HIST_2D:
		# Create bins of 1/4 magnitude for hax #
		__bin_x = np.linspace(np.min(plotHaxMag), np.max(plotHaxMag), math.ceil(4.0*(np.max(plotHaxMag) - np.min(plotHaxMag))))
		# Create bins of 1/20 magnitude for vax #
		__ylow = np.min(np.array(plotDeltaMag)) # -1*np.max(abs(np.array(__vax_color[i]))) ?
		__yhigh = np.max(np.array(plotDeltaMag)) # np.max(abs(np.array(__vax_color[i]))) ?
		__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))

		# Force symmetric vertical axis #
		__sym_ylim = np.mean([abs(__ylow), __yhigh])

		__lw = 0.9
		#FIXME force bin size here
                corner.hist2d(plotHaxMag, plotDeltaMag, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, levels=LVLS, color=PT_COLORS[band], contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidth':__lw})
		# Work-around for contour labels #
		for j in np.arange(0, len(LVLS)):
			plt.plot([0.5, 0.5], [0.5, 0.5], color=CLRS_LABEL[j], label='$P_{'+str(round(LVLS[j], 2))[2:]+'}$', linewidth=__lw)
		plt.legend().draggable()

		if MAG_YLOW is not None and MAG_YHIGH is not None:
			plt.ylim([MAG_YLOW, MAG_YHIGH])
		if MAG_YLOW is None and MAG_YHIGH is None:
			plt.ylim([-1*__sym_ylim, __sym_ylim])

	### Axes labels ###
        # Horizontal axis labels #
        if SWAP_HAX: plt.xlabel(str(mag_axlabel2))
        if SWAP_HAX is False: plt.xlabel(str(mag_axlabel1))
        # Vertical axis label #
     	plt.ylabel('('+ vax_label + ') / $\sigma$')

        plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)

        # Adjust vertical axes limits #
        if MAG_YLOW is not None and MAG_YHIGH is not None: plt.ylim([MAG_YLOW, MAG_YHIGH])


        ### Plot legend ###
        if PLOT_MAG_ERR and CORNER_HIST_2D is False: plt.legend(fontsize=10).draggable()


        ### Title for subplot ###
        plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 4)) + '%')

	if SUBPLOT is False:

                fn_plot = fn_plot.replace('griz', band)

                ### Title for  ###
                plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

                ### Save plot ###
		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

                if SHOW_PLOT: plt.show()

        #plt.gca().set_aspect('equal')


	return percentRecoveredFlagsIncluded 




#TODO error-->mag_err, 
def magnitude_difference_plotter(mag_hdr1, mag_hdr2, cbar_data, mag_err1, mag_err2, band, clean_mag1, full_mag1, mag_axlabel1, clean_mag2, mag_axlabel2, plot_title, realization, tile, cbar_label, fn_plot, fn_main_log, fn_mag_err_log, vax_label, fn_mag_diff_outliers_log):
	"""Produce a single plot of 'mag1' versus 'mag1'-'mag2'. 

	Parameters
	----------
	mag_hdr1, mag_hdr2 (str) -- Repeated


	cbar_data (list of floats) -- Data used to create a colorbar. Can be `None` in which case no colorbar is added.

	cbar_label (str) -- Label for the colorbar.

	vax_label (str) -- Label for the vertical axis. Contains LaTeX \bf{} formatting.

	mag_axlabel1, mag_axlabel2 (str) -- Axes label for magnitude; will be used to label horizontal axis depending on `SWAP_HAX`. Contains LaTeX \bf{} formatting.

	clean_mag1, clean_mag2 (list of floats) -- Magnitudes in matched catalog of `MATCH_CAT1`, `MATCH_CAT2` with flagged objects removed.

	full_mag1 (list of floats) -- Refers to the matched (join=1and2) catalog. No flagged objects removed.

	plot_title (str) -- Title for main plot (as opposed to individual subplots). 

	fn_plot (str) -- Complete filename for plot save name.

	bins (list of floats) -- Bins used to calculate error.

	fd* FIXME repeat

	band (str)

	realization (str) 

	tile (str)

        Returns
	-------
	"""

	### Values to plot ###
	__mag_diff = np.array(clean_mag1) - np.array(clean_mag2) 

	if SWAP_HAX:
		__hax_mag = clean_mag2
	if SWAP_HAX is False:
		__hax_mag = clean_mag1


	### 1sigma_mag curve ###
	if PLOT_MAG_ERR:
		haxBinMedian, vaxBinMedian, magErrBinMedians, initialBins, haxBins, vaxBins, outlierCleanedHaxMag, outlierCleanedVaxMag = bin_and_cut_measured_magnitude_error(mag_err1=mag_err1, mag_err2=mag_err2, clean_mag1=clean_mag1, clean_mag2=clean_mag2, band=band, tile=tile, realization=realization, fn_mag_err_log=fn_mag_err_log, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)

		### Remove zeros from x, y, and err (zeros were placeholders for instances in which there were no objects in a particular magnitude bin) ###
		err = [temp for temp in magErrBinMedians if temp is not None]
		hax = [temp for temp in haxBinMedian if temp is not None]
		vax = [temp for temp in vaxBinMedian if temp is not None]

		if CORNER_HIST_2D is False:
			### Plot 1sigma_mag curve ###
			if CENTER_ERR_ABT_ZERO:
				plt.plot(hax, np.array(err), color='red', linestyle='-', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
				plt.plot(hax, -1*np.array(err), color='red', linestyle='-', linewidth=0.7)
			if CENTER_ERR_ABT_ZERO is False:
				plt.plot(hax, np.array(vax) + np.array(err), color='red', linestyle='-', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
				plt.plot(hax, np.array(vax) - np.array(err), color='red', linestyle='-', linewidth=0.7)


	### Write to log files ### 
	percent_1sig, percentRecoveredFlagsIncluded, percentRecoveredFlagsNotIncluded = main_logger(vax_mag=outlierCleanedVaxMag, band=band, clean_mag1=clean_mag1, full_mag1=full_mag1, realization=realization, tile=tile, mag_bins_for_mag_err=initialBins, hax_mag=outlierCleanedHaxMag, fn_main_log=fn_main_log, mag_err_bin_medians=magErrBinMedians, vax_mag_bin_medians=vaxBinMedian)


	if PRINTOUTS:
                print 'Plotting ', len(clean_mag1), ' objects ... \n'

	### Plot ###
	# One colorbar at a time. This error is caught at beginning of script #
	if SCATTER:
		plt.scatter(__hax_mag, __mag_diff, color=PT_COLORS[band], alpha=0.25, s=0.25)

	
	if CM_T_ERR_COLORBAR or CM_T_COLORBAR:
		plt.scatter(__hax_mag, __mag_diff, c=cbar_data, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
		plt.colorbar(label=cbar_label)


	if HEXBIN:
		grid = 500
		plt.hexbin(__hax_mag, __mag_diff, gridsize=grid, cmap=CMAPS[band], bins='log')
		plt.colorbar(label='log(counts)')


	if HIST_2D:
		# 1/10 the bin size of that used in error calculation #
		bin_x = np.arange(min(__hax_mag), max(__hax_mag), 0.5/10)
		if MAG_YLOW is not None and MAG_YHIGH is not None:
			# Somewhat using reported 1% error in magnitude #
			bin_y = np.arange(MAG_YLOW, MAG_YHIGH, (MAG_YHIGH-MAG_YLOW)*0.01)
		if MAG_YLOW is None and MAG_YHIGH is None:
			bin_y = np.arange(min(__mag_diff), max(__mag_diff), (max(__mag_diff)-min(__mag_diff))*0.01) 
		plt.hist2d(__hax_mag, __mag_diff, bins=[bin_x, bin_y], cmap=CMAPS[band], norm=matplotlib.colors.LogNorm())
		plt.colorbar()


	if CORNER_HIST_2D: 
		__lw = 0.9

		# Create bins of 1/4 magnitude for hax #
		__bin_x = np.linspace(np.min(__hax_mag), np.max(__hax_mag), math.ceil(4.0*(np.max(__hax_mag) - np.min(__hax_mag))))
		# Create bins of 1/20 magnitude for vax #
		__ylow = np.min(np.array(__mag_diff)) # -1*np.max(abs(np.array(__vax_color[i]))) ?
		__yhigh = np.max(np.array(__mag_diff)) # np.max(abs(np.array(__vax_color[i]))) ?
		__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))

		# Force symmetric vertical axis #
		__sym_ylim = np.mean([abs(__ylow), __yhigh])

		# Plot points and not density bins #
		#corner.hist2d(__hax_mag, __mag_diff, plot_density=False, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, color=get_color(band=band)[0], levels=LVLS, contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw}, data_kwargs={'alpha':0.25, 'ms':1.75})

		# Only the densest regions of the plot are binned so increase bin size of plt.hist2d() #
		# SLACK channel corner.hist2d "draws 1- and 2-sigma contours automatically."Correct 1sigma levels: http://corner.readthedocs.io/en/latest/pages/sigmas.html #
		corner.hist2d(__hax_mag, __mag_diff, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, levels=LVLS, color=PT_COLORS[band], contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw}) 

		if MAG_YLOW is None and MAG_YHIGH is None:
			plt.ylim([-1*__sym_ylim, __sym_ylim])
		if MAG_YLOW is not None and MAG_YHIGH is not None:
			plt.ylim([MAG_YLOW, MAG_YHIGH])

		# Work-around for contour labels #
		for j in np.arange(0, len(LVLS)):
			plt.plot([0.5, 0.5], [0.5, 0.5], color=CLRS_LABEL[j], label='$P_{'+str(round(LVLS[j], 2))[2:]+'}$', linewidth=__lw)
		plt.legend().draggable()

	### Axes labels ###
	# Horizontal axis labels #
	if SWAP_HAX: plt.xlabel(str(mag_axlabel2))
	if SWAP_HAX is False: plt.xlabel(str(mag_axlabel1))
	# Vertical axis label #
        plt.ylabel(vax_label)

	plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)

	# Adjust vertical axes limits #
        if MAG_YLOW is not None and MAG_YHIGH is not None: plt.ylim([MAG_YLOW, MAG_YHIGH])


	### Plot legend ###
	if PLOT_MAG_ERR and CORNER_HIST_2D is False: plt.legend(fontsize=10).draggable()


	### Title for subplot ###
	plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 4)) + '%')


	if SUBPLOT is False:

		fn_plot = fn_plot.replace('placehold', band)

		### Title for  ###
		plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

		### Save plot ###
		if SAVE_PLOT: print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

		if SHOW_PLOT: plt.show()

        #plt.gca().set_aspect('equal')


        return percentRecoveredFlagsIncluded 




def magnitude_difference_subplotter(df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, tile, fn_flag_log, fn_mag_err_log,  fn_mag_diff_outliers_log, fn_main_log):
	"""Combine four subplots into a single plot with four panels (2-by-2). 

	Parameters
	----------
	df_1and2 (pandas DataFrame) -- DataFrame for the matched catalog. 

	mag_hdr1, mag_hdr2 (str) -- Headers for the magnitude of `MATCH_CAT1` and `MATCH_CAT2` respectively after the catalogs have been matched.

	mag_err_hdr1, mag_err_hdr2 (str) -- Headers for the magnitude error of `MATCH_CAT1` and `MATCH_CAT2` respectively after matching. Can be `None`.

	fn_plot (str) -- Complete name for plot.

	plot_title (str) -- Title for 2-by-2 plot. 

	realization (str) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.

	tile (str)

        Returns
	-------
	0
	"""



	# Counter for flag type() printout #
	counter_flag_type_printout = 0


        ### Create 4-by-4 subplot ###
	counter_subplot = 1
	# Figure size units: inches #
	plt.figure(figsize=(12, 10))


        ### Create one subplot for each griz band ###
	for b in ALL_BANDS:


		### Define variables ###
		magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel = get_magnitude_plot_variables(band=b, df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot)


		if CM_T_COLORBAR or CM_T_ERR_COLORBAR:
			cbarData, cbarLabel = get_colorbar_for_magnitude_plot_properties(df_1and2=df_1and2, cm_t_hdr=CM_T_HDR2, cm_t_err_hdr=CM_T_ERR_HDR2, idx_good=idxGood, clean_mag1=cleanMag1, clean_mag2=cleanMag2, meas_or_true_cat=AXLABEL2, inj=INJ2, inj_percent=INJ2_PERCENT) #FIXME use same suf
		else: 
			cbarData, cbarLabel = None, None



		### Subplot ###
		if SUBPLOT: 
			plt.subplot(2, 2, counter_subplot)

		if NORMALIZE is False:
			percentRecoveredWithFlags = magnitude_difference_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_data=cbarData, plot_title=plot_title, mag_err1=magErr1, mag_err2=magErr2, band=b, full_mag1=fullMag1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, mag_axlabel1=magAxLabel1, mag_axlabel2=magAxLabel2, realization=realization, tile=tile, cbar_label=cbarLabel, fn_plot=fn_plot, fn_main_log=fn_main_log, fn_mag_err_log=fn_mag_err_log, vax_label=magVaxLabel, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)

		if NORMALIZE:
			percentRecoveredWithFlags = normalized_magnitude_difference_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_data=cbarData, plot_title=plot_title, mag_err1=magErr1, mag_err2=magErr2, band=b, full_mag1=fullMag1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, mag_axlabel1=magAxLabel1, mag_axlabel2=magAxLabel2, realization=realization, tile=tile, cbar_label=cbarLabel, fn_plot=fn_plot, fn_main_log=fn_main_log, fn_mag_err_log=fn_mag_err_log, vax_label=magVaxLabel, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)

		counter_subplot += 1

	# Plot after all panels have been filled #
	if SUBPLOT:

		### Show or save the plot once all four subplots have been filled ###
		plt.subplots_adjust(hspace=0.4)
		#plt.subplots_adjust(wspace=0.3)
		#plt.tight_layout(pad=3, h_pad=2.5)


		### Title ###
		if percentRecoveredWithFlags is not None:
			plot_title = ' '.join([plot_title, 'Recovered (with flags):', str(round(percentRecoveredWithFlags, 4))+'%'])
		plt.suptitle(plot_title, fontweight='bold')

		### Save plot ###
		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

		### Show plot ###
		if SHOW_PLOT: plt.show()

	
	return 0 




def get_plot_suptitle(realization, tile, number_of_stacked_realizations, number_of_stacked_tiles):
	"""Generate plot title.

	Parameters
	----------
	realization (str)


	tile (str)

	number_of_stacked_realizations (int) -- Number of catalogs in stacked realization catalog. Can be `None`.

	number_of_stacked_tiles (int) -- Number of catalogs in stacked tile catalog. Can be `None`.

	Returns
	-------
	__plot_title (str) -- Ex: '10% Inj MOF Cat & 10% Inj Truth Cat' 
	"""

	if STACK_REALIZATIONS:
		realization = 'stacked '+str(number_of_stacked_realizations)
	if STACK_TILES:
		tile = 'stacked ' + str(number_of_stacked_tiles)

	__plot_title = str(TITLE_PIECE1) + ' & ' + str(TITLE_PIECE2) +'. Tile: ' + str(tile) 

	if cmd_line_realizations[0] != 'None': 
		__plot_title = __plot_title + '. Realization: ' + str(realization) + '.'	
	
	if RUN_TYPE == 'ok': 
		__plot_title = __plot_title + ' Unchanged FOF groups.'
	if RUN_TYPE == 'rerun':
		__plot_title = __plot_title + ' Changed FOF groups.'

	if NORMALIZE:
		__plot_title = 'Normalized. ' + __plot_title 

	return __plot_title 




def get_plot_save_filename(realization, tile):
        """Generate name of the plot that will be used in plt.savefig() if `SAVE_PLOT=True`.
	Relies on directory structure: /`OUTDIR`/plots/`BALROG_RUN`/`match_type`/{tile}/{realization}/plots/`plot_obs`/ `plot_obs` can be: 'color' 'magnitude' 'flux'.
	The FOF analysis plots are saved differently: /`OUTDIR`/plots/`BALROG_RUN`/`match_type`/{tile}/{realization}/plots/{plot_obs}/'fof_analysis'/

        Parameters
	----------
	realization (str) 

	tile (str)

        Returns
	-------
	fn (str)
        """

	### `plot_obs` describes the observable ### 
        if PLOT_COLOR: plot_obs = 'color'
        if PLOT_MAG: plot_obs = 'magnitude'
        if PLOT_FLUX: plot_obs = 'flux'
	#if PLOT_DIFF_ON_VAX: plot_type1 = 'diff' + plot_obs #TODO are we interested in x vs y plots?

	# By construction, completeness plots contain both 10% and 20% injected catalogs. Remove this distinction #
        if PLOT_COMPLETENESS: match_type = MATCH_TYPE.replace('10%_', ''); match_type = match_type.replace('20%_', '')
        if PLOT_COMPLETENESS is False: match_type = MATCH_TYPE

        plot_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, match_type, tile, realization, 'plots', plot_obs)
	if RUN_TYPE is not None: plot_dir = os.path.join(plot_dir, 'fof_analysis')


	# Get plot type #
	if PLOT_FLUX: plot_type = 'norm_flux_hist_'
        if CORNER_HIST_2D: plot_type = 'cornerhist2d_'
	if HIST_2D: plot_type = 'hist2d_'
        if PLOT_COMPLETENESS: plot_type = 'completeness_'
        if SCATTER: plot_type = 'scatter_'
        if CM_T_COLORBAR: plot_type = 'cbar_cm_t_'
        if CM_T_ERR_COLORBAR: plot_type = 'cbar_cm_t_err_'
        if HEXBIN: plot_type = 'hexbin_'
        if 'plot_type' not in locals(): plot_type = '' 
	

	### End of filename after lowest-level directory: {tile}_{realization}_{band(s)/color}_{match_type}_{vertical_axis_limits}.png ###
	if MAG_YLOW is None and MAG_YHIGH is None:
                # Default scale for the vertical axis (vax) #
                ylim = 'defaultvax'
        if MAG_YLOW is not None and MAG_YHIGH is not None:
                ylim = str(MAG_YLOW)+'y'+str(MAG_YHIGH)

	# `placehold` will be replaced with {band(s)/color} when `savefig()` is called within various functions #
	if RUN_TYPE is None:	
		endname = '_'.join([tile, realization, 'placehold', match_type, ylim+'.png'])
	if RUN_TYPE is not None:
		endname = '_'.join([tile, realization, 'placehold', match_type, RUN_TYPE, ylim+'.png'])
	if PLOT_COMPLETENESS:
		# By construction, completeness plots contain 10% and 20% injected catalogs. Remove this distinction #
                endname = endname.replace('10%_', ''); endname = endname.replace('20%_', '')
	
	# Filename after lowest-level directory /*/*/*/`fn`/ #
	fn = plot_type + endname


	#if PLOT_DIFF_ON_VAX: outname = 'vax_diff' + outname


	### Check for directory existence ###
	if os.path.isdir(plot_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(plot_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_plot_save_filename() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', plot_dir, '...\n'
			os.makedirs(plot_dir)


	### Get complete filename (including path) ###
        if NORMALIZE:
		__fn_plot = os.path.join(plot_dir, 'norm_' + str(fn))
        if NORMALIZE is False:
		__fn_plot = os.path.join(plot_dir, fn)


        return __fn_plot




def get_coadd_catalog_magnitude_and_magnitude_error(fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z, mag_hdr, mag_err_hdr):
	"""Creates a list of magnitudes and magnitude errors of form '(mag_g, mag_r, mag_i, mag_z)' from four catalogs. Solely for use with coadd catalogs.

	Parameters
	----------
	fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z (str)
		Complete filename for the coadd catalogs for the g, r, i, z bands. Must be FITS files.

	mag_hdr (str) -- Header for magnitude. Headers refer to columns in the matched catalog.

	err_hdr (str) -- Header for magnitude error.

	Returns
	-------
	__coadd_mag_griz (list of str) -- Stores magnitude of each band in form '(mag_g, mag_r, mag_i, mag_z)'

	__coadd_mag_err_griz (list of str) -- Stores error in magnitude of each band in form '(mag_g, mag_r, mag_i, mag_z)'
	"""

	# Files have not yet been matched, and do not have hdr_1 #
	mag_hdr = mag_hdr[:-2]
	mag_err_hdr = mag_err_hdr[:-2]

	# Open FITS files #
	hdu_g = fits.open(fn_coadd_cat_g); hdu_r = fits.open(fn_coadd_cat_r); hdu_i = fits.open(fn_coadd_cat_i); hdu_z = fits.open(fn_coadd_cat_z)
	
	# Read data #
	data_g = hdu_g[1].data; data_r = hdu_r[1].data; data_i = hdu_i[1].data; data_z = hdu_z[1].data

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	err_g = data_g[mag_err_hdr]; err_r = data_r[mag_err_hdr]; err_i = data_i[mag_err_hdr]; err_z = data_z[mag_err_hdr]

	__coadd_mag_griz, __coadd_mag_err_griz = [], []

        for i in np.arange(0, len(m_g)):
                __coadd_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")
		__coadd_mag_err_griz.append("'("+ str(err_g[i])+ ', ' + str(err_r[i])+ ', ' + str(err_i[i]) + ', ' + str(err_z[i]) + ")'")

        return __coadd_mag_griz, __coadd_mag_err_griz




def get_star_truth_catalog_magnitude(df_1and2, suf):
	"""Computes and creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'. Solely for use with star truth catalogs.

	Parameters
	----------
	df_1and2 (pandas DataFrame)

	suf (str)

	Returns
	-------
	__star_truth_mag_griz (list of str) 
	"""

	m_g = df_1and2['g_Corr'+suf]

	m_r = df_1and2['g_Corr'+suf] - df_1and2['gr_Corr'+suf]

	m_i = df_1and2['g_Corr'+suf] - df_1and2['gr_Corr'+suf] - df_1and2['ri_Corr'+suf]

	m_z = df_1and2['g_Corr'+suf] - df_1and2['gr_Corr'+suf] - df_1and2['ri_Corr'+suf] - df_1and2['iz_Corr'+suf]

	__star_truth_mag_griz = []

	for i in np.arange(0, len(m_g)):
		__star_truth_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")	
        
	return __star_truth_mag_griz 




#TODO mag_hdr --> y3_gold_mag_hdr? Similarly for star_truth and coadd
def get_y3_gold_catalog_magnitude(df_1and2, mag_hdr):
	"""Creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'. Solely for use with Y3 Gold catalogs.

	Parameters
	----------
	df_1and2 (pandas DataFrame)


	mag_hdr (str)

        Returns
	-------
	__y3_gold_mag_griz (list of str) 
        """

	# Get headers, which are dependent on band #
	hdr_g = mag_hdr[:-2] + '_G' + mag_hdr[-2:]; hdr_r = mag_hdr[:-2] + '_R' + mag_hdr[-2:]
	hdr_i = mag_hdr[:-2] + '_I' + mag_hdr[-2:]; hdr_z = mag_hdr[:-2] + '_Z' + mag_hdr[-2:]
	
	# Read magnitudes from DataFrame #
	m_g = df_1and2[hdr_g]; m_r = df_1and2[hdr_r]; m_i = df_1and2[hdr_i]; m_z = df_1and2[hdr_z]

	__y3_gold_mag_griz = []

        for i in np.arange(0, len(m_g)):
                __y3_gold_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")

        return __y3_gold_mag_griz 




def get_catalog_filename(cat_type, inj, inj_percent, realization, tile, band):
        """Get catalog to analyze.
	
        Parameters
	----------
	cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd', 'y3_gold'. Set by `MATCH_CAT1` or `MATCH_CAT2`.

	inj_percent (int)
		If `inj` is False this will not be used nor referenced.

	inj (bool)

	realization (str) -- 

	tile -- 

	band (str) -- Only used with coadd catalogs. Ignored if `cat_type` is not 'coadd'. 

        Returns
	-------
	__fn_cat (str) -- Complete catalog filename.
        """

	if (cat_type == 'gal_truth' or cat_type == 'star_truth') and inj is False:
		sys.exit('No non-injected truth catalog exists.')


	if BALROG_RUN == 'test':
		if cat_type == 'gal_truth' and inj_percent == 20 and inj:
                        __fn_cat = os.path.join(BASEPATH, tile, tile+'_0_balrog_truth_cat_gals.fits')
		if cat_type == 'sof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, tile, 'real_0_'+tile+'_sof.fits')
		if cat_type == 'mof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, tile, 'real_0_'+tile+'_mof.fits')


	# TODO GoogleDoc link to Balrog log to see which catalogs can exist within ea `BALROG_RUN` 

	if BALROG_RUN == 'COSMOS': #Jun 18
                if cat_type == 'gal_truth' and inj_percent == 20 and inj:
                        __fn_cat = os.path.join(BASEPATH, tile, tile+'_0_balrog_truth_cat_gals.fits')
                if cat_type == 'sof' and inj_percent == 20 and inj:
                        __fn_cat = os.path.join(BASEPATH, tile, 'real_0_'+tile+'_sof.fits')
                if cat_type == 'mof' and inj_percent == 20 and inj:
                        __fn_cat = os.path.join(BASEPATH, tile, 'real_0_'+tile+'_mof.fits')

	if BALROG_RUN == 'DES0236+0001_COSMOS':
		if cat_type == 'gal_truth' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, tile+'_0_balrog_truth_cat_gals.fits')
		if cat_type == 'sof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, 'real_0_'+tile+'_sof.fits')
		if cat_type == 'mof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, 'real_0_'+tile+'_mof.fits')

	# TODO As of May 2018 only TAMU_Balrog catalogs have 20% injections
	if BALROG_RUN == 'TAMU_Balrog':
		__fn_cat = get_tamu_catalog_filename(cat_type=cat_type, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=inj)


	if BALROG_RUN != 'TAMU_Balrog' and 'COSMOS' not in BALROG_RUN and BALROG_RUN != 'test':

		# $/data/des61.a/data/kuropat/blank_test/DES0102-4914/ #
		if cat_type == 'gal_truth' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_gals.fits')
		if cat_type == 'mof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'mof', tile+'_mof.fits')			
		if cat_type == 'sof' and inj_percent == 20 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'sof', tile+'_sof.fits')


		if cat_type == 'gal_truth' and inj_percent == 10 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_gals.fits')

		if cat_type == 'star_truth' and inj_percent == 10 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_stars.fits')

		if cat_type == 'sof' and inj_percent == 10 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'sof', tile+'_sof.fits')
		if cat_type == 'sof' and inj is False:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', tile, 'sof', tile+'_sof.fits')

		if cat_type == 'mof' and inj_percent == 10 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'mof', tile+'_mof.fits')
		if cat_type == 'mof' and inj is False:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', tile, 'mof', tile+'_mof.fits')

		if cat_type == 'coadd' and inj_percent == 10 and inj:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'coadd', tile+'_'+band+'_cat.fits')
		if cat_type == 'coadd' and inj is False:
			__fn_cat = os.path.join(BASEPATH, 'y3v02', tile, 'coadd', tile+'_'+band+'_cat.fits')


	# Y3 catalogs cannot be injected #
	#TODO 'y3_gold_2_0' and 'y3_gold_2_2'
	if cat_type == 'y3_gold_2_0':
		# !!!!! User may need to alter path to Y3 Gold catalog #
		__fn_cat = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile+'_y3_gold_2_0.fits')
	if cat_type == 'y3_gold_2_2':
		__fn_cat = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile+'_y3_gold_2_2.fits')

        return __fn_cat 




def get_tamu_catalog_filename(cat_type, inj, inj_percent, realization, tile, band):
	"""Get catalog for TAMU runs.

	Parameters
	----------
	inj_percent (int) -- Set by `INJ1_PERCENT` or `INJ2_PERCENT` typically.

	inj (bool) -- 

	cat_type (str)
		Catalog type. This is set by `MATCH_CAT1` or `MATCH_CAT2`.

	realization (str) --

	band (str)

	Returns
	-------
	__fn_tamu_cat (str) -- Complete catalog filename.
	"""

	# 20% injection directory #	
	if inj_percent == 20:
		if cat_type == 'mof' and inj:
                        __fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', 'real_' + realization + '_' + tile + '_mof.fits')
		if cat_type == 'mof' and inj is False:
			__fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', 'base_' + tile + '_mof.fits')

                if cat_type == 'sof' and inj:
                        __fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', 'real_' + realization + '_' + tile + '_sof.fits')
		if cat_type == 'sof' and inj is False:
			__fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', 'base_' + tile + '_sof.fits')

                if cat_type == 'gal_truth':
                        __fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', tile + '_' + realization + '_balrog_truth_cat_gals.fits')

		if cat_type == 'star_truth':
                        __fn_tamu_cat = os.path.join(BASEPATH, tile + '_20', tile + '_' + realization + '_balrog_truth_cat_stars.fits')

	# 10% injection directory #
	if inj_percent == 10:
		if cat_type == 'mof' and inj:
			__fn_tamu_cat = os.path.join(BASEPATH, tile, 'real_' + realization + '_' + tile + '_mof.fits')
		if cat_type == 'mof' and inj is False:
			__fn_tamu_cat = os.path.join(BASEPATH, tile, 'base_' + tile + '_mof.fits')

		if cat_type == 'sof' and inj:
			__fn_tamu_cat = os.path.join(BASEPATH, tile, 'real_' + realization + '_' + tile + '_sof.fits') 
		if cat_type == 'sof' and inj is False:
			__fn_tamu_cat = os.path.join(BASEPATH, tile, 'base_' + tile + '_sof.fits')

		if cat_type == 'gal_truth':
			__fn_tamu_cat = os.path.join(BASEPATH, tile, tile + '_' + realization + '_balrog_truth_cat_gals.fits')
		if cat_type == 'star_truth':
			__fn_tamu_cat = os.path.join(BASEPATH, tile, tile + '_' + realization + '_balrog_truth_cat_stars.fits')

	return __fn_tamu_cat 




def matcher(realization, tile, band, inj1, inj2, inj1_percent, inj2_percent):
        """Match two catalogs on RA and DEC with a tolerance of 1 arcsecond via STILTS.

        Parameters
	----------
	realization (str or int) -- Currently allowed values: 0 1 2 3 4 5 6 7 8 9 depending on the basepath.


	tile (str) -- Currently allowed values: DES0347-5540  DES2329-5622  DES2357-6456 DES2247-4414 depending on the basepath.

	band (str)

        Returns
	-------

	__fn_match_1and2 (str) -- Name of catalog matched via join=1and2. Note that headers of matched catalogs will have '_1' appended or '_2' appended.

	__fn_match_1not2 (str) -- Name of catalog matched via join=1not2.

	__fn_match_2not1 (str) -- Name of catalof matched via join=2not1.
        """

        ### Get arguments to pass to ms_matcher ###

        # Input catalogs for STILTS #
	if MATCH_CAT1 is not 'coadd':
		in1 = get_catalog_filename(cat_type=MATCH_CAT1, inj=inj1, inj_percent=inj1_percent, realization=realization, tile=tile, band=band)
	if MATCH_CAT1 == 'coadd':
		in1 =  get_coadd_catalog_for_matcher(cat_type=MATCH_CAT1, inj_percent=inj1_percent,  realization=realization, tile=tile, mag_hdr=M_HDR1, mag_err_hdr=M_ERR_HDR1)

	if MATCH_CAT2 is not 'coadd':
		in2 = get_catalog_filename(cat_type=MATCH_CAT2, inj=inj2, inj_percent=inj2_percent, realization=realization, tile=tile, band=band)	
	if MATCH_CAT2 == 'coadd':
		in2 =  get_coadd_catalog_for_matcher(cat_type=MATCH_CAT2, inj_percent=inj2_percent,  realization=realization, tile=tile, mag_hdr=M_HDR2, mag_err_hdr=M_ERR_HDR2, inj=INJ2)

	print 'Matching or already matched:'
	print '', in1
	print '', in2

	if PLOT_COMPLETENESS: 
		### Rewrite `match_type` ###
		if RUN_TYPE is not None:
			__class1 = get_class(cat_type=MATCH_CAT1, inj_percent=inj1_percent, inj=inj1, suf='')
		if RUN_TYPE is None:
			__class1 = get_class(cat_type=MATCH_CAT1, inj_percent=inj1_percent, inj=inj1, suf='_1')
		__class2 = get_class(cat_type=MATCH_CAT2, inj_percent=inj2_percent, inj=inj2, suf='_2')
		__title1, __title2 = __class1.title_piece, __class2.title_piece,
		match_type = get_match_type(title_piece1=__title1, title_piece2=__title2)

	if PLOT_COMPLETENESS is False: match_type = MATCH_TYPE

        # !!!!! User may wish to edit directory structure. Output catalog name for STILTS #
	match_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, match_type, tile, realization, 'catalog_compare')	

	# Check for directory existence #
	if os.path.isdir(match_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(match_dir) + ' does not exist. \n Change directory structure in ms_plotter.matcher() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', match_dir, '...\n'
			os.makedirs(match_dir)

	# Filenames #
        __fn_match_1and2 = os.path.join(match_dir, tile+'_'+realization+'_'+str(match_type)+'_match1and2.csv')
	__fn_match_1not2 = os.path.join(match_dir, tile+'_'+realization+'_'+str(match_type)+'_match1not2.csv')
	__fn_match_2not1 = os.path.join(match_dir, tile+'_'+realization+'_'+str(match_type)+'_match2not1.csv')

	# Overwrite matched catalogs if one already exists? # 
        overwrite = False 

	# Check existence #
	if os.path.isfile(__fn_match_2not1) is False or (os.path.isfile(__fn_match_2not1) and overwrite):


		### Matching done in ms_matcher. Parameters in1, in2, out, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2, overwrite ###
		# !!!!! Ensure that path to ms_matcher is correct #
		subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher', in1, in2, __fn_match_1and2, __fn_match_1not2, __fn_match_2not1, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2])

	print ' Matched file: ', __fn_match_1and2, '\n'

        return __fn_match_1and2, __fn_match_1not2, __fn_match_2not1




def fof_matcher(realization, tile):
        """Get catalogs to analyze. Return FOF-analysed catalogs.

        Parameters
	----------
	realization (str) 

	tile (str) 

        Returns
	-------
	__fn_ok_1and2 OR __fn_rerun_1and2 (str) -- Complete filename for catalogs of type join=1and2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned. 

	__fn_ok_1not2 OR __fn_rerun_1not2 (str) -- Complete filename for catalogs of type join=1not2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.

	__fn_ok_2not1 OR __fn_rerun_2not1 (str) -- Complete filename for catalogs of type join=2not1. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.
        """

        ### Filenames for input catalogs used in ms_fof_matcher ###
        # MOF or SOF #
        if MOF:
                mof = os.path.join(BASEPATH, 'y3v02', tile, 'mof', tile+'_mof.fits')
                inj_mof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'mof', tile+'_mof.fits')
                fof = os.path.join(BASEPATH, 'y3v02', tile, 'mof', tile+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'mof', tile+'_fofslist.fits')
        if MOF is False:
                mof = os.path.join(BASEPATH, 'y3v02', tile, 'sof', tile+'_sof.fits')
                inj_mof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'sof', tile+'_sof.fits')
                fof = os.path.join(BASEPATH, 'y3v02', tile, 'sof', tile+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'sof', tile+'_fofslist.fits')
        # Coadds. Using i-band #
        coadd = os.path.join(BASEPATH, 'y3v02', tile, 'coadd', tile+'_i_cat.fits')
        inj_coadd = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization, tile, 'coadd', tile+'_i_cat.fits')


        ### Filenames for outputs of fof_matcher ###
	os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, 'fof_analysis_catalog_compare')
        # Repeated #
        inj_outdir = os.path.join(outdir, tile, realization)
        inj_outname = tile + '_' + realization


        ### Check directory existence. Make directories if not present. ###
        if os.path.isdir(inj_outdir) is False:
                os.makedirs(inj_outdir)
        if os.path.isdir(os.path.join(outdir, tile)) is False:
                os.makedirs(os.path.join(outdir, tile))


	### Create filenames for output catalogs created by ms_fof_matcher ###
        fofcoadd = os.path.join(outdir, tile, tile+ '_num-match_fof_coadd.csv')
        fofgroups = os.path.join(outdir, tile, tile+ 'fofgroups.csv')
	inj_fofcoadd = os.path.join(inj_outdir, inj_outname + '_num-match_inj_fof_inj_coadd.csv')
        inj_fofgroups = os.path.join(inj_outdir, inj_outname + '_inj_fofgroups.csv')
        origfof_injfof = os.path.join(inj_outdir, inj_outname + '_inj_fofgroup_fofgroup_match1and2.csv')
        ok = os.path.join(inj_outdir, inj_outname + '.ok')
        rerun = os.path.join(inj_outdir, inj_outname + '.rerun')
        ok_inj_mof = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof.csv')
        rerun_inj_mof = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof.csv')
        ok_mof = os.path.join(inj_outdir, inj_outname + '_ok_mof.csv')
        rerun_mof = os.path.join(inj_outdir, inj_outname + '_rerun_mof.csv')

        __fn_ok_1and2 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1and2.csv')
	__fn_ok_1not2 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1not2.csv')
        __fn_ok_2not1 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match2not1.csv')

        __fn_rerun_1and2 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1and2.csv')
	__fn_rerun_1not2 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1not2.csv')
	__fn_rerun_2not1 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match2not1.csv')

        # Output directory for files made in par.py #
        parpy_outdir = os.path.join(inj_outdir, inj_outname)


        # WARNING: May need to overwrite if matching was interupted #
        overwrite = False 

        ### Check file existence of last file made in fof_matcher ###
        if os.path.isfile(__fn_rerun_1and2) is False or (os.path.isfile(__fn_rerun_1and2) and overwrite):

                ### Run fof_matcher ###
                subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher', fof, inj_fof, mof, inj_mof, coadd, inj_coadd, parpy_outdir, fofcoadd, fofgroups, inj_fofcoadd, inj_fofgroups, origfof_injfof, ok, rerun, ok_inj_mof, rerun_inj_mof, ok_mof, rerun_mof, __fn_ok_1and2, __fn_rerun_1and2, __fn_ok_1not2, __fn_rerun_1not2, __fn_ok_2not1, __fn_ok_2not1])


	if RUN_TYPE == 'ok':
		return __fn_ok_1and2, __fn_ok_1not2, __fn_ok_2not1
	if RUN_TYPE == 'rerun':
		return __fn_rerun_1and2, __fn_rerun_1not2, __fn_rerun_2not1




def stack_tiles(realization):
	"""Concatenate catalogs with multiple tiles and fixed realization.
	
	Parameters
	----------
	realization (str)

	Returns
	-------
	fn_stack_1and2 (str) 
		Complete filename for stacked catalog of type join=1and2.

	fn_stack_1not2 (str) 
		Complete filename for stacked catalog of type join=1not2.

	fn_stack_2not1 (str) 
		Complete filename for stacked catalog of type join=2not1.

	"""

	# Directory for stacked catalog #
	stack_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, 'stack', realization)

	# Check directory existence and handle nonexistence #
	if os.path.isdir(stack_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(stack_dir) + ' does not exist. \n Change directory structure in ms_plotter. or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', stack_dir, '...\n'
			os.makedirs(stack_dir)

	# Filename for stacked catalogs #
	fn_end_1and2 = '_'.join(['stacked', realization, MATCH_TYPE, 'match1and2.csv'])
	__fn_stack_tiles_1and2 = os.path.join(stack_dir, fn_end_1and2)

	fn_end_1not2 = '_'.join(['stacked', realization, MATCH_TYPE, 'match1not2.csv'])
	__fn_stack_tiles_1not2 = os.path.join(stack_dir, fn_end_1not2)

	fn_end_2not1 = '_'.join(['stacked', realization, MATCH_TYPE, 'match2not1.csv'])
	__fn_stack_tiles_2not1 = os.path.join(stack_dir, fn_end_2not1)

	# Check if stacked tile catalog already exists #
	overwrite = False

	if os.path.isfile(fn_stack_2not1) and overwrite is False:
		print 'Stacked tile catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(__fn_stack_tiles_1and2)
		df1not2 = pd.read_csv(__fn_stack_tiles_1not2)
		df2not1 = pd.read_csv(__fn_stack_tiles_2not1)

	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(__fn_stack_tiles_2not1) is False or overwrite:

		all_fn_1and2, all_fn_1not2, all_fn_2not1 = [], [], []

		for t in ALL_TILES:

			if RUN_TYPE is None:
				fn_1and2, fn_1not2, fn_2not1 = matcher(realization=realization, tile=t, band=None, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)

			if RUN_TYPE is not None:
				fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=realization, tile=t) #TODO update with inj?

			all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		print 'Stacking tiles. ', len(all_fn_1and2), 'files ...'	
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		print 'Stacking complete ... \n'


		# Save stacked catalogs as DataFrame #
		df1and2.to_csv(__fn_stack_tiles_1and2, sep=','); df1not2.to_csv(__fn_stack_tiles_1not2, sep=','); df2not1.to_csv(__fn_stack_tiles_2not1, sep=',')
		print '-----> Saving stacked tile catalogs as ', __fn_stack_tiles_1and2 
		print '----->', __fn_stack_tiles_1not2 
		print '----->', __fn_stack_tiles_2not1 

	__number_of_stacked_tiles = len(ALL_TILES)

	return __fn_stack_tiles_1and2, __fn_stack_tiles_1not2, __fn_stack_tiles_2not1, __number_of_stacked_tiles




def stack_realizations(tile):
	"""Concatenate catalogs with multiple realizations and fixed tile.

	Parameters
	----------
	tile (str) -- One stacked realization catalog created per tile.

	Returns
	-------
	fn_stack_1and2 (str) -- Complete filename for stacked catalog of join=1and2.

	fn_stack_1not2 (str) -- Complete filename for stacked catalog of type join=1not2.

	fn_stack_2not1 (str) -- Complete filename for stacked catalog of type join=2not1. 

	len(ALL_REALIZATIONS) (int) -- Number of catalogs stacked.
	"""

	# Directory for stacked catalog #
	stack_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile, 'stack')

	# Check directory existence and handle nonexistence #
	if os.path.isdir(stack_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(stack_dir) + ' does not exist. \n Change directory structure in ms_plotter. or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', stack_dir, '...\n'
			os.makedirs(stack_dir)

	# Filename for stacked catalogs #
	__fn_stack_reals_1and2 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match1and2.csv')
	__fn_stack_reals_1not2 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match1not2.csv')
	__fn_stack_reals_2not1 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match2not1.csv')

	# Check if stacked realization catalog already exists #
	overwrite = False

	if os.path.isfile(__fn_stack_reals_2not1) and overwrite is False:
		print 'Stacked realization catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(__fn_stack_reals_1and2)
		df1not2 = pd.read_csv(__fn_stack_reals_1not2)
		df2not1 = pd.read_csv(__fn_stack_reals_2not1)


	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(__fn_stack_reals_2not1) is False or overwrite:

		all_fn_1and2, all_fn_1not2, all_fn_2not1 = [], [], []

		for r in ALL_REALIZATIONS:

			if RUN_TYPE is None:
				fn_1and2, fn_1not2, fn_2not1 = matcher(realization=r, tile=tile, band=None, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)

			if RUN_TYPE is not None:
				fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=r, tile=tile) #FIXME

			all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		print 'Stacking realizations. ', len(all_fn_1and2), 'files ...'
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		print 'Stacking complete ... \n'


		# Save stacked catalog as DataFrame #
		df1and2.to_csv(__fn_stack_reals_1and2, sep=','); df1not2.to_csv(__fn_stack_reals_1not2, sep=','); df2not1.to_csv(__fn_stack_reals_2not1, sep=',')
		print '-----> Saving stacked realization catalogs as ', __fn_stack_reals_1and2 
	        print '----->', __fn_stack_reals_1not2 
                print '----->', __fn_stack_reals_2not1 

	__number_of_stacked_realizations = len(ALL_REALIZATIONS)

	return __fn_stack_reals_1and2, __fn_stack_reals_1not2, __fn_stack_reals_2not1, __number_of_stacked_realizations




def get_percent_recovered(full_data, clean_data, inj_percent, tile, band, realization):
	"""TODO"""

	# Number of injections in a single truth catalog #
        if inj_percent == 10:
                __inj_single = 5000.0
        if inj_percent == 20:
                __inj_single = 10000.0
	#TODO mixed injection percent? __inj_single = 15000.0/2

	constant = 1.0
	if STACK_TILES: constant = constant * len(ALL_TILES)
	if STACK_REALIZATIONS: constant = constant * len(ALL_REALIZATIONS)

        # Total number of injections # 
        __inj_total = __inj_single*constant


	### Remove flags* from truth catalog ###
	# Get filename for truth catalog #
	if 'truth' in MATCH_CAT1:
		fn_truth_catalog = get_catalog_filename(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=INJ1_PERCENT, realization=realization, tile=tile, band=band) 
		__mag_hdr = M_HDR1[:-2]
		__flag_hdr = FLAGS_HDR1[:-2]
		__cm_flag_hdr = CM_FLAGS_HDR1[:-2] 
	if 'truth' in MATCH_CAT2:
		fn_truth_catalog = get_catalog_filename(cat_type=MATCH_CAT2, inj=INJ2, inj_percent=INJ2_PERCENT, realization=realization, tile=tile, band=band)
		__mag_hdr = M_HDR2[:-2]
		__flag_hdr = FLAGS_HDR2[:-2]
                __cm_flag_hdr = CM_FLAGS_HDR2[:-2]

	# Read catalog #
	hdul = fits.open(fn_truth_catalog)
	data = hdul[1].data
	truth_magnitude_griz = data[__mag_hdr]	
	# For ONE band #
	truth_magnitude = []
	for mag in truth_magnitude_griz:
		truth_magnitude.append(mag[BAND_INDEX[band]])
	# List --> array #
	truth_magnitude = np.array(truth_magnitude)
	flag = data[__flag_hdr]
        cm_flag = data[__cm_flag_hdr]
	# Remove flags* #
	__idx_good_true = np.where( (abs(truth_magnitude) != 9999.0) & (abs(truth_magnitude) != 99.0) & (abs(truth_magnitude) != 37.5) & (flag == 0) & (cm_flag == 0) )[0]
	__number_of_flags_in_truth_cat = len(truth_magnitude) - len(__idx_good_true)

	if 'truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2:
                __percent_recovered_flags_in = 100.0*len(full_data)/__inj_total
                # Not including flags #
                ###__percent_recovered_flags_rm = 100.0*len(clean_data)/__inj_total
		__percent_recovered_flags_rm = 100.0*len(clean_data)/(__inj_total - __number_of_flags_in_truth_cat)

	return __percent_recovered_flags_in, __percent_recovered_flags_rm, __number_of_flags_in_truth_cat 




def get_dataframe(realization, tile, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, inj1, inj2, inj1_percent, inj2_percent):
	"""Get pandas DataFrame and read it.

	Parameters
	----------
	inj1, inj2 (bool)

	inj1_percent, inj2_percent (int)

	mag_err_hdr1, mag_err_hdr2 (str)

	mag_hdr1, mag_hdr2 (str)

	realization (str)

	tile (str)

	Returns
	-------
		df1and2, df1not2, df2not1, __mag_hdr1, __mag_hdr2, __mag_err_hdr1, __mag_err_hdr2, number_of_stacked_tiles, number_of_stacked_realizations
	"""

	print 'Getting and reading DataFrame ... \n'

	__mag_hdr1, __mag_hdr2, __mag_err_hdr1, __mag_err_hdr2 = mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2

	#FIXME if PLOT_COMPLTENESS and STACK_TILES

	#TODO make one general fn_1and2

	if PLOT_COMPLETENESS is False:
		### Stack tiles ###
		for r in ALL_REALIZATIONS:
			if STACK_TILES and STACK_REALIZATIONS is False: 
				fn_stack_1and2, fn_stack_1not2, fn_stack_2not1, __number_of_stacked_tiles = stack_tiles(realization=r)
				__number_of_stacked_realizations = None


		### Stack realizations ###
		for t in ALL_TILES:
			if STACK_REALIZATIONS and STACK_TILES is False:
				fn_stack_1and2, fn_stack_1not2, fn_stack_2not1, __number_of_stacked_realizations = stack_realizations(tile=t)
				__number_of_stacked_tiles = None

		### Read stacked catalog ###
		if STACK_TILES or STACK_REALIZATIONS:
			# Get DataFrame for stacked catalogs #
			if STACK_REALIZATIONS or STACK_TILES:
				__df1and2 = pd.read_csv(fn_stack_1and2)
				__df1not2 = pd.read_csv(fn_stack_1not2)
				__df2not1 = pd.read_csv(fn_stack_2not1)

	
	if (STACK_REALIZATIONS is False and STACK_TILES is False) or (PLOT_COMPLETENESS and STACK_TILES):
		__number_of_stacked_tiles, __number_of_stacked_realizations = None, None
		# Filenames for catalogs #
		if RUN_TYPE is None:
			fn_1and2, fn_1not2, fn_2not1 = matcher(realization=realization, tile=tile, band=None, inj1=inj1, inj2=inj2, inj1_percent=inj1_percent, inj2_percent=inj2_percent)
		if RUN_TYPE is not None:
			fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=realization, tile=tile) #FIXME

		# Get DataFrame #
		__df1and2 = pd.read_csv(fn_1and2)
		__df1not2 = pd.read_csv(fn_1not2)
		__df2not1 = pd.read_csv(fn_2not1)


	### Handle star truth catalogs ###
	# Star truth catalogs matched then combined #
	if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
		print 'Adding new column to matched csv ...\n'
		# mag_star #
		new_hdr = 'mag_s'

		if MATCH_CAT1 == 'star_truth':
			suf = '_1'
			__mag_hdr1 = new_hdr
			#__mag_hdr1 set at the beginning of function
		if MATCH_CAT2 == 'star_truth':
			suf = '_2'
			__mag_hdr2 = new_hdr

		star_mag = get_star_truth_catalog_magnitude(df_1and2=__df1and2, suf=suf)
		# New header must be of the form {base}_x where x is a single character because of the way m_axlabel is created from m_hdr #
		__df1and2.insert(len(__df1and2.columns), new_hdr, star_mag)


	### Handle Y3 Gold catalogs ###
	# Y3 catalogs are matched then combined #
	if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
		print 'Adding new column to matched csv ...\n'

		# New header name #
		if 'star' in MATCH_CAT1 or 'star' in MATCH_CAT2:
			new_hdr = 'psf_mag_y'
		if 'star' not in MATCH_CAT1 and 'star' not in MATCH_CAT2:
			new_hdr = 'cm_mag_y'

		if 'y3_gold' in MATCH_CAT1:
			hdr = M_HDR1
			__mag_hdr1 = new_hdr
		if 'y3_gold' in MATCH_CAT2:
			hdr = M_HDR2
			__mag_hdr2 = new_hdr

		y3_gold_mag = get_y3_gold_catalog_magnitude(df_1and2=__df1and2, mag_hdr=hdr)
		# Add new column to df #
		__df1and2.insert(len(__df1and2.columns), new_hdr, y3_gold_mag)


	### Handle coadd catalogs. New column has been added with name 'mag_c'. Catalog combined then matched so has suffix (unlike star) #
	if MATCH_CAT1 == 'coadd':
		__mag_hdr1 = 'mag_c_1'
		__mag_err_hdr1 = 'mag_err_c_1'
	if MATCH_CAT2 == 'coadd':
		__mag_hdr2 = 'mag_c_2'
		__mag_err_hdr2 = 'mag_err_c_2'

	#TODO numberOfStacked...
	return __df1and2, __df1not2, __df2not1, __mag_hdr1, __mag_hdr2, __mag_err_hdr1, __mag_err_hdr2, __number_of_stacked_tiles, __number_of_stacked_realizations




def make_plots(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2):
	"""Makes plots.

	Parameters
	----------
	mag_hdr1, mag_hdr2 (str)
		Headers for magnitude. May be altered, hence passed as parameters.

	mag_err_hdr1, mag_err_hdr2 (str)
		Headers for magnitude error. May be altered if new columns are added (#FIXME), hence passed as parameters.

	Returns
	-------
		0
	"""

	if STACK_TILES is False: list_of_tiles_for_loop = ALL_TILES
        if STACK_REALIZATIONS is False: list_of_realizations_for_loop = ALL_REALIZATIONS

	if STACK_TILES and STACK_REALIZATIONS is False: list_of_tiles_for_loop = ['stack']
	if STACK_REALIZATIONS and STACK_TILES is False: list_of_realizations_for_loop = ['stack']

	for t in list_of_tiles_for_loop:
                for r in list_of_realizations_for_loop:

			# Stacked completeness plots do NOT use stacked tile catalog ... #
			if PLOT_COMPLETENESS is False: 
				df1and2, df1not2, df2not1, magHdr1, magHdr2, magErrHdr1, magErrHdr2, numStackTile, numStackReal = get_dataframe(realization=r, tile=t, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)

			#if PLOT_COMPLETENESS:
			#TODO remove last two outputs from get_df
			numStackTile, numStackReal = len(ALL_TILES), len(ALL_REALIZATIONS)

			fnFlagLog, fnMagErrLog, fnMainLog, fnColorLog, fnMagDiffOutliersLog, fnMagCompletenessLog, fnGaussAperLog = get_log_filenames(tile=t, realization=r)
			#fn_flag_log, fn_mag_err_log, fn_main_log, __fn_color_log, __fn_magL_outliers, fnCompletenessLog, fnAper = get_log_filenames(tile=t, realization=r)
			
			fnFlagLog, fnMagErrLog, fnMainLog, fnColorLog, fnMagDiffOutliersLog, fnMagCompletenessLog = write_log_file_headers(fn_main_log=fnMainLog, fn_mag_err_log=fnMagErrLog, fn_flag_log=fnFlagLog, fn_color_log=fnColorLog, fn_mag_diff_outliers_log=fnMagDiffOutliersLog, fn_mag_completeness_log=fnMagCompletenessLog)

			# Name for plt.savefig() #
                        fnPlot = get_plot_save_filename(realization=r, tile=t)

			title = get_plot_suptitle(realization=r, tile=t, number_of_stacked_realizations=numStackReal, number_of_stacked_tiles=numStackTile)

			### Region files ###
			#TODO PLOT_COMPLETENESS df1and2 created in another function --> cannot call write_to_region_files()
			if PLOT_COMPLETENESS is False:
				if MAKE_REG: write_to_region_files(df_1and2=df1and2, df_1not2=df1not2, df_2not1=df2not1, realization=r, tile=t)	


			### Call plotting functions ###	
			if PLOT_MAG and PLOT_COMPLETENESS is False:
				magnitude_difference_subplotter(df_1and2=df1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, fn_plot=fnPlot, plot_title=title, realization=r, tile=t, fn_mag_err_log=fnMagErrLog, fn_main_log=fnMainLog, fn_flag_log=fnFlagLog, fn_mag_diff_outliers_log=fnMagDiffOutliersLog) 

			if PLOT_COLOR and PLOT_COMPLETENESS is False:
				for b in ALL_BANDS[:-1]:
					color_subplotter(band=b, df_1and2=df1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, fn_flag_log=fnFlagLog, plot_title=title, fn_plot=fnPlot, fn_color_log=fnColorLog)

			### Completeness plots (for truth catalogs vs injected {mof/sof/coadd/...}) ###
			if ('truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2) and (INJ1 and INJ2) and PLOT_COMPLETENESS:
				if PLOT_MAG and STACK_TILES is False:
					magnitude_completeness_subplotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, plot_title=title, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog)
				if PLOT_MAG and STACK_TILES:
					stacked_magnitude_completeness_subplotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, plot_title=title, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog)

				if PLOT_COLOR:
					for b in ALL_BANDS[:-1]: 
						color_completeness_subplotter(band=b, df_1and2=df1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, plot_title=title, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog)


			if PLOT_FLUX:
				normalized_flux_histogram_subplotter(df_1and2=df1and2, flux_hdr1=CM_FLUX_HDR1, flux_hdr2=CM_FLUX_HDR2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=title, fn_plot=fnPlot, tile=t, realization=r, fn_gauss_aper_log=fnGaussAperLog)


	return 0




def get_coadd_catalog_for_matcher(cat_type, inj_percent, inj, realization, mag_hdr, mag_err_hdr, tile):
	"""Make FITS file that includes a column of form '(m_g, m_r, m_i, m_z)' where m is magnitude. Column will be added to '..._i_cat.fits'. This will be used in matcher(). Relies on directory structure /`OUTDIR`/outputs/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/catalog_compare/

	Parameters
	----------
	cat_type (str)
		Catalog type. This is set by `MATCH_CAT1` or `MATCH_CAT2`. 

	inj (bool)
		Is the catalog Balrog-injected? Determined by `INJ1` and `INJ2`.

	inj_percent (int)
		Injection percent? Currently '10' and '20' are allowed values referring to 5,000 and 10,000 injected objects, respectively.

	realization (str)

	mag_hdr (str)
		Header for magnitude. Headers refer to columns in the matched catalog.

	mag_err_hdr (str)
		Header for magnitude #FIXME error. Headers refer to columns in the matched catalog.

	tile (str)

	Returns
	-------
	__fn_coadd_for_matcher (str)
		Complete filename for catalog with added column. Is a FITS file.
	"""

	dir_new = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, 'catalog_compare')
	if os.path.isdir(dir_new) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(dir_new) + ' does not exist...')
		if NO_DIR_MAKE:
			os.makedirs(dir_new)

	__fn_coadd_for_matcher = os.path.join(dir_new, str(tile) + '_i_cat_combo.fits')

	# Check if new coadd catalog has already been created #
	if os.path.isfile(__fn_coadd_for_matcher):
		print 'New coadd catalog already exists ...\n'


	if os.path.isfile(__fn_coadd_for_matcher) is False:	
		print 'Adding a column to i-band coadd catalog. Will take a moment ...\n'

		# Get list of filenames #
		fn_griz = []
		for b in ALL_BANDS:
			fn_griz.append(cat_type=cat_type, inj=inj, inj_percent=inj_percent, realization=realization, tile=tile, band=b)
		fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z = fn_griz

		# Get coadd magnitude (mag_c) and magnitude error to be of form '(m_g, m_r, m_i, m_z)'. Recall that this is a string #
		mag_c, mag_err_c = get_coadd_catalog_magnitude_and_magnitude_error(fn_coadd_cat_g=fn_coadd_cat_g, fn_coadd_cat_r=fn_coadd_cat_r, fn_coadd_cat_i=fn_coadd_cat_i, fn_coadd_cat_z=fn_coadd_cat_z, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr)
 
	       # Create new table #
		mag_c = Column(mag_c, name='mag_c')
		mag_err_c = Column(mag_err_c, name='mag_err_c')

		# Add new table to i-band coadd catalog #
		table = Table.read(fn_coadd_cat_i)
		table.add_column(mag_c, index=0)
       		table.add_column(mag_err_c, index=1)

		# Save new table as FITS #
		table.write(__fn_coadd_for_matcher)

	return __fn_coadd_for_matcher 




#TODO accept idx_good as input param? Will only be used for df_match. 
def write_to_region_files(df_1and2, df_1not2, df_2not1, realization, tile):
	"""Make DS9 region files for catalogs matched via join=1and2, join=1not2, and join=2not1 (join type is STILTS parameter set in ms_matcher or ms_fof_matcher).

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		Catalog for matches via join=1and2 between `MATCH_CAT1` and `MATCH_CAT2`.

	df_1not2 (pandas DataFrame)
		Catalog for matches via join=1not2 between `MATCH_CAT1` and `MATCH_CAT2`.

	df_2not1 (pandas DataFrame)
		Catalog for matches via join=2not1 between `MATCH_CAT1` and `MATCH_CAT2`.

	realization (str)

	tile (str)

	Returns
	-------
	fnRegion1and2 (str)
		Complete filename for region file containing regions pertaining to objects in `df_1and2`.

	fnRegion1not2 (str)
		Complete filename for region file containing regions pertaining to objects in `df_1not2`.

	fn_reg_2not1 (str)
		Complete filename for region file containing regions pertaining to objects in `df_2not1`.
	"""

	### Get filenames and open files ###
	fnRegion1and2, fnRegion1not2, fnRegion2not1 = get_region_filenames(tile=tile, realization=realization)

	overwrite = False
	
	if os.path.isfile(fnRegion2not1) and overwrite is False:
		print 'Region files already exist. Not overwriting ...'
	

	if os.path.isfile(fnRegion2not1) is False or (os.path.isfile(fnRegion2not1) and overwrite):
		fd_match = open(fnRegion1and2, 'w'); fd_1not2 = open(fnRegion1not2, 'w'); fd_2not1 = open(fnRegion2not1, 'w')
		# Write coordinate system #
		fd_match.write('J2000 \n'); fd_1not2.write('J2000 \n'), fd_2not1.write('J2000 \n')

		# Handle matched catalog #
		if RUN_TYPE is None:
			ra1 = RA_HDR1 + str('_1'); dec1 = DEC_HDR1 + str('_1')
			ra2 = RA_HDR2 + str('_2'); dec2 = DEC_HDR2 + str('_2')
		if RUN_TYPE is not None:
			# MOF or SOF catalogs #
			ra1 = 'ra'; dec1 = 'dec'
			ra2 = 'ra_2'; dec2 = 'dec_2' 

		### Get position. Arbitrarily using MATCH_CAT1 for RA and DEC ###
		ra_match, dec_match = df_1and2[ra2], df_1and2[dec2] 
		ra_1not2, dec_1not2 = df_1not2[ra1], df_1not2[dec1]
		ra_2not1, dec_2not1 = df_2not1[ra2], df_2not1[dec2]

		### Write to region file for matched catalog. Units are arcseconds. ###
		# Coadds allow for elliptical regions #
		if 'coadd' in (MATCH_CAT1, MATCH_CAT2) or 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
			### Get semimajor and semiminor axes (a and b, respectively) and orientation. Coadds and Y3 Gold have these values. ###
			a_match, b_match = df_1and2[MAJOR_AX_HDR1], df_1and2[MINOR_AX_HDR1]
			a_1not2, b_1not2 = df_1not2[MAJOR_AX_HDR1], df_1not2[MINOR_AX_HDR1]
			a_2not1, b_2not2 = df_2not1[MAJOR_AX_HDR2], df_2not1[MINOR_AX_HDR2]
			orientation_match, orientation_1not2, orientation_2not1 = df_1and2[ANGLE1], df_1not2[ANGLE1], df_2not1[ANGLE1]
			
			for i in np.arange(0, len(ra_match)):
				fd_match.write('ellipse ' + str(ra_match[i]) + ' ' + str(dec_match[i]) + ' ' + str(a_match[i]) + '" ' + str(b_match[i]) + '" ' + str(90+orientation[i]) + ' #color=green width=3 \n')
			for i in np.arange(0, len(ra_1not2)):
				fd_1not2.write('ellipse ' + str(ra_1not2[i]) + ' ' + str(dec_1not2[i]) + ' ' + str(a_1not2[i]) + '" ' + str(b_1not2[i]) + '" #color=yellow width=3 \n')
			for i in np.arange(0, len(ra_2not1)):
				fd_2not1.write('ellipse ' + str(ra_2not1[i]) + ' ' + str(dec_2not1[i]) + ' ' + str(a_2not1[i]) + '" ' + str(b_2not1[i]) + '" #color=blue width=3 \n')


		# Non-coadd catalogs allow for circular regions #
		if MATCH_CAT1 != 'coadd' and MATCH_CAT2 != 'coadd':
			size_sq_match = df_1and2[CM_T_HDR1]
			size_sq_1not2 = df_1not2[CM_T_HDR1]
			size_sq_2not1 = df_2not1[CM_T_HDR2]

			# Note: Can use a typical radius of 2 arcsec #
			for i in np.arange(0, len(ra_match)):
				if size_sq_match[i] > 0:# and np.isnan(size_sq_match[i]) is False:
					fd_match.write('circle ' + str(ra_match[i]) + ' ' + str(dec_match[i]) + ' ' + str(size_sq_match[i]**0.5) + '" #color=green width=3 \n')
			for i in np.arange(0, len(ra_1not2)):
				if size_sq_1not2[i] > 0: # and np.isnan(size_sq_1not2[i]) is False:
					fd_1not2.write('circle ' + str(ra_1not2[i]) + ' ' + str(dec_1not2[i]) + ' ' + str(size_sq_1not2[i]**0.5) + '" #color=yellow width=3 \n')
			for i in np.arange(0, len(ra_2not1)):
				if size_sq_2not1[i] > 0: # and np.isnan(size_sq_2not1[i]) is False:
					fd_2not1.write('circle ' + str(ra_2not1[i]) + ' ' + str(dec_2not1[i]) + ' ' + str(size_sq_2not1[i]**0.5) + '" #color=blue width=3 \n')

		# Close files #
		fd_match.close(); fd_1not2.close(); fd_2not1.close()

	print '-----> Saving region files as: ', fnRegion1and2
	print ' -----> ', fnRegion1not2 
	print ' ----->', fnRegion2not1 

	return fnRegion1and2, fnRegion1not2, fnRegion2not1 




################################################################### Run script. 0 returned when complete. ###################################################################


YLOOP = False 
### !!!!! Run once. Log files are closed once 0 is returned. ###
if YLOOP is False:
	print make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)

### Loop over vertical axis limits. Suggestions: for normalized plot with star truth catalog use y=[3, 10], for normalized plot with galaxy truth catalog use y=[3, 20]. For non-normalized plot with star truth catalog or galaxy truth catalog use y=[0.5, None]. ###
# !!!!! Loop over vertical axis limits? #

if NORMALIZE:
	ylist = [10, 20, None]
if NORMALIZE is False:
	ylist = [0.5, 1, 3, None]
if YLOOP:
	for y in ylist: 
		if y is None:
			MAG_YLOW, MAG_YHIGH = None, None
		if y is not None:
			MAG_YLOW, MAG_YHIGH = -1*y, y

		# Must pass constants as parameters here because used for plot name and manipulated. #
		print make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



### !!!!! Loop over possible colorbars? ###
CBAR_LOOP = False
if CBAR_LOOP:
	# Possible combinations for colorbars. Only one True allowed at a time and NORMALIZE and HEXBIN must both be False. #
	cbar_bools_list =[[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False]]
	for cbar_bools in cbar_bools_list:
		# Reset constants #
		#FIXME no longer sure if this works and renames plots correctly
		CM_T_COLORBAR, CM_T_ERR_COLORBAR, NORMALIZE, HEXBIN = cbar_bools
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



RUN_TYPE_LOOP = False 
if RUN_TYPE_LOOP:
	run_type_list = [None, 'ok', 'rerun']
	for run_type in run_type_list:
		RUN_TYPE = run_type
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)

