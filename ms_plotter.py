"""Creates various plots for Balrog validation testing.


To run: $python ms_plotter.py base_path_to_catalogs output_directory realization tile
Example: $python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /Users/mspletts/BalVal/ all DES2247-4414

Relies on ms_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` with the correct path to ms_matcher.
FOF analysis relies on ms_fof_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher` with the correct path to ms_fof_matcher.

Plot attributes are specified with constants (many of them booleans) at the beginning of this script. See also the table in README.md. 
Constants that the user may wish to change are indicated by: # !!!!! {description and/or warnings} #. For example, user may wish to set `PRINTOUTS = False` or comment out `NOTICE`.

# Comments are ABOVE the code they correspond to (with the exception of FIXMEs and TODOs). In general, function parameters are named_with_underscores, variables defined within functions (and that should only be accessed by said function) are __named_with_leading_underscores, and variable namesWithoutSpaces are defined by a function via nameWithoutSpace = function(param1=x). # 

Megan Splettstoesser mspletts@fnal.gov"""

#TODO Update docstrings

# `astropy` is needed only if analyzing a coadd catalog or making a completeness plot #
from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
from scipy import stats
# For Python environment with `corner`, run: $source activate des17a or $source activate des18a #
import corner
import csv
import fileinput
import math
import matplotlib.pyplot as plt
import matplotlib #FIXME
import numpy as np
import os
import pandas as pd
import subprocess
import sys




### Command line args ###
# Catch error from inadequate number of command line args #
if len(sys.argv) != 5:
        sys.exit("Args: basepath (location of catalogs), output directory, realizations (can be 'all', None (non-injected catalog), a list of form: real1,real2,...), tiles (can be 'all', a file, or a list of form: tile1,tile2,...) \n")
# Convert `realizations` and `tiles` into lists (may be one-element list). Note 'None' entered at command line is interpreted as a str #
BASEPATH, OUTDIR, realizations, tiles = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4].split(',')


### For directory structure ###
# '/data/des71.a/data/kuropat/des2247-4414_sof/' --> 'des2247-4414_sof' #
BALROG_RUN = BASEPATH[BASEPATH[:-1].rfind('/')+1:-1]
# Rewrite #
if BALROG_RUN == 'Balrog': BALROG_RUN = 'TAMU_Balrog'


### Constants needed to loop over filters, realizations, and tiles ###
# Filters #
ALL_FILTERS = [ 'g', 'r', 'i', 'z' ]

# Realizations #
if realizations[0] != 'all': ALL_REALIZATIONS = realizations
# !!!!! Number of realizations depends on the tile. User may need to manually set ALL_REALIZATIONS or supply a list at command line. #
if realizations[0] == 'all': ALL_REALIZATIONS = [ '0', '1', '2' ] # !!!!!

# Tiles #
if tiles[0] != 'all': ALL_TILES = tiles
# !!!!! User may need to manually set ALL_TILES or supply a list at command line. #
if tiles[0] == 'all': ALL_TILES = ['DES0347-5540', 'DES2329-5622', 'DES2357-6456'] # !!!!!
# Accept .dat file of tile names at command line #
if '.dat' in tiles[0]:
	ALL_TILES = []
	for line in fileinput.input(tiles):
		# Get rid of newline character \n #
		ALL_TILES.append(line[:-1])





################################################################### Specify plot and catalog attributes ###################################################################

### Catalog attributes ###
# !!!!! Allowed values: y3_gold_2_0, y3_gold_2_2, sof, mof, star_truth, gal_truth, coadd. Both can be 'sof' and both can be 'mof' if INJ1_10PERCENT and INJ2_10PERCENT are different. Note that truth catalogs always have INJ=True. #
MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'sof'
# !!!!! Booleans. Examine injected catalogs? `INJ1_10PERCENT` `INJ2_10PERCENT` refer to 10% injections #
# Are catalogs injected? #
INJ1, INJ2 = True, True
# If so, what % injection? If not, which directory to get base catalog? #
INJ1_10PERCENT, INJ2_10PERCENT = True, True 
INJ1_20PERCENT, INJ2_20PERCENT = False, False

# !!!!! What to do with the plot? #
SAVE_PLOT = False
SHOW_PLOT = True

### Plot type ###
HIST_2D = False
CORNER_HIST_2D = True
HEXBIN = False
SCATTER = False
# `_COLORBAR` cannot be True if NORMALIZE is True #
CM_T_S2N_COLORBAR = False
CM_T_ERR_COLORBAR = False
CM_T_COLORBAR = False
BIN_CM_T_S2N = False
# If True, normalizes plot to 1sigma_mag. If True, PLOT_1SIG must be True else errors will not be computed and normalization cannot be performed #
NORMALIZE = False


# !!!!! Plot 1sigma_mag curve? Must be True if NORMALIZE is True. If NORMALIZE is True will also plot the 68th percentile of the data in each error bin. #
PLOT_1SIG = True

# Plot colors not magnitudes?
PLOT_COLOR = False 

PLOT_COMPLETENESS = False

# 1D histogram #
PLOT_FLUX_HIST = True

# !!!!! If `True` plots x1 versus (x1-x2). If `False` plots x1 versus x2. #
PLOT_DELTA_VAX = True


# !!!!! Limits for the vertical axis. 'None' is an allowed value and will result in default scaling #
YLOW, YHIGH = None, None

CENTER_ERR_ABT_ZERO = True 

# Swap horizontal axis? Default is magnitude1. Matching script ms_matcher determines which catalog is 1 and which is 2. Generally SWAP_HAX does not need to be changed unless the truth catalog values are not on the horizontal axis. #
SWAP_HAX = False



### Handle nonsensical combinations ###
# Always non-injections #
if 'y3_gold' in MATCH_CAT1:
	INJ1_10PERCENT, INJ1_20PERCENT = False, False
if 'y3_gold' in MATCH_CAT2:
        INJ2_10PERCENT, INJ2_20PERCENT = False, False
if realizations[0] == 'None':
	INJ1_10PERCENT, INJ2_10PERCENT, INJ1_20PERCENT, INJ2_20PERCENT = False, False, False, False


# !!!!! Only used if MATCH_CAT1 or MATCH_CAT2 is 'y3_gold'. If False, SOF catalogs exists in subdirectories of BASEPATH #
Y3_MOF = None
if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
	if 'mof' in (MATCH_CAT1, MATCH_CAT2):
		Y3_MOF = True
	if 'sof' in (MATCH_CAT1, MATCH_CAT2):
		Y3_MOF = False
 

# !!!!! Must be used with realization=all at command line #
STACK_REALIZATIONS = False
STACK_TILES = False


# !!!!! Make 2x2 subplots of each griz filter? Or make individual plots? #
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
        INJ1_10PERCENT, INJ2_10PERCENT = True, False


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
SHOW_FLAG_TYPE = False
# Use quality cuts introduced by Eric Huff? Link: https://github.com/sweverett/Balrog-GalSim/blob/master/plots/balrog_catalog_tests.py. Can only be performed if catalog has all the necessary headers: cm_s2n_r, cm_T, cm_T_err, and psfrec_T. #
EH_CUTS = False



# Catch errors from plot attributes and command line arguments #
def catch_error():
	"""Identify errors."""

	msg = None

	if STACK_REALIZATIONS and realizations[0] != 'all': msg = 'STACK_REALIZATIONS is True must be used with realization = all'
	if STACK_TILES and ('.dat' not in tiles[0] and len(tiles) == 1): msg = 'STACK_TILES must be used with multiple tiles'
	if YLOW is not None and YHIGH is not None and YHIGH == YLOW: msg = 'YLOW and YHIGH cannot be equal'
	if NORMALIZE and PLOT_1SIG is False: msg = 'If NORMALIZE is True so must be PLOT_1SIG'

	# Colorbar errors #
	cbar_counter = 0
	if HEXBIN: cbar_counter += 1
	if CM_T_S2N_COLORBAR: cbar_counter += 1
	if CM_T_ERR_COLORBAR: cbar_counter += 1
	if HIST_2D: cbar_counter += 1
	if CORNER_HIST_2D: cbar_counter += 1
	if CM_T_COLORBAR: cbar_counter += 1 
	if BIN_CM_T_S2N: cbar_counter += 1
	#TODO add NORMALIZE
	if cbar_counter > 1: msg = 'Only one colorbar can be used. Edit HEXBIN, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, HIST_2D, CM_T_COLORBAR, BIN_CM_T_S2N'

	if INJ1_10PERCENT is False and INJ2_10PERCENT is False and INJ1_20PERCENT is False and INJ2_20PERCENT is False and realizations[0] != 'None': msg = 'If INJ1_10PERCENT and INJ2_10PERCENT are False realizations must be None at cmd line'

	if ('truth' in MATCH_CAT1 and INJ1_10PERCENT is False and INJ1_20PERCENT is False) or ('truth' in MATCH_CAT2 and INJ2_10PERCENT is False and INJ2_20PERCENT is False): msg = 'Truth catalogs are injected'

	plt_type_counter = 0
	if PLOT_COMPLETENESS: plt_type_counter += 1
	if PLOT_FLUX_HIST: plt_type_counter += 1
	if plt_type_counter > 1: msg = 'Pick one plot type...' #TODO

	#TODO not normalize and plot_color simult.
	return msg

if catch_error() is not None:
        sys.exit('Error: ' + catch_error())



# !!!!! Check that plots will not be overwritten, etc # #FIXME update this
NOTICE = raw_input(' \n !! CHECK BEFORE RUNNING !! \n Using catalogs -- ' + str(MATCH_CAT1) + ' & ' + str(MATCH_CAT2) + '\n Injected catalogs -- ' + str(INJ1_10PERCENT) + ' & ' + str(INJ2_10PERCENT) + ' \n Plotting color: ' + str(PLOT_COLOR) + ' \n Save plot(s) -- ' + str(SAVE_PLOT) + '\n Showing plot(s) -- ' + str(SHOW_PLOT) + '\n Normalize plot(s) -- ' + str(NORMALIZE) + '\n Stacking realizations -- ' + str(STACK_REALIZATIONS) + '\n Hexbin -- ' + str(HEXBIN) + '\n cm_T colorbar -- ' + str(CM_T_COLORBAR) + '\n cm_T_err colorbar -- ' + str(CM_T_ERR_COLORBAR) + '\n cm_T_s2n -- ' + str(CM_T_S2N_COLORBAR) + '\n Plot limits -- ' + str(YLOW) + ', ' + str(YHIGH) + '\n Plotting 1sigma curve -- ' + str(PLOT_1SIG) +'\n Plotting flagged objects -- ' + str(PLOT_FLAGGED_OBJS) + '\n Print flags and flag types -- ' + str(SHOW_FLAG_TYPE) + '\n Logging flags -- ' + str(LOG_FLAGS) + '\n --> Press enter to proceed, control+c to stop...\n')









################################################################### Store catalog information ###################################################################

class CoaddCat():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{filter}_cat.fits. There is a separate catalog for each filter."""

        # Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_10percent, inj_20percent, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj_10percent (bool)
			If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog. 

		inj_20percent (bool)
			If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog. 

		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# For plot title #
		if inj_10percent:
			self.title_piece = '10% Inj Coadd Cat'
		if inj_20percent:
			self.title_piece = '20% Inj Coadd Cat'
		if inj_10percent is False and inj_20percent is False:
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
	def __init__(self, inj_10percent, inj_20percent, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj_10percent (bool)
			If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

		inj_20percent (bool)
			If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# `inj` forced True for truth catalogs #
		if inj_10percent:
			self.title_piece = '10% Inj Gal Truth Cat'
		if inj_20percent:
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
        def __init__(self, inj_10percent, inj_20percent, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_10percent (bool)
                        If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                inj_20percent (bool)
                        If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

                if inj_10percent:
                        self.title_piece = '10% Inj SOF Cat'
		if inj_20percent:
			self.title_piece = '20% Inj SOF Cat'
                if inj_10percent is False and inj_20percent is False:
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
        def __init__(self, inj_10percent, inj_20percent, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_10percent (bool)
                        If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                inj_20percent (bool)
                        If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

		# For plot title #
		if inj_10percent:
			self.title_piece = '10% Inj MOF Cat'
		if inj_20percent:
			self.title_piece = '20% Inj MOF Cat'
		if inj_10percent is False and inj_20percent is False:
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
	def __init__(self, inj_10percent, inj_20percent, suf):
                """Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_10percent (bool)
                        If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                inj_20percent (bool)
                        If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """	

		# `inj` forced True for truth catalogs #	
		if inj_10percent:
			self.title_piece = '10% Inj Star Truth Cat'
		if inj_20percent:
			self.title_piece = '20% Inj Star Truth Cat'
		self.axlabel = 'true'
		# Magnitude #
		self.mag_hdr = 'g_Corr' + str(suf) 
		self.mag_axlabel = 'mag_true' 
		# For error calculation #
                # Is of form 'PSF_MAG_ERR_{filter}' + str(suf) #
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
	def __init__(self, inj_10percent, inj_20percent, suf):
                """Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

                Parameters
                ----------
                inj_10percent (bool)
                        If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                inj_20percent (bool)
                        If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

                suf (str)
                        Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
                """

		inj = False

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

		# Note: filter dependent headers -- MAG, MAG_ERR, FLUX, FLUX_COV, PSF_FLAGS #

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









def get_class(cat_type, inj_10percent, inj_20percent, suf):
        """Get the appropriate class for the catalog type.

        Parameters
	----------
	cat_type (str)
                Catalog type. Set by `MATCH_CAT1` or `MATCH_CAT2`.

	inj_10percent (bool)
		If True refers to 10% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

	inj_20percent (bool)
		If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog.

	suf (str)
		Refers to the order in which `cat_type` was matched (via join=1and2) in ms_matcher (order set by `in1` and `in2` in STILTS script). Allowed values: '_1' '_2'.

        Returns
	-------
	_class_ (class)
		Points to the appropriate class and class constants.
        """

        if cat_type == 'gal_truth':
                __class = GalTruthCat(inj_10percent=inj_10percent, inj_20percent=inj_20percent, suf=suf)

        if cat_type == 'mof':
                __class = MOFCat(inj_10percent=inj_10percent, inj_20percent=inj_20percent, suf=suf)

        if cat_type == 'star_truth':
                __class = StarTruthCat(inj_10percent=inj_10percent, inj_20percent=inj_20percent, suf=suf)

        if cat_type == 'sof':
                __class = SOFCat(inj_10percent=inj_10percent, inj_20percent=inj_20percent, suf=suf)

        if cat_type == 'coadd':
                __class = CoaddCat(inj_10percent=inj_10percent, inj_20percent=inj_20percent, suf=suf)

	if 'y3_gold' in cat_type:
		__class = Y3Gold(inj_10percent=inj_10percent, suf=suf)


        return __class









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

        return '_'.join([title_piece1.lower().replace(' ', '_'), title_piece2.lower().replace(' ', '_')]) 









def get_log_file_names(tile_name, realization_number):
        """Generate names for log files. Relies on directory structure: /`OUTDIR`/log_files/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/log_files/.

        Parameters
	----------
	tile_name (str)

	realization_number (str)

        Returns
	-------
	fn1 (str)
		Complete filename for flag log file.

	fn2 (str)
		Complete filename for error calculation log file.

	fn3 (str)
		Complete filename for log file containting number of objects matched, number of objects flagged, number of objects within 1sigma_mag, etc.

	fn4 (str)

	fn5 (str)

	fn6 (str)
        """

        # !!!!! User may wish to edit directory structure #
	log_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'log_files') 
        if RUN_TYPE is not None:
		log_dir = os.path.join(log_dir, 'fof_analysis')

	### Check for directory existence ###
        if os.path.isdir(log_dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory ' + str(log_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_log_file_names() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', log_dir, '...\n'
                        os.makedirs(log_dir)

	# Repeated #
	fn_rep = '_'.join([tile_name, realization_number, MATCH_TYPE])
	fn1 = os.path.join(log_dir, 'flag_log_'+fn_rep+'.log')
	fn2 = os.path.join(log_dir, 'err_calc_'+fn_rep+'.log')
	fn3 = os.path.join(log_dir, 'num_flags_num_1sig_'+fn_rep+'.log')
	fn4 = os.path.join(log_dir, 'color_plot_'+fn_rep+'.log')
	fn5 = os.path.join(log_dir, 'outlier_delta_magnitudes_'+fn_rep+'.log')
	fn6 = os.path.join(log_dir, 'completeness_'+fn_rep+'.log')

        if RUN_TYPE is not None:
		fn1 = '_'.join([fn1[:-4], RUN_TYPE, fn1[-4:]])
		fn2 = '_'.join([fn2[:-4], RUN_TYPE, fn2[-4:]])
		fn3 = '_'.join([fn3[:-4], RUN_TYPE, fn3[-4:]])
		fn4 = '_'.join([fn4[:-4], RUN_TYPE, fn4[-4:]])
		fn5 = '_'.join([fn5[:-4], RUN_TYPE, fn5[-4:]])
		fn6 = '_'.join([fn6[:-4], RUN_TYPE, fn6[-4:]])

        print '-----> Saving log file for flags as: ', fn1, '\n'
        print '-----> Saving log file for magnitude and error bins as: ', fn2, '\n'
        print '-----> Saving log file for number of flags, number of objects within 1sigma_mag as: ', fn3, '\n'
	print '-----> Saving color plot log as: ', fn4, '\n'
	print '-----> Saving outlier DeltaMagnitude log as: ', fn5, '\n'
	print '-----> Saving completeness plot log as: ', fn6, '\n'
	
	return fn1, fn2, fn3, fn4, fn5, fn6










def get_reg_names(tile_name, realization_number):
	"""Generate names for region files of different join types (join types specified in ms_matcher or ms_fof_matcher). Relies on directory structure `/OUTDIR/outputs/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/region_files/`

        Parameters
	----------
	tile_name (str)

	realization_number (str)

        Returns
	-------
	fn1 (str)
		Complete filename for region file containing regions with join=1and2.

	fn2 (str)
		Complete filename for region file containing regions with join=1not2.

	fn3 (str)
		Complete filename for region file containing regions with join=2not1.
        """

        # !!!!! User may wish to edit directory structure #
	reg_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'region_files')
        if RUN_TYPE is not None:
		reg_dir = os.path.join(reg_dir, 'fof_analysis')

	### Check for directory existence ###
        if os.path.isdir(reg_dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory ' + str(reg_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_reg_names() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', reg_dir, '...\n'
                        os.makedirs(reg_dir)

	# Repeated #
        fn_rep = '_'.join([tile_name, realization_number, MATCH_TYPE])
	fn1 = os.path.join(reg_dir, fn_rep+'_match1and2.reg')
	fn2 = os.path.join(reg_dir, fn_rep+'_match1not2.reg')
	fn3 = os.path.join(reg_dir, fn_rep+'_match2not1.reg')


	if RUN_TYPE is not None:
		fn1 = fn1[:-15] + '_' + str(RUN_TYPE) + fn1[-15:]
		fn2 = fn2[:-15] + '_' + str(RUN_TYPE) + fn2[-15:]
		fn3 = fn3[:-15] + '_' + str(RUN_TYPE) + fn3[-15:]


	return fn1, fn2, fn3







################################################################### Declare necessary constants ###################################################################

### For data analysis ###
# CLASS1 refers to in1 in ms_matcher. in1 appends _1 to all the headers, hence suf=1. ms_fof_matcher is done such that injected catalogs have no suffix #
if RUN_TYPE is not None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, suf='')
if RUN_TYPE is None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, suf='_1')
CLASS2 = get_class(cat_type=MATCH_CAT2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT, suf='_2')


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

# Used if LOG_FLAGS is True #
flag_idx = []


### For plot names, plot titles, log file names ###
TITLE_PIECE1, TITLE_PIECE2 = CLASS1.title_piece, CLASS2.title_piece
MATCH_TYPE = get_match_type(title_piece1=TITLE_PIECE1, title_piece2=TITLE_PIECE2)

'''
fn_full_log = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, 'all_flags_all_1sig.csv') 
fd_full_log = open(fn_full_log, 'w')
fd_full_log.write('TILE\tREALIZATION\tFILTER\tRUN_TYPE\tTOTAL_MATCHES\tTOTAL_FLAGS\tPERCENT_FLAGS\tTOTAL_1SIGMA_MAG\tPERCENT_1SIGMA_MAG\n')

fn_full_recovered_log = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, 'all_recovered.csv')
fd_full_recovered_log = open(fn_full_recovered_log, 'w')
fd_full_recovered_log.write('TILE\tREALIZATION\tFILTER\tPERCENT_RECOVERED\n')
'''
# Make pandas-readable csv #
log_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE)
if os.path.isdir(log_dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory ' + str(log_dir) + ' does not exist. \n Change directory structure.')
                if NO_DIR_MAKE:
                        print 'Making directory ', log_dir, '...\n'
                        os.makedirs(log_dir)
FN_LOG = os.path.join(log_dir, 'log_'+MATCH_TYPE+'.csv')
with open(FN_LOG, 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter=',')
	# Write headers #
	writer.writerow(['TILE', 'REALIZATION', 'FILTER', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJECTS', 'PERCENT_FLAGGED_OBJECTS', 'TOTAL_1SIGMA_MAG', 'PERCENT_1SIGMA_MAG', 'PERCENT_RECOVERED_FLAGS_INCLUDED', 'PERCENT_RECOVERED_FLAGS_REMOVED'])








def fd_first_write(fn_main_log, fn_mag_bins, fn_flag, fn_color_plot_log, fn_delta_mag_outliers, fn_completeness_log):
	"""Write headers to log files.

	Parameters
		fn_main_log (str) -- Complete filename for log file containing number of objects matched, number of objects flagged, number of objects within 1sigma_mag, etc.
		fn_mag_bins (str) -- Complete filename for log file containing details of error analysis such as bins used, median of the error in each bin, etc.
		fn_flag (str) -- Complete filename for log file that examines ALL flags declared above. File is empty if `LOG_FLAGS=False`. 
	Returns
		fd* (file descriptors) -- File descriptor for each fn*.
	"""

	### Open files ###
	fd_main_log = open(fn_main_log, 'w'); fd_mag_bins = open(fn_mag_bins, 'w'); fd_flag = open(fn_flag, 'w'); fd_color_plot_log = open(fn_color_plot_log, 'w')
	fd_delta_mag_outliers = open(fn_delta_mag_outliers, 'w'); fd_completeness_log = open(fn_completeness_log, 'w')

	### Write headers ###
	fd_main_log.write('TILE\tREALIZATION\tFILTER\tRUN_TYPE\tTOTAL_MATCHES\tTOTAL_FLAGS\tPERCENT_FLAGS\tTOTAL_1SIGMA_MAG\tPERCENT_1SIGMA_MAG\n')

	fd_mag_bins.write('TILE\tREALIZATION\tFILTER\tNUM_OBJS_IN_BIN\tBIN_LHS\tBIN_RHS\tMEDIAN_HAXIS_MAG\tMEDIAN_ERROR\n')

	fd_flag.write('TILE\tREALIZATION\tFILTER\tRUN_TYPE\tFLAG1_HEADER\tFLAG2_HEADER\tFLAG1_VALUE\tFLAG2_VALUE\tMAG1\tMAG2\n')

	if LOG_FLAGS is False:
		fd_flag.write('Flags not logged because LOG_FLAGS is False.')

	fd_color_plot_log.write('TILE\tREALIZATION\tCOLOR\tMAG_BIN\tOBJS_IN_MAG_BIN\tTOTAL_OBJECTS_NO_FLAGS\tRUN_TYPE\n')

	fd_delta_mag_outliers.write('TILE\tREALIZATION\tFILTER\tMAG1\tMAG2\tDELTA_MAG\n')

	fd_completeness_log.write('TILE\tREALIZATION\tINJ_PERCENT\tFILTER\tMAG_BIN_L\tMAG_BIN_R\tMATCH_CAT_OBJS_IN_BIN\tTRUTH_CAT_OBJS_IN_BIN\n')

	return fd_main_log, fd_mag_bins, fd_flag, fd_color_plot_log, fd_delta_mag_outliers, fd_completeness_log






################################################################### Analysis ###################################################################
def get_floats_from_string(df, filter_name, hdr):
	"""Transform a list of strings of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]' to a list of floats of form '[1,1,...]' (if filter_name="g"), '[2,2,...]' ("r"), '[3,3,...]' ("i"), or '[4,4,...]' ("z"). This is necessary for CSVs created from ms_matcher or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings. 

	Parameters
	----------
        df (pandas DataFrame)
        
	filter_name (str)
		Allowed values: 'g' 'r' 'i' 'z'.
        
	hdr (str)
		Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '(1,2,3,4)'.
        
	Returns
	-------
        list_a (list of floats)
		Collection of the numbers corresponding to a particular index in a list of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]. 
	"""

	strings = df[hdr]; list_a = []

	# Each element (elmt) is of form '(1, 2, 3, 4)' #
	for elmt in strings:

		if filter_name == 'g':
			i = 1
			idx1 = elmt.find('(') + i
			idx2 = elmt.find(',')

		if filter_name == 'r':
			i = 2
			idx1 = elmt.replace(',', ';', 0).find(',') + i
			idx2 = elmt.replace(',', ';', 1).find(',')

		if filter_name == 'i':
			i = 2
			idx1 = elmt.replace(',', ';', 1,).find(',') + i
			idx2 = elmt.replace(',', ';', 2).find(',')

		if filter_name == 'z':
			i = 2
			idx1 = elmt.replace(',', ';', 2).find(',') + i
			idx2 = elmt.find(')')

		list_a.append(float(elmt[idx1:idx2]))

	if PRINTOUTS_MINOR:
		print 'Got ', hdr, ' for filter ', filter_name, '...'
		print ' Check: ', strings[0], ' & ', list_a[0], '\n'


	return list_a









def get_matrix_diagonal_element(df, filter_name, hdr):
	"""Transforms a list of 4x4 matrices where each element is a string of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' into a list of either the 11 (if filter_name is "g"), 22 ("r"), 33 ("i"), or 44 ("z") matrix elements. This is necessary for CSVs created from ms_matcher or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings.

	Parameters
	----------
	df (pandas DataFrame)

	filter_name (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	hdr (str)
		Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))'.

        Returns
	-------
	list_aa (list of floats)
		Collection of the numbers corresponding to a particular diagonal element in a list of 4-by-4 matrices.
	"""

	matrices = df[hdr]; list_aa = []

	# Each element in `matrices` is a matrix of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' #
	for matrix in matrices:

		if filter_name == 'g':
 			i, j = 2, 0
			idx1 = 0
			idx2 = matrix.find(',')

		if filter_name == 'r':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 4).find(',')
			idx2 = matrix.replace(',', ';', 5).find(',')

		if filter_name == 'i':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 9).find(',')
			idx2 = matrix.replace(',', ';', 10).find(',')

		if filter_name == 'z':
			i, j = 2, -1
			idx1 = matrix.replace(',', ';', 14).find(',')
			idx2 = matrix.replace(',', ';', 15).find(',')

		list_aa.append(float(matrix[idx1+i:idx2+j]))

	if PRINTOUTS_MINOR:
		print 'Got ', hdr, ' for filter ', filter_name
		print ' Check: ', matrices[0], ' & ', list_aa[0], '\n'


	return list_aa









def get_good_index_using_primary_flags(df, full_magnitude1, full_magnitude2, cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2, filter_name):
	"""Get indices of objects without flags* where flags* used are indicated in README.md. Store flagged indices if PLOT_FLAGGED_OBJS is True.

	Parameters
	----------
	df (pandas DataFrame)

	filter_name (str)
		Used if analysing Y3 Gold catalog. Can be `None`.

	full_magnitude1, full_magnitude2 (list of floats)
		Uncleaned lists containing magnitudes. 

	Returns
	-------
	idx_good (list of ints)
	
	idx_bad (list of ints)
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
        flag1, flag2 = df[flag_hdr1], df[flag_hdr2]
        cm_flag1, cm_flag2 = df[cm_flag_hdr1], df[cm_flag_hdr2]


	# Make arrays to take absolute value in next step #
	full_magnitude1, full_magnitude2 = np.array(full_magnitude1), np.array(full_magnitude2)

	# Get rid of these objects; 37.5 corresponds to a negative flux #
	if 'y3_gold' not in MATCH_CAT1 and 'y3_gold' not in MATCH_CAT2:
		idx_good = np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) )[0]

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

		idx_good = np.where( (df['SEXTRACTOR_FLAGS_'+filter_name.upper()+suf] == 0) & (df['IMAFLAGS_ISO_'+filter_name.upper()+suf] == 0) & (df[hdr] == 0) & (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) )[0]

		
	if PLOT_FLAGGED_OBJS:
		idx_bad = np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & ((flag2 != 0) | (flag1 != 0) | (cm_flag1 != 0) | (cm_flag2 != 0)) )[0]

	if PLOT_FLAGGED_OBJS is False:
		idx_bad = None

		
        if PRINTOUTS:
		if 'y3_gold' not in MATCH_CAT1 and 'y3_gold' not in MATCH_CAT2:
			print 'Eliminated ', len(full_magnitude1) - len(idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ' ... \n'
		# For Y3 #
		if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
			print 'Eliminated ', len(full_magnitude1) - len(idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ',', hdr, ', ', 'SEXTRACTOR_FLAGS_'+filter_name.upper(), ', ',  'IMAFLAGS_ISO_'+filter_name.upper(),  ' ... \n'

	return idx_good, idx_bad	









def get_good_index_using_quality_cuts(df, full_magnitude1, full_magnitude2):
	"""Get indices of objects that satisfy quality cuts introduced by Eric Huff. Also get indices of objects without flags* as described in README.md. Store the flagged indices if PLOT_FLAGGED_OBJS is True.

	Parameters
	----------
	df (pandas DataFrame)

	*_hdr (str)
		Headers refer to column names in the matched catalog.

	full_magnitude1, full_magnitude2 (list of floats)
		Values read directly from pandas DataFrame or passed through `get_floats_from_string()`; no flags removed. 		

        Returns
	-------
	idx_good (list of ints)
		Indices of objects without flags and objects which met criteria for quality cuts.

	idx_bad (list of ints)
		Is empty if PLOT_FLAGGED_OBJS is False.
	"""

	if 'true' in AXLABEL1 and 'true' in AXLABEL2:
		sys.exit('ERROR. Cuts should be performed on measured catalog, not truth catalog.')

	# Ignore truth catalog if one is present. Preserve ability to check both catalogs. #
	if 'true' in AXLABEL1 and 'meas' in AXLABEL2:
		CM_T_HDR1, CM_T_ERR_HDR1, CM_S2N_R_HDR1, PSFREC_T_HDR1 = CM_T_HDR2, CM_T_ERR_HDR2, CM_S2N_R_HDR2, PSFREC_T_HDR2

	if 'true' in AXLABEL2 and 'meas' in AXLABEL1:
		CM_T_HDR2, CM_T_ERR_HDR2, CM_S2N_R_HDR2, PSFREC_T_HDR2 = CM_T_HDR1, CM_T_ERR_HDR1, CM_S2N_R_HDR1, PSFREC_T_HDR1

	idx_good, idx_bad = [], []

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if cm_flag_hdr2 is None:
		cm_flag_hdr2 = cm_flag_hdr1

	if cm_flag_hdr1 is None:
		cm_flag_hdr1 = cm_flag_hdr2


	### Define flags ###
	flag1, flag2 = df[flag_hdr1], df[flag_hdr2]
	cm_flag1, cm_flag2 = df[cm_flag_hdr1], df[cm_flag_hdr2]

	### Define parameters needed for quality cuts ###
	# Size squared of object #
	cm_t1, cm_t2 = df[CM_T_HDR1], df[CM_T_HDR2]
	# Size error #
	cm_t_err1, cm_t_err2 = df[CM_T_ERR_HDR1], df[CM_T_ERR_HDR2]
	# Signal to noise #
	cm_s2n_r1, cm_s2n_r2 = df[CM_S2N_R_HDR1], df[CM_S2N_R_HDR2]
	# PSF size #
	psfrec_t1, psfrec_t2 = df[PSFREC_T_HDR1], df[PSFREC_T_HDR2]

	# Cast into array to take absolute value in next step #
	full_magnitude1 = np.array(full_magnitude1)
	full_magnitude2 = np.array(full_magnitude2)

	idx_good = np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) & (cm_s2n_r1 > 10) & (cm_s2n_r2 > 10) & (cm_t1/cm_t_err1 > 0.5) & (cm_t2/cm_t_err2 > 0.5) & (cm_t1/psfrec_t1 > 0.5) & (cm_t1/psfrec_t2 > 0.5) )[0]


	idx_bad = []


	if PLOT_FLAGGED_OBJS:
		counter_bad_mag = 0
		for i in np.arange(0, len(full_magnitude1)):

			# Get rid of these objects #
			if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99:
				counter_bad_mag += 1
				if flag1[i] != 0 or flag2[i] != 0 or cm_flag1[i] != 0 or cm_flag2[i] != 0 or cm_s2n_r1[i] < 10 or cm_s2n_r2[i] < 10 or cm_t1[i]/cm_t_err1[i] < 0.5 or cm_t2[i]/cm_t_err2[i] < 0.5 or cm_t1[i]/psfrec_t1[i] < 0.5 or cm_t2[i]/psfrec_t2[i] < 0.5:
					idx_bad.append(i)

	if PRINTOUTS:
		print 'Eliminated objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', CM_FLAGS_HDR1, ', ', cm_flag_hdr2, ' ...'
		print ' Eliminated objects with signal-to-noise ratio < 10 ...'
		print ' Eliminated objects with cm_T/cm_T_err < 0.5 ...'
		print ' Eliminated objects with cm_T/psfrec_T < 0.5 ...'

	if PLOT_FLAGGED_OBJS is False:
		idx_bad = None


	return idx_good, idx_bad









def handle_flags(df, flag_hdr1, flag_hdr2, filter_name, full_magnitude1, full_magnitude2, realization_number, tile_name, fd_flag):
	"""Examine a particular flag and write to log file. Can also be used to check all flags in a list of flags.

	Parameters
	----------
	df (pandas DataFrame)

	filter_name (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	full_magnitude1, full_magnitude2 (numpy.ndarray if directly from `df[hdr]` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.

	realization_number (str)
		Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.

	tile_name (str) 

        Returns
	-------
	idx_good (list of ints)
		Indices of objects with flags values of zero.

	idx_bad (list of ints)
		Indices of objects with nonzero flag values.
	"""

	idx_good, idx_bad = [], []; counter_idx_bad = 0

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if flag_hdr1 is None:
		flag_hdr1 = flag_hdr2

	if flag_hdr2 is None:
		flag_hdr2 = flag_hdr1


	### psf_flags are strings of form '(0,0,0,0)' ###
	if flag_hdr1 is not None and flag_hdr2 is not None and 'psf' not in flag_hdr1 and 'psf' not in flag_hdr2:
		flag1 = df[flag_hdr1]
		flag2 = df[flag_hdr2]

	if 'psf' in flag_hdr1 and 'psf' in flag_hdr2:
		flag1 = get_floats_from_string(df=df, hdr=flag_hdr1, filter_name=filter_name)
		flag2 = get_floats_from_string(df=df, hdr=flag_hdr2, filter_name=filter_name)


	### Check for flags ###
	for i in np.arange(0, len(full_magnitude1)):

		if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99:
			
			if flag1[i] == 0 and flag2[i] == 0:
				idx_good.append(i)
	
			if flag1[i] != 0 or flag2[i] != 0:
				idx_bad.append(i)
                                counter_idx_bad += 1

				### Write flags to file with headers TILE, REALIZATION, FILTER, RUN_TYPE, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2 ### 
				fd_flag.write(str(tile_name) + '\t' + str(realization_number) + '\t' + str(filter_name) + '\t' + str(RUN_TYPE) + '\t' + str(flag_hdr1) + '\t' + str(flag_hdr2) + '\t' + str(flag1[i]) + '\t' + str(flag2[i]) + '\t' + str(full_magnitude1[i]) + '\t' + str(full_magnitude2[i]) +'\n')



	if PRINTOUTS:
		print 'For tile: ', str(tile_name), ' and filter: ', str(filter_name), ', checked flags: ', flag_hdr1, ' & ', flag_hdr2, '...'

	### Check if flags were found ###
	if counter_idx_bad > 0 and PRINTOUTS:
		print ' Number of flags for magnitudes values 9999, 99, 37.5 and flags ', str(flag_hdr1), ' and ', str(flag_hdr2), ' : ', counter_idx_bad, '\n'


	return idx_good, idx_bad









def calculate_total_fractional_magnitude_error(cov_hdr, df, filter_name, flux_hdr, idx_good):
	"""Calculate the magnitude error via 1.08 * (flux_cov[i][i])^0.5 / flux[i]. Ignore flagged objects in error calculation.

	Parameters
	----------
	cov_hdr (str) 
		Header for flux covariance matrix in the matched catalog.

	df (pandas DataFrame)

	filter_name (str)
		Allowed values: 'g' 'r' 'i' 'z'.

	flux_hdr (str)
		Headers refer to column names in the matched catalog.

	idx_good (list of ints)
		Indices with flag values equal to zero.

	Returns
	-------
	error (list of floats)
		The magnitude error corresponding to each object.
	"""

	# Uncleaned lists for flux and flux covariance #
	full_flux = get_floats_from_string(df=df, hdr=flux_hdr, filter_name=filter_name)
	full_flux_cov = get_matrix_diagonal_element(df=df, hdr=cov_hdr, filter_name=filter_name)

	error, flux, fluxcov = [], [], []; counter_neg = 0

	# 'Safe' indices #
	for i in idx_good:
		flux.append(full_flux[i])
		fluxcov.append(full_flux_cov[i])

	# Calculations #
	for i in np.arange(0, len(flux)):

		# Throw out negative fluxcov (error calculation involves taking the square root of a negative) #
		if fluxcov[i] < 0:
			error.append(0)
			counter_neg += 1

		if fluxcov[i] == 0:
			print 'cm_flux_cov is 0'

		if fluxcov[i] > 0:
			# Pogsons number = 1.08 #
			err = 1.08 * fluxcov[i]**0.5 / flux[i] 
			error.append(err)

	if PRINTOUTS:
		print 'Calculated the magnitude error for filter: ', filter_name
		print ' Number of negative cm_flux_cov: ', counter_neg, ' / ', len(flux), '\n'

	#np.array(error[idx_good])

	return error









def calculate_and_bin_cm_T_signal_to_noise(cm_t_hdr, cm_t_err_hdr, df, idx_good, clean_magnitude1, clean_magnitude2):
	"""Calculate measured signal-to-noise ratio via cm_T/cm_T_err. Cuts on these values are performed on truth catalogs so consider measured catalogs only.

	Parameters
	----------
	cm_t_hdr (str)
		Header for the size squared of object. Headers refer to column names in the matched catalog.

	cm_t_err_hdr (str)
		Header for error on the size squared of object. Headers refer to column names in the matched catalog.

	df (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.

	idx_good (list of ints)
		Indices where no flags exist. 

	clean_magnitude1, clean_magnitude2 (list of floats)
		Magnitudes with flags removed.

	Returns
	----------
	s2n (list of floats) #FIXME cm_t_s2n
		cm_T signal-to-noise at each index without flags*.
	"""

	cm_t = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
	cm_t_err = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)

	# cm_T signal to noise (s2n) #
	cm_t_s2n = []
	for i in np.arange(0, len(cm_t)):
		cm_t_s2n.append(abs(cm_t[i]/cm_t_err[i]))

	# Bin signal to noise #
	# Bins suggested by Spencer Everett #
	bins = [0, 1, 9, 20, max(cm_t_s2n)]
	if PRINTOUTS:
		print 'Binning cm_T_s1n with bins: ', bins, ' and headers/axlabels:', cm_t_hdr, ', ', cm_t_err_hdr, '...'
		print ' Min and max absolute value of cm_T signal-to-noise: ', min(cm_t_s2n), ' and ', max(cm_t_s2n), '...'

	# idx_bins is a list of lists to preserve bin structure #	
	binned_s2n, binned_hax_mag, binned_vax_mag, idx_bins = [], [], [], []

	for j in np.arange(0, len(bins)-1):
		idx_temp = []
		for i in np.arange(0, len(cm_t_s2n)):
			if cm_t_s2n[i] > bins[j] and cm_t_s2n[i] < bins[j+1]:
				idx_temp.append(i)	
		if PRINTOUTS:
			print ' For cm_T_s2n, number of objects in bin ', bins[j], '-', bins[j+1], ': ', len(idx_temp)
		idx_bins.append(idx_temp)
		#idx_temp = np.where(s2n > bins[j] & (s2n < bins[j+1]))

	if PRINTOUTS:
		print ' '


	return idx_bins, bins, cm_t_s2n









def get_68percentile_from_normalized_data(norm_dm_bins, bins, hax_mag_bins):
	"""Calculate the point on the normalized vertical axis corresponding to the 68th percentile of each bin used in the error calculation. This percentile can be calculated in several ways: centered about zero, centered about the median of each bin, calculated using the absolute value of the data, calculated by examining the 34th percentile of the positive and negative portions of the data separately. 

	Parameters
	----------
	norm_dm_bins (list of list of floats)
		Normalized delta magnitudes. Bin structure preserved.

	bins (list of floats)
		Bins used in error calculation.

	hax_mag_bins (list of list of floats)
		Magnitudes on the horizontal axis. Bin structure preserved.

	Returns
	-------
	vax_68percentile (list of floats)
		Point on the vertical axis (vax) corresponding to 68 percentile. Each element in the list corresponds to a different bin.

	bins (list of floats)
		Bins used in error calculation.

	neg_vax_34percentile (list of floats)

	pos_vax_34percentile (list of floats)

	"""

	vax_68percentile, neg_vax_34percentile, pos_vax_34percentile = [], [], []

	PLOT_HIST = False

	# Loop through bins (b) #
	for b in np.arange(0, len(norm_dm_bins)):

		if norm_dm_bins[b] is None:
			vax_68percentile.append(None)
			neg_vax_34percentile.append(None)
			pos_vax_34percentile.append(None)

		if norm_dm_bins[b] is not None:

			### Find 68th percentile about zero ### 	
			# Values in current bin (icb) #
			vax_mag_bins_icb = norm_dm_bins[b]
			# Take absolute value of each point in bin #
			abs_vax_mag_bins_icb = [abs(elmt) for elmt in vax_mag_bins_icb]
			# Percentile sorts the data #
			vax_68percentile.append(np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'))	

			if PRINTOUTS_MINOR:
				# Check the percentile because interpolation='lower' was used #
				num = 0
				for j in np.arange(0, len(norm_dm_bins[b])):
					if abs(norm_dm_bins[b][j]) <= np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'):
						num += 1
				print 'Number of objects within 68 percentile via np.percentile(interpolation=lower): ', float(num)/len(norm_dm_bins[b]), '...\n'


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
				neg_vax_34percentile.append(np.percentile(neg_vax, 34, interpolation='lower'))
			if counter_pos > 0:
				pos_vax_34percentile.append(np.percentile(pos_vax, 34, interpolation='lower'))
			if counter_neg  == 0:
				neg_vax_34percentile.append(None)
			if counter_pos  == 0:
				pos_vax_34percentile.append(None)


			# Plot histogram to see distrubtion of data (data is not normally distributed) #
                        if PLOT_HIST:
                                plt.figure()
                                norm_dm = [abs(elmt) for elmt in norm_dm_bins[b]]
                                plt.hist(norm_dm)
                                plt.title('Bin LHS: ' + str(bins[b]))
                                plt.xlabel(r'$\Delta M$')
                                plt.axvline(x=0, color='black', linestyle=':', linewidth=0.5)
                                plt.show()


	return vax_68percentile, bins, neg_vax_34percentile, pos_vax_34percentile









def bin_and_cut_measured_magnitude_error(clean_magnitude1, clean_magnitude2, error1, error2, filter_name, tile_name, realization_number, fd_mag_bins, fd_delta_mag_outliers):
        """Clean error. Bin error according to horizontal axis of plot. Remove error values corresponding to objects with |DeltaMagnitude|>3. Do not consider error corresponding to empty bins nor bins with a small number of objects.

        Parameters
	----------
	clean_magnitude1, clean_magnitude2 (list of floats)
		Objects with flag values of zero and/or quality cuts performed.
                
	error1, error2 (list of floats)
		1 and 2 refer to the matched catalogs. 
        
	Returns
	-------
	binned_hax_mag_median (list of floats)
		List of medians of the horizontal axis magnitude in each bin.
              
	binned_vax_mag_median (list of floats)
		List of medians of the vertical axis magnitude in each bin. Vertical axis is computed via clean_magnitude1 - clean_magnitude2.

	binned_err_median (list of floats)
		Median of the error in each bin.

	bins (list of floats)
		Bins used. Binned according to horizontal axis.
        """

	### !!!!! Comment this block out if errors are to be computed using both catalogs regardless of origin (measured catalog or truth catalog) ###
	if 'meas' in AXLABEL1 and 'meas' not in AXLABEL2:
		error2 = np.zeros(len(error1))
		if PRINTOUTS:
			print 'Using measured catalog (catalog1) for error calculation ... '

	if 'meas' in AXLABEL2 and 'meas' not in AXLABEL1:
		error1 = np.zeros(len(error2))
		if PRINTOUTS:
			print 'Using measured catalog (catalog2) for error calculation ... '

	if 'meas' in AXLABEL1 and 'meas' in AXLABEL2:
		if PRINTOUTS:
			print 'Using measured catalogs (catalog1 AND catalog2) for error calculation ... '

	if 'true' in AXLABEL1 and 'true' in AXLABEL2:
		sys.exit('Errors are to be computed using the measured catalog(s), not the truth catalog(s).')


	### Define bins ###
        __step = 0.5
        # Find the absolute min and max of the magnitudes in the matched catalog #
        __lim_low1, __lim_low2 = min(clean_magnitude1), min(clean_magnitude2)
        __lim_high1, __lim_high2 = max(clean_magnitude1), max(clean_magnitude2)
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
        __bins = np.arange(__lim_low, __lim_high+(__step*0.1), __step)


	# Stores median of values in each bin #
	__hax_mag_bin_median, __vax_mag_bin_median, __err_mag_bin_median = [], [], []
	# List of lists. Stores all values in each bin #
	__hax_mag_bins, __vax_mag_bins, __err_mag_bins = [], [], []
	# Counter for empty bins and bins with a small number of objects #
        __counter_empty_bin = 0


        # Bin magnitude errors according to the magnitude on the horizontal axis #
        if SWAP_HAX:
                __hax_mag = clean_magnitude2
        if SWAP_HAX is False:
                __hax_mag = clean_magnitude1

	# Magnitude on the vertical axis (vax) #
	__vax_mag = np.array(clean_magnitude1) - np.array(clean_magnitude2)

	__clean_vax_mag, __clean_hax_mag = [], []


	__cutoff_delta_mag = 3
	__log_cutoff_delta_mag = 1
	__counter_delta_mag_geq1 = 0

	### Populate each bin ###
	for b in __bins: 

		__hax_mag_in_bin, __vax_mag_in_bin, __err_mag_in_bin = [], [], []
		__counter = 0

                for i in np.arange(0, len(clean_magnitude1)):

			if abs(__vax_mag[i]) > __log_cutoff_delta_mag:
				__counter_delta_mag_geq1 += 1
				# Write to log file with headers TILE REALIZATION FILTER MAG1 MAG2 DELTA_MAG #
				fd_delta_mag_outliers.write(str(tile_name) + '\t' + str(realization_number) + '\t' + str(filter_name) + '\t' + str(clean_magnitude1[i]) + '\t' + str(clean_magnitude2[i]) + '\t' + str(__vax_mag[i]) + '\n')

                        # Do not calculate errors using outlier magnitudes (chosen to be |Delta-M| > 3). Bin magnitude errors according to the magnitude on the horizontal axis of the plot #
                        if __hax_mag[i] >= b and __hax_mag[i] < b+__step and abs(__vax_mag[i]) < __cutoff_delta_mag: 
				__err_mag_in_bin.append((error1[i]**2 + error2[i]**2)**0.5)
				__vax_mag_in_bin.append(__vax_mag[i])
				__hax_mag_in_bin.append(__hax_mag[i])
				__clean_vax_mag.append(__vax_mag[i])
				__clean_hax_mag.append(__hax_mag[i])
                                __counter += 1

		# Written in log file, hence 'minor' #
                if PRINTOUTS_MINOR:
                        print ' For magnitude, number of objects in bin ', round(j, 2), '-', round(j+__step, 2), ': ', __counter, '...'


		### Write to log file ###
		if __counter == 0: write_median, write_err = None, None
		if __counter > 0: write_median, write_err = np.median(__hax_mag_in_bin), np.median(__err_mag_in_bin)
		# TILE, REALIZATION, FILTER, NUM_OBJS_IN_BIN, BIN_LHS, BIN_RHS, MEDIAN_HAXIS_MAG, MEDIAN_ERROR #
		fd_mag_bins.write( str(tile_name) + ' \t ' + str(realization_number) + ' \t ' + str(filter_name) + ' \t ' + str(__counter) + ' \t ' + str(round(b, 2)) + ' \t ' + str(round(b+__step, 2)) + ' \t ' + str(write_median)+ ' \t ' + str(write_err) + '\n')


                ### Tame error calculation and normalization by adding `None` to empty bins and bins with a small number of points ###
		# Define 'small' #
		if STACK_REALIZATIONS: CONST = 30
		if STACK_REALIZATIONS is False: CONST = 10 

		if __counter <= CONST:
                        __counter_empty_bin += 1
                        __err_mag_bin_median.append(None)
                        __hax_mag_bin_median.append(None)
                        __vax_mag_bin_median.append(None)
			# Add to list of lists to keep bin structure #
			__err_mag_bins.append(None)
                        __hax_mag_bins.append(None)
                        __vax_mag_bins.append(None)		

                if __counter > CONST:
                        __err_mag_bin_median.append(np.median(__err_mag_in_bin))
                        __hax_mag_bin_median.append(np.median(__hax_mag_in_bin))
                        __vax_mag_bin_median.append(np.median(__vax_mag_in_bin))
			# Add to list of lists to keep bin structure #	
			__err_mag_bins.append(__err_mag_in_bin)
                        __hax_mag_bins.append(__hax_mag_in_bin)
                        __vax_mag_bins.append(__vax_mag_in_bin)


	if PRINTOUTS:
                if SWAP_HAX:
                        print ' Binned clean_magnitude2 (from', MATCH_CAT2, ') with step size: ', __step, ', and minimum: ', __lim_low, ', and maximum: ', __lim_high, '...'
		if SWAP_HAX is False:
                        print ' Binned clean_magnitude1 (from', MATCH_CAT1, ') with step size: ', __step, ', and minimum: ', __lim_low, ', and maximum: ', __lim_high, '...'
		print ' Calculated errors using objects where |DeltaM| <', __cutoff_delta_mag, ' ... '
		print ' Excluded ', __counter_empty_bin, ' bins with less than ', CONST, ' objects ... \n'
		print ' Number of objects with |DeltaM| > ', __log_cutoff_delta_mag, ' : ', __counter_delta_mag_geq1, ' / ', str(len(clean_magnitude1)), ' ...\n'

        return __hax_mag_bin_median, __vax_mag_bin_median, __err_mag_bin_median, __bins, __hax_mag_bins, __vax_mag_bins, __err_mag_bins, __clean_hax_mag, __clean_vax_mag









def normalize_plot_maintain_bin_structure(clean_magnitude1, clean_magnitude2, error1, error2, filter_name, tile_name, realization_number, fd_mag_bins, fd_delta_mag_outliers):
	"""Normalize the vertical axis using 1sigma_mag and preserve the bin structure of the horizontal axis. 

	Parameters
	----------
	clean_magnitude1, clean_magnitude2 (list of floats)

	error1, error2 (list of floats)

	Returns
	-------
	norm_dm_bins (list of list of floats)

	bins (list of floats)
		Bins used in error calculation. 
	"""

	# List of lists. Stores all values in each bin #
	__norm_vax_bins, __hax_bins = [], []
	# Stores the median of each bin #
	__norm_vax_bin_median = []

	haxBinMedian, vaxBinMedian, errorBinMedian, bins, haxBins, vaxBins, errorBins, hax, vax = bin_and_cut_measured_magnitude_error(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, filter_name=filter_name, tile_name=tile_name, realization_number=realization_number, fd_mag_bins=fd_mag_bins, fd_delta_mag_outliers=fd_delta_mag_outliers)
 

	# Loop through bins (b) #
	for b in np.arange(0, len(vaxBins)):

		# Normalized Delta-Magnitudes (dm) in current bin (icb) #
		__norm_vax_in_bin, __hax_in_bin = [], []

		# 0 is a placeholder for empty bins and bins with few objects #
		if errorBinMedian[b] is None:
			__norm_vax_bins.append(None)
			__hax_bins.append(None)
			__norm_vax_bin_median.append(None)

		#if vax_mag_icb != 0:
		if errorBinMedian[b] is not None:
			for i in np.arange(0, len(vaxBins[b])):
				__norm_vax_in_bin.append(vaxBins[b][i]/errorBinMedian[b])
				__hax_in_bin.append(haxBins[b][i])

			# List of lists to keep bin structure #
			__norm_vax_bins.append(__norm_vax_in_bin)
			__hax_bins.append(__hax_in_bin)
			__norm_vax_bin_median.append(np.median(__norm_vax_in_bin))

	return __norm_vax_bins, bins, __hax_bins, errorBinMedian, __norm_vax_bin_median









def normalize_plot(norm_delta_mag_bins, bins, hax_mag_bins):
	"""Normalize plot to 1sigma_mag curve using tame magnitude errors

	Parameters
	----------
	norm_dm_bins (list of list of floats)
		Normalized delta magnitudes in each bin. 

	bins (list of floats)
		Bins used in error calculation.

	hax_mag_bins (list of list of floats)
		Magnitudes on the horizontal axis. Bin structure preserved.

        Returns
	-------
	norm_dm (list of floats)
		Delta-Magnitude normalized by error. Delta-Magnitude computed via magnitude1 - magnitude2. 

        hax_mag (list of floats)
		Magnitude to be plotted on the horizontal axis.
	"""

	### Remove `None` so that lists can be flattened. `None` is a placeholder for missing lists due to empty or small bin. ###
	__idx_good = []
        for i in np.arange(0, len(norm_delta_mag_bins)):
                if norm_delta_mag_bins[i] is not None:
                        __idx_good.append(i)
	#bins, norm_delta_mag_bins, hax_mag_bins = np.array(bins)[__idx_good], np.array(norm_delta_mag_bins)[__idx_good], np.array(hax_mag_bins)[__idx_good]
	__bins, __norm_vax_bins, __hax_bins = np.array(bins)[__idx_good], np.array(norm_delta_mag_bins)[__idx_good], np.array(hax_mag_bins)[__idx_good] 

	### Flatten lists ###
	__hax = [item for sublist in __hax_bins for item in sublist]
	__norm_vax = [item for sublist in __norm_vax_bins for item in sublist]

	# Get __idx_relevant #	
	__idx_relevant = np.where((__hax >= min(bins)) & (__hax < max(bins)))[0]
	__hax, __norm_vax = np.array(__hax)[__idx_relevant], np.array(__norm_vax)[__idx_relevant]

	return __norm_vax, __hax, __bins
	#return norm_dm, hax_mag, bins










def one_sigma_counter(delta_mag, clean_magnitude1, bins, hax_mag, error_bins, vax_bins_median):
	"""Find the number of objects within 1sigma_mag. This function is called if `NORMALIZE` is False.

	Parameters
	----------
	delta_mag (list of floats)
		NON-normalized Delta-Magnitude.

	error (list of floats)
		Is median error IN EACH BIN

	Returns
	-------
	counter_1sig (int)
		Number of objects within 1sigma_mag curve.
        """

	if len(bins) != len(error_bins):
                sys.exit('len(bins) not equal to len(error)')

	tot = len(delta_mag)
	counter_1sig = 0; counter_objs_considered = 0

	
	# Cutoffs were introduced in error calculation. Consider only points not cutoff #
	# Get rid of `None` placeholders #
	idx_good = []
	for i in np.arange(0, len(error_bins)):
		if error_bins[i] is not None:
			idx_good.append(i)
	bins, error_bins, vax_bins_median = np.array(bins)[idx_good], np.array(error_bins)[idx_good], np.array(vax_bins_median)[idx_good]

	# Examine objects within the relevant bin bounds #
	__idx_relevant = np.where((hax_mag >= min(bins)) & (hax_mag < max(bins)))[0]
	hax_mag, delta_mag = np.array(hax_mag)[__idx_relevant], np.array(delta_mag)[__idx_relevant]


	if CENTER_ERR_ABT_ZERO:
		for b in np.arange(0, len(bins)-1):
			if error_bins[b] is not None: 
				for i in np.arange(0, len(hax_mag)):
					if hax_mag[i] >= bins[b] and hax_mag[i] < bins[b+1]:
						counter_objs_considered += 1
						if abs(delta_mag[i]) < error_bins[b]:
							counter_1sig += 1


	# Center normalization about the median (of vertical axis) of each bin #
	if CENTER_ERR_ABT_ZERO is False:
		print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(bins)-1):
                        if error_bins[b] is not None:
				print ' Considering objects with magnitudes (on the horizontal axis) in [', bins[b], ',', bins[b+1], ')'
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= bins[b] and hax_mag[i] < bins[b+1]:
						counter_objs_considered += 1
						if delta_mag[i] < error_bins[b] + vax_bins_median[b] and delta_mag[i] >= -1.0*error_bins[b] + vax_bins_median[b]: 
							counter_1sig += 1

	print ' NOT Normalized '
	print ' Total objects: ', tot
	print ' Number of objects considered:', counter_objs_considered
	print ' Number of objects within 1sigma_mag: ', counter_1sig
	print ' Fraction within 1sigma_mag: ', float(counter_1sig)/counter_objs_considered, '\n'

	return counter_1sig		









def norm_one_sigma_counter(norm_delta_mag, clean_magnitude1, bins, hax_mag, error_bins, norm_vax_bins):
	"""Find the number of objects within 1sigma_mag. This function is only called if `NORMALIZE` is True.

	Parameters
	----------
	norm_delta_mag (list of floats)
		Normalized Delta-Magnitude. #FIXME what about colors? 

        Returns
	-------
        counter_1sig (int)
		Number of objects within 1sigma_mag curve.
	"""

	tot = len(norm_delta_mag)
	counter_1sig, counter_objs_considered = 0, 0


        # Cutoffs were introduced in error calculation. Consider only points not cutoff #
        # Get rid of `None` placeholders #
	'''
	norm_vax_median = []
        for l in norm_vax_bins:
                if l is not None:
			print l
                        norm_vax_median.append(np.median(l))
	'''
        idx_good = []
        for i in np.arange(0, len(error_bins)):
                if error_bins[i] is not None:
                        idx_good.append(i)
        # Binned values #
        error_bins, bins = np.array(error_bins)[idx_good], np.array(bins)[idx_good]
	norm_vax_median = np.array(norm_vax_bins)[idx_good]

        # Examine objects within the relevant bin bounds #
        __idx_relevant = np.where((hax_mag >= min(bins)) & (hax_mag < max(bins)))[0]
        hax_mag, norm_delta_mag = np.array(hax_mag)[__idx_relevant], np.array(norm_delta_mag)[__idx_relevant]


	if CENTER_ERR_ABT_ZERO:
                for k in norm_delta_mag:
			counter_objs_considered += 1
                        if abs(k) < 1.0:
                                counter_1sig += 1
		#counter_objs_considered = len(norm_delta_mag)


	if CENTER_ERR_ABT_ZERO is False:
		print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(bins)-1):
			if error_bins[b] is not None:
				print ' Considering objects with magnitudes (on the horizontal axis) in [', bins[b], ',', bins[b+1], ')'
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= bins[b] and hax_mag[i] < bins[b+1]:
                                                counter_objs_considered += 1
                                                if norm_delta_mag[i] < 1 + norm_vax_median[b] and norm_delta_mag[i] >= -1.0 + norm_vax_median[b]: 
							counter_1sig += 1

	print ' Normalized '
	print ' Total objects: ', tot
        print ' Number of objects considered:', counter_objs_considered
        print ' Number of objects within 1sigma_mag: ', counter_1sig
        print ' Fraction within 1sigma_mag: ', float(counter_1sig)/counter_objs_considered, '\n'

	return counter_1sig









def get_flag_type(df, k):
	"""Print the flag type() once.

	Parameters
	----------
        df (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.

        k (int)
		Counter implemented so printout is not repeated.

        Returns
	-------
		0
	"""

	if k == 0:
		for flag_hdr in FLAG_HDR_LIST:
			print 'HEADER:', str(flag_hdr), ' -- EXAMPLE:', df[flag_hdr][0], ' -- TYPE:', type(df[flag_hdr][0])
			k += 1
        return 0









def get_color(filter_name):
	"""Color code plot such that each griz band is a different color.

	Parameters
	----------
	filter_name (str)
		Allowed values: 'g' 'r' 'i' 'z'

	Returns
	-------
	color (str)
		Color for plotting data points.

	cmap (str)
		Colormap for plotting density. 
	"""

	if filter_name == 'g': color, cmap = 'green', 'Greens'

	if filter_name == 'r': color, cmap = 'purple', 'Purples'

	if filter_name == 'i': color, cmap = 'darkgrey', 'Greys'

	if filter_name == 'z': color, cmap = 'navy', 'Blues'

	return color, cmap









def logger(vax_mag_bin_median, tile_name, filter_name, realization_number, clean_magnitude1, full_magnitude1, bins, hax_mag, fd_main_log, error, vax_mag, fraction_recovered):
	"""Write to various log files. Records the number of objects matched, the number of objects flagged, the number of objects within 1sigma_mag, etc.

	Parameters
	----------
	vax_mag (list of floats)
		Is normalized to 1sigma_mag if `NORMALIZE` is True.

	vax_mag_bin_median
		Is normalized to 1sigma_mag if `NORMALIZE` is True.

	fd_main_log (file descriptor)

	clean_magnitude1 (list of floats)
		Magnitude of catalog1 in matched (join=1and2) catalog with ~flagged~ objects removed.

	full_magnitude1 (list of floats)
		Magnitude read directly from pandas Dataframe of matched (join=1and2) catalog.

	bins

	hax_mag

	error

	fraction_recovered
	
        filter_name (str) 

        realization_number (str) 

	tile_name

        Returns
	-------
        percent_in_1sig (float)
		Percent of objects within 1sigma_mag.
	"""

	if NORMALIZE:
		num_1sig = norm_one_sigma_counter(norm_delta_mag=vax_mag, clean_magnitude1=clean_magnitude1, bins=bins, hax_mag=hax_mag, error_bins=error, norm_vax_bins=vax_mag_bin_median)
	if NORMALIZE is False:
		num_1sig = one_sigma_counter(delta_mag=vax_mag, clean_magnitude1=clean_magnitude1, bins=bins, hax_mag=hax_mag, error_bins=error, vax_bins_median=vax_mag_bin_median) #was median 

	num_flags = len(full_magnitude1)-len(clean_magnitude1)

	# TILE, REALIZATION, FILTER, RUN_TYPE, TOTAL_MATCHES, TOTAL_FLAGS, PERCENT_FLAGS, TOTAL_1SIGMA_MAG, PERCENT_1SIGMA_MAG #
	fd_main_log.write( str(tile_name) + ' \t ' + str(realization_number) + ' \t ' + str(filter_name) + ' \t ' + str(RUN_TYPE) + ' \t ' + str(len(full_magnitude1)) + ' \t ' + str(num_flags) + ' \t ' + str(float(num_flags)/len(full_magnitude1)*100) + ' \t ' + str(num_1sig) + ' \t ' + str(float(num_1sig)/len(clean_magnitude1)*100) + '\n')

	#fd_full_log.write( str(tile_name) + ' \t ' + str(realization_number) + ' \t ' + str(filter_name) + ' \t ' + str(RUN_TYPE) + ' \t ' + str(len(full_magnitude1)) + ' \t ' + str(num_flags) + ' \t ' + str(float(num_flags)/len(full_magnitude1)*100) + ' \t ' + str(num_1sig) + ' \t ' + str(float(num_1sig)/len(clean_magnitude1)*100) + '\n')

	percent_in_1sig = float(num_1sig)/len(clean_magnitude1)


	# Write to log.csv with headers: 'TILE', 'REALIZATION', 'FILTER', 'TOTAL_MATCHES', 'TOTAL_FLAGGED_OBJECTS', 'PERCENT_FLAGGED_OBJECTS', 'TOTAL_1SIGMA_MAG', 'PERCENT_1SIGMA_MAG', 'PERCENT_RECOVERED_FLAGS_INCLUDED', 'PERCENT_RECOVERED_FLAGS_REMOVED' #
	if 'truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2:
		if INJ1_20PERCENT and INJ2_20PERCENT: const1 = 10000.0
		if INJ1_20PERCENT is False and INJ2_20PERCENT is False: const1 = 5000.0
		#FIXME mixed % inj
		if (INJ1_20PERCENT and INJ2_20PERCENT) is False and (INJ1_10PERCENT and INJ2_10PERCENT) is False: const1=1; const2=15000.0/2
		if STACK_TILES: const2 = len(ALL_TILES)
		if STACK_REALIZATIONS: const2 = len(ALL_REALIZATIONS)
		if STACK_TILES is False and STACK_REALIZATIONS is False: const2 = 1
		#FIXME this is duplicate for including flags
		# Including flags #
		percent_recovered_flags = 100.0*len(full_magnitude1)/(const1*const2)
		# Not including flags #
		percent_recovered_flags_rm = 100.0*len(clean_magnitude1)/(const1*const2)

	if 'truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2:
		percent_recovered_flags = None
		percent_recovered_flags_rm = None

	# Append to csv #	
	with open(FN_LOG, 'a') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow([str(tile_name), str(realization_number), str(filter_name), str(len(full_magnitude1)), str(num_flags), str(100.0*num_flags/len(full_magnitude1)), str(num_1sig), str(100.0*num_1sig/len(clean_magnitude1)), str(percent_recovered_flags), str(percent_recovered_flags_rm)])


        return percent_in_1sig, percent_recovered_flags









def get_colorbar_for_magnitude_plot_properties(df, cm_t_hdr, cm_t_err_hdr, idx_good, clean_magnitude1, clean_magnitude2, axlabel, inj_10percent, inj_20percent):
	"""Get data that will be used for the colorbar of plot.

	Parameters
	----------
	df (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.
	
	cm_t_hdr, cm_t_err_hdr (str)
		Matched (join=1and2) catalog headers for cm_T (size) and cm_T_err. Can be `None`.

	idx_good (list of ints)
		Indices of objects without flags.		
		
	clean_magnitude1, clean_magnitude2 (list of floats)
		Magnitudes ..FIXME rep

	axlabel (str)
		Allowed values: 'true' 'meas'	

	Returns
	-------
	cbar_val (list of floats)
		Data used to produce colorbar. Can be `None` if not colorbar is to be plotted.

	cbar_idx_bins ???

	cbar_idx_bins ???

	cbar_label (str)
		Label for colorbar. Includes LaTeX \bf{} formatting. Can be `None`.
	"""

	if 'true' in axlabel:
                sys.exit('ERROR. Colorbars should describe measured catalog values, not truth catalog values.')


	### Colorbar label ###
	# Prefix to labels #
        if inj_10percent or inj_20percent:
                pref = 'inj_'
        if inj_10percent is False and inj_20percent is False:
		if 'y3_gold' not in match_cat:
                        pref = 'base_'
		if 'y3_gold' in match_cat:
                        pref = 'Y3_'

        if CM_T_COLORBAR: cbar_label = pref + cm_t_hdr[:-2] + '_' + str(axlabel)
        if CM_T_ERR_COLORBAR: cbar_label = pref + cm_t_err_hdr + '_' + str(axlabel)
        if CM_T_S2N_COLORBAR: cbar_label = 'cm_T_s2n_' + str(axlabel)


	### Colorbar value ###
        if CM_T_S2N_COLORBAR:
                cbar_idx_bins, cbar_bins, cbar_val = calculate_and_bin_cm_T_signal_to_noise(cm_t_hdr=cm_t_hdr, cm_t_err_hdr=cm_t_err_hdr, df=df, idx_good=idx_good, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2)

        if CM_T_ERR_COLORBAR:
                # For measured catalog, cuts performed on truth catalogs #
                cbar_val = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
		cbar_idx_bins, cbar_bins = None, None

        if CM_T_COLORBAR:
                cbar_val = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
                cbar_idx_bins, cbar_bins = None, None

	if CM_T_S2N_COLORBAR is False and CM_T_ERR_COLORBAR is False and CM_T_COLORBAR is False:
		cbar_val, cbar_idx_bins, cbar_bins, cbar_label = None, None, None, None


	return cbar_val, cbar_idx_bins, cbar_bins, cbar_label









def get_magnitude_error(mag_err_hdr, flux_hdr, cov_hdr, df, filter_name, idx_good, match_cat):
	"""Get errors for plot data.

	Parameters
	----------
	mag_err_hdr (str)
		Matched catalog header for the magnitude error. Can be `None` if error is to be calculated using flux and flux covariance matrix.

	flux_hdr (str)
		Matched catalog header refering to the flux. FIXME is flux_g, flux_r ?? Can be `None` (??). Used if mag_err_hdr is `None`.

	cov_hdr (str)
		Matched catalog header refering to the flux covariance matrix. Can be `None` (??). Used if mag_err_hdr is `None`.

	df (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.

	filter_name (str)

	idx_good (list of ints)
		Indices of objects without flags.

	match_cat (str)
		Catalog containing the data. set by `MATCH_CAT1` or `MATCH_CAT2`.

	Returns
	-------
	error (list of floats)
		Error in magnitude or color (set by `PLOT_COLOR`). `error` will be `None` if `PLOT_1SIG` is False.
	"""

        if PLOT_1SIG:
                if mag_err_hdr is None:
                        magError = calculate_total_fractional_magnitude_error(df=df, flux_hdr=flux_hdr, cov_hdr=cov_hdr, filter_name=filter_name, idx_good=idx_good)
                if mag_err_hdr is not None:
			if match_cat == 'coadd':
				magError = get_floats_from_string(df=df, hdr=mag_err_hdr, filter_name=filter_name)
			if match_cat == 'star_truth' or 'y3_gold' in match_cat:
                                magError = df[str(mag_err_hdr[:-2]) + '_' + filter_name.upper() + str(mag_err_hdr[-2:])]
			# Pass good indices #
			magError = np.array(magError)[idx_good]

        if PLOT_1SIG is False:
		magError = None

	return magError 





def get_color_plot_error(mag_err_hdr, flux_hdr, cov_hdr, df, filter_name, idx_good, match_cat):
	""" TODO. See @get_magnitude_plot_error docstring """

	magErr = get_magnitude_error(mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, cov_hdr=cov_hdr, df=df, filter_name=filter_name, idx_good=idx_good, match_cat=match_cat)

	if filter_name == 'g': filterFollow = 'r'
	if filter_name == 'r': filterFollow = 'i'
	if filter_name == 'i': filterFollow = 'z'

	if mag_err_hdr is None:
		magErrorFollow = calculate_total_fractional_magnitude_error(df=df, flux_hdr=flux_hdr, cov_hdr=cov_hdr, filter_name=filterFollow, idx_good=idx_good)
	if mag_err_hdr is None:
		if match_cat == 'coadd':
			magErrorFollow = get_floats_from_string(df=df, hdr=mag_err_hdr, filter_name=filterFollow)
			if match_cat == 'star_truth' or 'y3_gold' in match_cat:
				magErrorFollow = df[str(mag_err_hdr[:-2]) + '_' + filterFollow.upper() + str(mag_err_hdr[-2:])]

	__color_error = (np.array(magErr)**2 + np.array(magErrorFollow)**2)**0.5

	return __color_error





def get_good_data(df, hdr, idx_good, magnitude, filter_name):
	"""Get the data corresponding to good indices (no flags or post quality cuts).

	Parameters
	----------
	df (pandas DataFrame)
		DataFrame for the matched catalog.
	
	hdr (str)
		Header refers to the matched catalog.

	idx_good (list of floats)
		Indices of objects without flags.
 
	magnitude (bool)
		If True the data is a str of form '(data_g, data_r, data_i, data_z)'.

	filter_name (str) 
		Can be `None`. Used if `magnitude=True`.

	Returns
	-------
	clean_data (list of floats)
		Data with flagged objects removed.
	"""

	if magnitude:
		full_data = get_floats_from_string(df=df, hdr=hdr, filter_name=filter_name)
	if magnitude is False:
                full_data = df[hdr]

	clean_data = np.array(full_data)[idx_good]

	return clean_data







def get_color_from_binned_magnitude(df, hdr1, hdr2, clean_magnitude_1a, clean_magnitude_2a, filter_name, idx_good):
	"""Get color with magnitudes binned of form '(color_[21-22), color_[22-23), color_[23-24), color_[24-25) where [low-high) refer to magnitudes)'
        Parameters
        ----------
        df (pandas DataFrame)
                DataFrame for the matched catalog.

        hdr (str)
                Header for magnitude.

        clean_magnitude_a (list of floats)
                List of magnitudes with flagged objects removed.

        filter_name (str)

        idx_good (list of ints)
                Indices of objects without flags. #FIXME point to flags used somewhere in README.md.

        Returns
        -------
	"""

	if filter_name == 'g': __filter_follow = 'r'
        if filter_name == 'r': __filter_follow = 'i'
        if filter_name == 'i': __filter_follow = 'z'
	if filter_name == 'z': __filter_follow = 'g'
	# `From MATCH_CAT1` #
	__clean_magnitude_1b = np.array(get_floats_from_string(df=df, hdr=hdr1, filter_name=__filter_follow))[idx_good]	
	# `From MATCH_CAT2` #
	__clean_magnitude_2b = np.array(get_floats_from_string(df=df, hdr=hdr2, filter_name=__filter_follow))[idx_good] 

	__color1 = np.array(clean_magnitude_1a) - np.array(__clean_magnitude_1b)
	__color2 = np.array(clean_magnitude_2a) - np.array(__clean_magnitude_2b)

	# Indices in each bin #
	__idx1, __idx2, __idx3, __idx4 = [], [], [], []

	__bins = [20, 21, 22, 23, 24]

	# Bin magnitude according to `clean_magnitude_1a` TODO #
	for i in np.arange(0, len(clean_magnitude_1a)):
		if clean_magnitude_1a[i] >= __bins[0] and clean_magnitude_1a[i] < __bins[1]:
			__idx1.append(i) 
		if clean_magnitude_1a[i] >= __bins[1] and clean_magnitude_1a[i] < __bins[2]:
			__idx2.append(i)
		if clean_magnitude_1a[i] >= __bins[2] and clean_magnitude_1a[i] < __bins[3]:
			__idx3.append(i)
		if clean_magnitude_1a[i] >= __bins[3] and clean_magnitude_1a[i] < __bins[4]:
			__idx4.append(i)

	#__color_bin1 = np.array(clean_magnitude_a)[__idx1] - np.array(clean_magnitude_b)[__idx1]
	__color1_bin1, __color1_bin2, __color1_bin3, __color1_bin4 = __color1[__idx1], __color1[__idx2], __color1[__idx3], __color1[__idx4]
	__color2_bin1, __color2_bin2, __color2_bin3, __color2_bin4 = __color2[__idx1], __color2[__idx2], __color2[__idx3], __color2[__idx4]

	__color1_bins, __color2_bins = [], []
	# Append. Want list of lists #
	__color1_bins.append(__color1[__idx1]); __color1_bins.append(__color1[__idx2]); __color1_bins.append(__color1[__idx3]); __color1_bins.append(__color1[__idx4])
	__color2_bins.append(__color2[__idx1]); __color2_bins.append(__color2[__idx2]); __color2_bins.append(__color2[__idx3]); __color2_bins.append(__color2[__idx4])

	return __color1_bins, __color2_bins, __bins



def get_magnitude_axlabel(inj_10percent, mag_hdr, axlabel, match_cat, filter_name, inj_20percent):
	"""Get labels for the horizontal axis. 
	
	Parameters
	----------
	inj_10percent (bool)
		If `inj_10percent=True` refers to 10% Balrog-injected catalog. Set by `INJ1_10PERCENT` or `INJ2_10PERCENT`.
		
	inj_20percent (bool)
		If `inj_20percent=True` refers to 20% Balrog-injected catalog. Set by `INJ1_20PERCENT` or `INJ2_20PERCENT`.
	
	mag_hdr, cm_t_hdr, cm_t_err (str) 
		Headers in the matched catalog that refer to the magnitude, cm_T (size), and cm_T error.
	
	match_cat (str)
		Catalog containing the data the hax_label describes. Set by `MATCH_CAT1` or `MATCH_CAT2`. 

	filter_name (str)

	Returns
	-------
		hax_label (str) -- Label for the horizontal axis. Includes LaTeX \bf{} formatting. 
	"""

	### Prefix to labels ###
	if inj_10percent:
		__pref = '10%_inj_'
	if inj_20percent:
		__pref = '20%_inj_'

	if inj_10percent is False and inj_20percent is False:
		if 'y3_gold' not in match_cat:
			__pref = 'base_'
		if 'y3_gold' in match_cat:
			__pref = 'Y3_'

	### Horizontal axis label ###

	### Magnitude label ###
	# Transform 'mag_2' to 'mag_true' or 'mag_meas' #
	__mag_axlabel = str(mag_hdr[:-2]) + '_' + str(axlabel)
	# 'mag_true' --> 'mag_{filter}_true' with {filter} bolded #
	__mag_axlabel = __mag_axlabel[:-4] + '$\\bf{' + str(filter_name) + '}$_' + __mag_axlabel[-4:]	
	__mag_axlabel = __pref + __mag_axlabel 
        # Coadd catalogs. Combined to get '(m_g, m_r, m_i, m_z)' then matched. #
	if match_cat == 'coadd':
		__mag_axlabel = 'MAG_AUTO_' + '$\\bf{' + str(filter_name) + '}$_' + str(axlabel) 

	return __mag_axlabel



 
def get_color_axlabel(inj_10percent, axlabel, match_cat, filter_name, inj_20percent):
	"""Get axes labels for color plots.

	Parameters
	----------
	axlabel (str)
		Allowed values: 'meas' 'true'. Set by `AXLABEL1` or `AXLABEL2`.	
	
	Returns
	-------
	__color_axlabel (str)
		Contains LaTeX formatting. Ex: 'inj_(g-r)_true'.
	"""

	#TODO add 10%_ or 20%_ pref
	if inj_10percent:
                __pref = '10%_inj_'
        if inj_20percent:
                __pref = '20%_inj_'
	#if inj_10percent or inj_20percent:
                #__pref = 'inj_'
        if inj_10percent is False and inj_20percent is False:
		if 'y3_gold' not in match_cat:
                        __pref = 'base_'
		if 'y3_gold' in match_cat:
                        __pref = 'Y3_'

	if filter_name == 'g': __color_axlabel = '$\\bf{(g-r)}$_' + str(axlabel)
        if filter_name == 'r': __color_axlabel = '$\\bf{(r-i)}$_' + str(axlabel)
        if filter_name == 'i': __color_axlabel = '$\\bf{(i-z)}$_' + str(axlabel)
        __color_axlabel = __pref + __color_axlabel	

	return __color_axlabel



def completeness_magnitude_subplotter(filter_name, df, flag_idx, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_name, plot_title, realization_number, tile_name, fd_flag, fd_completeness_log):
	"""Creates two completeness plots with two panels of completeness for 10% injections (first subplot panel) and 20% injections (second subplot panel).
	This function is only to be called when comparing Balrog-injected catalogs and truth catalogs.

	Parameters
	----------
	
	Returns
	-------
	"""

	#FIXME make function parameters for: INJ_20P and INJ_10P
	global INJ1_10PERCENT; global INJ2_10PERCENT; global INJ1_20PERCENT; global INJ2_20PERCENT

	# Force 10% inj to be plotted in first panel # 
	INJ1_10PERCENT = True; INJ2_10PERCENT = True; INJ1_20PERCENT = False; INJ2_20PERCENT = False

	if INJ1_10PERCENT: print 'Plotting 10% magnitude completeness ... \n'

	err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, haxLabel1, haxLabel2, vaxLabel = get_magnitude_plot_variables(filter_name=filter_name, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=realization_number, tile_name=tile_name, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag=fd_flag, plot_title=plot_title, plot_name=plot_name)

	completenessFraction1, binMedian1 = get_magnitude_completeness(clean_magnitude1=cleanMag1, clean_magnitude2=cleanMag2, full_magnitude1=fullMag1, full_magnitude2=fullMag2, tile_name=tile_name, realization_number=realization_number, filter_name=filter_name, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, error1=err1, error2=err2, fd_completeness_log=fd_completeness_log)#, inj10_percent=INJ, inj_20percent=False)

	# Get horizontal axis label #
        if 'truth' in MATCH_CAT1: __hax_label1 = haxLabel1
        if 'truth' in MATCH_CAT2: __hax_label1 = haxLabel2


	# Reverse bool value and recalculate #
	INJ1_10PERCENT = not INJ1_10PERCENT; INJ2_10PERCENT = not INJ2_10PERCENT; INJ1_20PERCENT = not INJ1_20PERCENT; INJ2_20PERCENT = not INJ2_20PERCENT

	if INJ1_20PERCENT: print 'Plotting 20% magnitude completeness ... \n'

	err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, haxLabel1, haxLabel2, vaxLabel = get_magnitude_plot_variables(filter_name=filter_name, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=realization_number, tile_name=tile_name, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag=fd_flag, plot_title=plot_title, plot_name=plot_name)

	completenessFraction2, binMedian2 = get_magnitude_completeness(clean_magnitude1=cleanMag1, clean_magnitude2=cleanMag2, full_magnitude1=fullMag1, full_magnitude2=fullMag2, tile_name=tile_name, realization_number=realization_number, filter_name=filter_name, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, error1=err1, error2=err2, fd_completeness_log=fd_completeness_log)


	if 'truth' in MATCH_CAT1: __hax_label2 = haxLabel1
        if 'truth' in MATCH_CAT2: __hax_label2 = haxLabel2


	plt.figure(figsize=(12, 10))

	plt.subplot(1, 2, 1)
	plt.plot(binMedian1, completenessFraction1, color='blue')
	plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
	plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
	plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
	plt.title(__hax_label1[:3]+' Injection')
	plt.ylabel('Magnitude Completeness')
	plt.xlabel(__hax_label1)
	plt.grid(linestyle='dotted')

	plt.subplot(1, 2, 2)
	plt.plot(binMedian2, completenessFraction2, color='green')
        plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
        plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
        plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
        plt.title(__hax_label2[:3]+' Injection')
        plt.ylabel('Magnitude Completeness')
        plt.xlabel(__hax_label2)
	plt.grid(linestyle='dotted')

	if '10% ' in plot_title:
		plot_title = plot_title.replace('10% ', '')
	if '20% ' in plot_title:
                plot_title = plot_title.replace('20% ', '')
	plt.suptitle(plot_title, fontweight='bold')
	plt.show()

	return 0



def get_magnitude_completeness(clean_magnitude1, clean_magnitude2, full_magnitude1, full_magnitude2, tile_name, realization_number, filter_name, mag_hdr1, mag_hdr2, error1, error2, fd_completeness_log):#, inj_10percent, inj_20percent):
	"""TODO after function is done"""

	if INJ1_10PERCENT: 
		print 'Calculating completeness of 10% injected catalogs ... \n'
		__write_inj_percent = '10%'
        if INJ2_20PERCENT: 
		print 'Calculating completeness of 20% injected catalogs ... \n'
		__write_inj_percent = '20%'

	__completeness_bins, __plot_bins = [], []


	### Photometry cuts on matched catalog (|DeltaM| < 3sigma) ###
        match_mag1, match_mag2 = [], []
        for k in np.arange(0, len(clean_magnitude1)):
                #print '3sigma_mag:', 3.0*(error1[k]**2 + error2[k]**2)**0.5
                #if abs(clean_magnitude1[k] - clean_magnitude2[k]) < 3.0 * (error1[k]**2 + error2[k]**2)**0.5:
                if True is True:
                        match_mag1.append(clean_magnitude1[k])
                        match_mag2.append(clean_magnitude2[k])

	
	### Read truth catalog which has magnitudes in form (m_g, m_r, m_i, m_z) ###
	if 'truth' in MATCH_CAT1:
		fn_truth_cat = get_catalog(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name, inj=INJ1)	
		mag_hdr = mag_hdr1
		# Use error in measured catalogs only #
		error1 = np.zeros(len(error2))
		truth_match_mag = match_mag1
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
		flag_hdr = FLAGS_HDR1[:-2]
		cm_flag_hdr = CM_FLAGS_HDR1[:-2] 

	if 'truth' in MATCH_CAT2:
		fn_truth_cat = get_catalog(cat_type=MATCH_CAT2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name, inj=INJ2)
		mag_hdr = mag_hdr2
		# Use error in measured catalogs only #
		error2 = np.zeros(len(error1))
		truth_match_mag = match_mag2
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
                flag_hdr = FLAGS_HDR2[:-2]
                cm_flag_hdr = CM_FLAGS_HDR2[:-2]


	if filter_name == 'g': idx = 0
	if filter_name == 'r': idx = 1
	if filter_name == 'i': idx = 2
	if filter_name == 'z': idx = 3
	hdu = fits.open(fn_truth_cat)
	data = hdu[1].data
	truth_mag_griz = data[mag_hdr[:-2]]
	truth_mag = []
	for mag_griz in truth_mag_griz:
		truth_mag.append(mag_griz[idx])
	truth_mag = np.array(truth_mag)

	flag = data[flag_hdr]
	cm_flag = data[cm_flag_hdr]

	# Primary flag cuts to truth catalog #
	__idx_good = np.where( (abs(truth_mag) != 9999.0) & (abs(truth_mag) != 99.0) & (abs(truth_mag) != 37.5) & (flag == 0) & (cm_flag == 0) )[0]
	truth_mag = truth_mag[__idx_good]


	### Bin according to unmatched truth catalog (bin sizes used for denominator) ###
        __step = 1
        __bins = np.arange(int(np.min(truth_mag)), int(np.max(truth_mag)), __step)


	for b in np.arange(0, len(__bins)-1):

		# Count completeness in each bin #
		__counter_truth, __counter_match = 0, 0

		for i in np.arange(0, len(truth_mag)):
			if truth_mag[i] >= __bins[b] and truth_mag[i] < __bins[b+1]:
				__counter_truth += 1

		for j in np.arange(0, len(truth_match_mag)):
			if truth_match_mag[j] >= __bins[b] and truth_match_mag[j] < __bins[b+1]:
				__counter_match += 1 

		if __counter_truth != 0:
			__completeness_bins.append((1.0*__counter_match)/__counter_truth)
			__plot_bins.append(np.median([__bins[b], __bins[b+1]]))

		# Write to log file: 'TILE\tREALIZATION\t\INJ_PERCENT\tFILTER\tMAG_BIN_L\tMAG_BIN_R\tMATCH_CAT_OBJS_IN_BIN\tTRUTH_CAT_OBJS_IN_BIN #
		fd_completeness_log.write(str(tile_name) + '\t' + str(realization_number) + '\t' + str(__write_inj_percent) + '\t' + str(filter_name) + '\t' + str(__bins[b]) + '\t' + str(__bins[b+1]) + '\t' + str(__counter_match) + '\t' + str(__counter_truth) + '\n')

	'''
	# Check with alternative method, is consistent #
	__completeness_bins = []
	numer = 1.0*np.histogram(truth_match_mag, __bins)[0]
	denom = np.histogram(truth_mag, __bins)[0]
	denom = denom[denom != 0]
	numer = numer[denom != 0]	
	__plot_bins = np.array(__plot_bins)[denom != 0]
	__completeness_bins = numer/denom
	'''
	return __completeness_bins, __plot_bins 



def gaussian(x, mu, sig):
	"""Normalized Gaussian."""
	return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2.0)





def norm_flux_histogram_plotter(filter_name, df, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title):
	"""Create a normalized 1D histogram of DeltaFlux/SigmaFlux.

	Parameters
	----------

	Returns
	-------
	"""

	### Plot DeltaFlux/SigmaFlux as normalized histogram ###
	plotData = get_flux_plot_variables(filter_name=filter_name, df=df, flux_hdr1=flux_hdr1, flux_hdr2=flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)
	__bins = 5*10**2
	# Normalize by counts via `density=True` (area under histogram equals 1) #
	plt.hist(plotData, __bins, density=True, histtype='step', color=get_color(filter_name=filter_name)[0], label='data', rwidth=1)


	### Plot standard normalized Gaussian and vertical line through mean to guide eye ###
        gx = np.linspace(np.min(plotData), np.max(plotData), 1000)
        gy = gaussian(x=gx, mu=0.0, sig=1.0)
        plt.plot(gx, gy, linestyle='--', color='maroon', label=r'$\mu=0 \, ; \, \sigma=1$', linewidth=0.9)
        plt.ylabel('Normalized Count')
        plt.xlabel(r'$\Delta F/\sigma_F$')
        plt.axvline(x=0, linestyle=':', color='black')


	### Fit Gaussian to data http://danielhnyk.cz/fitting-distribution-histogram-using-python/ ###
	#FIXME is this fit normalized?
	mean, stddev = stats.norm.fit(plotData)
	fit_data = stats.norm.pdf(gx, mean, stddev)
	plt.plot(gx, fit_data, color='cyan', label=r'fit: $\mu=$'+str(round(mean, 2))+' ; $\sigma=$'+str(round(stddev, 2)))
	haxLabel = get_flux_histogram_hax_label(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, filter_name=filter_name)
	plt.xlabel(haxLabel)

	#FIXME
	plt.xlim([-25, 25])
	plt.legend().draggable()
        plt.show()

	return 0





def get_flux_histogram_hax_label(mag_hdr1, mag_hdr2, filter_name):
	"""Get the horizontal axis label for the normalized flux histogram.

	Parameters
	----------

	Returns
	-------
	"""

	if 'cm_mag' in mag_hdr1 and 'cm_mag' in mag_hdr2:
		label = 'cm_flux_$\\bf{'+str(filter_name)+'}$ true$-$meas) / $\sigma_{flux\_meas}$'

	#TODO will we deal with psf mag ever?

	if INJ1_10PERCENT and INJ2_10PERCENT and (INJ1_20PERCENT and INJ2_20PERCENT) is False:
		label = '(10%_inj_'+label
	if INJ1_20PERCENT and INJ2_20PERCENT and (INJ1_10PERCENT and INJ2_10PERCENT) is False:
                label = '(20%_inj_'+label

	if (INJ1_10PERCENT and INJ2_10PERCENT and INJ1_10PERCENT and INJ2_10PERCENT) is False:
		label = '(base_'+label

	return label


def get_flux_from_magnitude(m):
	"""Compute the flux from magnitude.

	Parameters
        ----------

	Returns
        -------
	"""

	return 10**((30-m)/2.5)



def get_flux_error_from_magnitude(m, m_err):
	"""Compute the error in flux.

	Parameters
	----------
	m (list of floats)
		Magnitude

	m_err (list of floats)
		Magnitude error

	Returns
	-------
		Error in flux
	"""

	# flux --> magnitude: f = 10**((30-m)/2.5) #
	# Note: np.log is ln #
	return m_err * np.log(10.0) * 10**((30.0-m)/2.5) * (-1.0/2.5)




def color_subplotter(filter_name, df, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization_number, tile_name, fd_flag, plot_title, plot_name, fd_color_plot_log):
	"""Plot. Get titles.

	Parameters
        ----------

	Returns
	-------
	"""	

	# For log file #
	if filter_name == 'g': __write_color = 'g-r'
	if filter_name == 'r': __write_color = 'r-i'	
	if filter_name == 'i': __write_color = 'i-z'

	axLabel1 = get_color_axlabel(inj_10percent=INJ1_10PERCENT, axlabel=AXLABEL1, match_cat=MATCH_CAT1, filter_name=filter_name, inj_20percent=INJ1_20PERCENT)
	axLabel2 = get_color_axlabel(inj_10percent=INJ2_10PERCENT, axlabel=AXLABEL2, match_cat=MATCH_CAT2, filter_name=filter_name, inj_20percent=INJ2_20PERCENT)

	vaxDeltaColorLabel = get_delta_color_axlabel(label1=axLabel1, label2=axLabel2, filter_name=filter_name)

	### Get plot data ###
	errMag1, errMag2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, haxMagLabel1, haxMagLabel2, vaxMagLabel = get_magnitude_plot_variables(filter_name=filter_name, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=realization_number, tile_name=tile_name, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag=fd_flag, plot_title=plot_title, plot_name=plot_name)

	#TODO get_magnitude_error is called twice
	err1 = get_color_plot_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, cov_hdr=CM_FLUX_COV_HDR1, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT1)
        err2 = get_color_plot_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, cov_hdr=CM_FLUX_COV_HDR2, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT2)


	if filter_name != 'z':
		cleanColor1, cleanColor2, magBins = get_color_from_binned_magnitude(df=df, hdr1=mag_hdr1, hdr2=mag_hdr2, filter_name=filter_name, idx_good=idxGood, clean_magnitude_1a=cleanMag1, clean_magnitude_2a=cleanMag2)


	if PLOT_DELTA_VAX is False:
		__vax_color = cleanColor2

	### Plot ###
        if SWAP_HAX is False:
                __hax_label = axLabel1 
		__vax_label = axLabel2
		__subtitle_pref = haxMagLabel1
        if SWAP_HAX:
                __hax_label = axLabel2 
		__vax_label = axLabel1
		__subtitle_pref = haxMagLabel2

	if PLOT_DELTA_VAX:
                __vax_color = np.array(cleanColor1) - np.array(cleanColor2)
                __vax_label = vaxDeltaColorLabel 


	# Percentile levels for contours. Default is 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2) via https://github.com/dfm/corner.py/blob/master/corner/corner.py #
	# Correct 1sigma given by http://corner.readthedocs.io/en/latest/pages/sigmas.html #
	__lvls = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5]) ** 2)
	__clrs = ['red', 'blue', 'yellow']
	# FIXME `colors` passed to contour_kwargs are passed/applied in reversed order so reverse order for labels #
	__clrs_label = __clrs[::-1]
	__lw = 1.1

        ### Create 4-by-4 subplot. Figure size units: inches ###
	if PLOT_DELTA_VAX is False:
		plt.figure(figsize=(8, 12))

	if PLOT_DELTA_VAX:
		plt.figure(figsize=(12, 10))


        for i in np.arange(0, len(cleanColor1)):

		# Write to log file with headers 'TILE REALIZATION COLOR MAG_BIN OBJS_IN_MAG_BIN TOTAL_OBJS_FLAGS_INC RUN_TYPE #
		__write_bin = '['+str(magBins[i])+','+str(magBins[i+1])+')'
		fd_color_plot_log.write(str(tile_name) + '\t' + str(realization_number) + '\t' + str(__write_color) + '\t' + str(__write_bin) + '\t' + str(len(cleanColor1[i])) + '\t' + str(len(cleanMag1)) + '\t' + str(RUN_TYPE) + '\n')

		print 'Plotting ', len(cleanColor1[i]), ' objects ... \n'
                plt.subplot(2, 2, i+1)

		if CORNER_HIST_2D:
			# Create bins of 1/4 magnitude for hax and 1/20 magnitude for vax # 
			__ylow = np.min(np.array(__vax_color[i])) # -1*np.max(abs(np.array(__vax_color[i]))) ?
			__yhigh = np.max(np.array(__vax_color[i])) # np.max(abs(np.array(__vax_color[i]))) ?
			__bin_x = np.linspace(np.min(cleanColor1[i]), np.max(cleanColor1[i]), math.ceil(4.0*(np.max(cleanColor1[i]) - np.min(cleanColor1[i]))))
			__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))
			__ylims = np.max(abs(np.array(__vax_color[i]))) #? 
			#__ylims = np.min(np.array(__vax_color[i]))

			# Pass `levels` to `corner.hist2d` (interpreted as percentiles). `bins` are centered about zero along the vertical axis. #
                        corner.hist2d(cleanColor1[i], __vax_color[i], bins=np.array([__bin_x, __bin_y]), no_fill_contours=False, color=get_color(filter_name=filter_name)[0], levels=__lvls, contour_kwargs={'colors':__clrs, 'cmap':None, 'linewidths':__lw})

			# Plot points and not density bins #
			#corner.hist2d(cleanColor1[i], __vax_color[i], plot_density=False, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, color=get_color(filter_name=filter_name)[0], levels=__lvls, contour_kwargs={'colors':__clrs, 'cmap':None, 'linewidths':__lw}, data_kwargs={'alpha':0.25, 'ms':1.75})

			if YLOW is None and YHIGH is None:
				# Force symmetric vertical axis #
				plt.ylim([-1*__ylims, __ylims])
			if YLOW is not None and YHIGH is not None:
				plt.ylim([YLOW, YHIGH])

			if PLOT_DELTA_VAX is False:
				plt.axis('scaled') #plt.gca().set_aspect('equal')
				# Plot x=y line to guide eye #
				lims = [np.min([plt.xlim(), plt.ylim()]), np.max([plt.xlim(), plt.ylim()])] 
				plt.plot(lims, lims, color='k', linestyle=':', linewidth=0.7)
			if PLOT_DELTA_VAX:
				plt.axhline(y=0, color='black', linestyle=':', linewidth=0.7)
				plt.subplots_adjust(hspace=0.6)

			# Work-around for contour label #
			for j in np.arange(0, len(__lvls)):
				plt.plot([0.5, 0.5], [0.5, 0.5], color=__clrs_label[j], label='$P_{'+str(round(__lvls[j], 2))[2:]+'}$', linewidth=__lw)
			plt.legend().draggable()

                plt.ylabel(__vax_label)
                plt.xlabel(__hax_label)
                plt.title(__subtitle_pref+' bins: ' + __write_bin) 

        plt.subplots_adjust(hspace=0.4)
        plt.suptitle(plot_title, fontweight='bold')

        if SAVE_PLOT:
                print '-----> Saving plot as: ', plot_name.replace('placehold', __write_color)
                plt.savefig(plot_name.replace('placehold', __write_color))

        if SHOW_PLOT:
                plt.show()


	return 0 



def get_flux_plot_variables(filter_name, df, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2): 
	"""Get variables needed for flux histogram.
	
	Parameters
        ----------

        Returns
        -------
	"""

	# Gives same answer as that from using flux headers directly
        #cleanFlux1 = get_flux_from_magnitude(m=np.array(fullMag1)[idxGood])
        #cleanFlux2 = get_flux_from_magnitude(m=np.array(fullMag2)[idxGood])

	fullMag1 = get_floats_from_string(df=df, hdr=mag_hdr1, filter_name=filter_name)
        fullMag2 = get_floats_from_string(df=df, hdr=mag_hdr2, filter_name=filter_name)

        idxGood = get_good_index_using_primary_flags(df=df, full_magnitude1=fullMag1, full_magnitude2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, filter_name=filter_name)[0]


	### Get flux ###
	fullFlux1 = get_floats_from_string(df=df, hdr=flux_hdr1, filter_name=filter_name)
        fullFlux2 = get_floats_from_string(df=df, hdr=flux_hdr2, filter_name=filter_name)
        cleanFlux1 = np.array(fullFlux1)[idxGood]; cleanFlux2 = np.array(fullFlux2)[idxGood]


        ### Flux error from measured catalog only ###
        if 'truth' in MATCH_CAT1:

                magErr = get_magnitude_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, cov_hdr=CM_FLUX_COV_HDR1, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT1)
                fluxErr = get_flux_error_from_magnitude(m=np.array(fullMag2)[idxGood], m_err=np.array(magErr))

        if 'truth' in MATCH_CAT2:
                magErr = get_magnitude_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, cov_hdr=CM_FLUX_COV_HDR2, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT2)
                fluxErr = get_flux_error_from_magnitude(m=np.array(fullMag1)[idxGood], m_err=np.array(magErr))

	if 'truth' in MATCH_CAT1:
                __plot_data = (cleanFlux1-cleanFlux2)/fluxErr
        if 'truth' in MATCH_CAT2:
                __plot_data = (cleanFlux2-cleanFlux1)/fluxErr

	return __plot_data





def get_magnitude_plot_variables(filter_name, df, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization_number, tile_name, fd_flag, plot_title, plot_name, mag_axlabel1, mag_axlabel2):
	"""Get variables needed for magnitude plots.

	Parameters
	----------
	
	Returns
	-------
	"""

	### Axes labels ###
        haxLabel1 = get_magnitude_axlabel(inj_20percent=INJ1_20PERCENT, inj_10percent=INJ1_10PERCENT, mag_hdr=mag_hdr1, axlabel=AXLABEL1, match_cat=MATCH_CAT1, filter_name=filter_name)
        haxLabel2 = get_magnitude_axlabel(inj_20percent=INJ2_20PERCENT, inj_10percent=INJ2_10PERCENT, mag_hdr=mag_hdr2, axlabel=AXLABEL2, match_cat=MATCH_CAT2, filter_name=filter_name)
        vaxLabel = get_delta_magnitude_axlabel(label1=haxLabel1, label2=haxLabel2, filter_name=filter_name)


        # Magnitudes in full #
        fullMag1 = get_floats_from_string(df=df, hdr=mag_hdr1, filter_name=filter_name)
        fullMag2 = get_floats_from_string(df=df, hdr=mag_hdr2, filter_name=filter_name)

        ### Clean the data: removed flags and/or perform quality cuts ###
        if EH_CUTS:
                idxGood = get_good_index_using_quality_cuts(df, full_magnitude1=fullMag1, full_magnitude2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2)[0]

        if EH_CUTS is False:
                idxGood = get_good_index_using_primary_flags(df=df, full_magnitude1=fullMag1, full_magnitude2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, filter_name=filter_name)[0]

        # Magnitudes cleaned #
        cleanMag1 = get_good_data(df=df, hdr=mag_hdr1, idx_good=idxGood, magnitude=True, filter_name=filter_name)
        cleanMag2 = get_good_data(df=df, hdr=mag_hdr2, idx_good=idxGood, magnitude=True, filter_name=filter_name)


        ### Calculate errors. get_error() handles the case when PLOT_COLOR=True. ###
        err1 = get_magnitude_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, cov_hdr=CM_FLUX_COV_HDR1, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT1)
        err2 = get_magnitude_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, cov_hdr=CM_FLUX_COV_HDR2, df=df, filter_name=filter_name, idx_good=idxGood, match_cat=MATCH_CAT2)



        ### Write flags to file ###
        if LOG_FLAGS:
                for i in np.arange(0, len(FLAG_HDR_LIST), 2):
                        # Bad index #
                        temp_idx = handle_flags(df=df, filter_name=f, realization_number=realization_number, flag_hdr1=FLAG_HDR_LIST[i], flag_hdr2=FLAG_HDR_LIST[i+1], full_magnitude1=fullMag1, full_magnitude2=fullMag2, tile_name=tile_name, fd_flag=fd_flag)[1]
                #flag_idx.append(temp_idx)
                flag_idx.extend(temp_idx)


        ### Print out the type() for each flag ###
        if SHOW_FLAG_TYPE:
                get_flag_type(df=df, k=counter_flag_type_printout)
                counter_flag_type_printout += 1



	return err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2,  haxLabel1, haxLabel2, vaxLabel





def get_delta_color_axlabel(label1, label2, filter_name):
	"""Shorten vertical axis label for ease of readability.

	Parameters
        ----------
        label1, label2 (str) -- Labels from the horizontal axis. Contain LaTeX \bf{} formatting.

        filter_name (str)

        Returns
        -------
        short_label (str) -- Shortened label for vertical axis. Contains LaTeX \bf{} formatting.
        """

	short_label = ''

	if '10%_inj' in label1 and '10%_inj' in label2: short_label += '10%_inj_'
        if '20%_inj' in label1 and '20%_inj' in label2: short_label += '20%_inj_'
        if '10%_inj' in label1 and '20%_inj' in label2: short_label += '(10%$-$20%)_inj_'
        if '20%_inj' in label1 and '10%_inj' in label2: short_label += '(20%$-$10%)_inj_'
        if 'base' in label1 and 'base' in label2: short_label += 'base_'

	if filter_name == 'g': short_label += '$\\bf{(g-r)}$_'
        if filter_name == 'r': short_label += '$\\bf{(r-i)}$_' 
        if filter_name == 'i': short_label += '$\\bf{(i-z)}$_'

        if 'true' in label1 and 'true' in label2: short_label += 'true'
        if 'meas' in label1 and 'meas' in label2: short_label += 'meas'

        ### Label prefix ###
        if 'inj' not in short_label and 'base' not in short_label:
                if 'inj' in label1: pref1 = 'inj'
                if 'inj' in label2: pref2 = 'inj'

                if 'base' in label1: pref1 = 'base'
                if 'base' in label2: pref2 = 'base'

                if 'Y3' in label1: pref1 = 'Y3'
                if 'Y3' in label2: pref2 = 'Y3'

        if 'pref1' not in locals(): pref = None
        if 'pref1' in locals(): pref = pref1 + '$-$' + pref2

        ### Label suffix ###
        if 'meas' not in short_label and 'true' not in short_label: suf = label1[-4:] + '$-$' + label2[-4:]
        if 'meas' in short_label or 'true' in short_label: suf = None

        if pref is not None: short_label = pref + '  ' + short_label
        if suf is not None: short_label = short_label[:-1] + '  ' + suf


	return short_label





def get_delta_magnitude_axlabel(label1, label2, filter_name):
	"""Shorten vertical axis label for ease of readability. Ex: '10%_inj_cm_mag_g_true - 10%_inj_cm_mag_g_meas' --> '10%_inj_cm_mag_g true-meas'

	Parameters
	----------
	label1, label2 (str) -- Labels from the horizontal axis. Contain LaTeX \bf{} formatting.

	filter_name (str)

	Returns
	-------
	short_label (str) -- Shortened label for vertical axis. Contains LaTeX \bf{} formatting.
	"""

	short_label = ''

	### Look for similarities in vertical axis labels and build a shared label ###
	# Inj #
        if '10%_inj' in label1 and '10%_inj' in label2: short_label += '10%_inj_'
	if '20%_inj' in label1 and '20%_inj' in label2: short_label += '20%_inj_'
	if '10%_inj' in label1 and '20%_inj' in label2: short_label += '(10%$-$20%)_inj_'
	if '20%_inj' in label1 and '10%_inj' in label2: short_label += '(20%$-$10%)_inj_'
	if 'base' in label1 and 'base' in label2: short_label += 'base_'

	# Magnitude type #
        if 'cm_mag' in label1 and 'cm_mag' in label2: short_label += 'cm_mag_$\\bf{' + filter_name + '}$_'

	# Catalog type #
        if 'true' in label1 and 'true' in label2: short_label += 'true'
        if 'meas' in label1 and 'meas' in label2: short_label += 'meas'

	### Label prefix ###
        if 'inj' not in short_label:
		if '10%_inj' in label1: pref1 = '10%_inj'
		if '20%_inj' in label1: pref1 = '20%_inj'
		if '10%_inj' in label2: pref2 = '10%_inj'
		if '20%_inj' in label2: pref2 = '20%_inj'
		
		if 'base' in label1: pref1 = 'base'
		if 'base' in label2: pref2 = 'base'

		if 'Y3' in label1: pref1 = 'Y3'
		if 'Y3' in label2: pref2 = 'Y3'

	if 'pref1' not in locals(): pref = None
	if 'pref1' in locals(): pref = pref1 + '$-$' + pref2

	### Label suffix ###
        if 'meas' not in short_label and 'true' not in short_label: suf = label1[-4:] + '$-$' + label2[-4:]
        if 'meas' in short_label or 'true' in short_label: suf = None

        if pref is not None: short_label = pref + '  ' + short_label
        if suf is not None: short_label = short_label[:-1] + '  ' + suf

	return short_label






def normalized_delta_magnitude_plotter(mag_hdr1, mag_hdr2, cbar_val, error1, error2, filter_name, clean_magnitude1, full_magnitude1, mag_axlabel1, clean_magnitude2, mag_axlabel2, plot_title, realization_number, tile_name, idx_bins, bins, cbar_axlabel, plot_name, fd_main_log, fd_mag_bins, vax_label, fraction_recovered):
	"""Produce a plot normalized to 1sigma_mag.

	Parameters
        ----------

	Returns
	-------
	"""

	#TODO PLOT_1SIG --> 1SIG_MAG

	# Args needed to call normalize_plot() #
	normVaxBins, initialBins, haxBins, errorBinMedian, vaxBinMedian = normalize_plot_maintain_bin_structure(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, filter_name=filter_name, tile_name=tile_name, realization_number=realization_number, fd_mag_bins=fd_mag_bins, fd_delta_mag_outliers=fd_delta_mag_outliers)


	# Percentiles #
	PLOT_68P, PLOT_34P_SPLIT = True, True

	if PLOT_1SIG and CORNER_HIST_2D is False:

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

		### Plot +/-34 percentile of each bin ###

		# Line width for top and sides of bins #
		lwt = 1.1; lws = 0.7

		if PLOT_34P_SPLIT:
			### Plot the 68th percentile calculated from np.percentile() ###
			vax_68percentile_bins, percentileBins, neg_vax_34percentile, pos_vax_34percentile = get_68percentile_from_normalized_data(norm_dm_bins=normVaxBins, bins=initialBins, hax_mag_bins=haxBins)
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
	plotDeltaMag, plotHaxMag, cleanBins = normalize_plot(norm_delta_mag_bins=normVaxBins, bins=initialBins, hax_mag_bins=haxBins)

	percent_1sig, fractionRecoveredFlagsIncluded = logger(vax_mag=plotDeltaMag, filter_name=filter_name, clean_magnitude1=clean_magnitude1, full_magnitude1=full_magnitude1, realization_number=realization_number, tile_name=tile_name, bins=initialBins, hax_mag=plotHaxMag, fd_main_log=fd_main_log, error=errorBinMedian, vax_mag_bin_median=vaxBinMedian, fraction_recovered=fraction_recovered)

	# One colorbar at a time. This error is caught at beginning of script #
        if SCATTER:
                plt.scatter(plotHaxMag, plotDeltaMag, color=get_color(filter_name=filter_name)[0], alpha=0.25, s=0.25)


        if CM_T_S2N_COLORBAR or CM_T_ERR_COLORBAR or CM_T_COLORBAR:
                plt.scatter(plotHaxMag, plotDeltaMag, c=cbar_val, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
                plt.colorbar(label=cbar_axlabel)


        if BIN_CM_T_S2N:
                colors = ['green', 'purple', 'cyan', 'orange', 'pink', 'yellow', 'black', 'blue']
                for i in np.arange(0, len(idx_bins)):
                        plt.scatter(np.array(plotHaxMag)[idx_bins[i]], np.array(plotDeltaMag)[idx_bins[i]], color=colors[i], alpha=0.25, s=0.25, label='%1.f'%bins[i]+'<cm_T_s2n<%1.f'%bins[i+1])


        if HEXBIN:
		grid = (100, 1000)
		print ' Normalized hexbin has a large number of grid cells. Will take a moment to plot ... \n'
                plt.hexbin(plotHaxMag, plotDeltaMag, gridsize=grid, cmap=get_color(filter_name=filter_name)[1], bins='log')
                plt.colorbar(label='log(counts)')



	if HIST_2D:
                # 1/10 the bin size of that used in error calculation #
                bin_x = np.arange(min(plotHaxMag), max(plotHaxMag), 0.5/10)
                if YLOW is not None and YHIGH is not None:
                        # Somewhat using reported 1% error in magnitude #
                        bin_y = np.arange(YLOW, YHIGH, (YHIGH-YLOW)*0.01)
                if YLOW is None and YHIGH is None:
                        bin_y = np.arange(min(plotDeltaMag), max(plotDeltaMag), (max(plotDeltaMag)-min(plotDeltaMag))*0.01)
                plt.hist2d(plotHaxMag, plotDeltaMag, bins=[bin_x, bin_y], cmap=get_color(filter_name=filter_name)[1], norm=matplotlib.colors.LogNorm())
                plt.colorbar()


        if CORNER_HIST_2D:
                # Only the densest regions of the plot are binned so increase bin size of plt.hist2d() #
                const = 5
                bin_x = np.arange(min(plotHaxMag), max(plotHaxMag), const*0.5/10)
                if YLOW is not None and YHIGH is not None:
                        bin_y = np.arange(YLOW, YHIGH, (YHIGH-YLOW)*const*0.01)
                if YLOW is None and YHIGH is None:
                        bin_y = np.arange(min(plotDeltaMag), max(plotDeltaMag), (max(plotDeltaMag)-min(plotDeltaMag))*0.01)

		# Sigmas #
		__lvls = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5]) ** 2)
		__clrs = ['red', 'cyan', 'yellow']
		__clrs_label = __clrs[::-1]
		__lw = 0.9

                corner.hist2d(plotHaxMag, plotDeltaMag, no_fill_contours=True, levels=__lvls, color=get_color(filter_name=filter_name)[0], contour_kwargs={'colors':__clrs, 'cmap':None, 'linewidths':__lvls})
		# Work-around for contour labels #
		for j in np.arange(0, len(__lvls)):
			plt.plot([0.5, 0.5], [0.5, 0.5], color=__clrs_label[j], label='$P_{'+str(round(__lvls[j], 2))[2:]+'}$', linewidth=__lw)
		plt.legend().draggable()

	### Axes labels ###
        # Horizontal axis labels #
        if SWAP_HAX: plt.xlabel(str(mag_axlabel2))
        if SWAP_HAX is False: plt.xlabel(str(mag_axlabel1))
        # Vertical axis label #
     	plt.ylabel('('+ vax_label + ') / $\sigma$')

        plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)

        # Adjust vertical axes limits #
        if YLOW is not None and YHIGH is not None: plt.ylim([YLOW, YHIGH])


        ### Plot legend ###
        if PLOT_1SIG and BIN_CM_T_S2N is False and CORNER_HIST_2D is False: plt.legend(fontsize=10).draggable()
        if BIN_CM_T_S2N:
                # Increase marker size and opacity in legend #
                lgnd = plt.legend(markerscale=4, fontsize=8)
                for l in lgnd.legendHandles:
                        l.set_alpha(1)



        ### Title for subplot ###
        plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 4)*100) + '%')

	if SUBPLOT is False:

                plot_name = plot_name.replace('griz', filter_name)

                ### Title for  ###
                plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

                ### Save plot ###
                if SAVE_PLOT:
                        print '-----> Saving plot as: ', plot_name
                        plt.savefig(plot_name)

                if SHOW_PLOT:
                        plt.show()

        #plt.gca().set_aspect('equal')


	return 0







def delta_magnitude_plotter(mag_hdr1, mag_hdr2, cbar_val, error1, error2, filter_name, clean_magnitude1, full_magnitude1, mag_axlabel1, clean_magnitude2, mag_axlabel2, plot_title, realization_number, tile_name, idx_bins, bins, cbar_axlabel, plot_name, fd_main_log, fd_mag_bins, vax_label, fraction_recovered, fd_delta_mag_outliers):
	"""Produce a single plot of 'mag1' versus 'mag1'-'mag2'. 

	Parameters
	----------
	mag_hdr1, mag_hdr2 (str) -- Repeated


	cbar_val (list of floats) -- Data used to create a colorbar. Can be `None` in which case no colorbar is added.

	cbar_axlabel (str) -- Label for the colorbar.

	vax_label (str) -- Label for the vertical axis. Contains LaTeX \bf{} formatting.

	mag_axlabel1, mag_axlabel2 (str) -- Axes label for magnitude; will be used to label horizontal axis depending on `SWAP_HAX`. Contains LaTeX \bf{} formatting.

	error1, error2 (list of floats) -- Error in `MATCH_CAT1`, `MATCH_CAT2`. Error refers to the error in magnitude or color depending on the `PLOT_COLOR`.

	clean_magnitude1, clean_magnitude2 (list of floats) -- Magnitudes in matched catalog of `MATCH_CAT1`, `MATCH_CAT2` with flagged objects removed.

	full_magnitude1 (list of floats) -- Refers to the matched (join=1and2) catalog. No flagged objects removed.

	plot_title (str) -- Title for main plot (as opposed to individual subplots). 

	plot_name (str) -- Complete filename for plot save name.

	idx_bins -- ??

	bins (list of floats) -- Bins used to calculate error.

	fd* FIXME repeat

	filter_name (str)

	realization_number (str) 

	tile_name (str)

        Returns
	-------
		0
	"""

	### Values to plot ###
	__delta_mag = np.array(clean_magnitude1) - np.array(clean_magnitude2) 

	if SWAP_HAX:
		__hax_mag = clean_magnitude2
	if SWAP_HAX is False:
		__hax_mag = clean_magnitude1


	### 1sigma_mag curve ###
	if PLOT_1SIG:
		haxBinMedian, vaxBinMedian, errorBinMedian, initialBins, haxBins, vaxBins, errorBins, outlierCleanedHaxMag, outlierCleanedVaxMag = bin_and_cut_measured_magnitude_error(error1=error1, error2=error2, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, filter_name=filter_name, tile_name=tile_name, realization_number=realization_number, fd_mag_bins=fd_mag_bins, fd_delta_mag_outliers=fd_delta_mag_outliers)

		### Remove zeros from x, y, and err (zeros were placeholders for instances in which there were no objects in a particular magnitude bin) ###
		err = [temp for temp in errorBinMedian if temp is not None]
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
	percent_1sig, fractionRecoveredFlagsIncluded = logger(vax_mag=outlierCleanedVaxMag, filter_name=filter_name, clean_magnitude1=clean_magnitude1, full_magnitude1=full_magnitude1, realization_number=realization_number, tile_name=tile_name, bins=initialBins, hax_mag=outlierCleanedHaxMag, fd_main_log=fd_main_log, error=errorBinMedian, vax_mag_bin_median=vaxBinMedian, fraction_recovered=fraction_recovered)


	if PRINTOUTS:
                print 'Plotting ', len(clean_magnitude1), ' objects ... \n'

	### Plot ###
	# One colorbar at a time. This error is caught at beginning of script #
	if SCATTER:
		plt.scatter(__hax_mag, __delta_mag, color=get_color(filter_name=filter_name)[0], alpha=0.25, s=0.25)

	
	if CM_T_S2N_COLORBAR or CM_T_ERR_COLORBAR or CM_T_COLORBAR:
		'''To plot only the worst (smallest) s2n ratio:
		plt.scatter(np.array(__hax_mag)[idx_bins[0]], np.array(__delta_mag)[idx_bins[0]], color='purple', alpha=1, s=1, label='%1.f'%bins[0]+'<cm_T_s2n<%1.f'%bins[1])
		'''
		plt.scatter(__hax_mag, __delta_mag, c=cbar_val, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
		plt.colorbar(label=cbar_axlabel)


	if BIN_CM_T_S2N:
		colors = ['green', 'purple', 'cyan', 'orange', 'pink', 'yellow', 'black', 'blue']
		for i in np.arange(0, len(idx_bins)):
			plt.scatter(np.array(__hax_mag)[idx_bins[i]], np.array(__delta_mag)[idx_bins[i]], color=colors[i], alpha=0.25, s=0.25, label='%1.f'%bins[i]+'<cm_T_s2n<%1.f'%bins[i+1])


	if HEXBIN:
		grid = 500
		plt.hexbin(__hax_mag, __delta_mag, gridsize=grid, cmap=get_color(filter_name=filter_name)[1], bins='log')
		plt.colorbar(label='log(counts)')


	if HIST_2D:
		# 1/10 the bin size of that used in error calculation #
		bin_x = np.arange(min(__hax_mag), max(__hax_mag), 0.5/10)
		if YLOW is not None and YHIGH is not None:
			# Somewhat using reported 1% error in magnitude #
			bin_y = np.arange(YLOW, YHIGH, (YHIGH-YLOW)*0.01)
		if YLOW is None and YHIGH is None:
			bin_y = np.arange(min(__delta_mag), max(__delta_mag), (max(__delta_mag)-min(__delta_mag))*0.01) 
		plt.hist2d(__hax_mag, __delta_mag, bins=[bin_x, bin_y], cmap=get_color(filter_name=filter_name)[1], norm=matplotlib.colors.LogNorm())
		plt.colorbar()


	if CORNER_HIST_2D: 
		__lvls = 1.0 - np.exp(-0.5 * np.array([0.5, 1.0, 1.5]) ** 2)
		__clrs = ['red', 'cyan', 'yellow']
		__clrs_label = __clrs[::-1]
		__lw = 0.9

		# Only the densest regions of the plot are binned so increase bin size of plt.hist2d() #
		# SLACK channel corner.hist2d "draws 1- and 2-sigma contours automatically."Correct 1sigma levels: http://corner.readthedocs.io/en/latest/pages/sigmas.html #
		corner.hist2d(__hax_mag, __delta_mag, no_fill_contours=True, levels=__lvls, color=get_color(filter_name=filter_name)[0], contour_kwargs={'colors':__clrs, 'cmap':None, 'linewidths':__lw}) 
		# Work-around for contour labels #
		for j in np.arange(0, len(__lvls)):
			plt.plot([0.5, 0.5], [0.5, 0.5], color=__clrs_label[j], label='$P_{'+str(round(__lvls[j], 2))[2:]+'}$', linewidth=__lw)
		plt.legend().draggable()

		'''
		# Bin magnitudes. Using bins for error but not necessary #
		__counter_ax = 0
		for b in np.arange(0, len(initialBins)-1):
			__hax_mag_bin, __delta_mag_bin = [], []
			__counter_bin = 0
			for i in np.arange(0, len(__hax_mag)):
				if __hax_mag[i] >= initialBins[b] and __hax_mag[i] < initialBins[b+1]:
					__hax_mag_bin.append(__hax_mag[i])
					__delta_mag_bin.append(__delta_mag[i])
					__counter_bin += 1
			if __counter_bin > 0:
				#FIXME not being stacked. `ax` is NoneType 
				corner.hist2d(np.array(__hax_mag_bin), np.array(__delta_mag_bin), color=get_color(filter_name=filter_name)[0], contour_kwargs={'colors':'red', 'cmap':None, 'linewidths':0.7})
		'''
	### Axes labels ###
	# Horizontal axis labels #
	if SWAP_HAX: plt.xlabel(str(mag_axlabel2))
	if SWAP_HAX is False: plt.xlabel(str(mag_axlabel1))
	# Vertical axis label #
        plt.ylabel(vax_label)

	plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)

	# Adjust vertical axes limits #
        if YLOW is not None and YHIGH is not None: plt.ylim([YLOW, YHIGH])


	### Plot legend ###
	if PLOT_1SIG and BIN_CM_T_S2N is False and CORNER_HIST_2D is False: plt.legend(fontsize=10).draggable()
	if BIN_CM_T_S2N:
		# Increase marker size and opacity in legend #
		lgnd = plt.legend(markerscale=4, fontsize=8)
		for l in lgnd.legendHandles:
			l.set_alpha(1)



	### Title for subplot ###
	plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 4)*100) + '%')


	if SUBPLOT is False:

		plot_name = plot_name.replace('griz', filter_name)

		### Title for  ###
		plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

		### Save plot ###
		if SAVE_PLOT:
			print '-----> Saving plot as: ', plot_name
			plt.savefig(plot_name)

		if SHOW_PLOT:
			plt.show()

        #plt.gca().set_aspect('equal')


        return 0









def delta_magnitude_subplotter(df, flag_idx, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_name, plot_title, realization_number, tile_name, fd_flag, fd_main_log, fd_mag_bins, fraction_recovered, fd_delta_mag_outliers):
	"""Combine four subplots into a single plot with four panels (2-by-2). 

	Parameters
	----------
	df (pandas DataFrame) -- DataFrame for the matched catalog. 

	mag_hdr1, mag_hdr2 (str) -- Headers for the magnitude of `MATCH_CAT1` and `MATCH_CAT2` respectively after the catalogs have been matched.

	mag_err_hdr1, mag_err_hdr2 (str) -- Headers for the magnitude error of `MATCH_CAT1` and `MATCH_CAT2` respectively after matching. Can be `None`.

	plot_name (str) -- Complete name for plot.

	plot_title (str) -- Title for 2-by-2 plot. 

	fraction_recovered (float) -- Only applied to truth catalogs.

	realization_number (str) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.

	tile_name (str)

	fd_flag (file descriptor) -- Log file that collects flags.

	fd_main_log (file descriptor) -- Log file for the number of matched objects, number of flagged objects, number of objects recovered (if using truth catalog), etc.

	fd_mag_bins (file descriptor) -- Log file for error calculation including bins used, median error in each bin, etc.

	flag_idx (list of ints) -- Stores indices with nonzero flag values if `LOG_FLAGS` is True.

        Returns
	-------
	flag_idx (list of ints) -- If `LOG_FLAGS` is True, will check for all nonzero flag values in `FLAG_HDR_LIST` and `flag_idx` will contain indices that have nonzero flag values. Will be empty if `LOG_FLAGS` is False.
	"""



	# Counter for flag type() printout #
	counter_flag_type_printout = 0


        ### Create 4-by-4 subplot ###
	counter_subplot = 1
	if PLOT_COLOR is False:
	# Figure size units: inches #
		plt.figure(figsize=(12, 10))


        ### Create one subplot for each griz filter ###
	for f in ALL_FILTERS:


		### Define variables ###
		err1, err2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, haxLabel1, haxLabel2, vaxLabel = get_magnitude_plot_variables(filter_name=f, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=realization_number, tile_name=tile_name, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2, fd_flag=fd_flag, plot_title=plot_title, plot_name=plot_name)


		if CM_T_COLORBAR or CM_T_ERR_COLORBAR or CM_T_S2N_COLORBAR:
			cbarData, cbarIdxBins, cbarBins, cbarLabel = get_colorbar_for_magnitude_plot_properties(df=df, cm_t_hdr=CM_T_HDR2, cm_t_err_hdr=CM_T_ERR_HDR2, idx_good=idxGood, clean_magnitude1=cleanMag1, clean_magnitude2=cleanMag2, axlabel=AXLABEL2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT)
		else: 
			cbarData, cbarIdxBins, cbarBins, cbarLabel = None, None, None, None



		### Subplot ###
		if SUBPLOT: 
			plt.subplot(2, 2, counter_subplot)

		if NORMALIZE is False:
			delta_magnitude_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_val=cbarData, plot_title=plot_title, error1=err1, error2=err2, filter_name=f, full_magnitude1=fullMag1, clean_magnitude1=cleanMag1, clean_magnitude2=cleanMag2, mag_axlabel1=haxLabel1, mag_axlabel2=haxLabel2, realization_number=realization_number, tile_name=tile_name, idx_bins=cbarIdxBins, bins=cbarBins, cbar_axlabel=cbarLabel, plot_name=plot_name, fd_main_log=fd_main_log, fd_mag_bins=fd_mag_bins, vax_label=vaxLabel, fraction_recovered=fraction_recovered, fd_delta_mag_outliers=fd_delta_mag_outliers)

		if NORMALIZE:
			normalized_delta_magnitude_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_val=cbarData, plot_title=plot_title, error1=err1, error2=err2, filter_name=f, full_magnitude1=fullMag1, clean_magnitude1=cleanMag1, clean_magnitude2=cleanMag2, mag_axlabel1=haxLabel1, mag_axlabel2=haxLabel2, realization_number=realization_number, tile_name=tile_name, idx_bins=cbarIdxBins, bins=cbarBins, cbar_axlabel=cbarLabel, plot_name=plot_name, fd_main_log=fd_main_log, fd_mag_bins=fd_mag_bins, vax_label=vaxLabel, fraction_recovered=fraction_recovered, fd_delta_mag_outliers=fd_delta_mag_outliers)

		counter_subplot += 1

	# Plot after all panels have been filled #
	if SUBPLOT:

		### Show or save the plot once all four subplots have been filled ###
		plt.subplots_adjust(hspace=0.4)
		#plt.subplots_adjust(wspace=0.3)
		#plt.tight_layout(pad=3, h_pad=2.5)


		### Title ###
		if fraction_recovered is not None:
			plot_title = ' '.join([plot_title, 'Recovered (with flags):', str(round(fraction_recovered, 4)*100)+'%'])
		plt.suptitle(plot_title, fontweight='bold')

		### Save plot ###
		if SAVE_PLOT:
			print '-----> Saving plot as: ', plot_name
			plt.savefig(plot_name)

		### Show plot ###
		if SHOW_PLOT:
			plt.show()

	
	return flag_idx











def get_plot_suptitle(realization_number, tile_name, num_stack_real, num_stack_tile):
	"""Generate plot title.

	Parameters
	----------
	realization_number (str)


	tile_name (str)

	num_stack_real (int) -- Number of catalogs in stacked realization catalog. Can be `None`.

	num_stack_tile (int) -- Number of catalogs in stacked tile catalog. Can be `None`.

	Returns
	-------
	title (str) -- Ex: '10% Inj MOF Cat & 10% Inj Truth Cat' 
	"""

	if STACK_REALIZATIONS:
		realization_number = 'stacked '+str(num_stack_real)
	if STACK_TILES:
		tile_name = 'stacked ' + str(num_stack_tile)

	title = str(TITLE_PIECE1) + ' & ' + str(TITLE_PIECE2) +'. Tile: ' + str(tile_name) 

	if realizations[0] != 'None': 
		title = title + '. Realization: ' + str(realization_number) + '.'	
	
	if RUN_TYPE == 'ok': 
		title = title + ' Unchanged FOF groups.'
	if RUN_TYPE == 'rerun':
		title = title + ' Changed FOF groups.'

	if NORMALIZE:
		title = 'Normalized. ' + title

	return title









def get_plot_save_name(realization_number, tile_name):
        """Generate name of the plot that will be used in plt.savefig().
	Relies on directory structure: /`OUTDIR`/plots/`BALROG_RUN`/`MATCH_TYPE`/{tile}/realization}/plots/{plot_type}/ where allowed values for plot_type are: 'normalized' 'scatter' 'fof_analysis/normalized' 'fof_analysis/scatter'. 

        Parameters
	----------
	realization_number (str) 

	tile_name (str)

        Returns
	-------
	fn (str) -- Complete filename for plot. 
        """

	### Get filename ###
	if YLOW is None and YHIGH is None:
		# Default scale for the vertical axis (vax) #
		ylim = 'defaultvax'
	if YLOW is not None and YHIGH is not None:
		ylim = str(YLOW)+'y'+str(YHIGH)

	if RUN_TYPE is None:	
		#FIXME use join? endname = '_'.join([tile_name, realization_number, 'griz', MATCH_TYPE, ylim+'.png'])
		endname = str(tile_name) + '_' + str(realization_number) + '_griz_' + str(MATCH_TYPE) + '_' + str(ylim) + '.png'
	if RUN_TYPE is not None:
		endname = str(tile_name) + '_' + str(realization_number) + '_griz_' + str(MATCH_TYPE) + '_' + str(RUN_TYPE) + '_' + str(ylim) + '.png'

	# dm = delta magnitude #	
	if PLOT_COLOR is False:
		if CM_T_S2N_COLORBAR:
			outname = 'm_vs_dm_cm_t_s2n_' + endname
		if CM_T_COLORBAR:
			outname = 'm_vs_dm_cm_t_' + endname
		if CM_T_ERR_COLORBAR:
			outname = 'm_vs_dm_cm_t_err_' + endname
		if HEXBIN:
			outname = 'm_vs_dm_hexbin_' + endname
		if CM_T_S2N_COLORBAR is False and CM_T_COLORBAR is False and CM_T_ERR_COLORBAR is False and HEXBIN is False:
			outname = 'm_vs_dm_' + endname

	if PLOT_COLOR:
		if CORNER_HIST_2D:
			outname = 'color_placehold_'+endname

	# !!!!! User may wish to edit directory structure #
	plot_dir_pref = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'plots')

	if RUN_TYPE is not None:
		plot_dir_pref = os.path.join(plot_dir_pref, 'fof_analysis')

	if NORMALIZE:
		plot_dir = os.path.join(plot_dir_pref, 'normalized')
	if NORMALIZE is False:
		plot_dir = os.path.join(plot_dir_pref, 'scatter')

	
	### Check for directory existence ###
	if os.path.isdir(plot_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(plot_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_plot_save_name() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', plot_dir, '...\n'
			os.makedirs(plot_dir)


	### Get filename and path ###
        if NORMALIZE:
		fn = os.path.join(plot_dir, 'norm_' + str(outname))
        if NORMALIZE is False:
		fn = os.path.join(plot_dir, outname)


        return fn









def get_coadd_mag_and_mag_err(fn_g, fn_r, fn_i, fn_z, mag_hdr, mag_err_hdr):
	"""Creates a list of magnitudes and magnitude errors of form '(mag_g, mag_r, mag_i, mag_z)' from four catalogs. Solely for use with coadd catalogs.

	Parameters
	----------
	fn_g, fn_r, fn_i, fn_z -- Catalog filenames for each filter. Must be FITS files.


	mag_hdr (str) -- Header for magnitude. Headers refer to columns in the matched catalog.

	err_hdr (str) -- Header for magnitude error.

	Returns
	-------
	m_griz (list of str) -- Stores magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'

	m_err_griz (list of str) -- Stores error in magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'
	"""

	# Files have not yet been matched, and do not have hdr_1 #
	mag_hdr = mag_hdr[:-2]
	err_hdr = err_hdr[:-2]

	# Open FITS files #
	hdu_g = fits.open(fn_g); hdu_r = fits.open(fn_r); hdu_i = fits.open(fn_i); hdu_z = fits.open(fn_z)
	
	# Read data #
	data_g = hdu_g[1].data; data_r = hdu_r[1].data; data_i = hdu_i[1].data; data_z = hdu_z[1].data

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	err_g = data_g[err_hdr]; err_r = data_r[err_hdr]; err_i = data_i[err_hdr]; err_z = data_z[err_hdr]

	m_griz, m_err_griz = [], []

        for i in np.arange(0, len(m_g)):
                m_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")
		m_err_griz.append("'("+ str(err_g[i])+ ', ' + str(err_r[i])+ ', ' + str(err_i[i]) + ', ' + str(err_z[i]) + ")'")

        return m_griz, m_err_griz









def get_star_mag(df, suf):
	"""Computes and creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'. Solely for use with star truth catalogs.

	Parameters
	----------
	df (pandas DataFrame)

	suf (str) -- Allowed values: '_1' '_2'

	Returns
	-------
	m_griz (list of str) -- Stores magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'.
	"""

	m_g = df['g_Corr'+suf]

	m_r = df['g_Corr'+suf] - df['gr_Corr'+suf]

	m_i = df['g_Corr'+suf] - df['gr_Corr'+suf] - df['ri_Corr'+suf]

	m_z = df['g_Corr'+suf] - df['gr_Corr'+suf] - df['ri_Corr'+suf] - df['iz_Corr'+suf]

	m_griz = []

	for i in np.arange(0, len(m_g)):
		m_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")	
        
	return m_griz









def get_y3_gold_mag(df, mag_hdr):
	"""Creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'. Solely for use with Y3 Gold catalogs.

	Parameters
	----------
	df (pandas DataFrame)


	mag_hdr (str) -- General header for magnitude.

        Returns
	-------
	m_griz (list of str) -- Stores magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'.
        """

	# Get headers, which are dependent on filter #
	hdr_g = mag_hdr[:-2] + '_G' + mag_hdr[-2:]; hdr_r = mag_hdr[:-2] + '_R' + mag_hdr[-2:]
	hdr_i = mag_hdr[:-2] + '_I' + mag_hdr[-2:]; hdr_z = mag_hdr[:-2] + '_Z' + mag_hdr[-2:]
	
	# Read magnitudes from DataFrame #
	m_g = df[hdr_g]; m_r = df[hdr_r]; m_i = df[hdr_i]; m_z = df[hdr_z]

	m_griz = []

        for i in np.arange(0, len(m_g)):
                m_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")

        return m_griz









def get_catalog(cat_type, inj_10percent, inj_20percent, realization_number, tile_name, filter_name, inj):
        """Get catalog to analyze.
	
        Parameters
	----------
	cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd', 'y3_gold'. Set by `MATCH_CAT1` or `MATCH_CAT2`.

	inj_10percent (bool)

	inj_20percent (bool)

	realization_number (str) -- 

	tile_name -- 

	filter_name (str) -- Only used with coadd catalogs. Ignored if `cat_type` is not 'coadd'. 

        Returns
	-------
	fn (str) -- Complete catalog filename.
        """

	if (cat_type == 'gal_truth' or cat_type == 'star_truth') and inj_10percent is False and inj_20percent is False:
		sys.exit('No non-injected truth catalog exists.')

	# TODO As of May 2018 only TAMU_Balrog catalogs have 20% injections
	if BALROG_RUN == 'TAMU_Balrog':
		fn = get_tamu_catalog(cat_type=cat_type, inj_10percent=inj_10percent, inj_20percent=inj_20percent, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name, inj=inj)


	if BALROG_RUN != 'TAMU_Balrog':
		if cat_type == 'gal_truth' and inj_10percent:
			fn = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_gals.fits')
		#if cat_type == 'gal_truth' and inj_20percent

		if cat_type == 'star_truth' and inj_10percent:
			fn = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_stars.fits')

		if cat_type == 'sof' and inj_10percent:
			fn = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_sof.fits')
		if cat_type == 'sof' and inj_10percent is False:
			fn = os.path.join(BASEPATH, 'y3v02', tile_name, 'sof', tile_name+'_sof.fits')

		if cat_type == 'mof' and inj_10percent:
			fn = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_mof.fits')
		if cat_type == 'mof' and inj_10percent is False:
			fn = os.path.join(BASEPATH, 'y3v02', tile_name, 'mof', tile_name+'_mof.fits')

		if cat_type == 'coadd' and inj_10percent:
			fn = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'coadd', tile_name+'_'+filter_name+'_cat.fits')
		if cat_type == 'coadd' and inj_10percent is False:
			fn = os.path.join(BASEPATH, 'y3v02', tile_name, 'coadd', tile_name+'_'+filter_name+'_cat.fits')


	# Y3 catalogs cannot be injected #
	#TODO 'y3_gold_2_0' and 'y3_gold_2_2'
	if cat_type == 'y3_gold_2_0':
		# !!!!! User may need to alter path to Y3 Gold catalog #
		fn = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile_name+'_y3_gold_2_0.fits')
	if cat_type == 'y3_gold_2_2':
		fn = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile_name+'_y3_gold_2_2.fits')

        return fn









def get_tamu_catalog(cat_type, inj_10percent, inj_20percent, realization_number, tile_name, filter_name, inj):
	"""Get catalog for TAMU runs.

	Parameters
	----------
	inj_10percent (bool) -- Set by `INJ1_10PERCENT` or `INJ2_10PERCENT`.


	inj_20percent (bool) -- Set by `INJ1_20PERCENT` or `INJ2_20PERCENT`.

	cat_type (str) -- Allowed values: mof, sof, gal_truth, star_truth, coadd. Set by `MATCH_CAT1` or `MATCH_CAT2`.

	realization_number (str) --

	filter_name (str)

	Returns
	-------
	fn (str) -- Complete catalog filename.
	"""

	# 20% injection #	
	if inj_20percent:
		# TODO add general `inj` param to access base catalogs
		if cat_type == 'mof' and inj:
                        fn = os.path.join(BASEPATH, tile_name + '_20', 'real_' + realization_number + '_' + tile_name + '_mof.fits')
		if cat_type == 'mof' and inj is False:
			fn = os.path.join(BASEPATH, tile_name + '_20', 'base_' + tile_name + '_mof.fits')

                if cat_type == 'sof' and inj:
                        fn = os.path.join(BASEPATH, tile_name + '_20', 'real_' + realization_number + '_' + tile_name + '_sof.fits')
		if cat_type == 'sof' and inj is False:
			fn = os.path.join(BASEPATH, tile_name + '_20', 'base_' + tile_name + '_sof.fits')

                if cat_type == 'gal_truth':
                        fn = os.path.join(BASEPATH, tile_name + '_20', tile_name + '_' + realization_number + '_balrog_truth_cat_gals.fits')

		if cat_type == 'star_truth':
                        fn = os.path.join(BASEPATH, tile_name + '_20', tile_name + '_' + realization_number + '_balrog_truth_cat_stars.fits')

	# 10% injection #
	if inj_10percent:

		if cat_type == 'mof' and inj:
			fn = os.path.join(BASEPATH, tile_name, 'real_' + realization_number + '_' + tile_name + '_mof.fits')
		if cat_type == 'mof' and inj is False:
			fn = os.path.join(BASEPATH, tile_name, 'base_' + tile_name + '_mof.fits')

		if cat_type == 'sof' and inj:
			fn = os.path.join(BASEPATH, tile_name, 'real_' + realization_number + '_' + tile_name + '_sof.fits') 
		if cat_type == 'sof' and inj is False:
			fn = os.path.join(BASEPATH, tile_name, 'base_' + tile_name + '_sof.fits')

		if cat_type == 'gal_truth':
			fn = os.path.join(BASEPATH, tile_name, tile_name + '_' + realization_number + '_balrog_truth_cat_gals.fits')
		if cat_type == 'star_truth':
			fn = os.path.join(BASEPATH, tile_name, tile_name + '_' + realization_number + '_balrog_truth_cat_stars.fits')

	return fn









def matcher(realization_number, tile_name, filter_name):
        """Match two catalogs on RA and DEC with a tolerance of 1 arcsecond via STILTS.

        Parameters
	----------
	realization_number (str or int) -- Currently allowed values: 0 1 2 3 4 5 6 7 8 9 depending on the basepath.


	tile_name (str) -- Currently allowed values: DES0347-5540  DES2329-5622  DES2357-6456 DES2247-4414 depending on the basepath.

	filter_name (str)

        Returns
	-------

	outname_match (str) -- Name of catalog matched via join=1and2. Note that headers of matched catalogs will have '_1' appended or '_2' appended.

	outname_1not2 (str) -- Name of catalog matched via join=1not2.

	outname_2not1 (str) -- Name of catalof matched via join=2not1.
        """

        ### Get arguments to pass to ms_matcher ###

        # Input catalogs for STILTS #
	if MATCH_CAT1 is not 'coadd':
		in1 = get_catalog(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name, inj=INJ1)
	if MATCH_CAT1 == 'coadd':
		in1 =  get_coadd_matcher_catalog(cat_type=MATCH_CAT1, inj_10percent=INJ1_10PERCENT, inj_20percent=INJ1_20PERCENT, realization_number=realization_number, tile_name=tile_name, mag_hdr=M_HDR1, err_hdr=M_ERR_HDR1, inj=INJ1)

	if MATCH_CAT2 is not 'coadd':
		in2 = get_catalog(cat_type=MATCH_CAT2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name, inj=INJ2)
	if MATCH_CAT2 == 'coadd':
		in2 =  get_coadd_matcher_catalog(cat_type=MATCH_CAT2, inj_10percent=INJ2_10PERCENT, inj_20percent=INJ2_20PERCENT, realization_number=realization_number, tile_name=tile_name, mag_hdr=M_HDR2, err_hdr=M_ERR_HDR2, inj=INJ2)

        # !!!!! User may wish to edit directory structure. Output catalog name for STILTS #
	match_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'catalog_compare')	

	# Check for directory existence #
	if os.path.isdir(match_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(match_dir) + ' does not exist. \n Change directory structure in ms_plotter.matcher() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', match_dir, '...\n'
			os.makedirs(match_dir)

	# Filenames #
        outname_match = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match1and2.csv')
	outname_1not2 = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match1not2.csv')
	outname_2not1 = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match2not1.csv')

	# Overwrite matched catalogs if one already exists? # 
        overwrite = False 

	# Check existence #
	if os.path.isfile(outname_2not1) is False or (os.path.isfile(outname_2not1) and overwrite):


		### Matching done in ms_matcher. Parameters in1, in2, out, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2, overwrite ###
		# !!!!! Ensure that path to ms_matcher is correct #
		subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher', in1, in2, outname_match, outname_1not2, outname_2not1, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2])

        return outname_match, outname_1not2, outname_2not1









def fof_matcher(realization_number, tile_name):
        """Get catalogs to analyze. Return FOF-analysed catalogs.

        Parameters
	----------
	realization_number (str) 

	tile_name (str) 

        Returns
	-------
	ok_match OR rerun_match (str) -- Complete filename for catalogs of type join=1and2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned. 

	ok_1not2 OR rerun_1not2 (str) -- Complete filename for catalogs of type join=1not2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.

	ok_2not1 OR rerun_2not1 (str) -- Complete filename for catalogs of type join=2not1. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.
        """

        ### Filenames for input catalogs used in ms_fof_matcher ###
        # MOF or SOF #
        if MOF:
                mof = os.path.join(BASEPATH, 'y3v02', tile_name, 'mof', tile_name+'_mof.fits')
                inj_mof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_mof.fits')
                fof = os.path.join(BASEPATH, 'y3v02', tile_name, 'mof', tile_name+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_fofslist.fits')
        if MOF is False:
                mof = os.path.join(BASEPATH, 'y3v02', tile_name, 'sof', tile_name+'_sof.fits')
                inj_mof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_sof.fits')
                fof = os.path.join(BASEPATH, 'y3v02', tile_name, 'sof', tile_name+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_fofslist.fits')
        # Coadds. Using i-band #
        coadd = os.path.join(BASEPATH, 'y3v02', tile_name, 'coadd', tile_name+'_i_cat.fits')
        inj_coadd = os.path.join(BASEPATH, 'y3v02', 'balrog_images', realization_number, tile_name, 'coadd', tile_name+'_i_cat.fits')


        ### Filenames for outputs of fof_matcher ###
	os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'fof_analysis_catalog_compare')
        # Repeated #
        inj_outdir = os.path.join(outdir, tile_name, realization_number)
        inj_outname = tile_name + '_' + realization_number


        ### Check directory existence. Make directories if not present. ###
        if os.path.isdir(inj_outdir) is False:
                os.makedirs(inj_outdir)
        if os.path.isdir(os.path.join(outdir, tile_name)) is False:
                os.makedirs(os.path.join(outdir, tile_name))


	### Create filenames for output catalogs created by ms_fof_matcher ###
        fofcoadd = os.path.join(outdir, tile_name, tile_name+ '_num-match_fof_coadd.csv')
        fofgroups = os.path.join(outdir, tile_name, tile_name+ 'fofgroups.csv')
	inj_fofcoadd = os.path.join(inj_outdir, inj_outname + '_num-match_inj_fof_inj_coadd.csv')
        inj_fofgroups = os.path.join(inj_outdir, inj_outname + '_inj_fofgroups.csv')
        origfof_injfof = os.path.join(inj_outdir, inj_outname + '_inj_fofgroup_fofgroup_match1and2.csv')
        ok = os.path.join(inj_outdir, inj_outname + '.ok')
        rerun = os.path.join(inj_outdir, inj_outname + '.rerun')
        ok_inj_mof = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof.csv')
        rerun_inj_mof = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof.csv')
        ok_mof = os.path.join(inj_outdir, inj_outname + '_ok_mof.csv')
        rerun_mof = os.path.join(inj_outdir, inj_outname + '_rerun_mof.csv')

        ok_match = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1and2.csv')
	ok_1not2 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1not2.csv')
        ok_2not1 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match2not1.csv')

        rerun_match = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1and2.csv')
	rerun_1not2 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1not2.csv')
	rerun_2not1 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match2not1.csv')

        # Output directory for files made in par.py #
        parpy_outdir = os.path.join(inj_outdir, inj_outname)


        # WARNING: May need to overwrite if matching was interupted #
        overwrite = False 

        ### Check file existence of last file made in fof_matcher ###
        if os.path.isfile(rerun_match) is False or (os.path.isfile(rerun_match) and overwrite):

                ### Run fof_matcher ###
                subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher', fof, inj_fof, mof, inj_mof, coadd, inj_coadd, parpy_outdir, fofcoadd, fofgroups, inj_fofcoadd, inj_fofgroups, origfof_injfof, ok, rerun, ok_inj_mof, rerun_inj_mof, ok_mof, rerun_mof, ok_match, rerun_match, ok_1not2, rerun_1not2, ok_2not1, ok_2not1])


	if RUN_TYPE == 'ok':
		return ok_match, ok_1not2, ok_2not1
	if RUN_TYPE == 'rerun':
		return rerun_match, rerun_1not2, rerun_2not1










def stack_tiles(realization_number):
	"""Concatenate catalogs with multiple tiles and fixed realization.
	
	Parameters
	----------
	realization_number (str)

	Returns
	-------
	fn_stack_match (str) -- Complete filename for stacked catalog of type join=1and2.

	fn_stack_1not2 (str) -- Complete filename for stacked catalog of type join=1not2.

	fn_stack_2not1 (str) -- Complete filename for stacked catalog of type join=2not1.

	len(ALL_TILES) (int) -- Number of catalogs stacked.	
	"""

	# Directory for stacked catalog #
	stack_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, 'stack', realization_number)

	# Check directory existence and handle nonexistence #
	if os.path.isdir(stack_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(stack_dir) + ' does not exist. \n Change directory structure in ms_plotter. or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', stack_dir, '...\n'
			os.makedirs(stack_dir)

	# Filename for stacked catalogs #
	fn_end_1and2 = '_'.join(['stacked', realization_number, MATCH_TYPE, 'match1and2.csv'])
	fn_stack_match = os.path.join(stack_dir, fn_end_1and2)

	fn_end_1not2 = '_'.join(['stacked', realization_number, MATCH_TYPE, 'match1not2.csv'])
	fn_stack_1not2 = os.path.join(stack_dir, fn_end_1not2)

	fn_end_2not1 = '_'.join(['stacked', realization_number, MATCH_TYPE, 'match2not1.csv'])
	fn_stack_2not1 = os.path.join(stack_dir, fn_end_2not1)

	# Check if stacked tile catalog already exists #
	overwrite = False

	if os.path.isfile(fn_stack_2not1) and overwrite is False:
		print 'Stacked tile catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(fn_stack_match)
		df1not2 = pd.read_csv(fn_stack_1not2)
		df2not1 = pd.read_csv(fn_stack_2not1)

	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(fn_stack_2not1) is False or overwrite:

		all_fn_match, all_fn_1not2, all_fn_2not1 = [], [], []

		for t in ALL_TILES:

			if RUN_TYPE is None:
				fn_match, fn_1not2, fn_2not1 = matcher(realization_number=realization_number, tile_name=t, filter_name=None)

			if RUN_TYPE is not None:
				fn_match, fn_1not2, fn_2not1 = fof_matcher(realization_number=realization_number, tile_name=t)

			all_fn_match.append(fn_match); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		print 'Stacking tiles. ', len(all_fn_match), 'files ...'	
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_match))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		print 'Stacking complete ... \n'


		# Save stacked catalogs as DataFrame #
		df1and2.to_csv(fn_stack_match, sep=','); df1not2.to_csv(fn_stack_1not2, sep=','); df2not1.to_csv(fn_stack_2not1, sep=',')
		print '-----> Saving stacked tile catalogs as ', fn_stack_match
		print '----->', fn_stack_1not2
		print '----->', fn_stack_2not1

	return fn_stack_match, fn_stack_1not2, fn_stack_2not1, len(ALL_TILES)











def stack_realizations(tile_name):
	"""Concatenate catalogs with multiple realizations and fixed tile.

	Parameters
	----------
	tile_name (str) -- One stacked realization catalog created per tile.

	Returns
	-------
	fn_stack_match (str) -- Complete filename for stacked catalog of join=1and2.

	fn_stack_1not2 (str) -- Complete filename for stacked catalog of type join=1not2.

	fn_stack_2not1 (str) -- Complete filename for stacked catalog of type join=2not1. 

	len(ALL_REALIZATIONS) (int) -- Number of catalogs stacked.
	"""

	# Directory for stacked catalog #
	stack_dir = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, 'stack')

	# Check directory existence and handle nonexistence #
	if os.path.isdir(stack_dir) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(stack_dir) + ' does not exist. \n Change directory structure in ms_plotter. or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', stack_dir, '...\n'
			os.makedirs(stack_dir)

	# Filename for stacked catalogs #
	fn_stack_match = os.path.join(stack_dir, tile_name+'_stacked_'+str(MATCH_TYPE)+'_match1and2.csv')
	fn_stack_1not2 = os.path.join(stack_dir, tile_name+'_stacked_'+str(MATCH_TYPE)+'_match1not2.csv')
	fn_stack_2not1 = os.path.join(stack_dir, tile_name+'_stacked_'+str(MATCH_TYPE)+'_match2not1.csv')

	# Check if stacked realization catalog already exists #
	overwrite = False

	if os.path.isfile(fn_stack_2not1) and overwrite is False:
		print 'Stacked realization catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(fn_stack_match)
		df1not2 = pd.read_csv(fn_stack_1not2)
		df2not1 = pd.read_csv(fn_stack_2not1)


	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(fn_stack_2not1) is False or overwrite:

		all_fn_match, all_fn_1not2, all_fn_2not1 = [], [], []

		for r in ALL_REALIZATIONS:

			if RUN_TYPE is None:
				fn_match, fn_1not2, fn_2not1 = matcher(realization_number=r, tile_name=tile_name, filter_name=None)

			if RUN_TYPE is not None:
				fn_match, fn_1not2, fn_2not1 = fof_matcher(realization_number=r, tile_name=tile_name)

			all_fn_match.append(fn_match); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		print 'Stacking realizations. ', len(all_fn_match), 'files ...'
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_match))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		print 'Stacking complete ... \n'


		# Save stacked catalog as DataFrame #
		df1and2.to_csv(fn_stack_match, sep=','); df1not2.to_csv(fn_stack_1not2, sep=','); df2not1.to_csv(fn_stack_2not1, sep=',')
		print '-----> Saving stacked realization catalogs as ', fn_stack_match
	        print '----->', fn_stack_1not2
                print '----->', fn_stack_2not1

        return fn_stack_match, fn_stack_1not2, fn_stack_2not1, len(ALL_REALIZATIONS)









def get_fraction_recovered_include_flags(inj_10percent, inj_20percent, constant, df):
	"""Get fraction of injected objects recovered after matching.

	Parameters
	----------
	inj_10percent (bool) -- If `inj_10percent=True` refers to 10% Balrog-injected catalog.


	inj_20percent (bool) -- If `inj_20percent=True` refers to 20% Balrog-injected catalog.

	constant (int) -- Refers to the number of catalogs. `constant=1` unless catalogs are stacked (by tile or by realization).

	Returns
	-------
	fracRecovered (float) -- Fraction of Balrog-injections recovered from the truth catalog present in the full injected image.	
	"""

	notRecovered = df.shape[0]

	# Number of injections in a single truth catalog #
	if inj_10percent:
		# 10% injection #
		injSingleCat = 5000.0
	if inj_20percent:
		injSingleCat = 10000.0

	# Total number of injections. Differs from `injSingleCat` if catalogs were stacked #
	injTotal = injSingleCat*constant	

	'''
	### Alternative method to calculate `injTotal` ###
	# Requires: get_fraction_recovered(cat_type, inj_10percent, inj_20percent, realization_number, tile_name, df, constant)
	fn = get_catalog(cat_type=cat_type, inj_10percent=inj_10percent, inj_20percent=inj_20percent, realization_number=realization_number, tile_name=tile_name, filter_name=None)
	# Number injected is the same for all 20% realization catalogs (hdul). df2not1 is stacked. #
	hdul = fits.open(fn)
	data = hdul[1].data
	injTotal = data.shape[0]*constant
	'''

	# Fraction of objects recovered #
	fracRecovered = float(injTotal-notRecovered)/injTotal
	print 'Recovered: ', injTotal-notRecovered, '/', injTotal, '\n'

	return fracRecovered








def make_plots(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2):
	"""Makes plots.

	Parameters
	----------
	mag_hdr1, mag_hdr2 (str) -- Headers for magnitude. May be altered, hence passed as parameters.

	mag_err_hdr1, mag_err_hdr2 (str) -- Headers for magnitude error. May be altered, hence passed as parameters.
	Returns
	-------
		0
	"""

	if STACK_TILES is False: list_of_tiles_for_loop = ALL_TILES
	if STACK_REALIZATIONS is False: list_of_realizations_for_loop = ALL_REALIZATIONS


	### Stack tiles ###
	for r in ALL_REALIZATIONS:
		if STACK_TILES and STACK_REALIZATIONS is False:
			fn_stack_match, fn_stack_1not2, fn_stack_2not1, num_stack_tile = stack_tiles(realization_number=r)
			# Rewrite #
			list_of_tiles_for_loop = ['stack']



	### Stack realizations ###
	for t in ALL_TILES:

		if STACK_REALIZATIONS and STACK_TILES is False:

			fn_stack_match, fn_stack_1not2, fn_stack_2not1, num_stack_real = stack_realizations(tile_name=t)

			# Rewrite #
			list_of_realizations_for_loop = ['stack']



	### Objects recovered from truth catalog (truth catalogs ONLY) ###
	if STACK_TILES or STACK_REALIZATIONS:

		# Get DataFrame for stacked catalogs #
		if STACK_REALIZATIONS or STACK_TILES:
			df1and2 = pd.read_csv(fn_stack_match)
			df1not2 = pd.read_csv(fn_stack_1not2)
			df2not1 = pd.read_csv(fn_stack_2not1)

		# Define number of catalogs stacked #
		if STACK_TILES: const = num_stack_tile
		if STACK_REALIZATIONS: const = num_stack_real



	for t in list_of_tiles_for_loop:

		for r in list_of_realizations_for_loop:

			# Filenames for log files #
			#TODO fnFlagLog, fnMagErrLog, fnMainLog, fnColorPlotLog..forwhat, fnDeltaMagOutlierLog, fnCompletenessLog 
			fn_flag, fn_mag_bins, fn_main_log, __fn_color_plot_log, __fn_mag_outliers, fnCompletenessLog = get_log_file_names(tile_name=t, realization_number=r)
			# Write headers #
			fd_main_log, fd_mag_bins, fd_flag, __fd_color_plot_log, fdDeltaMagOutliers, fdCompletenessLog = fd_first_write(fn_main_log=fn_main_log, fn_mag_bins=fn_mag_bins, fn_flag=fn_flag, fn_color_plot_log=__fn_color_plot_log, fn_delta_mag_outliers=__fn_mag_outliers, fn_completeness_log=fnCompletenessLog)



			if STACK_REALIZATIONS is False and STACK_TILES is False:
				const = 1
				print 'Not stacking realizations and not stacking tiles ...\n'

				# Filenames for catalogs #
				if RUN_TYPE is None:
					fn_match, fn_1not2, fn_2not1 = matcher(realization_number=r, tile_name=t, filter_name=None)
				if RUN_TYPE is not None:
					fn_match, fn_1not2, fn_2not1 = fof_matcher(realization_number=r, tile_name=t)

				# Get DataFrame #
				df1and2 = pd.read_csv(fn_match)
				df1not2 = pd.read_csv(fn_1not2)
				df2not1 = pd.read_csv(fn_2not1)




			### Objects recovered from truth catalog (truth catalogs ONLY) ###
			if 'truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2:

				if 'truth' in MATCH_CAT1:
					 fractionRecovered = get_fraction_recovered_include_flags(inj_10percent=True, inj_20percent=INJ1_20PERCENT, df=df1not2, constant=const)
				if 'truth' in MATCH_CAT2:
					fractionRecovered = get_fraction_recovered_include_flags(inj_10percent=True, inj_20percent=INJ2_20PERCENT, df=df2not1, constant=const)

				# TILE \t REALIZATION \t FILTER \t PERCENT_RECOVERED #
				#fd_full_recovered_log.write(str(t) + '\t' + str(r) + '\t' + 'griz\t' + str(fractionRecovered*100) + '\n')

			if 'truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2:
				# This will not be used, is a placeholder #
				fractionRecovered = None

			
			### Region files ####
			if MAKE_REG: make_region_files(df_match=df1and2, df_1not2=df1not2, df_2not1=df2not1, realization_number=r, tile_name=t)


                        # Name for plt.savefig() #
                        fn = get_plot_save_name(realization_number=r, tile_name=t)

                        # Title for plot #
			if STACK_REALIZATIONS is False: num_stack_real = None
			if STACK_TILES is False: num_stack_tile = None
                        title = get_plot_suptitle(realization_number=r, tile_name=t, num_stack_real=num_stack_real, num_stack_tile=num_stack_tile) 


			### Handle star truth catalogs ###
			# Star truth catalogs matched then combined #
			if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
                                print 'Adding new column to matched csv ...\n'
				# mag_star #
				new_hdr = 'mag_s'

				if MATCH_CAT1 == 'star_truth':
					suf = '_1'
					mag_hdr1 = new_hdr
				if MATCH_CAT2 == 'star_truth':
					suf = '_2'
					mag_hdr2 = new_hdr

                                star_mag = get_star_mag(df=df1and2, suf=suf)
                                # New header must be of the form {base}_x where x is a single character because of the way m_axlabel is created from m_hdr #
                                df1and2.insert(len(df1and2.columns), new_hdr, star_mag)



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
                                        mag_hdr1 = new_hdr
				if 'y3_gold' in MATCH_CAT2:
                                        hdr = M_HDR2
                                        mag_hdr2 = new_hdr

                                y3_gold_mag = get_y3_gold_mag(df=df1and2, mag_hdr=hdr)
                                # Add new column to df #
                                df1and2.insert(len(df1and2.columns), new_hdr, y3_gold_mag)


	

			### Handle coadd catalogs. New column has been added with name 'mag_c'. Catalog combined then matched so has suffix (unlike star) #
			if MATCH_CAT1 == 'coadd':
				mag_hdr1 = 'mag_c_1'
				mag_err_hdr1 = 'mag_err_c_1'
			if MATCH_CAT2 == 'coadd':
				mag_hdr2 = 'mag_c_2'
				mag_err_hdr2 = 'mag_err_c_2'

			if PLOT_COLOR is False and PLOT_COMPLETENESS is False and PLOT_FLUX_HIST is False:
				delta_magnitude_subplotter(df=df1and2, flag_idx=flag_idx, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_name=fn, plot_title=title, realization_number=r, tile_name=t, fd_mag_bins=fd_mag_bins, fd_main_log=fd_main_log, fd_flag=fd_flag, fraction_recovered=fractionRecovered, fd_delta_mag_outliers=fdDeltaMagOutliers) 
			#FIXME call to normalized_delta_magnitude_subplotter()

			if PLOT_COLOR and PLOT_COMPLETENESS is False:
				for f in ['g', 'r', 'i']:
					color_subplotter(filter_name=f, df=df1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=r, tile_name=t, fd_flag=fd_flag, plot_title=title, plot_name=fn, fd_color_plot_log=__fd_color_plot_log)

			### Completeness plots (for truth catalogs vs ...) ###
			if ('truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2) and ((INJ1_10PERCENT and INJ2_10PERCENT) or (INJ1_20PERCENT and INJ2_20PERCENT)) and PLOT_COMPLETENESS and PLOT_COLOR is False:
				if PLOT_COLOR is False:
					for f in ['g', 'r', 'i', 'z']:
						completeness_magnitude_subplotter(filter_name=f, df=df1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=r, tile_name=t, plot_title=title, plot_name=fn, fd_flag=fd_flag, flag_idx=[], fd_completeness_log=fdCompletenessLog)#TODO change flag_idx?
				if PLOT_COLOR:
					print 'TODO'


			if PLOT_FLUX_HIST:
				for f in ['g', 'r', 'i', 'z']:
					norm_flux_histogram_plotter(filter_name=f, df=df1and2, flux_hdr1=CM_FLUX_HDR1, flux_hdr2=CM_FLUX_HDR2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=title) 

			### Close log files after each iteration over a realization ###
			fd_flag.close(); fd_mag_bins.close(); fd_main_log.close() 

	return 0







def get_coadd_matcher_catalog(cat_type, inj_10percent, inj_20percent, realization_number, mag_hdr, err_hdr, tile_name, inj):
	"""Make FITS file that includes a column of form '(m_g, m_r, m_i, m_z)' where m is magnitude. Column will be added to '..._i_cat.fits'. This will be used in matcher(). Relies on directory structure /`OUTDIR`/outputs/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/catalog_compare/

	Parameters
	----------
	cat_type (str) -- Catalog type. Allowed values: mof, sof, coadd, gal_truth, star_truth, y3_gold.


	inj_10percent (bool) -- If `inj_10percent=True` refers to 10% Balrog-injected catalog.

	inj_20percent (bool) -- If `inj_20percent=True` refers to 20% Balrog-injected catalog.

	realization_number (str)

	mag_hdr (str) -- Header for magnitude. Headers refer to columns in the matched catalog.

	err_hdr (str) -- Header for error. Headers refer to columns in the matched catalog.

	tile_name (str)

	Returns
	-------
	fn_new (str) -- Filename for catalog with added column. Is a FITS file.
	"""

	dir_new = os.path.join(OUTDIR, 'outputs', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'catalog_compare')
	if os.path.isdir(dir_new) is False:
		if NO_DIR_MAKE is False:
			sys.exit('Directory ' + str(dir_new) + ' does not exist...')
		if NO_DIR_MAKE:
			os.makedirs(dir_new)

	fn_new = os.path.join(dir_new, str(tile_name) + '_i_cat_combo.fits')

	# Check if new coadd catalog has already been created #
	if os.path.isfile(fn_new):
		print 'New coadd catalog already exists ...\n'


	if os.path.isfile(fn_new) is False:	
		print 'Adding a column to i-band coadd catalog. Will take a moment ...\n'

		# Get list of filenames #
		fn_griz = []
		for f in ALL_FILTERS:
			fn_griz.append(get_catalog(inj=inj, cat_type=cat_type, inj_10percent=inj_10percent, inj_20percent=inj_20percent, realization_number=realization_number, tile_name=tile_name, filter_name=f))
		fn_g, fn_r, fn_i, fn_z = fn_griz

		# Get coadd magnitude (mag_c) and magnitude error to be of form '(m_g, m_r, m_i, m_z)'. Recall that this is a string #
		mag_c, mag_err_c = get_coadd_mag_and_mag_err(fn_g=fn_g, fn_r=fn_r, fn_i=fn_i, fn_z=fn_z, mag_hdr=mag_hdr, mag_err_hdr=err_hdr)
 
	       # Create new table #
		mag_c = Column(mag_c, name='mag_c')
		mag_err_c = Column(mag_err_c, name='mag_err_c')

		# Add new table to i-band coadd catalog #
		table = Table.read(fn_i)
		table.add_column(mag_c, index=0)
       		table.add_column(mag_err_c, index=1)

		# Save new table as FITS #
		table.write(fn_new)

	return fn_new 









#TODO accept idx_good as input param? Will only be used for df_match. 
def make_region_files(df_match, df_1not2, df_2not1, realization_number, tile_name):
	"""Make DS9 region files for catalogs matched via join=1and2, join=1not2, and join=2not1 (join type is STILTS parameter set in ms_matcher or ms_fof_matcher).

	Parameters
	----------
	df* (pandas DataFrame)

	realization_number (str)

	tile_name (str)

	Returns
	-------
	fn* -- Filenames for each of the join types.
	"""

	### Get filenames and open files ###
	fn_match, fn_1not2, fn_2not1 = get_reg_names(tile_name=tile_name, realization_number=realization_number)

	overwrite = False
	
	if os.path.isfile(fn_2not1) and overwrite is False:
		print 'Region files already exist. Not overwriting ...'
	

	if os.path.isfile(fn_2not1) is False or (os.path.isfile(fn_2not1) and overwrite):
		fd_match = open(fn_match, 'w'); fd_1not2 = open(fn_1not2, 'w'); fd_2not1 = open(fn_2not1, 'w')
		# Write coordinate system #
		fd_match.write('J20000 \n'); fd_1not2.write('J20000 \n'), fd_2not1.write('J20000 \n')

		# Handle matched catalog #
		if RUN_TYPE is None:
			ra1 = RA_HDR1 + str('_1'); dec1 = DEC_HDR1 + str('_1')
			ra2 = RA_HDR2 + str('_2'); dec2 = DEC_HDR2 + str('_2')
		if RUN_TYPE is not None:
			# MOF or SOF catalogs #
			ra1 = 'ra'; dec1 = 'dec'
			ra2 = 'ra_2'; dec2 = 'dec_2' 

		### Get position. Arbitrarily using MATCH_CAT1 for RA and DEC ###
		ra_match, dec_match = df_match[ra2], df_match[dec2] 
		ra_1not2, dec_1not2 = df_1not2[ra1], df_1not2[dec1]
		ra_2not1, dec_2not1 = df_2not1[ra2], df_2not1[dec2]

		### Write to region file for matched catalog. Units are arcseconds. ###
		# Coadds allow for elliptical regions #
		if 'coadd' in (MATCH_CAT1, MATCH_CAT2) or 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
			### Get semimajor and semiminor axes (a and b, respectively) and orientation. Coadds and Y3 Gold have these values. ###
			a_match, b_match = df_match[MAJOR_AX_HDR1], df_match[MINOR_AX_HDR1]
			a_1not2, b_1not2 = df_1not2[MAJOR_AX_HDR1], df_1not2[MINOR_AX_HDR1]
			a_2not1, b_2not2 = df_2not1[MAJOR_AX_HDR2], df_2not1[MINOR_AX_HDR2]
			orientation_match, orientation_1not2, orientation_2not1 = df_match[ANGLE1], df_1not2[ANGLE1], df_2not1[ANGLE1]
			
			for i in np.arange(0, len(ra_match)):
				fd_match.write('ellipse ' + str(ra_match[i]) + ' ' + str(dec_match[i]) + ' ' + str(a_match[i]) + '" ' + str(b_match[i]) + '" ' + str(90+orientation[i]) + ' #color=green width=3 \n')
			for i in np.arange(0, len(ra_1not2)):
				fd_1not2.write('ellipse ' + str(ra_1not2[i]) + ' ' + str(dec_1not2[i]) + ' ' + str(a_1not2[i]) + '" ' + str(b_1not2[i]) + '" #color=yellow width=3 \n')
			for i in np.arange(0, len(ra_2not1)):
				fd_2not1.write('ellipse ' + str(ra_2not1[i]) + ' ' + str(dec_2not1[i]) + ' ' + str(a_2not1[i]) + '" ' + str(b_2not1[i]) + '" #color=blue width=3 \n')


		# Non-coadd catalogs allow for circular regions #
		if MATCH_CAT1 != 'coadd' and MATCH_CAT2 != 'coadd':
			size_sq_match = df_match[CM_T_HDR1]
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

	print '-----> Saving region files as: ', fn_match
	print ' -----> ', fn_1not2
	print ' ----->', fn_2not1

	return 0 








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
			YLOW, YHIGH = None, None
		if y is not None:
			YLOW, YHIGH = -1*y, y

		# Must pass constants as parameters here because used for plot name and manipulated. #
		print make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



### !!!!! Loop over possible colorbars? ###
CBAR_LOOP = False
if CBAR_LOOP:
	# Possible combinations for colorbars. Only one True allowed at a time and NORMALIZE and HEXBIN must both be False. #
	cbar_bools_list =[[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False]]
	for cbar_bools in cbar_bools_list:
		# Reset constants #
		CM_T_COLORBAR, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, NORMALIZE, HEXBIN = cbar_bools
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



RUN_TYPE_LOOP = False 
if RUN_TYPE_LOOP:
	run_type_list = [None, 'ok', 'rerun']
	for run_type in run_type_list:
		RUN_TYPE = run_type
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)


# Close files #
#fd_full_log.close()
#fd_full_recovered_log.close()
