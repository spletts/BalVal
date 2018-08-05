"""
Creates various plots for Balrog validation testing.

To run: $python ms_plotter.py base_path_to_catalogs output_directory realization tile
Example: $python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /Users/mspletts/BalVal/ 0 DES2247-4414

Relies on `stilts_matcher`. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/stilts_matcher` with the correct path to `stilts_matcher`.
FOF analysis relies on `fof_stilts_matcher`. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher` with the correct path to `ms_fof_matcher`.

Docstrings describe function parameters and function returns, but a GoogleSheet is also availabe: https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing

# Comments are ABOVE the code they correspond to (with the exception of FIXMEs and TODOs) #
In general, function input parameters are `named_with_underscores`. 
Variables defined within functions (and that should only be accessed by said function) are `__named_with_leading_underscores`.
Variable `namesWithoutSpaces` are defined by another function via `nameWithoutSpace = function(param1=x)` 

Note that `None` is typically used as a flag value throughout this script.

Megan Splettstoesser mspletts@fnal.gov
"""

#TODO Update docstrings

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.table import Table, vstack, Column
from collections import OrderedDict
from ngmix import gmix
from scipy import stats
import corner
import csv
import decimal
import esutil
import fileinput
import fitsio
import math
import matplotlib
import matplotlib.pyplot as plt
import ngmix
### To use`ngmixer` ###
# Need `matplotlib v2.2.2`. @des70 and @des71 meet this requirement #
# Run this setup script: $source /home/s1/mspletts/setup_ngmixer_gaussap.sh # 
import ngmixer
import numpy as np
import os
import pandas as pd
import subprocess
import sys
import tarfile

# From BalVal #
import calculate_injection_percent
import manipulate_catalogs 
import outputs 
### Import list of constants ###
from catalog_headers import *
from set_constants import *
import analysis
import plot_labels



### Command line args ###
# Catch error from inadequate number of command line args #
if len(sys.argv) != 5:
        sys.exit("Args: basepath (location of Balrog catalogs), output directory, realizations (can be a list of form: real1,real2,...), tiles (can be a filename.dat, or a list of form: tile1,tile2,...) \n")
# Convert `realizations` and `tiles` into lists (may be one-element list). Note 'None' entered at command line is interpreted as a str #
BASE_PATH_TO_CATS, OUTPUT_DIRECTORY, CMD_LINE_REALIZATIONS, CMD_LINE_TILES = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4].split(',')


### For directory structure ###
if BASE_PATH_TO_CATS[-1] != '/': BASE_PATH_TO_CATS = BASE_PATH_TO_CATS + '/'
# '/data/des71.a/data/kuropat/des2247-4414_sof/' --> 'des2247-4414_sof' #
BALROG_RUN = BASE_PATH_TO_CATS[BASE_PATH_TO_CATS[:-1].rfind('/')+1:-1]
# Rewrite #
if BALROG_RUN == 'Balrog': BALROG_RUN = 'TAMU_Balrog'
if BALROG_RUN == 'DES0102-4914': BALROG_RUN = 'blank_test_DES0102-4914'
# /data/des61.a/data/severett/grid_empty_bkg_skysig/y3v02/balrog_images/0/DES0347-5540 #
if BALROG_RUN == 'DES0347-5540': BALROG_RUN = 'se_DES0347-5540'

### Constants needed to loop over bands, realizations, and tiles ###
# Bands #
ALL_BANDS = [ 'g', 'r', 'i', 'z' ]

# Realization(s) #
ALL_REALIZATIONS = CMD_LINE_REALIZATIONS

# Tile(s) #
if '.dat' not in CMD_LINE_TILES[0]: ALL_TILES = CMD_LINE_TILES
# Accept .dat file of tile names at command line #
if '.dat' in CMD_LINE_TILES[0]:
	ALL_TILES = []
	for line in fileinput.input(CMD_LINE_TILES):
		# Get rid of newline character \n #
		ALL_TILES.append(line[:-1])



# User must hardcode these if multiple injection density are in the same `BASE_PATH_TO_CATALOG`. Also edit `get_catalog_filename()` # 
try:
	INJ1_PERCENT = calculate_injection_percent.get_injection_percent(cat_types=[MATCH_CAT1, MATCH_CAT2], tile=CMD_LINE_TILES[0], realization=CMD_LINE_REALIZATIONS[0], balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)
except:
	INJ1_PERCENT = 20 

try:
	INJ2_PERCENT = calculate_injection_percent.get_injection_percent(cat_types=[MATCH_CAT1, MATCH_CAT2], tile=CMD_LINE_TILES[0], realization=CMD_LINE_REALIZATIONS[0], balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)
except:
        INJ2_PERCENT = 20 



# `catch_error()` from `set_constants` #
if catch_error(cmd_line_realizations=CMD_LINE_REALIZATIONS, cmd_line_tiles=CMD_LINE_TILES) is not None:
	sys.exit('Error: ' + catch_error())

if NOTICE:
	# Which observable is being plotted? #
	ALL_OBS = {'Magnitude':PLOT_MAG, 'Color':PLOT_COLOR, 'Flux':PLOT_FLUX}
	idxo = np.where(ALL_OBS.values())[0][0]
	obs = ALL_OBS.keys()[idxo]

	lines = [
	['\ncatalog1 (name, injected?, injection percent) :\t {}, {}, {}'.format(MATCH_CAT1, INJ1, INJ1_PERCENT)],
	['catalog2 (name, injected?, injection percent) :\t {}, {}, {}'.format(MATCH_CAT2, INJ2, INJ2_PERCENT)],
	['Observable: \t {}'.format(obs)],
	['Plot type: \t {}'.format(outputs.DISPLAY)], 
	['Stack matched catalogs by tile? \t '.format(STACK_TILES)],
	['Stack matched catalogs by tile? \t '.format(STACK_REALIZATIONS)],
	['Show plot? :\t {}'.format(SHOW_PLOT)],
	['Save plot? :\t {}'.format(SAVE_PLOT)],
	['Limits for the axis of plot? :\t ({}, {})'.format(MAG_YLOW, MAG_YHIGH)],
	['If plotting magnityde, plot 1sigma_meas curve? :\t {}'.format(PLOT_MAG_ERR)],
	['If plotting 1sigma_meas curve, center about zero? (else centered about medians) :\t {}'.format(CENTER_ERR_ABT_ZERO)],
	['Normalize magnitude plot to magnitude error? :\t {}'.format(NORMALIZE)],
	['If PLOT_MAG and NORMALIZE, plot percentiles? :\t {}, {}'.format(PLOT_68P, PLOT_34P_SPLIT)],
	['Make region files? :\t {}'.format(MAKE_REG)],
	['FOF analysis? :\t {}'.format(RUN_TYPE)]
	]

	for line in lines:
		print('{:>12}'.format(*line))

	raw_input('\n Check the above before running! \n --> Press enter to proceed, control+c to stop...\n')





################################################################### Analysis ###################################################################




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

	if VERBOSE_ING: print 'Calculating {}-band magnitude error using 1.08 * (flux_cov[i][i])^0.5 / flux[i]...'.format(band)

	# Uncleaned lists for flux and flux covariance #
	full_flux = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr, band=band)
	full_flux_cov = analysis.get_matrix_diagonal_element(df=df_1and2, sq_matrices_hdr=flux_cov_hdr, band=band)

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
			print 'WARNING: cm_flux_cov is 0'

		if fluxcov[i] > 0:
			# Pogsons number = 1.08 #
			__mag_err_from_flux.append(1.08*fluxcov[i]**0.5/flux[i])

	if VERBOSE_ED:
		print ' Number of negative cm_flux_cov: {}/{}'.format(counter_neg, len(flux))

	if 'truth' in match_cat: __mag_err_from_flux = None

	return __mag_err_from_flux 




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

	if VERBOSE_ING: print 'Calculating percentiles of normalized magnitude difference...'

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

			if VERBOSE_ED:
				# Check the percentile because interpolation='lower' was used #
				num = 0
				for j in np.arange(0, len(binned_norm_mag_diff[b])):
					if abs(binned_norm_mag_diff[b][j]) <= np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'):
						num += 1
				if VERBOSE_ED: print 'Percent of objects within 68 percentile via np.percentile(interpolation=LOWER): {}...'.format(float(num)/len(binned_norm_mag_diff[b]))


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
                                plt.xlabel(r'$\Delta Mag$')
                                plt.axvline(x=0, color='black', linestyle=':', linewidth=0.5)
                                plt.show()


	return __norm_mag_diff_68p, mag_bins_for_mag_err, __neg_norm_mag_diff_34p, __pos_norm_mag_diff_34p




def bin_and_cut_measured_magnitude_error(clean_mag1, clean_mag2, mag_err1, mag_err2, band, tile, realization, fn_mag_err_log, fn_mag_diff_outliers_log):
        """Clean error. Bin error according to horizontal axis of plot. Remove error values corresponding to objects with |MagnitudeDifference|>cutoff. Do not consider error corresponding to empty bins nor bins with a small number of objects.

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

	if VERBOSE_ING: print 'Calculating magnitude error for {}-band...'.format(band)

	if 'meas' in AXLABEL1 and 'meas' not in AXLABEL2 and VERBOSE_ING:
		print 'Using measured catalog ({}) for error calculation...'.format(MATCH_CAT1)

	if 'meas' in AXLABEL2 and 'meas' not in AXLABEL1 and VERBOSE_ING:
		print 'Using measured catalog ({}) for error calculation...'.format(MATCH_CAT2)

	if 'meas' in AXLABEL1 and 'meas' in AXLABEL2 and VERBOSE_ING:
		print 'Using measured catalogs ({} & {}) for error calculation...'.format(MATCH_CAT1, MATCH_CAT2)

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

        if VERBOSE_ING: print ' Forcing magnitudes to be binned with max ', __lim_high, '...'

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
	if SWAP_ORDER_OF_SUBTRACTION:
		__vax_mag = np.array(clean_mag2) - np.array(clean_mag1)

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
					writer.writerow([tile, realization, band, clean_mag1[i], clean_mag2[i], __vax_mag[i]])

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


                if VERBOSE_ED: print ' For {}-band magnitude, number of objects in bin [{}, {}): {}...'.format(band, round(b, 2), round(b+__step, 2), __counter)


		### Write to log file ###
		if __counter == 0: write_median, write_err = None, None
		if __counter > 0: write_median, write_err = np.median(__hax_mag_in_bin), np.median(__err_mag_in_bin)
		# Append to csv #
		with open(fn_mag_err_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			# TILE, REALIZATION, BAND, NUM_OBJS_IN_BIN, BIN_LHS, BIN_RHS, MEDIAN_HAXIS_MAG, MEDIAN_ERROR #
			writer.writerow([tile, realization, band, __counter, round(b, 2), round(b+__step, 2), write_median, write_err])


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


	if VERBOSE_ED:
                if SWAP_HAX:
			print ' Binned magnitudes from {} with bin size: {}, bin minimum: {}, and bin maximum: {}...'.format(MATCH_CAT2, __step, __lim_low, __lim_high)
		if SWAP_HAX is False:
			print ' Binned magnitudes from {} with bin size: {}, bin minimum: {}, and bin maximum: {}...'.format(MATCH_CAT1, __step, __lim_low, __lim_high)
		print ' Calculated magnitude error using objects where |MagnitudeDifference| < {}'.format(__cutoff_mag_diff)
		print ' Excluded {} bins with less than {} objects...'.format(__counter_empty_bin, __min_bin_pop)
		print ' Number of objects with magnitude difference > {}: {}/{}...\n'.format(__cutoff_mag_diff_for_log, __counter_mag_diff_geq1, str(len(clean_mag1)))

        return __hax_mag_bin_medians, __mag_diff_bin_medians, __mag_err_bin_medians, __mag_bins_for_mag_err, __binned_hax_mag, __binned_mag_diff, __clean_hax_mag, __clean_mag_diff




def normalize_magnitude_difference_to_error_and_maintain_bin_structure(clean_mag1, clean_mag2, mag_err1, mag_err2, band, tile, realization, fn_mag_err_log, fn_mag_diff_outliers_log):
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

	if VERBOSE_ING: print 'Normalizing magnitude difference...\n'

	# List of lists. Stores all values in each bin #
	__binned_norm_mag_diff, __binned_hax_mag = [], []
	# Stores the median of each bin #
	__norm_mag_diff_bin_medians = []

	#__mag_err_bin_medians
	haxBinMedian, vaxBinMedian, magErrBinMedians, magBinsForMagErr, haxBins, vaxBins, hax, vax = bin_and_cut_measured_magnitude_error(clean_mag1=clean_mag1, clean_mag2=clean_mag2, mag_err1=mag_err1, mag_err2=mag_err2, band=band, tile=tile, realization=realization, fn_mag_err_log=fn_mag_err_log, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)
 

	# Loop through bins (b) #
	for b in np.arange(0, len(vaxBins)):

		# Normalized magnitude difference in each bin #
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




def normalize_magnitude_difference_to_error(binned_norm_mag_diff, mag_bins_for_mag_err, binned_hax_mag):
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
		MagnitudeDifference normalized by error. MagnitudeDifference computed via magnitude1 - magnitude2. 

        hax_mag (list of floats)
		Magnitude to be plotted on the horizontal axis.
	"""

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




def one_sigma_magnitude_counter(mag_diff, clean_mag, full_mag, mag_bins_for_mag_err, hax_mag, binned_mag_err, vax_mag_bin_medians, tile, realization, band, fn_num_objs_in_1sig_mag_log):
	"""Find the number of objects within 1sigma_mag. This function is called if `NORMALIZE` is False.

	Parameters
	----------
	mag_diff (list of floats)
		NON-normalized MagnitudeDifference.

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
		if VERBOSE_ING: print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(mag_bins_for_mag_err)-1):
                        if binned_mag_err[b] is not None:
				if VERBOSE_ING: print ' Considering objects with magnitudes (on the horizontal axis) in [{}, {})...'.format(mag_bins_for_mag_err[b], mag_bins_for_mag_err[b+1])
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= mag_bins_for_mag_err[b] and hax_mag[i] < mag_bins_for_mag_err[b+1]:
						counter_objs_considered += 1
						if mag_diff[i] < binned_mag_err[b] + vax_mag_bin_medians[b] and mag_diff[i] >= -1.0*binned_mag_err[b] + vax_mag_bin_medians[b]: 
							__num_1sigma_mag += 1


	### Log ###
	with open(fn_num_objs_in_1sig_mag_log, 'a') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow([tile, realization, band, len(full_mag), len(full_mag)-len(clean_mag) , __num_1sigma_mag, str(NORMALIZE), str(RUN_TYPE)])
        csvfile.close()


	if VERBOSE_ED:
		print ' NOT Normalized '
		print ' Total objects: ', tot
		print ' Number of objects considered:', counter_objs_considered
		print ' Number of objects within 1sigma_mag: ', __num_1sigma_mag 
		print ' Fraction within 1sigma_mag: ', float(__num_1sigma_mag)/counter_objs_considered, '\n'

	# return a percent
	return 100.0*float(__num_1sigma_mag)/counter_objs_considered 




def one_sigma_magnitude_counter_for_normalized_magnitude_difference(full_mag, norm_mag_diff, clean_mag, mag_bins_for_mag_err, hax_mag, mag_err_bin_medians, norm_vax_mag_bin_medians, tile, realization, band, fn_num_objs_in_1sig_mag_log):
	"""Find the number of objects within 1sigma_mag using the normalized (to error) MagnitudeDifference. This function is only called if `NORMALIZE` is True.

	Parameters
	----------
	norm_mag_diff (list of floats)
		Normalized MagnitudeDifference. #FIXME what about colors? 

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
		if VERBOSE_ING: print 'Centering 1sigma_mag about vax median of each bin... \n'
		for b in np.arange(0, len(mag_bins_for_mag_err)-1):
			if mag_err_bin_medians[b] is not None:
				if VERBOSE_ING: print ' Considering objects with magnitudes (on the horizontal axis) in [', mag_bins_for_mag_err[b], ',', mag_bins_for_mag_err[b+1], ')'
                                for i in np.arange(0, len(hax_mag)):
                                        if hax_mag[i] >= mag_bins_for_mag_err[b] and hax_mag[i] < mag_bins_for_mag_err[b+1]:
                                                counter_objs_considered += 1
                                                if norm_mag_diff[i] < 1 + norm_vax_median[b] and norm_mag_diff[i] >= -1.0 + norm_vax_median[b]: 
							__num_1sigma_mag_norm += 1


	### Log ###
	with open(fn_num_objs_in_1sig_mag_log, 'a') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
                writer.writerow([tile, realization, band, len(full_mag), len(full_mag)-len(clean_mag) , __num_1sigma_mag_norm, str(NORMALIZE), str(RUN_TYPE)])
        csvfile.close()


	if VERBOSE_ED:
		print ' Normalized '
		print ' Total objects: ', tot
		print ' Number of objects considered:', counter_objs_considered
		print ' Number of objects within 1sigma_mag: ', __num_1sigma_mag_norm 
		print ' Fraction within 1sigma_mag: ', float(__num_1sigma_mag_norm)/counter_objs_considered, '\n'

	#Return a percent
	return 100.0* float(__num_1sigma_mag_norm)/counter_objs_considered 




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

        if CM_T_CBAR: __cbar_label = pref + cm_t_hdr[:-2] + '_' + str(meas_or_true_cat)
        if CM_T_ERR_CBAR: __cbar_label = pref + cm_t_err_hdr + '_' + str(meas_or_true_cat)


	### Colorbar value ###
        if CM_T_ERR_CBAR:
                # For measured catalog, cuts performed on truth catalogs #
                cbarData = analysis.get_good_data(df_1and2=df_1and2, hdr=cm_t_err_hdr, idx_good=idx_good, str_of_arr=False, band=None)

        if CM_T_CBAR:
                cbarData = analysis.get_good_data(df_1and2=df_1and2, hdr=cm_t_hdr, idx_good=idx_good, str_of_arr=False, band=None)


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
				__mag_error = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_err_hdr, band=band)
			if match_cat == 'star_truth' or 'y3_gold' in match_cat: 
                                __mag_error = df_1and2[str(mag_err_hdr[:-2]) + '_' + band.upper() + str(mag_err_hdr[-2:])]
			if match_cat == 'base':
				try:
					__mag_error = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=NEW_BASE_MAG_ERR_GRIZ_HDR+'_1', band=band)
				except:
					__mag_error = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=NEW_BASE_MAG_ERR_GRIZ_HDR+'_2', band=band)
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
			magErrorFollow = analysis.get_floats_from_string(df=df, four_elmt_arrs_hdr=mag_err_hdr, band=BANDS_FOR_COLOR[band])
			if match_cat == 'star_truth' or 'y3_gold' in match_cat:
				magErrorFollow = df[str(mag_err_hdr[:-2]) + '_' + BANDS_FOR_COLOR[band].upper() + str(mag_err_hdr[-2:])]

	__color_error = (np.array(magErr)**2 + np.array(magErrorFollow)**2)**0.5
	'''

	__color_error = (np.array(magErrBand1)**2 + np.array(magErrBand2)**2)**0.5

	return __color_error






def get_color_from_binned_magnitude(df_1and2, mag_hdr1, mag_hdr2, clean_mag1_band1, clean_mag2_band1, band, idx_good):
	"""Get color with magnitudes binned via [21-22), [22-23), [23-24), [24-25). TODO mag_true 

        Parameters
        ----------
        df_1and2 (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.

        mag_hdr1, mag_hdr2 (str) 
		Matched (via join=1and2) catalog headers for magnitude for `MATCH_CAT1`, `MATCH_CAT2`.

        clean_magnitude_a (list of floats)
                List of magnitudes with objects with flags* removed. 

        band (str)

        idx_good (list of ints)
                Indices of objects without flags. 

        Returns
        -------
	"""

	# From `MATCH_CAT1` #
	__clean_mag1_band2 = np.array(analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=BANDS_FOR_COLOR[band]))[idx_good]	
	# From `MATCH_CAT2` #
	__clean_mag2_band2 = np.array(analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=BANDS_FOR_COLOR[band]))[idx_good] 

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







#FIXME this function is not done
def get_stacked_magnitude_differences(realization, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2):
	"""TODO

	Parameters
	----------

	Returns
	-------
	"""

	__mag_diffs = np.empty(len(ALL_TILES), len(ALL_BANDS))

	# First axis: tiles #
	for i in np.arange(0, len(ALL_TILES)):
		t = ALL_TILEs[i]

		__df1and2 = get_dataframe_and_headers(realization=realization, tile=t, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)[0] 

		for j in np.arange(0, len(ALL_BANDS)):
			b = ALL_BANDS[j]

			magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=b, df_1and2=__df1, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=t, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot, fn_percent_recovered_log=fn_percent_recovered_log)

			#FIXME not sure if this will work bc __completeness has np.empty(x,y,z) and sets elements with [i][j]
			__mag_diffs[i][j] = np.array(cleanMag2) - np.array(cleanMag1)

		#FIXME need a way to plot error and bin error

	return 0	




def stacked_magnitude_completeness_subplotter(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, fn_flag_log, fn_mag_completeness_log, fn_percent_recovered_log):
	"""TODO

	Parameters
	----------

	Returns
	-------
	"""

	if STACK_TILES and VERBOSE_ING: print 'Stacking completeness for multiple tiles ...'
	if STACK_REALIZATIONS and VERBOSE_ING: print 'TODO' #print 'Stacking completeness for multiple realizations ...'

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
                        __df1 = get_dataframe_and_headers(realization=realization, tile=t, inj1=True, inj2=True, inj1_percent=10, inj2_percent=10, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

			# Fill second axis: band #
			for j in np.arange(0, len(ALL_BANDS)):
				b = ALL_BANDS[j]

                                magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1a, magAxLabel2a, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(fn_percent_recovered_log=fn_percent_recovered_log, band=b, df_1and2=__df1, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=t, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

                                completenessFraction1 = get_magnitude_completeness(inj_percent=10, idx_good=idxGood, df_1and2=__df1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=t, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]
			
				# Fill third axis: completeness (list) #
				__completeness1[i][j] = completenessFraction1


		### Get DataFrame for tile with 20% injection ###
		if VERBOSE_ING: print 'Getting DataFrame for tile: {}...'.format(t)

		#TODO df_b, completeness_b
		__df2 = get_dataframe_and_headers(realization=realization, tile=t, inj1=True, inj2=True, inj1_percent=20, inj2_percent=20, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

		# Fill second axis: band #
		for j in np.arange(0, len(ALL_BANDS)):
			b = ALL_BANDS[j]

			#TODO rename 21, 22, 21, et
                        magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1b, magAxLabel2b, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=b, fn_percent_recovered_log=fn_percent_recovered_log, df_1and2=__df2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=t, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

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
                        plt.ylabel('$\\bf{%s}$-band Magnitude Completeness' % ALL_BANDS[k])
			plt.xlabel(magAxLabel1a)
			if SWAP_HAX: plt.xlabel(magAxLabel2a)
                        plt.grid(linestyle='dotted')

                        plt.subplot(1, 2, 2)

                # Get completeness mean (amongst all tiles) for each magnitude bin --> `axis=0` #
                # Note: plot() ignores nans so there is no need to remove them from data in sep step #
		# `b` labels for secondary completeness #
                plt.plot(COMPLETENESS_PLOT_MAG_BINS, np.nanmean(__completeness2[:,k,:], axis=0), color='green')
		plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
                plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
                plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
                plt.title('20% Injection') #TODO automate this?
                plt.ylabel('$\\bf{%}$-band Magnitude Completeness' % ALL_BANDS[k])
		plt.xlabel(magAxLabel1b)
		if SWAP_HAX: plt.xlabel(magAxLabel2b)
                plt.grid(linestyle='dotted')

		plot_title = plot_title.replace('10% ', '')
		plot_title = plot_title.replace('20% ', '')

                plt.suptitle(plot_title, fontweight='bold') #FIXME check if there is a sigma cutoff

		if SHOW_PLOT: plt.show()

	return 0




def magnitude_completeness_subplotter(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, tile, fn_flag_log, fn_mag_completeness_log, fn_percent_recovered_log):
	"""Creates two completeness plots with two panels of completeness for 10% injections (first subplot panel) and 20% injections (second subplot panel).
	This function is only to be called when comparing Balrog-injected catalogs and truth catalogs.

	Parameters
	----------
	
	Returns
	-------
	0
	"""

	__completeness_griz1, __completeness_griz2 = [], []
	__bin_median_griz1, __bin_median_griz2 = [], []
	__mag_axlabel_griz_1a, __mag_axlabel_griz_2a, __mag_axlabel_griz_1b, __mag_axlabel_griz_2b = [], [], [], []


	# COSMOS `BALROG_RUN`s so far (Jun 21 2018) have only 20% injections #
	if 'COSMOS' not in BALROG_RUN: 

		# Get completeness for 10% inj via `inj_percent=10` #
		#FIXME instances where data_frame depends on band:
		__df1 = get_dataframe_and_headers(realization=realization, tile=tile, inj1=True, inj2=True, inj1_percent=10, inj2_percent=10, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

		for b in ALL_BANDS:
			#TODO this will be stuck at value for 'z'
			magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1a, magAxLabel2a, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=b, fn_percent_recovered_log=fn_percent_recovered_log, df_1and2=__df1, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot) #TODO this does not have inj_% parameter so axlabel is wrong

			__mag_axlabel_griz_1a.append(magAxLabel1a), __mag_axlabel_griz_2a.append(magAxLabel2a)

			completenessFraction1 = get_magnitude_completeness(inj_percent=10, idx_good=idxGood, df_1and2=__df1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]

			__completeness_griz1.append(completenessFraction1)


	# Get completeness for 20% injection via `inj_percent=20` #
	__df2 = get_dataframe_and_headers(realization=realization, tile=tile, inj1=True, inj2=True, inj1_percent=20, inj2_percent=20, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2)[0]

	for b in ALL_BANDS:
		magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1b, magAxLabel2b, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=b, fn_percent_recovered_log=fn_percent_recovered_log, df_1and2=__df2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot)

		__mag_axlabel_griz_1b.append(magAxLabel1b), __mag_axlabel_griz_2b.append(magAxLabel2b)

		completenessFraction2 = get_magnitude_completeness(inj_percent=20, df_1and2=__df2, idx_good=idxGood, clean_mag1=cleanMag1, clean_mag2=cleanMag2, full_mag1=fullMag1, full_mag2=fullMag2, tile=tile, realization=realization, band=b, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err1=magErr1, mag_err2=magErr2, fn_mag_completeness_log=fn_mag_completeness_log)[0]

		if VERBOSE_ING: print 'Getting magnitude completeness for {}-band...\n'.format(b)

		__completeness_griz2.append(completenessFraction2)



	if VERBOSE_ING: print 'Plotting magnitude completeness...\n'

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
			plt.ylabel('$\\bf{%s}$-band Magnitude Completeness' % ALL_BANDS[i])
			plt.xlabel(__mag_axlabel_griz_1a[i])
			if SWAP_HAX: plt.xlabel(__mag_axlabel_griz_1a[i])
			#plt.xlabel('10%_inj_cm_mag_$\\bf{'+str(ALL_BANDS[i])+'}$_true') #TODO automate this?
			plt.grid(linestyle='dotted')

			plt.subplot(1, 2, 2)

		plt.plot(COMPLETENESS_PLOT_MAG_BINS, __completeness_griz2[i], color='green')
		plt.axhline(y=1, color='black', linestyle='-', linewidth=0.7)
		plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
		plt.axhline(y=0.9, color='orange', linestyle='--', linewidth=0.7)
		plt.title('20% Injection') #TODO automate this?
		plt.ylabel('$\\bf{%s}$-band Magnitude Completeness' % ALL_BANDS[i])
		plt.xlabel(__mag_axlabel_griz_1b[i])
		if SWAP_HAX: plt.xlabel(__mag_axlabel_griz_2b[i])
		plt.grid(linestyle='dotted')

		plot_title = plot_title.replace('10% ', '')
		plot_title = plot_title.replace('20% ', '')

		plt.suptitle(plot_title, fontweight='bold') #FIXME check if there is a sigma cutoff

		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', ALL_BANDS[i]); print '-----> Saving plot as: ', fn_plot; plt.savefig(fn_plot)

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


	if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', WRITE_COLORS[band]); print '-----> Saving plot as: ', fn_plot; plt.savefig(fn_plot) 

	return 0
'''



def get_magnitude_completeness(idx_good, df_1and2, clean_mag1, clean_mag2, full_mag1, full_mag2, tile, realization, band, mag_hdr1, mag_hdr2, mag_err1, mag_err2, fn_mag_completeness_log, inj_percent): #FIXME inj?
	"""TODO after function is done
	Introduce cutoff for bin size.
	"""

	#TODO
	print '`df_1and2` and `full_mag` NOT BEING USED IN `get_magnitude_completeness()`'

	if VERBOSE_ING: print 'Calculating magnitude completeness of {}% injected catalogs ...'.format(inj_percent)


	__mag_completeness, __plot_bins = [], []
	mag_bin_l = 21

	### Read truth catalog which has magnitudes in form (m_g, m_r, m_i, m_z) ###
	if 'truth' in MATCH_CAT1:
		fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT1, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ1)	
		mag_hdr = mag_hdr1
		# Use error in measured catalogs only #
		mag_err1 = None 
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
		flag_hdr = FLAGS_HDR1[:-2]
		cm_flag_hdr = CM_FLAGS_HDR1[:-2] 

	if 'truth' in MATCH_CAT2:
		fn_truth_cat = get_catalog_filename(cat_type=MATCH_CAT2, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ2)
		mag_hdr = mag_hdr2
		ra_hdr, dec_hdr = RA_HDR2, DEC_HDR2
		# Use error in measured catalogs only #
		mag_err2 = None 
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
                flag_hdr = FLAGS_HDR2[:-2]
                cm_flag_hdr = CM_FLAGS_HDR2[:-2]

	data = fits.open(fn_truth_cat)[1].data
	truth_mag_griz = data[mag_hdr[:-2]]
	__truth_mag = []

	for mag_griz in truth_mag_griz:
		__truth_mag.append(mag_griz[BAND_INDEX[band]])
	__truth_mag = np.array(__truth_mag)

	flag = data[flag_hdr]
	cm_flag = data[cm_flag_hdr]

	# Primary flag cuts to truth catalog #
	__idx_good = np.where( (abs(__truth_mag) != 9999.0) & (abs(__truth_mag) != 99.0) & (abs(__truth_mag) != 37.5) & (flag == 0) & (cm_flag == 0) )[0]
	__truth_mag = __truth_mag[__idx_good]

	if VERBOSE_ING: print '{} flagged objects removed from {} for {}-band...'.format(len(__idx_good), fn_truth_cat[fn_truth_cat[:-1].rfind('/')+1:-1] , band)


	### Photometry cuts on matched catalog (|MagnitudeDifference| < 3sigma) where sigma refers to measured catalog only ###
        match_mag1, match_mag2 = [], []

	__cutoff_mag_diff = 2
	print 'Using objects with |MagnitudeDifference| < {}...'.format(__cutoff_mag_diff)
	
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

	if VERBOSE_ING: print ' Including bins with more than {} objects...'.format(__bin_size_cutoff)

	### Log ###
	for i in np.arange(0, len(COMPLETENESS_MAG_BINS)-1):
		if numer[i] != 0 and denom[i] != 0: #TODO better way to eliminate the large range of COMPLETENESS_MAG_BINS
			# Append to csv #
			with open(fn_mag_completeness_log, 'a') as csvfile:
				writer = csv.writer(csvfile, delimiter=',')
				# TILE, REALIZATION, BAND, INJ_PERCENT, TRUTH_MAG_BIN_L, TRUTH_MAG_BIN_R, MATCH_CAT_OBJS_IN_BIN, TRUTH_CAT_OBJS_IN_BIN #
				writer.writerow([tile, realization, band, inj_percent, COMPLETENESS_MAG_BINS[i], COMPLETENESS_MAG_BINS[i+1], numer[i], denom[i]])


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







def normalized_flux_difference_histogram_subplotter(fn_gauss_aper_log, fn_plot, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title, tile, realization, fn_flux_sigma_clip_log, gap_flux_hdr1, gap_flux_hdr2, fn_percent_recovered_log):
	"""Create 2x2 plot grid of normalized (area under the curve equals 1) 1D FluxDifference/SigmaFlux histograms. One subplot is created for each griz band. Each plot includes a standard Gaussian (mean=0, standard_deviation=1) and the best-fit Gaussian to FluxDifference/SigmaFlux. 

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

	if PLOT_GAUSS_APER_FLUX:
		### Write headers to csv outside of loop over bands ###
		with open(fn_gauss_aper_log, 'wb') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow(['TILE', 'REALIZATION', 'BAND', 'FLUX_MEAS', 'FLUX_TRUE', 'GAUSS_APER_FLUX_MEAS', 'GAUSS_APER_FLUX_TRUE', 'FLUX_ERR_MEAS', 'GAUSS_APER_FLUX_DIFF/FLUX_ERR_MEAS'])


	### Create 4-by-4 subplot ###
        counter_subplot = 1
        # Figure size units: inches #
        plt.figure(figsize=(12, 10))


        ### Create one subplot for each griz band ###
        for b in ALL_BANDS:
		plt.subplot(2, 2, counter_subplot)

		# Create one subplot #
		if PLOT_CM_FLUX:
			percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = normalized_flux_difference_histogram_plotter(tile=tile, realization=realization, band=b, df_1and2=df_1and2, flux_hdr1=flux_hdr1, flux_hdr2=flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=plot_title, fn_flux_sigma_clip_log=fn_flux_sigma_clip_log, fn_percent_recovered_log=fn_percent_recovered_log, label_hist_raw='cm', label_fit_raw='cm', label_hist_sig_clip='{}-sig_clip_cm'.format(int(N)), label_fit_sig_clip='sig_clip_cm', color_hist=FLUX_HIST[b], color_fit=FIT_TO_FLUX[b])

		if PLOT_GAUSS_APER_FLUX:
			percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = normalized_flux_difference_histogram_plotter(tile=tile, realization=realization, band=b, df_1and2=df_1and2, flux_hdr1=gap_flux_hdr1, flux_hdr2=gap_flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=plot_title, fn_flux_sigma_clip_log=fn_flux_sigma_clip_log, fn_percent_recovered_log=fn_percent_recovered_log, label_hist_raw='gap', label_fit_raw='gap', label_hist_sig_clip='{}-sig_clip_gap'.format(int(N)), label_fit_sig_clip='sig_clip_gap', color_hist=GAP_FLUX_HIST[b], color_fit=FIT_TO_GAP_FLUX[b])

		counter_subplot += 1

	# The percent recovered with flags included is the title because this is the same for all griz bands. However different bands have a different number of flags #
	if percentRecoveredFlagsIncluded is not None:
		plot_title = '{} Recovered: {}%'.format(plot_title, round(percentRecoveredFlagsIncluded, 2))
	plt.suptitle(plot_title, fontweight='bold')
	plt.subplots_adjust(hspace=0.3)
	#plt.tight_layout(pad=1, h_pad=1, w_pad=1)
	plt.text(0.1, 0.01, BASE_PATH_TO_CATS, fontsize=10, transform=plt.gcf().transFigure)

	if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as: ', fn_plot; plt.savefig(fn_plot)

	if SHOW_PLOT: plt.show()

	return 0




def normalized_flux_difference_histogram_plotter(tile, realization, band, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_title, fn_flux_sigma_clip_log, fn_percent_recovered_log, label_hist_raw, label_fit_raw, label_hist_sig_clip, label_fit_sig_clip, color_hist, color_fit):
	"""Create a normalized (area under curve equals 1) 1D histogram of FluxDifference/SigmaFlux for a particular `band`. Plot includes a standard Gaussian (mean=0, standard_deviation=1) and the best-fit Gaussian to FluxDifference/SigmaFlux.

	Parameters
	----------

	Returns
	-------
	0
	"""

	#__flux_hist_bins = 10**2
	# Equal bin size #
	__flux_hist_bins = np.linspace(FLUX_XLOW, FLUX_XHIGH, 10**2)

        normFluxDiff, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved, haxLabel = get_flux_plot_variables(tile=tile, realization=realization, band=band, df_1and2=df_1and2, flux_hdr1=flux_hdr1, flux_hdr2=flux_hdr2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, fn_percent_recovered_log=fn_percent_recovered_log)


        if VERBOSE_ING: print 'Plotting {} points for {}-band flux...'.format(len(normFluxDiff), band)


        if RAW_NORM_FLUX_DIFF:
                mean, stddev = stats.norm.fit(normFluxDiff)
                x_arr = np.linspace(np.min(normFluxDiff), np.max(normFluxDiff), len(normFluxDiff))
                fit_data = stats.norm.pdf(x_arr, mean, stddev)
                # Plot fit #
		if PLOT_GAUSSIAN_FIT:
			plt.plot(x_arr, fit_data, color=color_fit, label='{}:\n$\mu=${}, $\sigma$={}'.format(label_fit_raw, round(mean, 2), round(stddev, 2)))
                # Plot histogram #
                plt.hist(normFluxDiff, __flux_hist_bins, density=NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY, histtype='step', color=color_hist, label=label_hist_raw, rwidth=1)


	if SIGMA_CLIP_NORM_FLUX_DIFF:
		### N-sigma clip ###
		try:
			masked_clipped_norm_flux_diff = sigma_clip(normFluxDiff, sigma=N)
		except:
			masked_clipped_norm_flux_diff = sigma_clip(normFluxDiff, sig=N)
		# Remove masked objects that did not survive clipping #
		clipped_norm_flux_diff = masked_clipped_norm_flux_diff.data[~masked_clipped_norm_flux_diff.mask]

		__num_clipped = len(masked_clipped_norm_flux_diff.data[masked_clipped_norm_flux_diff.mask])

		with open(fn_flux_sigma_clip_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			# 'TILE', 'REALIZATION', 'BAND', 'GAUSSIAN_APER_APPLIED', 'NUM_OBJS_FLAGS_RM', 'NUM_OBJS_CLIPPED' #
			writer.writerow([tile, realization, band, PLOT_GAUSS_APER_FLUX, len(normFluxDiff), __num_clipped])

                        if VERBOSE_ED: print ' Number of objects {}-sigma clipped for {}-band: {}'.format(int(N), band, __num_clipped)

                        ### Fit N-sigma clipped objects to normalized Gaussian http://danielhnyk.cz/fitting-distribution-histogram-using-python/ ###
                        sig_clip_mean, sig_clip_stddev = stats.norm.fit(clipped_norm_flux_diff)
                        sig_clip_x_arr = np.linspace(np.min(clipped_norm_flux_diff), np.max(clipped_norm_flux_diff), len(clipped_norm_flux_diff))
                        sig_clip_fit_data = stats.norm.pdf(sig_clip_x_arr, sig_clip_mean, sig_clip_stddev)
                        # Plot fit #
			if PLOT_GAUSSIAN_FIT:
				plt.plot(sig_clip_x_arr, sig_clip_fit_data, color=color_fit, label='{}:\n$\mu=${}, $\sigma$={}'.format(label_fit_sig_clip, round(sig_clip_mean, 2), round(sig_clip_stddev, 2)))
                        # Plot histogram #
			plt.hist(clipped_norm_flux_diff, __flux_hist_bins, density=NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY, histtype='step', color=color_hist, label=label_hist_sig_clip)


	if TRIM_NORM_FLUX_DIFF:

		if VERBOSE_ING: print 'Trimming normalized {}-band flux difference...'.format(band)

		# Returns value (not index) of 2nd and 98th percentiles #
		__val_2p, __val_98p = np.percentile(normFluxDiff, [2, 98], interpolation='nearest')
		normFluxDiff= np.sort(normFluxDiff)
		idx1, idx2 = normFluxDiff.tolist().index(__val_2p), normFluxDiff.tolist().index(__val_98p)
		# Rewrite. Note this is done AFTER sigma clipping #
		normFluxDiff= normFluxDiff[idx1:idx2]

		# Check for non-unique idx #
		counter1, counter2 = 0, 0
		for elem in normFluxDiff:
			if elem == __val_2p: counter1 += 1
			if elem == __val_98p: counter2 += 1
		if counter1 > 1 or counter2 > 1: print 'ERROR: non-unique index.'

		mean, stddev = stats.norm.fit(normFluxDiff)
                x_arr = np.linspace(np.min(normFluxDiff), np.max(normFluxDiff), len(normFluxDiff))
                fit_data = stats.norm.pdf(x_arr, mean, stddev)
                # Plot fit #
		if PLOT_GAUSSIAN_FIT:
			plt.plot(x_arr, fit_data, color=color_fit, label='trim_{}:\n$\mu=${}, $\sigma$={}'.format(label_fit_raw, round(mean, 2), round(stddev, 2)))
                # Plot histogram #
                plt.hist(normFluxDiff, __flux_hist_bins, density=NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY, histtype='step', color=color_hist, label='trim_{}'.format(label_hist_raw), rwidth=1)



	# Note: this is the peak of the raw distribution #
	if PLOT_PEAK:
		### Get median of data. Get x-value corresponding to peak of data distribution ###
		__bin_size, __bin_edges = np.histogram(normFluxDiff, __flux_hist_bins, density=NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY)
		# Find bin with max counts #
		__max_val = max(__bin_size)
		for i in np.arange(0, len(__bin_size)):
			if __bin_size[i] == __max_val:
				__ymax_idx = i

		if VERBOSE_ED:
			print 'For {}-band:'.format(band)
			print 'Corresponding hax value of max vax value:', __bin_edges[__ymax_idx]
			print 'Data median: {}\n'.format(np.median(normFluxDiff))

		data_median = np.median(normFluxDiff)
		# Middle of bin if `bin` is given by an array #
		peak = __bin_edges[__ymax_idx] + (1.0*__flux_hist_bins[1]-__flux_hist_bins[0])/2.0
		plt.axvline(x=data_median, color=color_fit, linestyle='--', label='{}_median: {}'.format(label_fit_raw, round(data_median, 2)))
		plt.axvline(x=peak, color=color_fit, linestyle=':', label='{}_peak: {}'.format(label_fit_raw, round(peak, 2))) 


        ### Improve plot readability ###
        # Plot standard Gaussian #
        gx = np.linspace(plt.xlim()[0], plt.xlim()[1], 5*10**3)
        plt.plot(gx, gaussian(x=gx, mu=0.0, sig=1.0), linestyle=':', color='black', linewidth=0.7)#, label=r'$\mu=0 \,, \, \sigma=1$')
        # Plot guide for eye #
        plt.axvline(x=0, linestyle=':', color='black', linewidth=0.7)
        ### Labels ###
        plt.xlabel(haxLabel)
        if NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY:
                plt.ylabel('Normalized Count')
        if NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY is False:
                plt.ylabel('Count')
        plt.legend(fontsize=10).draggable()
        if FLUX_XLOW is not None and FLUX_XHIGH is not None:
                plt.xlim([FLUX_XLOW, FLUX_XHIGH])


	return percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved




def color_subplotter(band, df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization, tile, fn_flag_log, plot_title, fn_plot, fn_color_log, fn_percent_recovered_log):
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

	axLabel1 = plot_labels.get_color_axlabel(inj_percent=INJ1_PERCENT, inj=INJ1, meas_or_true_cat=AXLABEL1, match_cat=MATCH_CAT1, band=band)
	axLabel2 = plot_labels.get_color_axlabel(inj_percent=INJ2_PERCENT, inj=INJ2, meas_or_true_cat=AXLABEL2, match_cat=MATCH_CAT2, band=band)
	# Assumes `axlabel_a` minus `axlabel_b` #
	vaxColorDiffLabel = plot_labels.get_short_difference_axlabel(axlabel_a=axLabel1, axlabel_b=axLabel2, band=band)
	if SWAP_ORDER_OF_SUBTRACTION:
		# Assumes `axlabel_a` minus `axlabel_b` so set `axlabel_a=axLabel2` to reverse order of subtraction in label #
		vaxColorDiffLabel = plot_labels.get_short_difference_axlabel(axlabel_a=axLabel2, axlabel_b=axLabel1, band=band)

	### Get plot data ###
	errMag1, errMag2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=band, df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot, fn_percent_recovered_log=fn_percent_recovered_log)

	#TODO get_magnitude_error is called twice
	'''
	if 'truth' not in MATCH_CAT1:
		err1 = get_color_plot_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, flux_cov_hdr=CM_FLUX_COV_HDR1, df=df, band=band, idx_good=idxGood, match_cat=MATCH_CAT1)
	if 'truth' not in MATCH_CAT2:
		err2 = get_color_plot_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, flux_cov_hdr=CM_FLUX_COV_HDR2, df=df, band=band, idx_good=idxGood, match_cat=MATCH_CAT2)
	'''

	if band != 'z':
		cleanColor1, cleanColor2, magBinsForColor = get_color_from_binned_magnitude(df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, band=band, idx_good=idxGood, clean_mag1_band1=cleanMag1, clean_mag2_band1=cleanMag2)



	### Plot ###
	# Assuming `SWAP_HAX=False` and `SWAP_ORDER_OF_SUBTRACTION=False` #
	__hax_label = axLabel1
	__hax_color = np.array(cleanColor1)
	__subtitle_pref = magAxLabel1
	__vax_color = np.array(cleanColor1) - np.array(cleanColor2)

	if SWAP_HAX:
		__hax_label = axLabel2
		__hax_color = np.array(cleanColor2)
		__subtitle_pref = magAxLabel2
	if SWAP_ORDER_OF_SUBTRACTION:
		__vax_color = np.array(cleanColor2) - np.array(cleanColor1)


	__lw = 1.1

        ### Create 4-by-4 subplot. Figure size units: inches ###
	plt.figure(figsize=(12, 10))


        for i in np.arange(0, len(__hax_color)):

		# Write to log file with headers 'TILE REALIZATION COLOR MAG_BIN OBJS_IN_MAG_BIN TOTAL_OBJS_FLAGS_INC RUN_TYPE #
		__write_bin = '['+str(magBinsForColor[i])+','+str(magBinsForColor[i+1])+')'
		# Append to csv #
		with open(fn_color_log, 'a') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow([tile, realization, WRITE_COLORS[band], __write_bin, len(__hax_color[i]), len(cleanMag1), str(RUN_TYPE)])
		csvfile.close()

		if VERBOSE_ING: print 'Plotting {} objects for {}...'.format(len(__hax_color[i]), WRITE_COLORS[band])

                plt.subplot(2, 2, i+1)


		if SCATTER:
			plt.scatter(__hax_color[i], __vax_color[i], color=PT_COLORS[band], alpha=0.25, s=0.25) 
			#TODO symmetric y axis?


		if CORNER_HIST_2D:
			# Create bins of 1/4 magnitude for hax # 
			__bin_x = np.linspace(np.min(__hax_color[i]), np.max(__hax_color[i]), math.ceil(4.0*(np.max(__hax_color[i]) - np.min(__hax_color[i]))))
			# Create bins of 1/20 magnitude for vax #
			__ylow = np.min(np.array(__vax_color[i])) # -1*np.max(abs(np.array(__vax_color[i]))) ?
			__yhigh = np.max(np.array(__vax_color[i])) # np.max(abs(np.array(__vax_color[i]))) ?
			__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))

			# Force symmetric vertical axis #
			__sym_ylim = np.mean([abs(__ylow), __yhigh])
			if abs(abs(__ylow) - __yhigh) > 1:
				__sym_ylim = np.min([abs(__ylow), __yhigh])

			# Plot density bins #
                        #corner.hist2d(__hax_color[i], __vax_color[i], bins=np.array([__bin_x, __bin_y]), no_fill_contours=False, color=get_color(band=band)[0], levels=LVLS, contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw})

			# Plot points and not density bins #
			corner.hist2d(__hax_color[i], __vax_color[i], plot_density=False, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, color=PT_COLORS[band], levels=LVLS, contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidths':__lw}, data_kwargs={'alpha':0.35, 'ms':1.75})

			if COLOR_YLOW is None and COLOR_YHIGH is None:
				# Force symmetric vertical axis #
				plt.ylim([-1*__sym_ylim, __sym_ylim])

			# Work-around for contour label #
			for j in np.arange(0, len(LVLS)):
				plt.plot([0.5, 0.5], [0.5, 0.5], color=CLRS_LABEL[j], label='$P_{'+str(round(LVLS[j], 2))[2:]+'}$', linewidth=__lw)
			plt.legend().draggable()


		if COLOR_YLOW is not None and COLOR_YHIGH is not None:
			plt.ylim([COLOR_YLOW, COLOR_YHIGH])

		'''
		plt.axis('scaled') #plt.gca().set_aspect('equal')
		# Plot x=y line to guide eye #
		lims = [np.min([plt.xlim(), plt.ylim()]), np.max([plt.xlim(), plt.ylim()])]
		plt.plot(lims, lims, color='k', linestyle=':', linewidth=0.7)
		'''

		plt.axhline(y=0, color='black', linestyle=':', linewidth=0.7)
		plt.subplots_adjust(hspace=0.6)

                plt.ylabel(vaxColorDiffLabel)
                plt.xlabel(__hax_label)
                plt.title(__subtitle_pref+' bins: ' + __write_bin) 

        plt.subplots_adjust(hspace=0.4)

	if percentRecoveredFlagsIncluded is not None:
		plot_title = '{} Percent Recovered: {}%'.format(plot_title, round(percentRecoveredFlagsIncluded, 2))
        plt.suptitle(plot_title, fontweight='bold')

	plt.text(0.1, 0.01, BASE_PATH_TO_CATS, fontsize=10, transform=plt.gcf().transFigure)

        if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', WRITE_COLORS[band]); print '-----> Saving plot as: ', fn_plot; plt.savefig(fn_plot)

        if SHOW_PLOT: plt.show()

	return 0 




def get_flux_plot_variables(tile, realization, band, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_percent_recovered_log): 
	"""Get variables needed for flux histogram.
	
	Parameters
        ----------
	df_1and2 (pandas DataFrame)
		Contains FluxDifference/SigmaFlux where FluxDifference is FluxMeasured-FluxTrue

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
		Contains FluxDifference/SigmaFlux with flags* removed. FluxDifference is given by FluxMeasured-FluxTrue. SigmaFlux is given by 'cm_flux_cov_{band}_{band}'^0.5

	idxGood (list of ints)

	cleanFlux2-cleanFlux1 (list of floats)
		Contains FluxDifference where FluxDifference is FluxMeasured-FluxTrue.
	"""

	fullMag1 = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=band)
        fullMag2 = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=band)

        idxGood = analysis.get_good_indices_using_primary_flags(df_1and2=df_1and2, full_mag1=fullMag1, full_mag2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, band=band)[0]

	### Get flux ###
	fullFlux1 = np.array(analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr1, band=band))
        fullFlux2 = np.array(analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr2, band=band))
        cleanFlux1 = np.array(fullFlux1)[idxGood]; cleanFlux2 = np.array(fullFlux2)[idxGood]



	fluxErr = get_total_flux_error(df_1and2=df_1and2, idx_good=idxGood, band=band)

	# For label #
        magAxLabel1 = plot_labels.get_magnitude_axlabel(inj_percent=INJ1_PERCENT, inj=INJ1, mag_hdr=mag_hdr1, meas_or_true_cat=AXLABEL1, match_cat=MATCH_CAT1, band=band)
        magAxLabel2 = plot_labels.get_magnitude_axlabel(inj_percent=INJ2_PERCENT, inj=INJ2, mag_hdr=mag_hdr2, meas_or_true_cat=AXLABEL2, match_cat=MATCH_CAT2, band=band)
	
	# Subtract: {mof/sof/coadd} - truth if possible #
	if 'truth' in MATCH_CAT1:
                __norm_flux_diff = (cleanFlux2-cleanFlux1)/fluxErr
		__flux_diff = cleanFlux2-cleanFlux1
		# Assumes `axlabel_a` minus `axlabel_b`. Note that 'mag' is replaced with 'flux' if `PLOT_FLUX`. #
		haxLabel = plot_labels.get_short_difference_axlabel(axlabel_a=magAxLabel2, axlabel_b=magAxLabel1, band=band)
        if 'truth' in MATCH_CAT2 or ('truth' not in MATCH_CAT1 and 'truth' not in MATCH_CAT2):
                __norm_flux_diff = (cleanFlux1-cleanFlux2)/fluxErr
		__flux_diff = cleanFlux1-cleanFlux2
		# Assumes `axlabel_a` minus `axlabel_b` #
		haxLabel = plot_labels.get_short_difference_axlabel(axlabel_a=magAxLabel1, axlabel_b=magAxLabel2, band=band)

	### Percent recovered if applicable ###
        if MATCH_CAT1 in ('gal_truth', 'star_truth'):
                # Note that len(full_mag1) = len(full_mag2) & len(clean_mag1) = len(clean_mag2) so it does not matter which is passed to get_percent_recovered() #
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanFlux1, inj_percent=INJ1_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)
        if MATCH_CAT2 in ('gal_truth', 'star_truth'):
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanFlux1, inj_percent=INJ2_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)

	if MATCH_CAT1 not in ('gal_truth', 'star_truth') and MATCH_CAT2 not in ('gal_truth', 'star_truth'):
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = None, None


	return __norm_flux_diff, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved, haxLabel





def get_total_flux_error(df_1and2, idx_good, band):
	"""Compute total flux error. Total means taken from both `MATCH_CAT1` and `MATCH_CAT2`."""

	### Flux error from measured catalog only ###
        if MATCH_CAT1 == 'coadd':
                #__flux_err1 = df1and2[NEW_COADD_FLUX_ERR_GRIZ_HDR+'_1'][idx_good]
		__flux_err1 = np.array(analysis.get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_COADD_FLUX_ERR_GRIZ_HDR+'_1'))[idx_good]
        if MATCH_CAT2 == 'coadd':
                #__flux_err2 = df1and2[NEW_COADD_FLUX_ERR_GRIZ_HDR+'_2'][idx_good]
		__flux_err2 = np.array(analysis.get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_COADD_FLUX_ERR_GRIZ_HDR+'_2'))[idx_good]

	# Do not use error from truth catalogs #
        if 'truth' in MATCH_CAT1:
		__flux_err1 = None
        if 'truth' in MATCH_CAT2:
		__flux_err2 = None

	if 'y3_gold' in MATCH_CAT1:
		__flux_err1 = np.sqrt(df_1and2[NEW_Y3_GOLD_CM_FLUX_COV_HDR+'_'+band.upper()+'_'+band.upper()+'_1'][idx_good])
	if 'y3_gold' in MATCH_CAT2:
		__flux_err2 = np.sqrt(df_1and2[NEW_Y3_GOLD_CM_FLUX_COV_HDR+'_'+band.upper()+'_'+band.upper()+'_2'][idx_good])

	if 'y3_gold' not in MATCH_CAT1 and 'truth' not in MATCH_CAT1 and MATCH_CAT1 != 'coadd':
		__flux_err1 = np.sqrt(analysis.get_matrix_diagonal_element(df=df_1and2, band=band, sq_matrices_hdr=CM_FLUX_COV_HDR1)[idx_good])
	if 'y3_gold' not in MATCH_CAT2 and 'truth' not in MATCH_CAT2 and MATCH_CAT2 != 'coadd':
		__flux_err2 = np.sqrt(analysis.get_matrix_diagonal_element(df=df_1and2, band=band, sq_matrices_hdr=CM_FLUX_COV_HDR2)[idx_good])

	if __flux_err1 is not None and __flux_err2 is not None: __flux_err = np.sqrt(np.power(__flux_err2, 2)+np.power(__flux_err1, 2))
	if __flux_err1 is None: __flux_err = np.sqrt(np.power(__flux_err2, 2))
	if __flux_err2 is None: __flux_err = np.sqrt(np.power(__flux_err1, 2))
	
	return __flux_err




def get_magnitude_plot_variables(band, df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization, tile, fn_flag_log, plot_title, fn_plot, fn_percent_recovered_log):
	"""Get variables needed for magnitude plots.

	Parameters
	----------
	
	Returns
	-------
	"""

	### Axes labels ###
        magAxLabel1 = plot_labels.get_magnitude_axlabel(inj_percent=INJ1_PERCENT, inj=INJ1, mag_hdr=mag_hdr1, meas_or_true_cat=AXLABEL1, match_cat=MATCH_CAT1, band=band)
        magAxLabel2 = plot_labels.get_magnitude_axlabel(inj_percent=INJ2_PERCENT, inj=INJ2, mag_hdr=mag_hdr2, meas_or_true_cat=AXLABEL2, match_cat=MATCH_CAT2, band=band)
	# Assumes `axlabel_a` minus `axlabel_b` #
        magVaxLabel = plot_labels.get_short_difference_axlabel(axlabel_a=magAxLabel1, axlabel_b=magAxLabel2, band=band)
        if SWAP_ORDER_OF_SUBTRACTION:
		# Assumes `axlabel_a` minus `axlabel_b` so set `axlabel_a=magAxLabel2` to reverse order of subtraction in label #
                magVaxLabel = plot_labels.get_short_difference_axlabel(axlabel_a=magAxLabel2, axlabel_b=magAxLabel1, band=band)


        # Magnitudes with no flags* removed # 
        fullMag1 = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=band)
        fullMag2 = analysis.get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=band)

        ### Remove objects with flags* or perform quality cuts ###
	idxGood = analysis.get_good_indices_using_primary_flags(df_1and2=df_1and2, full_mag1=fullMag1, full_mag2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, band=band)[0]


        # Magnitudes with flags* removed #
        cleanMag1 = analysis.get_good_data(df_1and2=df_1and2, hdr=mag_hdr1, idx_good=idxGood, str_of_arr=True, band=band)
        cleanMag2 = analysis.get_good_data(df_1and2=df_1and2, hdr=mag_hdr2, idx_good=idxGood, str_of_arr=True, band=band)


        ### Calculate errors. get_magnitude_error() will return array of zeros for truth catalogs. ###
        magErr1 = get_magnitude_error(mag_err_hdr=mag_err_hdr1, flux_hdr=CM_FLUX_HDR1, flux_cov_hdr=CM_FLUX_COV_HDR1, df_1and2=df_1and2, band=band, idx_good=idxGood, match_cat=MATCH_CAT1)
        magErr2 = get_magnitude_error(mag_err_hdr=mag_err_hdr2, flux_hdr=CM_FLUX_HDR2, flux_cov_hdr=CM_FLUX_COV_HDR2, df_1and2=df_1and2, band=band, idx_good=idxGood, match_cat=MATCH_CAT2)


	### Percent recovered if applicable ### 
        if MATCH_CAT1 in ('gal_truth', 'star_truth'):
                # Note that len(full_mag1) = len(full_mag2) & len(clean_mag1) = len(clean_mag2) so it does not matter which is passed to get_percent_recovered() #
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanMag1, inj_percent=INJ1_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)
        #TODO % recovered with star_truth (no flags)
        if MATCH_CAT2 in ('gal_truth', 'star_truth'):
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanMag1, inj_percent=INJ2_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)

	if MATCH_CAT1 not in ('gal_truth', 'star_truth') and MATCH_CAT2 not in ('gal_truth', 'star_truth'):
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = None, None


	return magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved 



def normalized_magnitude_difference_plotter(mag_hdr1, mag_hdr2, cbar_data, mag_err1, mag_err2, band, clean_mag1, full_mag1, mag_axlabel1, clean_mag2, mag_axlabel2, plot_title, realization, tile, cbar_label, fn_plot, fn_mag_err_log, vax_label, fn_mag_diff_outliers_log, fn_num_objs_in_1sig_mag_log):
	"""Produce a mangitude plot normalized to 1sigma_mag. If specified with `PLOT_68P` and `PLOT_34P_SPLIT` the 68th percentiles of each magnitude bin (see `normalize_magnitude_difference_to_error_and_maintain_bin_structure()` for bins) are plotted.

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

	# Args needed to call normalize_magnitude_difference_to_error() #
	normVaxBins, initialBins, haxBins, magErrBinMedians, vaxBinMedian = normalize_magnitude_difference_to_error_and_maintain_bin_structure(clean_mag1=clean_mag1, clean_mag2=clean_mag2, mag_err1=mag_err1, mag_err2=mag_err2, band=band, tile=tile, realization=realization, fn_mag_err_log=fn_mag_err_log, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log)


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
	plotMagnitudeDifferenceag, plotHaxMag, cleanBins = normalize_magnitude_difference_to_error(binned_norm_mag_diff=normVaxBins, mag_bins_for_mag_err=initialBins, binned_hax_mag=haxBins)

	percent_1sig = one_sigma_magnitude_counter_for_normalized_magnitude_difference(norm_mag_diff=plotMagnitudeDifferenceag, full_mag=full_mag1, clean_mag=clean_mag1, mag_bins_for_mag_err=initialBins, hax_mag=plotHaxMag, mag_err_bin_medians=magErrBinMedians, norm_vax_mag_bin_medians=vaxBinMedian, tile=tile, realization=realization, band=band, fn_num_objs_in_1sig_mag_log=fn_num_objs_in_1sig_mag_log) 


	# One colorbar at a time. This error is caught at beginning of script #
        if SCATTER:
                plt.scatter(plotHaxMag, plotMagnitudeDifferenceag, color=PT_COLORS[band], alpha=0.25, s=0.25)


        if CM_T_ERR_CBAR or CM_T_CBAR:
                plt.scatter(plotHaxMag, plotMagnitudeDifferenceag, c=cbar_data, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
                plt.colorbar(label=cbar_label)


        if HEXBIN:
		grid = (100, 1000)
		if VERBOSE_ING: print ' Normalized hexbin has a large number of grid cells. Plotting will take a moment... \n'
                plt.hexbin(plotHaxMag, plotMagnitudeDifferenceag, gridsize=grid, cmap=CMAPS[band], bins='log')
                plt.colorbar(label='log(counts)')



	if HIST_2D:
                # 1/10 the bin size of that used in error calculation #
                bin_x = np.arange(min(plotHaxMag), max(plotHaxMag), 0.5/10)
                if MAG_YLOW is not None and MAG_YHIGH is not None:
                        # Somewhat using reported 1% error in magnitude #
                        bin_y = np.arange(MAG_YLOW, MAG_YHIGH, (MAG_YHIGH-MAG_YLOW)*0.01)
                if MAG_YLOW is None and MAG_YHIGH is None:
                        bin_y = np.arange(min(plotMagnitudeDifferenceag), max(plotMagnitudeDifferenceag), (max(plotMagnitudeDifferenceag)-min(plotMagnitudeDifferenceag))*0.01)
                plt.hist2d(plotHaxMag, plotMagnitudeDifferenceag, bins=[bin_x, bin_y], cmap=CMAPS[band], norm=matplotlib.colors.LogNorm())
                plt.colorbar()


        if CORNER_HIST_2D:
		# Create bins of 1/4 magnitude for hax #
		__bin_x = np.linspace(np.min(plotHaxMag), np.max(plotHaxMag), math.ceil(4.0*(np.max(plotHaxMag) - np.min(plotHaxMag))))
		# Create bins of 1/20 magnitude for vax #
		__ylow = np.min(np.array(plotMagnitudeDifferenceag)) # -1*np.max(abs(np.array(__vax_color[i]))) ?
		__yhigh = np.max(np.array(plotMagnitudeDifferenceag)) # np.max(abs(np.array(__vax_color[i]))) ?
		__bin_y = np.linspace(__ylow, __yhigh, 20.0*(__yhigh-__ylow))

		# Force symmetric vertical axis #
		__sym_ylim = np.mean([abs(__ylow), __yhigh])

		__lw = 0.9
                corner.hist2d(plotHaxMag, plotMagnitudeDifferenceag, bins=np.array([__bin_x, __bin_y]), no_fill_contours=True, levels=LVLS, color=PT_COLORS[band], contour_kwargs={'colors':CLRS, 'cmap':None, 'linewidth':__lw})
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
        plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 2)) + '%')

	if SUBPLOT is False:

                fn_plot = fn_plot.replace('griz', band)

                ### Title for  ###
                plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

                ### Save plot ###
		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

                if SHOW_PLOT: plt.show()

        #plt.gca().set_aspect('equal')


	#FIXME 
	prcentRecoveredFlagsIncluded = None
	return percentRecoveredFlagsIncluded 




def magnitude_difference_plotter(mag_hdr1, mag_hdr2, cbar_data, mag_err1, mag_err2, band, clean_mag1, full_mag1, mag_axlabel1, clean_mag2, mag_axlabel2, plot_title, realization, tile, cbar_label, fn_plot, fn_mag_err_log, vax_label, fn_mag_diff_outliers_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log):
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

	band (str)

	realization (str) 

	tile (str)

        Returns
	-------
	"""

	### Values to plot ###
	__mag_diff = np.array(clean_mag1) - np.array(clean_mag2) 
	if SWAP_ORDER_OF_SUBTRACTION:
		__mag_diff = np.array(clean_mag2) - np.array(clean_mag1)

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
	percent_1sig = one_sigma_magnitude_counter(full_mag=full_mag1, mag_diff=__mag_diff, binned_mag_err=magErrBinMedians, clean_mag=clean_mag1, mag_bins_for_mag_err=initialBins, hax_mag=outlierCleanedHaxMag, vax_mag_bin_medians=vaxBinMedian, tile=tile, realization=realization, band=band, fn_num_objs_in_1sig_mag_log=fn_num_objs_in_1sig_mag_log) 



	if VERBOSE_ING: print 'Plotting {} objects for {}-band...\n'.format(len(clean_mag1), band)

	### Plot ###
	# One colorbar at a time. This error is caught at beginning of script #
	if SCATTER:
		plt.scatter(__hax_mag, __mag_diff, color=PT_COLORS[band], alpha=0.25, s=0.25)

	
	if CM_T_ERR_CBAR or CM_T_CBAR:
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
	plt.title('Objects in 1$\sigma_{mag}$: ' + str(round(percent_1sig, 2)) + '%')


	if SUBPLOT is False:

		fn_plot = fn_plot.replace('placehold', band)

		### Title for  ###
		plt.title(plot_title + '\n% objs in 1$\sigma$: ' + str(percent_1sig))

		### Save plot ###
		if SAVE_PLOT: print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

		if SHOW_PLOT: plt.show()

        #plt.gca().set_aspect('equal')


        return 0 




def magnitude_difference_subplotter(df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_plot, plot_title, realization, tile, fn_flag_log, fn_mag_err_log,  fn_mag_diff_outliers_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log):
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
		magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_magnitude_plot_variables(band=b, df_1and2=df_1and2, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=realization, tile=tile, fn_flag_log=fn_flag_log, plot_title=plot_title, fn_plot=fn_plot, fn_percent_recovered_log=fn_percent_recovered_log)


		### Get colorbar from measured catalog ###
		if CM_T_CBAR or CM_T_ERR_CBAR:
			if 'truth' not in MATCH_CAT2:
				cbarData, cbarLabel = get_colorbar_for_magnitude_plot_properties(df_1and2=df_1and2, cm_t_hdr=CM_T_HDR2, cm_t_err_hdr=CM_T_ERR_HDR2, idx_good=idxGood, clean_mag1=cleanMag1, clean_mag2=cleanMag2, meas_or_true_cat=AXLABEL2, inj=INJ2, inj_percent=INJ2_PERCENT) 
			if 'truth' not in MATCH_CAT1:
				cbarData, cbarLabel = get_colorbar_for_magnitude_plot_properties(df_1and2=df_1and2, cm_t_hdr=CM_T_HDR1, cm_t_err_hdr=CM_T_ERR_HDR1, idx_good=idxGood, clean_mag1=cleanMag1, clean_mag2=cleanMag2, meas_or_true_cat=AXLABEL1, inj=INJ1, inj_percent=INJ1_PERCENT)
		else: 
			cbarData, cbarLabel = None, None



		### Subplot ###
		if SUBPLOT: 
			plt.subplot(2, 2, counter_subplot)

		if NORMALIZE is False:
			magnitude_difference_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_data=cbarData, plot_title=plot_title, mag_err1=magErr1, mag_err2=magErr2, band=b, full_mag1=fullMag1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, mag_axlabel1=magAxLabel1, mag_axlabel2=magAxLabel2, realization=realization, tile=tile, cbar_label=cbarLabel, fn_plot=fn_plot, fn_mag_err_log=fn_mag_err_log, vax_label=magVaxLabel, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log, fn_percent_recovered_log=fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log=fn_num_objs_in_1sig_mag_log)

		if NORMALIZE:
			normalized_magnitude_difference_plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_data=cbarData, plot_title=plot_title, mag_err1=magErr1, mag_err2=magErr2, band=b, full_mag1=fullMag1, clean_mag1=cleanMag1, clean_mag2=cleanMag2, mag_axlabel1=magAxLabel1, mag_axlabel2=magAxLabel2, realization=realization, tile=tile, cbar_label=cbarLabel, fn_plot=fn_plot, fn_mag_err_log=fn_mag_err_log, vax_label=magVaxLabel, fn_mag_diff_outliers_log=fn_mag_diff_outliers_log, fn_num_objs_in_1sig_mag_log=fn_num_objs_in_1sig_mag_log)

		counter_subplot += 1

	# Plot after all panels have been filled #
	if SUBPLOT:

		### Show or save the plot once all four subplots have been filled ###
		plt.subplots_adjust(hspace=0.4)
		#plt.subplots_adjust(wspace=0.3)
		#plt.tight_layout(pad=3, h_pad=2.5)


		### Title ###
		# `percentRecoveredFlagsIncluded` is the same for all bands. Difference bands have different numbers of flags so this is not included in the suptitle #
		if percentRecoveredFlagsIncluded is not None:
			plot_title = '{} Recovered: {}%'.format(plot_title, round(percentRecoveredFlagsIncluded, 2))
		plt.suptitle(plot_title, fontweight='bold')

		# Include `BASE_PATH_TO_CATS` on the plot (using figure coordinates #
		plt.text(0.1, 0.01, BASE_PATH_TO_CATS, fontsize=10, transform=plt.gcf().transFigure)

		### Save plot ###
		if SAVE_PLOT: fn_plot = fn_plot.replace('placehold', 'griz'); print '-----> Saving plot as, ', fn_plot; plt.savefig(fn_plot)

		### Show plot ###
		if SHOW_PLOT: plt.show()

	
	return 0 











#TODO split this into flux and mag functions because the flux is not always needed
def get_coadd_catalog_observables(fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z, mag_hdr, mag_err_hdr, flux_hdr, flux_err_hdr):
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

	print 'Getting g-, r-, i-, z-band magnitudes and fluxes and corresponding errors for coadd catalogs...'

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
	mag_err_g = data_g[mag_err_hdr]; mag_err_r = data_r[mag_err_hdr]; mag_err_i = data_i[mag_err_hdr]; mag_err_z = data_z[mag_err_hdr]

	# Get fluxes #
	f_g = data_g[flux_hdr]; f_r = data_r[flux_hdr]; f_i = data_i[flux_hdr]; f_z = data_z[flux_hdr]

	# Get flux error #
	flux_err_g = data_g[flux_err_hdr]; flux_err_r = data_r[flux_err_hdr]; flux_err_i = data_i[flux_err_hdr]; flux_err_z = data_z[flux_err_hdr]

	__coadd_flux_griz, __coadd_flux_err_griz, __coadd_mag_griz, __coadd_mag_err_griz = [], [], [], []

        for i in np.arange(0, len(m_g)):
		__coadd_flux_griz.append("'("+ str(f_g[i]) + ', ' + str(f_r[i]) + ', ' + str(f_i[i]) + ', ' + str(f_z[i]) + ")'")
		__coadd_flux_err_griz.append("'("+ str(flux_err_g[i])+ ', ' + str(flux_err_r[i])+ ', ' + str(flux_err_i[i]) + ', ' + str(flux_err_z[i]) + ")'")
                __coadd_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")
		__coadd_mag_err_griz.append("'("+ str(mag_err_g[i])+ ', ' + str(mag_err_r[i])+ ', ' + str(mag_err_i[i]) + ', ' + str(mag_err_z[i]) + ")'")

        return __coadd_flux_griz, __coadd_flux_err_griz, __coadd_mag_griz, __coadd_mag_err_griz 




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

	if VERBOSE_ING: print 'Calculating g-, r-, i-, z-band magnitudes for star truth catalog...'

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

	if VERBOSE_ING: print 'Getting Y3 Gold g-, r-, i-, z-band magnitudes for', mag_hdr

	# Get headers, which are dependent on band #
	hdr_g = mag_hdr[:-2] + '_G' + mag_hdr[-2:]; hdr_r = mag_hdr[:-2] + '_R' + mag_hdr[-2:]
	hdr_i = mag_hdr[:-2] + '_I' + mag_hdr[-2:]; hdr_z = mag_hdr[:-2] + '_Z' + mag_hdr[-2:]
	
	# Read magnitudes from DataFrame #
	m_g = df_1and2[hdr_g]; m_r = df_1and2[hdr_r]; m_i = df_1and2[hdr_i]; m_z = df_1and2[hdr_z]

	__y3_gold_mag_griz = []

        for i in np.arange(0, len(m_g)):
                __y3_gold_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")

        return __y3_gold_mag_griz 




def get_y3_gold_catalog_flux_griz(df_1and2, flux_hdr, flux_cov_hdr):
        """Creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'. Solely for use with Y3 Gold catalogs.

        Parameters
        ----------
        df_1and2 (pandas DataFrame)

        Returns
        -------
        __y3_gold_flux_griz (list of str)
        """

        if VERBOSE_ING: print 'Getting Y3 Gold g-, r-, i-, z-band fluxes for', flux_hdr

        # Get headers, which are dependent on band #
	hdr_g = flux_hdr[:-2] + '_G' + flux_hdr[-2:]; hdr_r = flux_hdr[:-2] + '_R' + flux_hdr[-2:]
	hdr_i = flux_hdr[:-2] + '_I' + flux_hdr[-2:]; hdr_z = flux_hdr[:-2] + '_Z' + flux_hdr[-2:]

	hdr_flux_cov_g = '{}_G_G{}'.format(flux_cov_hdr[:-2], flux_cov_hdr[-2:])
	hdr_flux_cov_r = '{}_R_R{}'.format(flux_cov_hdr[:-2], flux_cov_hdr[-2:])
	hdr_flux_cov_i = '{}_I_I{}'.format(flux_cov_hdr[:-2], flux_cov_hdr[-2:])
	hdr_flux_cov_z = '{}_Z_Z{}'.format(flux_cov_hdr[:-2], flux_cov_hdr[-2:])

        # Read magnitudes from DataFrame #
        f_g = df_1and2[hdr_g]; f_r = df_1and2[hdr_r]; f_i = df_1and2[hdr_i]; f_z = df_1and2[hdr_z]


	flux_cov_g = df_1and2[hdr_flux_cov_g]
	flux_cov_r = df_1and2[hdr_flux_cov_r]
	flux_cov_i = df_1and2[hdr_flux_cov_i]
	flux_cov_z = df_1and2[hdr_flux_cov_z]

        __y3_gold_flux_griz = []
	#__y3_gold_flux_cov_griz = np.empty(len(f_g))
	__y3_gold_flux_cov_griz = []

        for i in np.arange(0, len(f_g)):
		#TODO .format()
                __y3_gold_flux_griz.append("'("+ str(f_g[i]) + ', ' + str(f_r[i]) + ', ' + str(f_i[i]) + ', ' + str(f_z[i]) + ")'")
		__y3_gold_flux_cov_griz.append('({}, {}, {}, {})'.format(flux_cov_g[i], flux_cov_r[i], flux_cov_i[i], flux_cov_r[i])) 
		
	
        return __y3_gold_flux_griz, __y3_gold_flux_cov_griz




def catch_catalog_filename_error(cat_type, inj, inj_percent):
	"""TODO"""

	__err = None

	return __err





# !!!!! #
def get_directory(tile, realization, low_level_dir):
        """Get directory"""

	if isinstance(low_level_dir, str):
		__dir = os.path.join(OUTPUT_DIRECTORY, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, low_level_dir)

	if isinstance(low_level_dir, list):
		if low_level_dir[1] == 'magnitude' and NORMALIZE:
			__dir = os.path.join(OUTPUT_DIRECTORY, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, low_level_dir[0], low_level_dir[1], 'normalized')
		else:
			__dir = os.path.join(OUTPUT_DIRECTORY, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, low_level_dir[0], low_level_dir[1])


	if RUN_TYPE is not None:
		__dir = os.path.join(OUTPUT_DIRECTORY, 'outputs', BALROG_RUN, MATCH_TYPE, tile, realization, low_level_dir, 'fof_analysis')

        if os.path.isdir(__dir) is False:
                if NO_DIR_MAKE is False:
                        sys.exit('Directory '+str(__dir)+' does not exist. Change directory structure in `get_directory()`.\n')
                if NO_DIR_MAKE:
                        if VERBOSE_ING: print 'Making directory {}...'.format(__dir)
                        os.makedirs(__dir)

        return __dir




def get_zipped_coadd_magnitudes(fn_base_g, fn_base_r, fn_base_i, fn_base_z):
	"""TODO"""

	if VERBOSE_ING: print 'Getting g-, r-, i-, z-band magnitudes and flags for zipped coadd catalogs...\n'

	mag_hdr = 'MAG_AUTO'
	mag_err_hdr = 'MAGERR_AUTO'
	flag1_hdr = 'FLAGS'
	flag2_hdr = 'IMAFLAGS_ISO'

	# Open FITS files #
	hdu_g = fits.open(fn_base_g); hdu_r = fits.open(fn_base_r); hdu_i = fits.open(fn_base_i); hdu_z = fits.open(fn_base_z)

	# Read data #
	data_g = hdu_g[1].data; data_r = hdu_r[1].data; data_i = hdu_i[1].data; data_z = hdu_z[1].data

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	err_g = data_g[mag_err_hdr]; err_r = data_r[mag_err_hdr]; err_i = data_i[mag_err_hdr]; err_z = data_z[mag_err_hdr]

	# Get flags #
	flag1_g = data_g[flag1_hdr]; flag1_r = data_r[flag1_hdr]; flag1_i = data_i[flag1_hdr]; flag1_z = data_z[flag1_hdr]
	flag2_g = data_g[flag2_hdr]; flag2_r = data_r[flag2_hdr]; flag2_i = data_i[flag2_hdr]; flag2_z = data_z[flag2_hdr]


	__base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz = [], [], [], []

	for i in np.arange(0, len(m_g)):
		__base_mag_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")
		__base_mag_err_griz.append("'("+ str(err_g[i])+ ', ' + str(err_r[i])+ ', ' + str(err_i[i]) + ', ' + str(err_z[i]) + ")'")
		__base_flag1_griz.append("'("+ str(flag1_g[i])+ ', ' + str(flag1_r[i])+ ', ' + str(flag1_i[i]) + ', ' + str(flag1_z[i]) + ")'")
		__base_flag2_griz.append("'("+ str(flag2_g[i])+ ', ' + str(flag2_r[i])+ ', ' + str(flag2_i[i]) + ', ' + str(flag2_z[i]) + ")'")

	return __base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz


def get_and_reformat_base_catalog(tile, realization):
	"""Handle new format for base catalogs. Zipped and SExtractor?

	Parameters
	----------

	Returns
	-------
	"""

	if VERBOSE_ING: print 'Unzipping and reformatting base coadd catalogs...\n'

	__overwrite = False 
	if __overwrite: raw_input('`__overwrite=True` in `get_and_reformat_base_catalog()`. Press enter to proceed, ctrl+c to stop')


	matchCatDir = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare')

	__fn_base_cat_for_matcher = os.path.join(matchCatDir, tile+'_base_i.fits')

	if os.path.isfile(__fn_base_cat_for_matcher) and __overwrite is False:
                if VERBOSE_ING: print ' Reformatted base catalog already exists. Not overwriting...\n'

	if os.path.isfile(__fn_base_cat_for_matcher) is False or __overwrite:
		if VERBOSE_ING: print ' Reformatting base catalog...'

		# Base catalogs are zipped #
		__fn_cat_bundle = os.path.join(BASE_PATH_TO_CATS, tile, 'base_cat.tgz')
		__tar = tarfile.open(__fn_cat_bundle, 'r')

		# Base catalogs are zipped #
		__fn_base_griz = []
		__tar.extractall(path=matchCatDir)
		for b in ALL_BANDS:
			__fn_b = '_'.join([tile, b, 'cat.fits'])
			__fn_base_griz.append(os.path.join(matchCatDir, __fn_b))
		__fn_base_g, __fn_base_r, __fn_base_i, __fn_base_z = __fn_base_griz


		__base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz = get_zipped_coadd_magnitudes(fn_base_g=__fn_base_g, fn_base_r=__fn_base_r, fn_base_i=__fn_base_i, fn_base_z=__fn_base_z) 

	
	       # Create new table #
		base_mag_col = Column(__base_mag_griz, name=NEW_BASE_MAG_GRIZ_HDR)
		base_mag_err_col = Column(__base_mag_err_griz, name=NEW_BASE_MAG_ERR_GRIZ_HDR)
		base_flag1_col = Column(__base_flag1_griz, name=NEW_BASE_FLAG_GRIZ_HDR)
		base_flag2_col = Column(__base_flag2_griz, name=NEW_BASE_IMAFLAG_ISO_GRIZ_HDR)
	
		# Add new table to i-band base catalog #
		table = Table.read(fn_base_i)
		table.add_column(base_mag_col, index=0)
		table.add_column(base_mag_err_col, index=1)
		table.add_column(base_flag1_col, index=2)
		table.add_column(base_flag2_col, index=3)

		# Save new table as FITS #
		table.write(__fn_base_cat_for_matcher, overwrite=__overwrite)

	print ' ----->', __fn_base_cat_for_matcher

	return __fn_base_cat_for_matcher 








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
	stack_dir = get_directory(tile='stack', realization=realization, low_level_dir='catalog_compare') 

	# Filename for stacked catalogs #
	fn_end_1and2 = '_'.join(['stacked', realization, MATCH_TYPE, 'match1and2.csv'])
	__fn_stack_tiles_1and2 = os.path.join(stack_dir, fn_end_1and2)

	fn_end_1not2 = '_'.join(['stacked', realization, MATCH_TYPE, 'match1not2.csv'])
	__fn_stack_tiles_1not2 = os.path.join(stack_dir, fn_end_1not2)

	fn_end_2not1 = '_'.join(['stacked', realization, MATCH_TYPE, 'match2not1.csv'])
	__fn_stack_tiles_2not1 = os.path.join(stack_dir, fn_end_2not1)

	# Check if stacked tile catalog already exists #
	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `stack_tiles()`. Press enter to procees and ctrl+c to stop.')

	if os.path.isfile(__fn_stack_tiles_2not1) and __overwrite is False:
		if VERBOSE_ING: print 'Stacked tile catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(__fn_stack_tiles_1and2)
		df1not2 = pd.read_csv(__fn_stack_tiles_1not2)
		df2not1 = pd.read_csv(__fn_stack_tiles_2not1)

	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(__fn_stack_tiles_2not1) is False or __overwrite:

		all_fn_1and2, all_fn_1not2, all_fn_2not1 = [], [], []

		for t in ALL_TILES:

			if RUN_TYPE is None:
				fn_1and2, fn_1not2, fn_2not1 = matcher(realization=realization, tile=t, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)[:3]

			if RUN_TYPE is not None:
				fn_1and2, fn_1not2, fn_2not1 = manipulate_catalogs.fof_matcher(realization=realization, tile=t) #TODO update with inj?

			all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		if VERBOSE_ING: print 'Stacking tiles. Stacking {} catalogs...'.format(len(all_fn_1and2))
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		if VERBOSE_ING: print 'Stacking complete ... \n'


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
	stack_dir = get_directory(tile=tile, realization='stack', low_level_dir='catalog_compare') 

	# Filename for stacked catalogs #
	__fn_stack_reals_1and2 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match1and2.csv')
	__fn_stack_reals_1not2 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match1not2.csv')
	__fn_stack_reals_2not1 = os.path.join(stack_dir, tile+'_stacked_'+str(MATCH_TYPE)+'_match2not1.csv')

	# Check if stacked realization catalog already exists #
	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `stack_realizations()`. Press enter to procees and ctrl+c to stop.')

	if os.path.isfile(__fn_stack_reals_2not1) and __overwrite is False:
		if VERBOSE_ING: print 'Stacked realization catalog exists. Not overwriting ... \n'
		df1and2 = pd.read_csv(__fn_stack_reals_1and2)
		df1not2 = pd.read_csv(__fn_stack_reals_1not2)
		df2not1 = pd.read_csv(__fn_stack_reals_2not1)


	# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
	if os.path.isfile(__fn_stack_reals_2not1) is False or __overwrite:

		all_fn_1and2, all_fn_1not2, all_fn_2not1 = [], [], []

		for r in ALL_REALIZATIONS:

			if RUN_TYPE is None:
				fn_1and2, fn_1not2, fn_2not1 = matcher(realization=r, tile=tile, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)[:3]

			if RUN_TYPE is not None:
				fn_1and2, fn_1not2, fn_2not1 = manipulate_catalogs.fof_matcher(realization=r, tile=tile) #FIXME

			all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

		if VERBOSE_ING: print 'Stacking realizations. ', len(all_fn_1and2), 'files ...'
		df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
		df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
		df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
		if VERBOSE_ING: print 'Stacking complete ... \n'


		# Save stacked catalog as DataFrame #
		df1and2.to_csv(__fn_stack_reals_1and2, sep=','); df1not2.to_csv(__fn_stack_reals_1not2, sep=','); df2not1.to_csv(__fn_stack_reals_2not1, sep=',')
		print '-----> Saving stacked realization catalogs as ', __fn_stack_reals_1and2 
	        print '----->', __fn_stack_reals_1not2 
                print '----->', __fn_stack_reals_2not1 

	__number_of_stacked_realizations = len(ALL_REALIZATIONS)

	return __fn_stack_reals_1and2, __fn_stack_reals_1not2, __fn_stack_reals_2not1, __number_of_stacked_realizations






def get_dataframe_and_headers(realization, tile, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, inj1, inj2, inj1_percent, inj2_percent):
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

	if VERBOSE_ING: print 'Getting (reading) DataFrame and getting headers for magnitude, magnitude error, flux, and flux error... \n'

	__mag_hdr1, __mag_hdr2, __mag_err_hdr1, __mag_err_hdr2, __flux_hdr1, __flux_hdr2 = mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, CM_FLUX_HDR1, CM_FLUX_HDR2 

	#FIXME if PLOT_COMPLTENESS and STACK_TILES

	#TODO make one general fn_1and2

	if PLOT_COMPLETENESS is False:
		### Stack tiles ###
		for r in ALL_REALIZATIONS:
			if STACK_TILES and STACK_REALIZATIONS is False: 
				fn_1and2, fn_1not2, fn_2not1, __number_of_stacked_tiles = stack_tiles(realization=r)
				__number_of_stacked_realizations = None


		### Stack realizations ###
		for t in ALL_TILES:
			if STACK_REALIZATIONS and STACK_TILES is False:
				fn_1and2, fn_1not2, fn_2not1, __number_of_stacked_realizations = stack_realizations(tile=t)
				__number_of_stacked_tiles = None


	
	if (STACK_REALIZATIONS is False and STACK_TILES is False) or (PLOT_COMPLETENESS and STACK_TILES):
		__number_of_stacked_tiles, __number_of_stacked_realizations = None, None
		# Filenames for catalogs #
		if RUN_TYPE is None:
			#fn_1and2, fn_1not2, fn_2not1 = matcher(realization=realization, tile=tile, inj1=inj1, inj2=inj2, inj1_percent=inj1_percent, inj2_percent=inj2_percent)[:3]
			fn_1and2, fn_1not2, fn_2not1 = manipulate_catalogs.match_catalogs(realization=realization, tile=tile, inj1=inj1, inj2=inj2, inj1_percent=inj1_percent, inj2_percent=inj2_percent, output_directory=OUTPUT_DIRECTORY, balrog_run=BALROG_RUN, base_path_to_catalogs=BASE_PATH_TO_CATS)[:3]
		if RUN_TYPE is not None:
			fn_1and2, fn_1not2, fn_2not1 = manipulate_catalogs.fof_matcher(realization=realization, tile=tile) #FIXME


	# Get DataFrame #
	__df1and2 = pd.read_csv(fn_1and2)
	__df1not2 = pd.read_csv(fn_1not2)
	__df2not1 = pd.read_csv(fn_2not1)


	### Handle star truth catalogs ###
	# Star truth catalogs matched then combined so a suffix (_1 or _2) is necessary #
	if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
		if VERBOSE_ING: print 'Adding new column to matched CSV. Column contains star truth catalog g-, r-, i-, z-band magnitudes...\n'

		if MATCH_CAT1 == 'star_truth':
			suf = '_1'
			__mag_hdr1 = NEW_STAR_TRUTH_MAG_GRIZ_HDR
		if MATCH_CAT2 == 'star_truth':
			suf = '_2'
			__mag_hdr2 = NEW_STAR_TRUTH_MAG_GRIZ_HDR

		star_mag = get_star_truth_catalog_magnitude(df_1and2=__df1and2, suf=suf)
		# New header must be of the form {base}_x where x is a single character because of the way m_axlabel is created from m_hdr # #FIXME check this
		__df1and2.insert(len(__df1and2.columns), NEW_STAR_TRUTH_MAG_GRIZ_HDR, star_mag)

	### Handle Y3 Gold catalogs ###
	# Y3 catalogs are matched then combined so a suffix (_1 or _2) is necessary #
	if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
		if VERBOSE_ING: print 'Adding two or four  new columns to matched CSV. Columns contain Y3 Gold 1) g-, r-, i-, z-band magnitudes, 2) g-, r-, i-, z-band magnitude errors; if PLOT_FLUX=True: 3) g-, r-, i-, z-band fluxes, and 4) g-, r-, i-, z-band flux errors...\n'

		# New header name #
		'''
		if 'star_truth' in (MATCH_CAT1, MATCH_CAT2):
			new_y3_gold_hdr = 'psf_mag_griz'
		if 'star_truth' not in (MATCH_CAT1, MATCH_CAT2):
			new_y3_gold_hdr = 'cm_mag_griz'
		'''
		if 'y3_gold' in MATCH_CAT1:
			# Headers in Y3 Gold #
			mhdr = M_HDR1
			fhdr = CM_FLUX_HDR1 
			fchdr = CM_FLUX_COV_HDR1
			# Headers to be added to DataFrame #
			__mag_hdr1 = NEW_Y3_GOLD_MAG_GRIZ_HDR 
			__flux_hdr1 = NEW_Y3_GOLD_FLUX_GRIZ_HDR
			__flux_cov_hdr1 = NEW_Y3_GOLD_CM_FLUX_COV_HDR 
		if 'y3_gold' in MATCH_CAT2:
			# Headers in Y3 Gold #
			mhdr = M_HDR2
			fhdr = CM_FLUX_HDR2 
			fchdr = CM_FLUX_COV_HDR2
			# Headers to be added to DataFrame #
			__mag_hdr2 = NEW_Y3_GOLD_MAG_GRIZ_HDR 
			__flux_hdr2 = NEW_Y3_GOLD_FLUX_GRIZ_HDR
			__flux_cov_hdr2 = NEW_Y3_GOLD_CM_FLUX_COV_HDR
		# Note: need magnitudes for `__idx_good` to get rid of mag=37.5 #
		y3_gold_mag = get_y3_gold_catalog_magnitude(df_1and2=__df1and2, mag_hdr=mhdr)
		# Add new column to df #
		__df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_MAG_GRIZ_HDR, y3_gold_mag)
		if PLOT_FLUX:
			y3_gold_flux_griz, y3_gold_flux_cov_griz = get_y3_gold_catalog_flux_griz(df_1and2=__df1and2, flux_hdr=fhdr, flux_cov_hdr=fchdr)
			__df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_FLUX_GRIZ_HDR, y3_gold_flux_griz)
			__df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_CM_FLUX_COV_HDR, y3_gold_flux_cov_griz)
		

	if MATCH_CAT1 == 'base':
		__mag_hdr1, __mag_err_hdr1 = NEW_BASE_MAG_GRIZ_HDR+'_1', NEW_BASE_MAG_ERR_GRIZ_HDR+'_1'
	if MATCH_CAT2 == 'base':
		__mag_hdr2, __mag_err_hdr2 = NEW_BASE_MAG_GRIZ_HDR+'_2', NEW_BASE_MAG_ERR_GRIZ_HDR+'_2'


	### Handle coadd catalogs. Catalog combined then matched so has suffix (unlike star) TODO re-explain #
	# Headers added to DataFrame #
	if MATCH_CAT1 == 'coadd':
		__mag_hdr1, __mag_err_hdr1 = NEW_COADD_MAG_GRIZ_HDR+'_1', NEW_COADD_MAG_ERR_GRIZ_HDR+'_1'
		__flux_hdr1, __flux_err_hdr1 = NEW_COADD_FLUX_GRIZ_HDR+'_1', NEW_COADD_FLUX_ERR_GRIZ_HDR+'_1' 
	if MATCH_CAT2 == 'coadd':
		__mag_hdr2, __mag_err_hdr2 = NEW_COADD_MAG_GRIZ_HDR+'_2', NEW_COADD_MAG_ERR_GRIZ_HDR+'_2'
		__flux_hdr2, __flux_err_hdr2 = NEW_COADD_FLUX_GRIZ_HDR+'_2', NEW_COADD_FLUX_ERR_GRIZ_HDR+'_2'

	if '__flux_err_hdr1' not in locals(): __flux_err_hdr1 = None
	if '__flux_err_hdr2' not in locals(): __flux_err_hdr2 = None


	### Log matches ###
	log_dir = get_directory(tile=tile, realization=realization, low_level_dir='log_files')

	fn_match_log = os.path.join(log_dir, tile+'_'+realization+'_'+MATCH_TYPE+'_matched_catalogs.log') 
	print ' -----> Saving match log as: ', fn_match_log, '\n'
	with open(fn_match_log, 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow(['TILE', 'REALIZATION', 'IN1', 'IN2', 'TOTAL_MATCHES_1AND2', 'TOTAL_MATCHES_1NOT2', 'TOTAL_MATCHES_2NOT1'])
                writer.writerow([tile, realization, MATCH_CAT1, MATCH_CAT2, len(__df1and2.index), len(__df1not2.index), len(__df2not1.index)])


	return __df1and2, __df1not2, __df2not1, __mag_hdr1, __mag_hdr2, __flux_hdr1, __flux_hdr2, __mag_err_hdr1, __mag_err_hdr2, __flux_err_hdr1, __flux_err_hdr2, __number_of_stacked_tiles, __number_of_stacked_realizations




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

	# Will be overwritten if necessary #
	magHdr1 = mag_hdr1
	magHdr2 = mag_hdr2

	if STACK_TILES is False: list_of_tiles_for_loop = ALL_TILES
        if STACK_REALIZATIONS is False: list_of_realizations_for_loop = ALL_REALIZATIONS

	if STACK_TILES and STACK_REALIZATIONS is False: list_of_tiles_for_loop = ['stack']
	if STACK_REALIZATIONS and STACK_TILES is False: list_of_realizations_for_loop = ['stack']

	for t in list_of_tiles_for_loop:
                for r in list_of_realizations_for_loop:

			# Stacked completeness plots do NOT use stacked tile catalog ... #
			if PLOT_COMPLETENESS is False: 
				df1and2, df1not2, df2not1, magHdr1, magHdr2, fluxHdr1, fluxHdr2, magErrHdr1, magErrHdr2, fluxErrHdr1, fluxErrHdr2, numStackTile, numStackReal = get_dataframe_and_headers(realization=r, tile=t, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT)

			if PLOT_COMPLETENESS: numStackTile, numStackReal = len(ALL_TILES), len(ALL_REALIZATIONS)

			fnFlagLog, fnMagErrLog, fnColorLog, fnMagDiffOutliersLog, fnMagCompletenessLog, fnGaussAperLog, fnFluxSigmaClipLog, fnPercentRecoveredLog, fnObjsIn1SigMagLog, fnMatchLog = outputs.get_log_filenames(balrog_run=BALROG_RUN, match_type=MATCH_TYPE, output_directory=OUTPUT_DIRECTORY, tile=t, realization=r)

			fnFlagLog, fnMagErrLog, fnColorLog, fnMagDiffOutliersLog, fnMagCompletenessLog, fnFluxSigmaClipLog, fnPercentRecoveredLog, fnObjsIn1SigMagLog = outputs.write_log_file_headers(fn_mag_err_log=fnMagErrLog, fn_flag_log=fnFlagLog, fn_color_log=fnColorLog, fn_mag_diff_outliers_log=fnMagDiffOutliersLog, fn_mag_completeness_log=fnMagCompletenessLog, fn_flux_sigma_clip_log=fnFluxSigmaClipLog, fn_percent_recovered_log=fnPercentRecoveredLog, fn_num_objs_in_1sig_mag_log=fnObjsIn1SigMagLog)

			# Name for plt.savefig() #
                        fnPlot = outputs.get_plot_filename(balrog_run=BALROG_RUN, match_type=MATCH_TYPE, output_directory=OUTPUT_DIRECTORY, realization=r, tile=t)

			# Plot title #
			plotTitle = plot_labels.get_plot_suptitle(realization=r, tile=t, number_of_stacked_realizations=numStackReal, number_of_stacked_tiles=numStackTile)

			### Region files ###
			#TODO PLOT_COMPLETENESS df1and2 created in another function --> cannot call write_to_region_files()
			if PLOT_COMPLETENESS is False:
				if MAKE_REG: write_to_region_files(df_1and2=df1and2, df_1not2=df1not2, df_2not1=df2not1, realization=r, tile=t)	


			### Call plotting functions ###	
			if PLOT_MAG and PLOT_COMPLETENESS is False:
				magnitude_difference_subplotter(df_1and2=df1and2, mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, fn_plot=fnPlot, plot_title=plotTitle, realization=r, tile=t, fn_mag_err_log=fnMagErrLog, fn_flag_log=fnFlagLog, fn_mag_diff_outliers_log=fnMagDiffOutliersLog, fn_percent_recovered_log=fnPercentRecoveredLog, fn_num_objs_in_1sig_mag_log=fnObjsIn1SigMagLog)

			if PLOT_COLOR and PLOT_COMPLETENESS is False:
				for b in ALL_BANDS[:-1]:
					color_subplotter(band=b, df_1and2=df1and2, mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, fn_flag_log=fnFlagLog, plot_title=plotTitle, fn_plot=fnPlot, fn_color_log=fnColorLog, fn_percent_recovered_log=fnPercentRecoveredLog)

			### Completeness plots (for truth catalogs vs injected {mof/sof/coadd/...}) ###
			if ('truth' in MATCH_CAT1 or 'truth' in MATCH_CAT2) and (INJ1 and INJ2) and PLOT_COMPLETENESS:
				if PLOT_MAG and STACK_TILES is False:
					magnitude_completeness_subplotter(mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, plot_title=plotTitle, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog, fn_percent_recovered_log=fnPercentRecoveredLog)
				if PLOT_MAG and STACK_TILES:
					stacked_magnitude_completeness_subplotter(mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, plot_title=plotTitle, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog, fn_percent_recovered_log=fnPercentRecoveredLog)

				if PLOT_COLOR:
					for b in ALL_BANDS[:-1]: 
						color_completeness_subplotter(band=b, df_1and2=df1and2, mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization=r, tile=t, plot_title=plotTitle, fn_plot=fnPlot, fn_flag_log=fnFlagLog, fn_mag_completeness_log=fnMagCompletenessLog)


			if PLOT_FLUX:
				normalized_flux_difference_histogram_subplotter(df_1and2=df1and2, flux_hdr1=fluxHdr1, flux_hdr2=fluxHdr2, mag_hdr1=magHdr1, mag_hdr2=magHdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_title=plotTitle, fn_plot=fnPlot, tile=t, realization=r, fn_gauss_aper_log=fnGaussAperLog, fn_flux_sigma_clip_log=fnFluxSigmaClipLog, gap_flux_hdr1=GAUSS_APER_FLUX_GRIZ_HDR+'_1', gap_flux_hdr2=GAUSS_APER_FLUX_GRIZ_HDR+'_2', fn_percent_recovered_log=fnPercentRecoveredLog)


	return 0



#FIXME get rid of flux_hdr dependence if this is a constant
def get_coadd_catalog_for_matcher(cat_type, inj_percent, inj, realization, mag_hdr, mag_err_hdr, tile, flux_hdr, flux_err_hdr):
	"""Make FITS file that includes a column of form '(m_g, m_r, m_i, m_z)' where m is magnitude. Column will be added to '..._i_cat.fits'. This will be used in matcher(). Relies on directory structure /`OUTPUT_DIRECTORY`/outputs/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/catalog_compare/

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

	tile (str)

	Returns
	-------
	__fn_coadd_for_matcher (str)
		Complete filename for catalog with added column. Is a FITS file.
	"""

	__overwrite = False 
	if __overwrite: raw_input('`__overwrite=True` in `get_coadd_catalog_for_matcher()`. Press enter to proceed, ctrl+c to exit.') 

	dir_new = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare') 

	__fn_coadd_for_matcher = os.path.join(dir_new, str(tile) + '_i_cat_combo.fits')

	# Check if new coadd catalog has already been created #
	if os.path.isfile(__fn_coadd_for_matcher):
		print 'New coadd catalog already exists ...\n'


	if os.path.isfile(__fn_coadd_for_matcher) is False or __overwrite:	
		print 'Adding a column to i-band coadd catalog. Will take a moment ...\n'

		# Get list of filenames #
		fn_griz = []
		for b in ALL_BANDS:
			fn_griz.append(manipulate_catalogs.get_catalog_filename(cat_type=cat_type, inj=inj, inj_percent=inj_percent, realization=realization, tile=tile, band=b, base_path_to_catalogs=BASE_PATH_TO_CATS, balrog_run=BALROG_RUN))
		fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z = fn_griz

		# Get coadd magnitude (mag_c) and magnitude error to be of form '(m_g, m_r, m_i, m_z)'. Recall that this is a string #
		coadd_flux, coadd_flux_err, mag_c, mag_err_c = get_coadd_catalog_observables(fn_coadd_cat_g=fn_coadd_cat_g, fn_coadd_cat_r=fn_coadd_cat_r, fn_coadd_cat_i=fn_coadd_cat_i, fn_coadd_cat_z=fn_coadd_cat_z, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_err_hdr=flux_err_hdr)
 
	       # Create new table #
		coadd_mag_col = Column(mag_c, name=NEW_COADD_MAG_GRIZ_HDR)
		coadd_mag_err_col = Column(mag_err_c, name=NEW_COADD_MAG_ERR_GRIZ_HDR)
		coadd_flux_col = Column(coadd_flux, name=NEW_COADD_FLUX_GRIZ_HDR)
		coadd_flux_err_col = Column(coadd_flux_err, name=NEW_COADD_FLUX_ERR_GRIZ_HDR)

		# Add new table to i-band coadd catalog #
		table = Table.read(fn_coadd_cat_i)
		table.add_column(coadd_mag_col, index=0)
       		table.add_column(coadd_mag_err_col, index=1)
		table.add_column(coadd_flux_col, index=2)
		table.add_column(coadd_flux_err_col, index=3)

		# Save new table as FITS #
		table.write(__fn_coadd_for_matcher, overwrite=__overwrite)

	return __fn_coadd_for_matcher 




#TODO accept idx_good as input param? Will only be used for df_match. 
def write_to_region_files(df_1and2, df_1not2, df_2not1, realization, tile):
	"""Make DS9 region files for catalogs matched via join=1and2, join=1not2, and join=2not1 (join type is STILTS parameter set in stilts_matcher or fof_stilts_matcher).

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

	if VERBOSE_ING: print 'Writing region files...'

	### Get filenames and open files ###
	fnRegion1and2, fnRegion1not2, fnRegion2not1 = outputs.get_region_filenames(balrog_run=BALROG_RUN, match_type=MATCH_TYPE, output_directory=OUTPUT_DIRECTORY, realization=realization, tile=tile)

	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `write_to_region_files()`. Press enter to procees and ctrl+c to stop.')
	
	if os.path.isfile(fnRegion2not1) and __overwrite is False:
		print 'Region files already exist. Not overwriting ...'
	
	#TODO make write_reg(a, b, color{}, ra, dec, angle) and call this three times. put in region_files.py

 
	if os.path.isfile(fnRegion2not1) is False or __overwrite:
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
				#TODO label units needed here
				fd_match.write('ellipse {} {} {}" {}" {} #color=green width=3\n'.format(ra_match[i], dec_match[i], a_match[i], b_match[i], 90+orientation_match[i]))

			for i in np.arange(0, len(ra_1not2)):
				fd_1not2.write('ellipse {} {} {}" {}" {} #color=green width=3\n'.format(ra_1not2[i], dec_1not2[i], a_1not2[i], b_1not2[i], 90+orientation_1not2[i]))

			for i in np.arange(0, len(ra_2not1)):
				fd_2not1.write('ellipse {} {} {}" {}" {} #color=green width=3\n'.format(ra_2not1[i], dec_2not1[i], a_2not1[i], b_2not1[i], 90+orientation_2not1[i]))


		# Non-coadd catalogs allow for circular regions #
		# Use cm_T #
		if MATCH_CAT1 != 'coadd' and MATCH_CAT2 != 'coadd':
			size_sq_match = df_1and2[CM_T_HDR1]
			size_sq_1not2 = df_1not2[CM_T_HDR1]
			size_sq_2not1 = df_2not1[CM_T_HDR2]

			# Note: Can use a typical radius of 2 arcsec #
			for i in np.arange(0, len(ra_match)):
				if size_sq_match[i] > 0:# and np.isnan(size_sq_match[i]) is False:
					fd_match.write('cirlce {} {} {}" #color=green width=3\n'.format(ra_match[i], dec_match[i], size_sq_match[i]**0.5))
			for i in np.arange(0, len(ra_1not2)):
				if size_sq_1not2[i] > 0: # and np.isnan(size_sq_1not2[i]) is False:
					fd_1not2.write('cirlce {} {} {}" #color=green width=3\n'.format(ra_1not2[i], dec_1not2[i], size_sq_1not2[i]**0.5))
			for i in np.arange(0, len(ra_2not1)):
				if size_sq_2not1[i] > 0: # and np.isnan(size_sq_2not1[i]) is False:
					fd_1not2.write('cirlce {} {} {}" #color=green width=3\n'.format(ra_2not1[i], dec_2not1[i], size_sq_2not1[i]**0.5))

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
		CM_T_CBAR, CM_T_ERR_CBAR, NORMALIZE, HEXBIN = cbar_bools
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



RUN_TYPE_LOOP = False 
if RUN_TYPE_LOOP:
	run_type_list = [None, 'ok', 'rerun']
	for run_type in run_type_list:
		RUN_TYPE = run_type
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)

