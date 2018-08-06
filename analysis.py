"""
Analysis for matched catalog.
E.g. Handle STILTS conversion of arrays and matrices to strings,
Remove flagged objects from matched (via `join=1and2`) catalog,
Calculate percent recovered (if applicable),
Get colorbar data (if applicable),
...

Comments are ABOVE the code they refer to.
"""

import csv
import fitsio
import numpy as np
import sys

# From BalVal #
from set_constants import *
from catalog_headers import *
import manipulate_catalogs
import plot_labels




def get_floats_from_string(df, band, four_elmt_arrs_hdr):
	"""Transform a list of strings of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]' 
	to a list of floats of form '[1,1,...]' (if band="g"), '[2,2,...]' ("r"), '[3,3,...]' ("i"), or '[4,4,...]' ("z"). 
	This is necessary for CSVs created from `stilts_matcher` or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings. 

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

	strings = df[four_elmt_arrs_hdr]

	__elmts = np.empty([len(df.index)])

	# Each element (elmt) is of form '(1, 2, 3, 4)' #
	for k in np.arange(0, len(strings)):

		elmt = strings[k]

		# Get first number in '(1, 2, 3, 4)' #
		if band == 'g':
			i = 1
			idx1 = elmt.find('(') + i
			idx2 = elmt.find(',')

		# Get second number in '(1, 2, 3, 4)' #
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

		__elmts[k] = float(elmt[idx1:idx2])

	if VERBOSE_MINOR:
		print 'Got {} for {}-band...'.format(four_elmt_arrs_hdr, band)
		print ' Check: ', strings[0], ' & ', __elmts[0], '\n'


	return __elmts




def get_matrix_diagonal_element(df, band, sq_matrices_hdr):
	"""Transforms a list of 4x4 matrices where each element is a string of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' into a list of either the 11 (if band is "g"), 22 ("r"), 33 ("i"), or 44 ("z") matrix elements. 
	This is necessary for CSVs created from `stilts_matcher` or ms_fof_matcher because arrays in FITS files of form (m_g, m_r, m_i, m_z) are converted to strings.

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

	matrices = df[sq_matrices_hdr]

	__matrix_diagonals = np.empty([len(df.index)])

	# Each element in `matrices` is a matrix of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' #
	for k in np.arange(0, len(matrices)):

		matrix = matrices[k]

		if band == 'g':
 			i, j = 2, 0
			idx1 = 0 + i
			idx2 = matrix.find(',')

		if band == 'r':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 4).find(',') + i
			idx2 = matrix.replace(',', ';', 5).find(',')

		if band == 'i':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 9).find(',') + i
			idx2 = matrix.replace(',', ';', 10).find(',')

		if band == 'z':
			i, j = 2, -1
			idx1 = matrix.replace(',', ';', 14).find(',') + i
			idx2 = matrix.replace(',', ';', 15).find(',')

		__matrix_diagonals[k] = float(matrix[idx1:idx2+j])

	if VERBOSE_MINOR:
		print 'Got {} for {}-band'.format(sq_matrices_hdr, band)
		print ' Check: {} & {}'.format(matrices[0], __matrix_diagonals[0])


	return __matrix_diagonals




def get_coadd_catalog_flags(df_1and2, suf, band):
	"""Get flags for coadd catalogs

	Parameters
	----------
	suf (str)
		Allowed values: '_1', '_2'
	"""

	__flags = []
	flag1 = get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_BASE_FLAG_GRIZ_HDR+suf)
	flag2 = get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_BASE_IMAFLAG_ISO_GRIZ_HDR+suf)
	__flags.append(flag1)
	__flags.append(flag2)

	return __flags




def get_y3_gold_catalog_flags(df_1and2, suf, band):
	"""Get flags for Y3 Gold catalogs.

        Parameters
        ----------
        suf (str)
                Allowed values: '_1', '_2'
        """

	__flags = []

	try:
		flag1 = df_1and2['{}_CM_MOF_FLAGS{}'.format(Y3_FIT, suf)]
		__flags.append(flag1)
	except:
		pass

	flag2 = df_1and2['SEXTRACTOR_FLAGS_{}{}'.format(band.upper(), suf)]
	flag3 = df_1and2['IMAFLAGS_ISO_{}{}'.format(band.upper(), suf)]

	__flags.append(flag2)
	__flags.append(flag2)

	return __flags




def get_additional_flags(cat_type, suf, band):
	"""

	Parameters
	----------

	Returns
	-------
	__flags (list)
	"""

	if cat_type == 'coadd':
		flags = get_coadd_catalog_flags(suf=suf, band=band)

	if 'y3_gold' in cat_type:
		flags = get_y3_gold_flags(suf=suf, band=band)

	return flags




#TODO pass mag hdrs?
def get_good_indices_using_primary_flags(df_1and2, full_mag1, full_mag2, cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2, band=None):
	"""Get indices of objects without flags* where flags* used are indicated in README.md. Store flagged indices if PLOT_FLAGGED_OBJS is True.
	Magnitude headers are necessary (not simply explicit flag headers) to search for mag=37.5 which is a negative flux.

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

	if VERBOSE_ING: print 'Removing flagged* objects from {}-band data...'.format(band)

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

	__num_additional_flags = 0
	# Additional flags for certain catalogs #
        if PLOT_GAUSS_APER_FLUX:
		# Two more flags from Gaussian aperture measurements on each catalog #
		__num_additional_flags += 2

	if MATCH_CAT1 == 'coadd' or 'y3_gold' in MATCH_CAT1:
		flags_list1 = get_additional_flags(cat_type=MATCH_CAT1, suf='_1')
		__num_additional_flags += len(flags_list1)

	if MATCH_CAT2 == 'coadd' or 'y3_gold' in MATCH_CAT2:
                flags_list2 = get_additional_flags(cat_type=MATCH_CAT2, suf='_2')
		__num_additional_flags += len(flags_list2)

	__ad_flags = np.empty([__num_additional_flags, len(df_1and2.index)]) #TODO include this in other functions here 


	if PLOT_GAUSS_APER_FLUX:
                __ad_flags[-2] = np.array(df_1and2['{}_1'.format(GAUSS_APER_FLAGS_HDR)].values)
                __ad_flags[-1] = df_1and2['{}_2'.format(GAUSS_APER_FLAGS_HDR)]

	__counter1 = 0
	try:
		flags_list1
		for i in np.array(0, len(flags_list1)):
			__ad_flags[i] = flags_list1[i]
		__counter1 = len(flags_list1)
	except NameError:
		pass
			 
        try:
                flags_list2
                for j in np.array(0, len(flags_list2)):
                        __ad_flags[__counter1+j] = flags_list2[i]
        except NameError:
                pass


	# Make arrays to take absolute value in next step #
	#FIXME can eventually remove np.array
	full_mag1, full_mag2 = np.array(full_mag1), np.array(full_mag2)


	#  Get rid of these objects; 37.5 corresponds to a negative flux #
	__idx_good = np.where(
			(abs(full_mag1) != 9999.0) &
			(abs(full_mag1) != 99.0) &
			(abs(full_mag1) != 37.5) &
			(abs(full_mag2) != 9999.0) &
			(abs(full_mag2) != 99.0) &
			(abs(full_mag2) != 37.5) &
			(flag1 == 0) &
			(flag2 == 0) &
			(cm_flag1 == 0) &
			(cm_flag2 == 0))[0]

	# `__ad_flags` is not empty if... #
	if 'coadd' in (MATCH_CAT1, MATCH_CAT2) or 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2 or PLOT_GAUSS_APER_FLUX:
		# Each row in `__ad_flags` corresponds to a separate flag #
		# Each column refers to an object #
		# Want ALL the cells in a column to be zero #
		# `np.where()` along an axis? Work-around #
		# Cannot predict number of indices so use a list #
		__ad_idx_list = []

		for i in np.arange(0, __ad_flags.shape[0]):
			__ad_idx_list.append(np.where(__ad_flags[i] == 0)[0])

		# Indices are only good if they are good for every flag type in `_ad_flags`, so find repeats #
		# Flatten list #
		__ad_idxs = [item for sublist in __ad_idx_list for item in sublist]
		# TODO try `unique_indices=True` 
		__bad = np.unique(np.array(__ad_idxs), return_inverse=True)[1] 

		__ad_idx_good = set(__ad_idxs) - set(__bad)

		__idx_good_tot = list(set(__idx_good) - set(__ad_idx_good))
	else:
		__idx_good_tot = __idx_good


	if PLOT_FLAGGED_OBJS:
		__all_idx = np.arange(0, len(df_1and2.index))
		__idx_bad = list(set(__all_idx) - set(__idx_good_tot))
	if PLOT_FLAGGED_OBJS is False:
		__idx_bad = None

		
        if VERBOSE_ED: print 'For {}-band, eliminated {} objects with flags..'.format(band, len(full_mag1) - len(__idx_good_tot))


	#TODO try just __idx_good
	return __idx_good_tot, __idx_bad	




def get_good_data(df_1and2, hdr, idx_good, str_of_arr, band):
	"""Get the data corresponding to indices of objects without flags*.

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (via join=1and2) catalog.

	hdr (str)
		Matched (via join=1and2) catalog header for the desired column of the DataFrame. 
	idx_good (list of ints)
		Indices of objects without flags* or indices of objects that satistfy quality cuts (if `

	str_of_arr (bool)
		Set to `True` if the desired column of the DataFrame is a string of form '(data_g, data_r, data_i, data_z)'.
	band (str) 
		Can be `None`. Used if `magnitude=True`.

	Returns
	-------
	clean_data (list of floats)
		Contains column of the DataFrame with objects with flags* removed (or objects that do not pass quality cuts removed if `
	"""

	if str_of_arr:
		fullData = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=hdr, band=band)
	if str_of_arr is False:
		fullData = df_1and2[hdr]

	__clean_data = np.array(fullData)[idx_good]

	return __clean_data




def get_percent_recovered(full_data, clean_data, inj_percent, tile, band, realization, fn_percent_recovered_log, balrog_run, base_path_to_catalogs):
	"""Calculate the percent of Balrog-injected objects recovered in the measured (real image + Balrog injections) catalog.
	Write this to log file.

	Parameters
	----------

	Returns
	-------
	__percent_recovered_flags_in (float)
		Percent recovered as calculated without removing any flagged objects.
	__percent_recovered_flags_rm (float)
		Percent recovered as calculated after removing flags from both catalogs.
	"""

	if VERBOSE_ING: print 'Calculating percent recovered for {}-band...\n'.format(band)

	__objs_per_tile = 50000.0
	__num_inj_objs_per_tile = (1.0*inj_percent/100.0)*__objs_per_tile

	constant = 1.0
	if STACK_TILES: constant = 1.0 * len(ALL_TILES)
	if STACK_REALIZATIONS: constant = 1.0 * len(ALL_REALIZATIONS)

	__inj_total = __num_inj_objs_per_tile*constant

	### Remove flags* from truth catalog (if possible) ###
	if 'gal_truth' in (MATCH_CAT1, MATCH_CAT2):

		if MATCH_CAT1 == 'gal_truth':
			fn_truth_catalog = manipulate_catalogs.get_catalog_filename(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=INJ1_PERCENT, realization=realization, tile=tile, band=band, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
			__mag_hdr, __flag_hdr, __cm_flag_hdr = M_HDR1[:-2], FLAGS_HDR1[:-2], CM_FLAGS_HDR1[:-2]
		if MATCH_CAT2 == 'gal_truth':
			fn_truth_catalog = manipulate_catalogs.get_catalog_filename(cat_type=MATCH_CAT2, inj=INJ2, inj_percent=INJ2_PERCENT, realization=realization, tile=tile, band=band, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
			__mag_hdr, __flag_hdr, __cm_flag_hdr = M_HDR2[:-2], FLAGS_HDR2[:-2], CM_FLAGS_HDR2[:-2]

		# Read catalog #
		data = fitsio.read(fn_truth_catalog, hud=1)
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



	# Star truth catalogs do not contain flags #
	if 'star_truth' in (MATCH_CAT1, MATCH_CAT2):
		__number_of_flags_in_truth_cat = 0 


	# Including flags #
	__percent_recovered_flags_in = 100.0*len(full_data)/__inj_total


	# Not including flags. This means removing flags from matched catalog AND truth catalog (if possible) #
	###__percent_recovered_flags_rm = 100.0*len(clean_data)/__inj_total
	if 'gal_truth' in (MATCH_CAT1, MATCH_CAT2):
		__percent_recovered_flags_rm = 100.0*len(clean_data)/(__inj_total - __number_of_flags_in_truth_cat)

	### Log ###
	with open(fn_percent_recovered_log, 'a') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		# 'TILE', 'REALIZATION', 'BAND', 'TOTAL_OBJS_IN_MATCH1AND2', 'TOTAL_FLAGGED_OBJS_IN_MATCH1AND2', 'TOTAL_FLAGGED_OBJS_IN_TRUTH', '%_RECOVERED_FLAGS_IN', '%_RECOVERED_FLAGS_RM', 'RUN_TYPE' #
		writer.writerow([tile, realization, band, len(full_data), len(full_data)-len(clean_data), __number_of_flags_in_truth_cat, __percent_recovered_flags_in, __percent_recovered_flags_rm, str(RUN_TYPE)]) 
	csvfile.close()


	return __percent_recovered_flags_in, __percent_recovered_flags_rm 




def get_colorbar_data(df_1and2, cm_t_hdr, cm_t_err_hdr, idx_good):
	"""Get colorbar for magnitude plot.

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.
	
	cm_t_hdr, cm_t_err_hdr (str)
		Matched (join=1and2) catalog headers for cm_T (size squared) and cm_T_err. Can be `None`.
	idx_good (list of ints)
		Indices of objects without flags.		

	Returns
	-------
	cbarData (list of floats)
		Data used to produce colorbar. Can be `None` if not colorbar is to be plotted.
	"""

        if CM_T_ERR_CBAR:
                cbarData = get_good_data(df_1and2=df_1and2, hdr=cm_t_err_hdr, idx_good=idx_good, str_of_arr=False, band=None)

        if CM_T_CBAR:
                cbarData = get_good_data(df_1and2=df_1and2, hdr=cm_t_hdr, idx_good=idx_good, str_of_arr=False, band=None)


	return cbarData, __cbar_label



def calculate_measured_magnitude_error_from_flux(flux_cov_hdr, df_1and2, band, flux_hdr, idx_good, match_cat):
	"""Calculate the magnitude error via 1.08 * (flux_cov[i][i])^0.5 / flux[i].
	Ignore flagged objects in error calculation.
	This error is a FRACTIONAL error.

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
	flux = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr, band=band)[idx_good]
	flux_cov = get_matrix_diagonal_element(df=df_1and2, sq_matrices_hdr=flux_cov_hdr, band=band)[idx_good]

	counter_neg = 0
	__mag_err_from_flux = np.empty([len(idx_good)])

	
	# Calculations #
	for i in np.arange(0, len(flux)):

		# Throw out negative and zero flux covariance elements (error calculation involves taking the square root of a negative) #
		if flux_cov[i] <= 0:
			__mag_err_from_flux[i] = 0
			counter_neg += 1

		if flux_cov[i] > 0:
			# Pogsons number = 1.08 #
			__mag_err_from_flux[i] = (1.08*np.sqrt(flux_cov[i])/flux[i])

	if VERBOSE_ED: print ' Number of negative or zero cm_flux_cov: {}/{}'.format(counter_neg, len(flux))

	if 'truth' in match_cat:
		__mag_err_from_flux = None

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

	__norm_mag_diff_68p, __neg_norm_mag_diff_34p, __pos_norm_mag_diff_34p = np.empty(len(binned_norm_mag_diff)), np.empty(len(binned_norm_mag_diff)), np.empty(len(binned_norm_mag_diff))

	# Loop through bins (b) #
	for b in np.arange(0, len(binned_norm_mag_diff)):

		if binned_norm_mag_diff[b] is None:
			__norm_mag_diff_68p[b] = None
			__neg_norm_mag_diff_34p[b] = None
			__pos_norm_mag_diff_34p[b] = None

		if binned_norm_mag_diff[b] is not None:

			### Find 68th percentile about zero ### 	
			# Values in current bin (icb) #
			vax_mag_bins_icb = binned_norm_mag_diff[b]
			# Take absolute value of each point in bin #
			abs_vax_mag_bins_icb = [abs(elmt) for elmt in vax_mag_bins_icb]
			# Percentile sorts the data #
			__norm_mag_diff_68p[b] = np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower')

			if VERBOSE_ED:
				# Check the percentile because interpolation='lower' was used #
				num = 0
				for j in np.arange(0, len(binned_norm_mag_diff[b])):
					if abs(binned_norm_mag_diff[b][j]) <= np.percentile(abs_vax_mag_bins_icb, 68, interpolation='lower'):
						num += 1
				if VERBOSE_ED: print 'Percent of objects within 68 percentile via np.percentile(interpolation=LOWER): {}...'.format(float(num)/len(binned_norm_mag_diff[b]))


			### Find 34th percentile of positive and negative values separately ###
			# Cannot predict number of negative and positive values, so use list #
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
				__neg_norm_mag_diff_34p[b] = np.percentile(neg_vax, 34, interpolation='lower')
			if counter_pos > 0:
				__pos_norm_mag_diff_34p[b] = np.percentile(pos_vax, 34, interpolation='lower')
			if counter_neg  == 0:
				__neg_norm_mag_diff_34p[b] = None
			if counter_pos  == 0:
				__pos_norm_mag_diff_34p[b] = None


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
	# TODO __vax_mag --> __mag_diff
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
				csvfile.close()

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

	if VERBOSE_ING: print 'Normalizing magnitude difference...'

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
	"""Normalize plot to 1sigma_mag curve using tame magnitude errors.

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
			if match_cat == 'base':
				try:
					__mag_error = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=NEW_BASE_MAG_ERR_GRIZ_HDR+'_1', band=band)
				except:
					__mag_error = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=NEW_BASE_MAG_ERR_GRIZ_HDR+'_2', band=band)
			# Pass good indices #
			__mag_error = __mag_error[idx_good]

        if PLOT_MAG_ERR is False or 'truth' in match_cat:
		__mag_error = None

	return __mag_error 




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
	__clean_mag1_band2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=BANDS_FOR_COLOR[band])[idx_good]	
	# From `MATCH_CAT2` #
	__clean_mag2_band2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=BANDS_FOR_COLOR[band])[idx_good] 

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




def get_magnitude_completeness(idx_good, df_1and2, clean_mag1, clean_mag2, full_mag1, full_mag2, tile, realization, band, mag_hdr1, mag_hdr2, mag_err1, mag_err2, fn_mag_completeness_log, inj_percent, balrog_run, base_path_to_catalogs): #FIXME inj?
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
		fn_truth_cat = manipulate_catalogs.get_catalog_filename(cat_type=MATCH_CAT1, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ1, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)	
		mag_hdr = mag_hdr1
		# Use error in measured catalogs only #
		mag_err1 = None 
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
		flag_hdr = FLAGS_HDR1[:-2]
		cm_flag_hdr = CM_FLAGS_HDR1[:-2] 

	if 'truth' in MATCH_CAT2:
		fn_truth_cat = manipulate_catalogs.get_catalog_filename(cat_type=MATCH_CAT2, inj_percent=inj_percent, realization=realization, tile=tile, band=band, inj=INJ2, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
		mag_hdr = mag_hdr2
		ra_hdr, dec_hdr = RA_HDR2, DEC_HDR2
		# Use error in measured catalogs only #
		mag_err2 = None 
		# Note: [:-2] to get rid of the suffix '_1' or '_2' added by STILTS because truth catalogs have not been matched #
                flag_hdr = FLAGS_HDR2[:-2]
                cm_flag_hdr = CM_FLAGS_HDR2[:-2]

	data = fitsio.read(fn_truth_cat, hdu=1)
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




def get_flux_plot_variables(tile, realization, band, df_1and2, flux_hdr1, flux_hdr2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, fn_percent_recovered_log, base_path_to_catalogs, balrog_run): 
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

	fullMag1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr1, band=band)
        fullMag2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=mag_hdr2, band=band)

        idxGood = analysis.get_good_indices_using_primary_flags(df_1and2=df_1and2, full_mag1=fullMag1, full_mag2=fullMag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2, band=band)[0]

	### Get flux ###
	fullFlux1 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr1, band=band)
        fullFlux2 = get_floats_from_string(df=df_1and2, four_elmt_arrs_hdr=flux_hdr2, band=band)
        cleanFlux1, cleanFlux2 = fullFlux1[idxGood], fullFlux2[idxGood]


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
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanFlux1, inj_percent=INJ1_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)
        if MATCH_CAT2 in ('gal_truth', 'star_truth'):
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = analysis.get_percent_recovered(full_data=fullMag1, clean_data=cleanFlux1, inj_percent=INJ2_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)

	if MATCH_CAT1 not in ('gal_truth', 'star_truth') and MATCH_CAT2 not in ('gal_truth', 'star_truth'):
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = None, None


	return __norm_flux_diff, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved, haxLabel



def get_total_flux_error(df_1and2, idx_good, band):
	"""Compute total flux error. Total means taken from both `MATCH_CAT1` and `MATCH_CAT2`.

	Parameters
	----------

	Returns
	-------
	"""

	### Flux error from measured catalog only ###
        if MATCH_CAT1 == 'coadd':
                #__flux_err1 = df1and2[NEW_COADD_FLUX_ERR_GRIZ_HDR+'_1'][idx_good]
		__flux_err1 = get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_COADD_FLUX_ERR_GRIZ_HDR+'_1')[idx_good]
        if MATCH_CAT2 == 'coadd':
                #__flux_err2 = df1and2[NEW_COADD_FLUX_ERR_GRIZ_HDR+'_2'][idx_good]
		__flux_err2 = get_floats_from_string(df=df_1and2, band=band, four_elmt_arrs_hdr=NEW_COADD_FLUX_ERR_GRIZ_HDR+'_2')[idx_good]

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
		__flux_err1 = np.sqrt(get_matrix_diagonal_element(df=df_1and2, band=band, sq_matrices_hdr=CM_FLUX_COV_HDR1)[idx_good])
	if 'y3_gold' not in MATCH_CAT2 and 'truth' not in MATCH_CAT2 and MATCH_CAT2 != 'coadd':
		__flux_err2 = np.sqrt(get_matrix_diagonal_element(df=df_1and2, band=band, sq_matrices_hdr=CM_FLUX_COV_HDR2)[idx_good])

	if __flux_err1 is not None and __flux_err2 is not None: __flux_err = np.sqrt(np.power(__flux_err2, 2)+np.power(__flux_err1, 2))
	if __flux_err1 is None: __flux_err = np.sqrt(np.power(__flux_err2, 2))
	if __flux_err2 is None: __flux_err = np.sqrt(np.power(__flux_err1, 2))
	
	return __flux_err




def get_magnitude_plot_variables(band, df_1and2, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization, tile, fn_flag_log, plot_title, fn_plot, fn_percent_recovered_log, balrog_run, base_path_to_catalogs):
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


	### Percent recovered if applicable ### 
	if MATCH_CAT1 in ('gal_truth', 'star_truth'):
		# Note that len(full_mag1) = len(full_mag2) & len(clean_mag1) = len(clean_mag2) so it does not matter which is passed to get_percent_recovered() #
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_percent_recovered(full_data=fullMag1, clean_data=cleanMag1, inj_percent=INJ1_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)
        #TODO % recovered with star_truth (no flags)
	if MATCH_CAT2 in ('gal_truth', 'star_truth'):
                percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = get_percent_recovered(full_data=fullMag1, clean_data=cleanMag1, inj_percent=INJ2_PERCENT, band=band, tile=tile, realization=realization, fn_percent_recovered_log=fn_percent_recovered_log, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)

	if MATCH_CAT1 not in ('gal_truth', 'star_truth') and MATCH_CAT2 not in ('gal_truth', 'star_truth'):
		percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved = None, None


	return magErr1, magErr2, cleanMag1, cleanMag2, idxGood, fullMag1, fullMag2, magAxLabel1, magAxLabel2, magVaxLabel, percentRecoveredFlagsIncluded, percentRecoveredFlagsRemoved 





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
	__coadd_mag_griz (list of str)
		Stores magnitude of each band in form '(mag_g, mag_r, mag_i, mag_z)'
	__coadd_mag_err_griz (list of str)
		Stores error in magnitude of each band in form '(mag_g, mag_r, mag_i, mag_z)'
	"""

	print 'Getting g-, r-, i-, z-band magnitudes and fluxes and corresponding errors for coadd catalogs...'

	# Files have not yet been matched, and do not have hdr_1 #
	mag_hdr = mag_hdr[:-2]
	mag_err_hdr = mag_err_hdr[:-2]

	# Open and read FITS files #
	data_g = fitsio.read(fn_coadd_cat_g, hdu=1); data_r = fitsio.read(fn_coadd_cat_r, hdu=1); data_i = fitsio.read(fn_coadd_cat_i, hdu=1); data_z = fitsio.read(fn_coadd_cat_z, hdu=1)

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	mag_err_g = data_g[mag_err_hdr]; mag_err_r = data_r[mag_err_hdr]; mag_err_i = data_i[mag_err_hdr]; mag_err_z = data_z[mag_err_hdr]

	# Get fluxes #
	f_g = data_g[flux_hdr]; f_r = data_r[flux_hdr]; f_i = data_i[flux_hdr]; f_z = data_z[flux_hdr]

	# Get flux error #
	flux_err_g = data_g[flux_err_hdr]; flux_err_r = data_r[flux_err_hdr]; flux_err_i = data_i[flux_err_hdr]; flux_err_z = data_z[flux_err_hdr]

	__coadd_flux_griz, __coadd_flux_err_griz, __coadd_mag_griz, __coadd_mag_err_griz = np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str)

	for i in np.arange(0, len(m_g)):
		__coadd_flux_griz[i] = '({}, {}, {}, {})'.format(f_g[i], f_r[i], f_i[i], f_z[i])
		__coadd_flux_err_griz[i] = '({}, {}, {}, {})'.format(flux_err_g[i], flux_err_r[i], flux_err_i[i], flux_err_z[i])
		__coadd_mag_griz[i] = '({}, {}, {}, {})'.format(m_g[i], m_r[i], m_i[i], m_z[i]) 
		__coadd_mag_err_griz[i] = '({}, {}, {}, {})'.format(mag_err_g[i], mag_err_r[i], mag_err_i[i], mag_err_z[i])

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

	__star_truth_mag_griz = np.empty(len(m_g), dtype=str)

	for i in np.arange(0, len(m_g)):
		__star_truth_mag_griz[i] = '({}, {}, {}, {})'.format(m_g[i], m_r[i], m_i[i], m_z[i]) 

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
	hdr_g = '{}_G_{}'.format(mag_hdr[:-2], mag_hdr[-1])
        hdr_r = '{}_R_{}'.format(mag_hdr[:-2], mag_hdr[-1])
        hdr_i = '{}_I_{}'.format(mag_hdr[:-2], mag_hdr[-1])
        hdr_z = '{}_Z_{}'.format(mag_hdr[:-2], mag_hdr[-1])
	'''	
	hdr_g = mag_hdr[:-2] + '_G' + mag_hdr[-2:]; hdr_r = mag_hdr[:-2] + '_R' + mag_hdr[-2:]
	hdr_i = mag_hdr[:-2] + '_I' + mag_hdr[-2:]; hdr_z = mag_hdr[:-2] + '_Z' + mag_hdr[-2:]
	'''
	# Read magnitudes from DataFrame #
	m_g = df_1and2[hdr_g]; m_r = df_1and2[hdr_r]; m_i = df_1and2[hdr_i]; m_z = df_1and2[hdr_z]

	__y3_gold_mag_griz = np.empty(len(m_g), dtype=str)

        for i in np.arange(0, len(m_g)):
                __y3_gold_mag_griz[i] = '({}, {}, {}, {})'.format(m_g[i], m_r[i], m_i[i], m_z[i]) 

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
        hdr_g = '{}_G_{}'.format(flux_hdr[:-2], flux_hdr[-1])
        hdr_r = '{}_R_{}'.format(flux_hdr[:-2], flux_hdr[-1])
        hdr_i = '{}_I_{}'.format(flux_hdr[:-2], flux_hdr[-1])
        hdr_z = '{}_Z_{}'.format(flux_hdr[:-2], flux_hdr[-1])

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

	__y3_gold_flux_griz, __y3_gold_flux_cov_griz = np.empty(len(f_g), dtype=str), np.empty(len(f_g), dtype=str)

	for i in np.arange(0, len(f_g)):
	#TODO .format()
		#__y3_gold_flux_griz.append("'("+ str(f_g[i]) + ', ' + str(f_r[i]) + ', ' + str(f_i[i]) + ', ' + str(f_z[i]) + ")'")
		__y3_gold_flux_griz[i] = '({}, {}, {}, {})'.format(f_g[i], f_r[i], f_i[i], f_r[i])
		__y3_gold_flux_cov_griz[i] = '({}, {}, {}, {})'.format(flux_cov_g[i], flux_cov_r[i], flux_cov_i[i], flux_cov_z[i]) 


	return __y3_gold_flux_griz, __y3_gold_flux_cov_griz




def get_zipped_coadd_magnitudes(fn_base_g, fn_base_r, fn_base_i, fn_base_z):
	"""TODO

	Parameters
	----------

	Returns
	-------
	"""

	if VERBOSE_ING: print 'Getting g-, r-, i-, z-band magnitudes and flags for zipped coadd catalogs...\n'

	mag_hdr = 'MAG_AUTO'
	mag_err_hdr = 'MAGERR_AUTO'
	flag1_hdr = 'FLAGS'
	flag2_hdr = 'IMAFLAGS_ISO'

	# Open and read FITS files #
	data_g = fitsio.read(fn_base_g, hdu=1); data_r = fitsio.read(fn_base_r, hdu=1); data_i = fitsio.read(fn_base_i, hdu=1); data_z = fitsio.read(fn_base_z, hdu=1)

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	err_g = data_g[mag_err_hdr]; err_r = data_r[mag_err_hdr]; err_i = data_i[mag_err_hdr]; err_z = data_z[mag_err_hdr]

	# Get flags #
	flag1_g = data_g[flag1_hdr]; flag1_r = data_r[flag1_hdr]; flag1_i = data_i[flag1_hdr]; flag1_z = data_z[flag1_hdr]
	flag2_g = data_g[flag2_hdr]; flag2_r = data_r[flag2_hdr]; flag2_i = data_i[flag2_hdr]; flag2_z = data_z[flag2_hdr]


	__base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz = np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str), np.empty(len(m_g), dtype=str)

	for i in np.arange(0, len(m_g)):
		__base_mag_griz[i] = '({}, {}, {}, {})'.format(m_g[i], m_r[i], m_i[i], m_z[i])
		__base_mag_err_griz[i] = '({}, {}, {}, {})'.format(err_g[i], err_r[i], err_i[i], err_z[i])
		__base_flag1_griz[i] = '({}, {}, {}, {})'.format(flag1_g[i], flag1_r[i], flag1_i[i], flag1_z[i])
		__base_flag2_griz[i] = '({}, {}, {}, {})'.format(flag2_g[i], flag2_r[i], flag2_i[i], flag2_z[i])

	return __base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz





