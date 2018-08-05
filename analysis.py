"""
Analysis for MATCHED catalog.
"""

import csv
import fitsio
import numpy as np
import sys

# From BalVal #
from set_constants import *
from catalog_headers import *
import manipulate_catalogs




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
		# Each row in `__ad_flags` corresponds to a separate flag header's values #
		# Each column refers to an object #
		# We want all the columns to be zero for an object #
		# `np.where()` along an axis? Work-around #
		# Cannot predict number of idx --> use list #
		__ad_idx_list = []
		for i in np.arange(0, __ad_flags.shape[0]):
			__ad_idx_list.append(np.where(__ad_flags[i] == 0)[0])
		# Indices are only good if they are good for every flag. Find repeats #
		__ad_idxs = [item for sublist in __ad_idx_list for item in sublist]
		__bad = np.unique(np.array(__ad_idxs), return_inverse=True)[1] #FIXME once this is working sue otherbool=True

		__ad_idx_good = set(__ad_idxs) - set(__bad)

		__idx_good_tot = list(set(__idx_good) - set(__ad_idx_good))
	else:
		__idx_good_tot = __idx_good


	if PLOT_FLAGGED_OBJS:
		__all_idx = np.arange(0, len(df_1and2.index))
		__idx_bad = list(set(__all_idx) - set(__idx_good_tot))
	if PLOT_FLAGGED_OBJS is False:
		__idx_bad = None

		
        if VERBOSE_ED:
		__flag_list = ''
		if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
			__flag_list += ' SEXTRACTOR_FLAGS_'+band.upper() + ', IMAFLAGS_ISO_'+band.upper()
		if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2 and Y3_FIT.upper() == 'MOF':
			 __flag_list += ', MOF_CM_MOF_FLAGS'+suf
		if 'base' in (MATCH_CAT1, MATCH_CAT2):
			__flag_list += NEW_BASE_FLAG_GRIZ_HDR+suf + ', ' + NEW_BASE_IMAFLAG_ISO_GRIZ_HDR+suf

		print 'For {}-band, eliminated {} objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for the following headers: {}, {}, {}, {}, {}'.format(band, len(full_mag1) - len(__idx_good), flag_hdr1, flag_hdr2, cm_flag_hdr1, cm_flag_hdr2, __flag_list)


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

	Parameters
	----------

	Returns
	-------
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




def get_colorbar_data():
	"""Get colorbar for magnitude plot.

	Parameters
	----------

	Returns
	-------
	"""

	return 0
