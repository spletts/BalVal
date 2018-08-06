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
	"""Get directory for various outputs. 

	Parameters
	----------
	low_level_dir (str OR list of str)

	Returns
	-------
	__dir (str)
		Directory path.
	"""

	if isinstance(low_level_dir, str):
		__dir = os.path.join(output_directory, 'outputs', balrog_run, match_type, tile, realization, low_level_dir)

	if isinstance(low_level_dir, list):
		if low_level_dir[1] == 'magnitude' and NORMALIZE:
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
	"""Generate names for log files.
	Directory structure: /`OUTPUT_DIRECTORY`/log_files/`BALROG_RUN`/`match_type`/{tile}/{realization}/log_files/.

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
	__fn_color_log (str)
	__fn_mag_diff_outliers_log (str)
	__fn_mag_completeness_log (str)
	__fn_gauss_aper_log (str)
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
	"""

	Parameters
	----------

	Returns
	-------
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
        """

        Parameters
        ----------

        Returns
        -------
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
        """

        Parameters
        ----------

        Returns
        -------
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

	if NORMALIZE and PLOT_MAG:
		__fn_plot = os.path.join(plotDirectory, 'norm_{}'.format(__fn_plot))


	return __fn_plot





def get_plot_filename(balrog_run, match_type, output_directory, realization, tile):
        """Generate name of plot.
	Relies on directory structure: /`OUTPUT_DIRECTORY`/plots/`BALROG_RUN`/`match_type`/{tile}/{realization}/plots/`plot_obs`/ `plot_obs` can be: 'color' 'magnitude' 'flux'.
	The FOF analysis plots are saved differently: /`OUTPUT_DIRECTORY`/plots/`BALROG_RUN`/`match_type`/{tile}/{realization}/plots/{plot_obs}/'fof_analysis'/

        Parameters
	----------
	realization (str) 
	tile (str)

        Returns
	-------
	__fn_plot (str)
		Complete filename for the plot. This will be used in `plt.savefig()` if `SAVE_PLOT=True`.
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
	"""Generate names for region files of different join types (join types specified in ms_matcher or ms_fof_matcher). 
	Relies on directory structure `/OUTPUT_DIRECTORY/outputs/BALROG_RUN/match_type/{tile}/{realization}/region_files/`

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
	"""

	Parameters
	----------
	tile (str)
		Allowed value: 'stack'

	realization (str)
		Allowed value: 'stack'

	Returns
	-------
	"""

	matchedCatalogDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

	repeated_str = '{}_{}_{}'.format(tile, realization, match_type)

	__fn_match_1and2 = os.path.join(matchedCatalogDirectory, '{}_match1and2.csv'.format(repeated_str))
	__fn_match_1not2 = os.path.join(matchedCatalogDirectory, '{}_match1not2.csv'.format(repeated_str))
	__fn_match_2not1 = os.path.join(matchedCatalogDirectory, '{}_match2not1.csv'.format(repeated_str))

	return __fn_match_1and2, __fn_match_1not2, __fn_match_2not1
	



def get_ngmix_compatible_catalog_filename(balrog_run, fn_unmatched_cat, match_type, output_directory, realization, tile):
	"""

	Parameters
	----------

	Returns
	-------
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
	"""

        Parameters
        ----------

        Returns
        -------
        """

	catalogFileDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=match_type) 

	__fn_coadd_for_matcher = os.path.join(catalogFileDirectory, '{}_i_cat_reformatted.fits'.format(tile))

	return __fn_coadd_for_matcher 



def make_region_files():


	return 0


def write_log_file_headers(fn_mag_err_log, fn_flag_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log, fn_flux_sigma_clip_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log):
	"""Write headers to log files.

	Parameters
	----------
	Returns
	-------
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
		writer.writerow(['TILE', 'REALIZATION', 'BAND', 'GAUSSIAN_APER_APPLIED', 'NUM_OBJS_FLAGS_RM', 'N-SIGMA_CLIP' 'NUM_OBJS_CLIPPED']) 
	csvfile.close()


	return fn_flag_log, fn_mag_err_log, fn_color_log, fn_mag_diff_outliers_log, fn_mag_completeness_log, fn_flux_sigma_clip_log, fn_percent_recovered_log, fn_num_objs_in_1sig_mag_log
