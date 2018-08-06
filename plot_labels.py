"""
Get plot labels.
Plot labels will be overwritten if `OVERWRITE_AXLABELS=True` in `set_constants.py`.

Comments are ABOVE the code they refer to.
"""

from difflib import SequenceMatcher

from set_constants import *
from catalog_headers import TITLE_PIECE1, TITLE_PIECE2




def get_magnitude_axlabel(inj, mag_hdr, meas_or_true_cat, match_cat, band, inj_percent):
	"""Get labels for the horizontal axis. 
	Note that `'{}'.format()` does not render TeX reliably.

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
	hax_label (str)
		Label for the horizontal axis. Includes LaTeX \bf{} formatting.
	"""

	# 'mag_2' --> 'mag_{band}_true' with {band} bolded. Note that `format()` does not render TeX #
	__mag_axlabel = '%s_$\\bf{%s}$_%s' % (mag_hdr[:-2], band, meas_or_true_cat)

	if 'y3_gold' in match_cat:
		__mag_axlabel = 'Y3_%s' % __mag_axlabel

	# 'mag_{band}_true' --> 'base_mag_{band}_true' or 'xx%_inj_mag_{band}_true' #
	if inj:
		__mag_axlabel = '%s%%_inj_%s' % (inj_percent, __mag_axlabel)
	if inj is False:
		__mag_axlabel = 'base_%s' %  __mag_axlabel

	if OVERWRITE_AXLABELS:
		__mag_axlabel = 'mag_\\bf{%s}' % band

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
		Contains LaTeX formatting. Ex: 'inj_\bf{(g-r)}_true'.
	"""

	if inj:
		__color_axlabel = '{}%_inj_'.format(inj_percent) + '$\\bf{(%s)}$_%s' % (WRITE_COLORS[band], meas_or_true_cat) 
	if inj is False:
		__color_axlabel = 'base_$\\bf{(%s)}$_%s' % (WRITE_COLORS[band], meas_or_true_cat)

	if 'y3_gold' in match_cat:
		__color_axlabel = 'Y3_%s' % __color_axlabel


	if OVERWRITE_AXLABELS:
		__color_axlabel = '$\\bf{(%s)}$_%s' % (WRITE_COLORS[band], meas_or_true_cat)

	return __color_axlabel 





def get_short_difference_axlabel(axlabel_a, axlabel_b, band):
	"""Create a short axis label for the difference (`axlabel_a` minus `axlabel_b`) between two observables.
	If the order of subtraction is not `axlabel_a` - `axlabel_b`, call this function via `axlabel_a=label2`, `axlabel_b=label1`.

	Parameters
	----------
	axlabel_a, axlabel_b (str)
		Axis labels for an observable in `MATCH_CAT1` or `MATCH_CAT2`.  `axlabel_a` does not need to correspond to `MATCH_CAT1`.
	Returns
	-------
	__short_axlabel (str)
	"""

	### Find longest common substring ###
	match = SequenceMatcher(None, axlabel_a, axlabel_b).find_longest_match(0, len(axlabel_a), 0, len(axlabel_b))	
	# `match.a` and `match.b` are indices for `axlabel_a` and `axlabel_b` #
	idx1, idx2 = match.a, match.b
	size = match.size
	# Note: this is the same as `axlabel_b[idx2: idx2+size]` # 
	shared_label = axlabel_a[idx1: idx1+size]

	### Recover substrings missing from `shared_label` ###
	# `axlabel_a` #
	__spare_label_begin_a = axlabel_a[:idx1]
	__spare_label_end_a = axlabel_a[idx1+size:]
	# `axlabel_b` #
	__spare_label_begin_b = axlabel_b[:idx2]
        __spare_label_end_b = axlabel_b[idx2+size:]
	
	# Look for '10%_inj' and '20%_inj' which causes '0%_inj' in `shared_label` #
	if __spare_label_begin_a in ['1', '2'] or __spare_label_begin_b in ['1', '2']:
		# '0%_inj' --> '_inj' #
		shared_label = shared_label[2:]
		# '1' --> '10%' #
		__spare_label_begin_a = '{}0%'.format(__spare_label_begin_a)
		__spare_label_begin_b = '{}0%'.format(__spare_label_begin2)

	__short_axlabel = '({}$-${})_{}_({}$-${})'.format(__spare_label_begin_a, __spare_label_begin_b, shared_label, __spare_label_end_a, __spare_label_end_b)

	# Check for empty `__spare_label_begin_a`, `__spare_label_begin_b`, `__spare_label_end_a`, and `__spare_label_end_b` #
	__short_axlabel = __short_axlabel.replace('($-$)', '')
	# Check for double underscores #
	__short_axlabel = __short_axlabel.replace('__', '_')
	# Check for leading underscore #
	if __short_axlabel[0] == '_':
		__short_axlabel = __short_axlabel[1:]


	if PLOT_FLUX:
		__short_axlabel = '%s/$\sigma_{flux\_meas}$' % __short_axlabel.replace('mag', 'flux')
		


	# Default labels #
	if OVERWRITE_AXLABELS:
		if PLOT_COLOR: __short_axlabel = '$\Delta(\\bf{%s})$' % WRITE_COLORS[band]
		if PLOT_FLUX: __short_axlabel = '$\Delta$flux_$\\bf{%s}/$meas_flux_err' % band
		if PLOT_MAG: __short_axlabel = '$\Delta$mag_$\\bf{%s}$' % band 


	return __short_axlabel





def get_plot_suptitle(realization, tile, number_of_stacked_realizations, number_of_stacked_tiles):
	"""Generate plot title.

	Parameters
	----------
	realization (str)
	tile (str)
	number_of_stacked_realizations (int)
		Number of catalogs in stacked realization catalog. Can be `None`.
	number_of_stacked_tiles (int)
		Number of catalogs in stacked tile catalog. Can be `None`.

	Returns
	-------
	__plot_title (str)
		Ex: '10% Inj MOF Cat & 10% Inj Truth Cat' 
	"""

	if STACK_REALIZATIONS:
		realization = 'stacked '+str(number_of_stacked_realizations)
	if STACK_TILES:
		tile = 'stacked ' + str(number_of_stacked_tiles)

	__plot_title = '{} & {}. Tile: {}'.format(TITLE_PIECE1, TITLE_PIECE2, tile)

	if realization is not None: 
		__plot_title = '{}. Realization: {}.'.format(__plot_title, realization)
	
	if RUN_TYPE == 'ok': 
		__plot_title = '{}. Unchanged FOF groups.'.format(__plot_title)
	if RUN_TYPE == 'rerun':
		__plot_title = ' Changed FOF groups.'.format(__plot_title)

	if NORMALIZE:
		__plot_title = 'Normalized. ' + __plot_title 

	return __plot_title 



def get_colorbar_for_magnitude_plot_axlabel(df_1and2, cm_t_hdr, cm_t_err_hdr, idx_good, clean_mag1, clean_mag2, meas_or_true_cat, inj_percent, inj):
	"""Get the label for the colorbar of the magnitude plot.
	This function will return `None` if no colorbar is to be added to the plot.

	Parameters
	----------
	df_1and2 (pandas DataFrame)
		DataFrame for the matched (join=1and2) catalog.

	cm_t_hdr, cm_t_err_hdr (str)
		Matched (join=1and2) catalog headers for cm_T (size squared) and cm_T_err. Can be `None`.
		idx_good (list of ints)
		Indices of objects without flags.		

	clean_mag1, clean_mag2 (list of floats)
		meas_or_true_cat (str)
		Allowed values: 'true' 'meas'	

	Returns
	-------
	__cbar_label (str)
		Label for colorbar. Includes LaTeX \bf{} formatting. Can be `None`.
	"""

	if 'true' in meas_or_true_cat:
		sys.exit('ERROR. Colorbars should describe measured catalog values, not truth catalog values.')


	### Colorbar label ###
	# Prefix to labels #
	if inj:
		pref = '{}%_inj'.format(inj_percent)
	if inj is False:
		pref = 'base'
	if 'y3_gold' in match_cat:
		pref = 'Y3_{}'.format(pref)

	# `[:-2]` to remove the suffix '_1' or '_2' from the header #
	if CM_T_CBAR:
		__cbar_label = '{}_{}_{}'.format(pref, cm_t_hdr[:-2], meas_or_true_cat)
	if CM_T_ERR_CBAR:
		__cbar_label = '{}_{}_{}'.format(pref, cm_t_err_hdr[:-2], meas_or_true_cat)


	return __cbar_label

