"""
Get injection percent.
Note that star truth catalogs and galaxy truth catalogs in the same `base_path_to_catalogs` may not have the same number of injected objects.
"""

import fitsio
import os

# From BalVal #
import manipulate_catalogs




def get_injection_percent(cat_types, tile, realization, base_path_to_catalogs, balrog_run):
	"""Get the injection percent of the catalogs in `base_path_to_catalogs`.
	Parameters
	----------
	cat_types (list)
		Types of both catalogs to be analysed.
	Returns
	-------
	__inj_percent (float)
		Balrog-injection percent of `MATCH_CAT1` and `MATCH_CAT2`. 
	"""

	__objs_per_tile = 50000.0

	### Get truth catalog ###

	# Star truth #
	if 'star_truth' in cat_types:
		# `inj_percent` is only used for TAMU Balrog runs, in which case `INJ1_PERCENT` is force set and this function will not be called. `band` is only used for coadd catalogs. #
		__fn_cat = get_catalog.get_catalog_filename(cat_type='star_truth', inj=True, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

	if 'gal_truth' in cat_types:
		__fn_cat = get_catalog.get_catalog_filename(cat_type='gal_truth', inj=True, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)


	### Read truth catalog ###
	data = fitsio.read(__fn_cat, hud=1)
	__num_rows = len(data)

	# Approximate objects per tile: 50,000 #
	__inj_percent = 100*(1.0*__num_rows/__objs_per_tile)

	__inj_percent = int(__inj_percent)

	
	return __inj_percent
