"""
Match catalogs with a call to `ms_matcher`.
This requires retrieving the catalog filenames.
"""
#TODO rename match_catalogs.py get_and_match_cats?
import os
from set_plot_and_catalog_constants import * 




def get_catalog_filename(cat_type, inj, inj_percent, realization, tile, band, base_path_to_catalogs, balrog_run):
	"""Get the complete filename of a catalog to match.
	Note that different Balrog runs (`balrog_run`) may have different directory structures and filename structures; an (incomplete) log is at https://docs.google.com/document/d/1hb3OlmU02jS4KkvCxvw80zr1fd15LbgzOK7tps1uqBM/edit?usp=sharing

	Parameters
	----------
	cat_type
		Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd', 'y3_gold'. Set by `MATCH_CAT1` or `MATCH_CAT2`.

	inj_percent (int)

	inj (bool)

	realization (str)

	tile

	band (str)
		Only used with coadd catalogs. Ignored if `cat_type` is not 'coadd'.

	base_path_to_catalogs

	balrog_run

	Returns
	-------
	__fn_cat (str)
		Complete catalog filename.
	"""

	# Note: only need inj percent for TAMU cat #

	### Balrog-injected MOF/SOF ###
	if cat_type in ('mof', 'sof') and inj: #and inj_percent == 20:
		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, 'real_0_'+tile+'_{}.fits'.format(cat_type))): 
			__fn_cat = os.path.join(base_path_to_catalogs, tile, 'real_0_'+tile+'_{}.fits'.format(cat_type))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'real_0_'+tile+'_{}.fits'.format(cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'real_0_'+tile+'_{}.fits'.format(cat_type))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}'.format(cat_type), tile+'_{}.fits'.format(cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}'.format(cat_type), tile+'_{}.fits'.format(cat_type))	


	### Base MOF/SOF ###
	if cat_type in ('mof', 'sof') and inj is False:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', tile, '{}'.format(cat_type), tile+'_{}.fits'.format(cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', tile, '{}'.format(cat_type), tile+'_{}.fits'.format(cat_type))


	### Galaxy truth catalogs ###
	if cat_type == 'gal_truth' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_gals.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_gals.fits')

		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, tile+'_0_balrog_truth_cat_gals.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, tile, tile+'_0_balrog_truth_cat_gals.fits')

		if os.path.isfile(os.path.join(base_path_to_catalogs, tile+'_0_balrog_truth_cat_gals.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, tile+'_0_balrog_truth_cat_gals.fits')

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat.fits')


	### Star truth catalogs ###
	if cat_type == 'star_truth' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, tile+'_0_balrog_truth_cat_stars.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, tile, tile+'_0_balrog_truth_cat_stars.fits')

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_stars.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, tile+'_'+realization+'_balrog_truth_cat_stars.fits')


	### Balrog-injected coadd catalogs ###
	if cat_type == 'coadd' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, 'coadd', tile+'_'+band+'_cat.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, 'coadd', tile+'_'+band+'_cat.fits')


	### Base coadd catalogs ###
	if cat_type == 'coadd' and inj is False:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', tile+'_'+band+'_cat.fits')):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', tile+'_'+band+'_cat.fits')

		#FIXME
		if os.path.isfile(get_and_reformat_base_catalog(tile=tile, realization=realization)):
			__fn_cat = get_and_reformat_base_catalog(tile=tile, realization=realization)


	### Y3 Gold catalogs ###
	# !!!!! #
	if 'y3_gold' in cat_type:
		__fn_cat = os.path.join('data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'y3_gold_catalogs', '{}_{}.fits'.format(tile, cat_type))
	'''
	if cat_type == 'y3_gold_2_0' and inj is False:
		__fn_cat = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile+'_y3_gold_2_0.fits')
	if cat_type == 'y3_gold_2_2' and inj is False:
		__fn_cat = os.path.join('/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/', tile+'_y3_gold_2_2.fits')
	'''

	### Deep SN catalogs ###
	# !!!!! #
	if cat_type == 'deep_sn_mof' and tile == 'TODO':
		__fn_cat = os.path.join('data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'deep_sn_catalogss', 'SN-C3_C28_r3499p02-Y3A2_DEEP-mof-001.fits')

	if cat_type == 'deep_sn_sof' and tile == 'TODO':
		__fn_cat = os.path.join('data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'deep_sn_catalogss', 'SN-C3_C28_r3499p02-Y3A2_DEEP-sof-001.fits')


	if balrog_run == 'TAMU_Balrog' and 'y3_gold' not in cat_type:
		__fn_cat = get_tamu_catalog_filename(cat_type=cat_type, inj_percent=inj_percent, realization=realization, tile=tile, inj=inj)


	return __fn_cat




def get_tamu_catalog_filename(cat_type, inj, inj_percent, realization, tile):
	"""Get catalog for TAMU runs.

	Parameters
	----------
	inj_percent (int)
	inj (bool)
	cat_type (str)
		Catalog type. This is set by `MATCH_CAT1` or `MATCH_CAT2`.
	realization (str)

	Returns
	-------
	__fn_tamu_cat (str)
		Complete catalog filename.
	"""

	# 20% injection directory #	
	if inj_percent == 20:
		if cat_type in ('mof', 'sof') and inj:
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile + '_20', 'real_' + realization + '_' + tile + '_{}.fits'.format(cat_type))

		if cat_type in ('mof', 'sof') and inj is False:
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile + '_20', 'base_' + tile + '_{}.fits'.format(cat_type))

		if cat_type == 'gal_truth':
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile + '_20', tile + '_' + realization + '_balrog_truth_cat_gals.fits')

		if cat_type == 'star_truth':
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile + '_20', tile + '_' + realization + '_balrog_truth_cat_stars.fits')


	# 10% injection directory #
	if inj_percent == 10:
		if cat_type in ('mof', 'sof') and inj:
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile, 'real_' + realization + '_' + tile + '_{}.fits'.format(cat_type))

		if cat_type in ('mof', 'sof') and inj is False:
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile, 'base_' + tile + '_{}.fits'.format(cat_type))

		if cat_type == 'gal_truth':
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile, tile + '_' + realization + '_balrog_truth_cat_gals.fits')

		if cat_type == 'star_truth':
			__fn_tamu_cat = os.path.join(base_path_to_catalogs, tile, tile + '_' + realization + '_balrog_truth_cat_stars.fits')


	return __fn_tamu_cat 




'''
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
'''



#get_and_reformat()

#match_catalogs(balrog_run, base_path_to_catalogs, output_directory)

#TODO other reformatting catalog processes here? or reformat_catalogs.py that includes get_and_reformat as above?

#fof_matcher()
