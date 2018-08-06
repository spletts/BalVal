"""
Manipulate catalogs in the following ways:
  Match catalogs with a call to `stilts_matcher`.
  This requires retrieving the catalog filenames.
  If necessary, reformat the catalog(s) (prior to matching).
    Reformat e.g.:
    Add a column necessary to perform Gaussian aperture measured flux (via `ngmixer`)
    Add columns to coadd catalog.
    ...

Comments are ABOVE the code they refer to.
"""

import fitsio
import ngmixer
import numpy as np
import os
import pandas as pd
import subprocess
import sys
from astropy.io import fits
from astropy.table import Table, Column

from set_constants import * 
from catalog_headers import *
import outputs




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
		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, 'real_0_{}_{}.fits'.format(tile, cat_type))): 
			__fn_cat = os.path.join(base_path_to_catalogs, tile, 'real_0_{}_{}.fits'.format(tile, cat_type))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'real_0_{}_{}.fits'.format(tile, cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'real_0_{}_{}.fits'.format(tile, cat_type))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, cat_type, '{}_{}.fits'.format(tile, cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, cat_type, '{}_{}.fits'.format(tile, cat_type))	


	### Base MOF/SOF ###
	if cat_type in ('mof', 'sof') and inj is False:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', tile, cat_type, '{}_{}.fits'.format(tile, cat_type))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', tile, cat_type, '{}_{}.fits'.format(tile, cat_type))


	### Galaxy truth catalogs ###
	if cat_type == 'gal_truth' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat_gals.fits'.format(tile, realization))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat_gals.fits'.format(tile, realization))

		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, '{}_0_balrog_truth_cat_gals.fits'.format(tile))):
			__fn_cat = os.path.join(base_path_to_catalogs, tile, '{}_0_balrog_truth_cat_gals.fits'.format(tile))

		if os.path.isfile(os.path.join(base_path_to_catalogs, '{}_0_balrog_truth_cat_gals.fits'.format(tile))):
			__fn_cat = os.path.join(base_path_to_catalogs, '{}_0_balrog_truth_cat_gals.fits'.format(tile))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat.fits').format(tile, realization)):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat.fits'.format(tile, realization))


	### Star truth catalogs ###
	if cat_type == 'star_truth' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, tile, '{}_0_balrog_truth_cat_stars.fits'.format(tile))):
			__fn_cat = os.path.join(base_path_to_catalogs, tile, '{}_0_balrog_truth_cat_stars.fits'.format(tile))

		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat_stars.fits'.format(tile, realization))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, '{}_{}_balrog_truth_cat_stars.fits'.format(tile, realization))


	### Balrog-injected coadd catalogs ###
	if cat_type == 'coadd' and inj:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, 'coadd', '{}_{}_cat.fits'.format(tile, band))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, 'coadd', '{}_{}_cat.fits'.format(tile, band))


	### Base coadd catalogs ###
	if cat_type == 'coadd' and inj is False:
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', '{}_{}_cat.fits'.format(tile, band))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', '{}_{}_cat.fits'.format(tile, band))

		#FIXME
		if os.path.isfile(get_and_reformat_base_catalog(tile=tile, realization=realization)):
			__fn_cat = get_and_reformat_base_catalog(tile=tile, realization=realization)


	### Y3 Gold catalogs ###
	# !!!!! #
	if 'y3_gold' in cat_type:
		# Download Y3 Gold catalog if it has not already been downloaded #
		if os.path.isfile(os.path.join('/data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'y3_gold_catalogs', '{}_{}.fits'.format(tile, cat_type))) is False:
			subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/query_by_tile', tile])

		__fn_cat = os.path.join('/data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'y3_gold_catalogs', '{}_{}.fits'.format(tile, cat_type))


	### Deep SN catalogs ###
	# !!!!! #
	if cat_type == 'deep_sn_mof' and tile == 'DES0328-2749':
		__fn_cat = os.path.join('/data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'deep_sn_catalogs', 'SN-C3_C28_r3499p02-Y3A2_DEEP-mof-001.fits')

	if cat_type == 'deep_sn_sof' and tile == 'DES0328-2749':
		__fn_cat = os.path.join('/data', 'des71.a', 'data', 'mspletts', 'balrog_validation_tests', 'deep_sn_catalogs', 'SN-C3_C28_r3499p02-Y3A2_DEEP-sof-001.fits')


	### Balrog runs for TAMU had both 10% and 20% injections. The catalog filename is dependent on the injection percent. ###
	if balrog_run == 'TAMU_Balrog' and 'y3_gold' not in cat_type:
		__fn_cat = get_tamu_catalog_filename(cat_type=cat_type, inj_percent=inj_percent, realization=realization, tile=tile, inj=inj)


	if cat_type == 'fofslist':
		if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, FOF_FIT, '{}_fofslist.fits'.format(tile))):
			__fn_cat = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, FOF_FIT, '{}_fofslist.fits'.format(tile))


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




#def get_and_reformat_base_catalog(tile, realization):
	"""Handle new format for base catalogs. Zipped and SExtractor?

	Parameters
	----------

	Returns
	-------
	"""
	'''
	if VERBOSE_ING: print 'Unzipping and reformatting base coadd catalogs...\n'

	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `get_and_reformat_base_catalog()`. Press enter to proceed, ctrl+c to stop')


	matchedCatalogDirectory = get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare')

	__fn_base_cat_for_matcher = os.path.join(matchedCatalogDirectory, tile+'_base_i.fits')

	if os.path.isfile(__fn_base_cat_for_matcher) and __overwrite is False:
		if VERBOSE_ING: print ' Reformatted base catalog already exists. Not overwriting...\n'

	if os.path.isfile(__fn_base_cat_for_matcher) is False or __overwrite:
		if VERBOSE_ING: print ' Reformatting base catalog...'

		# Base catalogs are zipped #
		__fn_cat_bundle = os.path.join(base_path_to_catalogs, tile, 'base_cat.tgz')
		__tar = tarfile.open(__fn_cat_bundle, 'r')

	# Base catalogs are zipped #
	__fn_base_griz = []

	__tar.extractall(path=matchedCatalogDirectory)
	for b in ALL_BANDS:
		__fn_b = '_'.join([tile, b, 'cat.fits'])
		__fn_base_griz.append(os.path.join(matchedCatalogDirectory, __fn_b))
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


#TODO other reformatting catalog processes here? or reformat_catalogs.py that includes get_and_reformat as above?

#get_coadd_catalog_for_matcher()

#fof_matcher()




def reformat_catalog(cat_type, inj, inj_percent, realization, tile, fn_match, mag_hdr, mag_err_hdr, flux_hdr, flux_err_hdr, base_path_to_catalogs, balrog_run, match_type, output_directory):
	"""Reformat catalog PRIOR to matching

	Parameters
	----------

	Returns
	-------

	"""

	__overwrite = False
	if __overwrite: raw_input('To rewrite press Enter')

	if cat_type == 'coadd':
		__fn_cat = get_coadd_catalog_for_matcher(cat_type=cat_type, inj_percent=inj_percent,  realization=realization, tile=tile, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_err_hdr=flux_err_hdr) 


	### Make catalog compatible for Gaussian aperture catalog ###
        if PLOT_GAUSS_APER_FLUX and PLOT_FLUX:
                #TODO reformat_catalog(): calls gauss aper, reformat_base, all add cols added before matching. make post_match_catalog_manipulation?
                # Will need to rematch this catalog if Gaussian aperture columns are not present #
                if os.path.isfile(fn_match):
                        try:
                                df = pd.read_csv(fn_match)
				# Can also do `GAUSS_APER_FLUX_GRIZ_HDR+'_2'`. Both exist simultaneously #
                                temp = df[GAUSS_APER_FLUX_GRIZ_HDR+'_1']
                        except:
                                __overwrite = True
                                if VERBOSE_ING: print 'Forcing overwrite to add Gaussian aperture flux measurements...'


		fnUnmatchedCatalog = get_catalog_filename(cat_type=cat_type, inj=inj, inj_percent=inj_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

                if os.path.isfile(fn_match) is False or __overwrite:
                        # Delete these files after the matched catalog is created #
                        __fn_cat = add_gaussian_aperture_flux_measurements_to_catalog(fn_unmatched_cat=fnUnmatchedCatalog, tile=tile, realization=realization, cat_type=cat_type, balrog_run=balrog_run, match_type=match_type, output_directory=output_directory)


		# We still need filename #
		if os.path.isfile(fn_match) and __overwrite is False:
			__fn_cat = make_ngmixer_gaussap_compatible_catalog(fn_unmatched_cat=fnUnmatchedCatalog, tile=tile, realization=realization, cat_type=cat_type, balrog_run=balrog_run, match_type=match_type, output_directory=output_directory)

	return __fn_cat




def match_catalogs(realization, tile, inj1, inj2, inj1_percent, inj2_percent, output_directory, balrog_run, base_path_to_catalogs):
	"""Match two catalogs on RA and Dec with a tolerance of 1 arcsecond, via STILTS.

	Parameters
	----------
		realization (str) 
		tile (str)
		band (str)
	Returns
	-------
		__fn_match_1and2 (str) -- Name of catalog matched via join=1and2. Note that headers of matched catalogs will have '_1' appended or '_2' appended.
		__fn_match_1not2 (str) -- Name of catalog matched via join=1not2.
		__fn_match_2not1 (str) -- Name of catalof matched via join=2not1.
	"""

	### Get arguments to pass to `stilts_matcher` ###
	# Args: in1, in2, out1, out2, out3, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2, overwrite # 
	# Get output filenames (`out` which is used three times in `stilts_matcher`) #

	__fn_match_1and2, __fn_match_1not2, __fn_match_2not1 = outputs.get_matched_catalog_filenames(balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory, realization=realization, tile=tile)	

	# Get input catalog filenames (`in1`, `in2`) #
	if 'coadd' == MATCH_CAT1 or (PLOT_GAUSS_APER_FLUX and PLOT_FLUX):
		input_cat1 = reformat_catalog(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=inj1_percent, realization=realization, tile=tile, mag_hdr=M_HDR1, mag_err_hdr=M_ERR_HDR1, flux_hdr=FLUX_HDR1, flux_err_hdr=FLUX_ERR_HDR1, fn_match=__fn_match_1and2, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory) 

	if 'coadd' == MATCH_CAT2 or (PLOT_GAUSS_APER_FLUX and PLOT_FLUX): 
		input_cat2 = reformat_catalog(cat_type=MATCH_CAT2, inj=INJ2, inj_percent=inj2_percent, realization=realization, tile=tile, mag_hdr=M_HDR2, mag_err_hdr=M_ERR_HDR2, flux_hdr=FLUX_HDR2, flux_err_hdr=FLUX_ERR_HDR2, fn_match=__fn_match_1and2, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory) 

	if 'coadd' not in (MATCH_CAT1, MATCH_CAT2) and PLOT_FLUX is False:
		input_cat1 = get_catalog_filename(cat_type=MATCH_CAT1, inj=inj1, inj_percent=inj1_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
		input_cat2 = get_catalog_filename(cat_type=MATCH_CAT2, inj=inj2, inj_percent=inj2_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)


	if VERBOSE_ING: print 'Matching or already matched: \n {}\n {}'.format(input_cat1, input_cat2)


	'''
	### Rewrite `match_type` for completeness plots ###
	if PLOT_COMPLETENESS: 
		if RUN_TYPE is not None:
			__class1 = get_catalog_headers.get_class(cat_type=MATCH_CAT1, cat_type_pair=MATCH_CAT2, inj_percent=inj1_percent, inj=inj1, suf='')
		if RUN_TYPE is None:
			__class1 = get_catalog_headers.get_class(cat_type=MATCH_CAT1, cat_type_pair=MATCH_CAT2, inj_percent=inj1_percent, inj=inj1, suf='_1')

		__class2 = get_catalog_headers.get_class(cat_type=MATCH_CAT2, cat_type_pair=MATCH_CAT1, inj_percent=inj2_percent, inj=inj2, suf='_2')

		__title1, __title2 = __class1.title_piece, __class2.title_piece,
		__match_type = get_match_type(title_piece1=__title1, title_piece2=__title2)

	if PLOT_COMPLETENESS is False: 
		__match_type = MATCH_TYPE
	'''


	# Overwrite matched catalogs if one already exists? # 
	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `matcher()`. Press enter to proceed, control+c to stop')

	# Check if matched catalog constains Gaussian aperture measurements #
	# If it does not, must re-match catalogs #
	try:
		df = pd.read_csv(__fn_match_1and2)
		# Can also do `GAUSS_APER_FLUX_GRIZ_HDR+'_2'`. Both exist simultaneously #
		temp = df[GAUSS_APER_FLUX_GRIZ_HDR+'_1']
	except:
		__overwrite = True
		if VERBOSE_ING: print 'Forcing overwrite to add Gaussian aperture flux measurements...'

	# Check existence #
	if os.path.isfile(__fn_match_1and2) is False or __overwrite:

		### Matching done in `stilts_matcher` ###
		# !!!!! #
		subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/stilts_matcher', input_cat1, input_cat2, __fn_match_1and2, __fn_match_1not2, __fn_match_2not1, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2])


	print ' ----->', __fn_match_1and2, '\n'


	# Delete files created for `ngmixer` once their contents are in the matched catalog #
	if PLOT_GAUSS_APER_FLUX:
		if os.path.isfile(__fn_match_1and2) is False or __overwrite:
			if VERBOSE_ING: print ' Deleting {}\n {}'.format(input_cat1, input_cat2)
			os.remove(input_cat1); os.remove(input_cat2) 


	return __fn_match_1and2, __fn_match_1not2, __fn_match_2not1, input_cat1, input_cat2





def get_coadd_catalog_observables(fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z, mag_hdr, mag_err_hdr, flux_hdr, flux_err_hdr):
	"""Creates a list of observables of form '({obs}_g, {obs}_r, {obs}_i, {obs}_z)' from four catalogs. Solely for use with coadd catalogs.
	Observables are: magnitude and magnitude error. if `PLOT_FLUX` then flux and flux error are included.

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

	if VERBOSE_ING: print 'Getting g-, r-, i-, z-band magnitudes and magnitude errors for coadd catalog...'

	# Files have not yet been matched, thus do not have 'hdr' --> 'hdr_1'. We know these constants # 
	mag_hdr = mag_hdr[:-2]
        mag_err_hdr = mag_err_hdr[:-2]
	flux_hdr = flux_hdr[:-2]
	flux_err_hdr = flux_err_hdr[:-2]

	# Open FITS files and read data #
	data_g = fits.open(fn_coadd_cat_g)[1].data
	data_r = fits.open(fn_coadd_cat_r)[1].data
        data_i = fits.open(fn_coadd_cat_i)[1].data
        data_z = fits.open(fn_coadd_cat_z)[1].data

	# Magnitude, magnitude error #
	mag_g, mag_err_g = data_g[mag_hdr], data_g[mag_err_hdr]
	mag_r, mag_err_r = data_r[mag_hdr], data_r[mag_err_hdr]
        mag_i, mag_err_i = data_i[mag_hdr], data_i[mag_err_hdr]
        mag_z, mag_err_z = data_z[mag_hdr], data_z[mag_err_hdr]

	if PLOT_FLUX:
		# Flux, flux error #
		flux_g, flux_err_g = data_g[flux_hdr], data_g[flux_err_hdr]
		flux_r, flux_err_r  = data_r[flux_hdr], data_r[flux_err_hdr]
                flux_i, flux_err_i = data_i[flux_hdr], data_i[flux_err_hdr]
                flux_z, flux_err_z  = data_z[flux_hdr], data_z[flux_err_hdr]

	__coadd_mag_griz, __coadd_mag_err_griz, __coadd_flux_griz, __coadd_flux_err_griz = np.empty(len(mag_g), dtype=str), np.empty(len(mag_g), dtype=str), np.empty(len(mag_g), dtype=str), np.empty(len(mag_g), dtype=str) 


	for i in np.arange(0, len(mag_g)):
		__coadd_mag_griz[i] = '({}, {}, {}, {})'.format(mag_g[i], mag_r[i], mag_i[i], mag_z[i])
		__coadd_mag_err_griz[i] = '({}, {}, {}, {})'.format(mag_err_g[i], mag_err_r[i], mag_err_i[i], mag_err_z[i])

		if PLOT_FLUX:
			__coadd_flux_griz[i] = '({}, {}, {}, {})'.format(flux_g[i], flux_r[i], flux_i[i], flux_z[i])
			__coadd_flux_err_griz[i] = '({}, {}, {}, {})'.format(flux_err_g[i], flux_err_r[i], flux_err_i[i], flux_err_z[i])


	return __coadd_mag_griz, __coadd_mag_err_griz, __coadd_flux_griz, __coadd_flux_err_griz




def get_coadd_catalog_for_matcher(cat_type, inj_percent, inj, realization, mag_hdr, mag_err_hdr, tile, flux_hdr, flux_err_hdr, output_directory, balrog_run):
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

	catDir = outputs.get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)

	__fn_coadd_for_matcher = os.path.join(catDir, '{}_i_cat_combo.fits'.format(tile))

	# Check if new coadd catalog has already been created #
	if os.path.isfile(__fn_coadd_for_matcher):
		print 'New coadd catalog already exists ...\n'


	if os.path.isfile(__fn_coadd_for_matcher) is False or __overwrite:	
		print 'Adding a column to i-band coadd catalog. Will take a moment ...\n'

		# Get list of filenames #
		fn_coadd_griz = []
		for b in ALL_BANDS:
			fn_coadd_griz.append(get_catalog_filename(cat_type=cat_type, inj=inj, inj_percent=inj_percent, realization=realization, tile=tile, band=b, base_path_to_catalogs=BASE_PATH_TO_CATS, balrog_run=BALROG_RUN))
		fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z = fn_coadd_griz

		# Get coadd magnitude (mag_c) and magnitude error to be of form '(m_g, m_r, m_i, m_z)'. Recall that this is a string #
		coadd_mags, coadd_mag_errs, coadd_fluxes, coadd_flux_errs = get_coadd_catalog_observables(fn_coadd_cat_g=fn_coadd_cat_g, fn_coadd_cat_r=fn_coadd_cat_r, fn_coadd_cat_i=fn_coadd_cat_i, fn_coadd_cat_z=fn_coadd_cat_z, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_err_hdr=flux_err_hdr)
 
	       # Create new table #
		coadd_mag_col = Column(coadd_mags, name=NEW_COADD_MAG_GRIZ_HDR)
		coadd_mag_err_col = Column(coadd_mag_errs, name=NEW_COADD_MAG_ERR_GRIZ_HDR)
		coadd_flux_col = Column(coadd_fluxes, name=NEW_COADD_FLUX_GRIZ_HDR)
		coadd_flux_err_col = Column(coadd_flux_errs, name=NEW_COADD_FLUX_ERR_GRIZ_HDR)

		# Add new table to i-band coadd catalog #
		table = Table.read(fn_coadd_cat_i)
		table.add_column(coadd_mag_col, index=0)
       		table.add_column(coadd_mag_err_col, index=1)
		table.add_column(coadd_flux_col, index=2)
		table.add_column(coadd_flux_err_col, index=3)

		# Save new table as FITS #
		table.write(__fn_coadd_for_matcher, overwrite=__overwrite)

	return __fn_coadd_for_matcher 




def add_gaussian_aperture_flux_measurements_to_catalog(fn_unmatched_cat, tile, realization, cat_type, balrog_run, match_type, output_directory): 
	"""

	Parameters
	----------
	fn_fits (str)
		Complete filename directly from `base_path_to_catalogs`. 

	Returns
	-------
	"""

	if VERBOSE_ING: print 'Adding Gaussian aperture measured fluxes to {}...'.format(fn_unmatched_cat)

	### Make compatible catalog ###
	fnNgmixCat = make_ngmixer_gaussap_compatible_catalog(fn_unmatched_cat=fn_unmatched_cat, tile=tile, realization=realization, cat_type=cat_type, balrog_run=balrog_run, match_type=match_type, output_directory=output_directory)


	### Get Gaussian aperture measurements ###
	# Units -- `pixel_scale`: arcsec/pixel, `weight_fwhm`: arcsec #
	# Note this inserts column id too as first column Ex: (511139065, 0, [1187.84212219, 3835.43660733, 6112.35110035, 8653.79708668]) 
	__gap_flags_and_flux_griz = ngmixer.gaussap.get_gauss_aper_flux_cat(cat=fitsio.read(fnNgmixCat, hdu=1), model='cm', pixel_scale=0.263, weight_fwhm=2.5, verbose=False, fit=Y3_FIT)
	#TODO/FIXME add_gauss_aper_flux_cat(cat=fitsio.read(fnNgmixCat), model='cm', pixel_scale=0.263, weight_fwhm=2.5, verbose=False)


	# Separate fluxes and flags. 4 bands #
	__gap_flux_griz, __gap_flags = np.empty((len(__gap_flags_and_flux_griz), 4)), np.empty(len(__gap_flags_and_flux_griz))

	for i in np.arange(0, len(__gap_flags_and_flux_griz)):
		__gap_flags[i] = __gap_flags_and_flux_griz[i][1]
		__gap_flux_griz[i] = __gap_flags_and_flux_griz[i][2]

	### Add to catalog ###
	table = Table.read(fnNgmixCat)

	try:
		# Cannot overwrite table column #
		del table[GAUSS_APER_FLAGS_HDR]
		__gap_flags_col = Column(__gap_flags, name=GAUSS_APER_FLAGS_HDR)
	except:
		__gap_flags_col = Column(__gap_flags, name=GAUSS_APER_FLAGS_HDR)

	try:
		del table[GAUSS_APER_FLUX_GRIZ_HDR]
		__gap_flux_griz_col = Column(__gap_flux_griz, name=GAUSS_APER_FLUX_GRIZ_HDR)
	except:
		__gap_flux_griz_col = Column(__gap_flux_griz, name=GAUSS_APER_FLUX_GRIZ_HDR)

	table.add_column(__gap_flags_col, index=0)
	table.add_column(__gap_flux_griz_col, index=1)
	table.write(fnNgmixCat, overwrite=True)


	return fnNgmixCat





def make_ngmixer_gaussap_compatible_catalog(fn_unmatched_cat, tile, realization, cat_type, balrog_run, match_type, output_directory):
	"""Add 'cm_pars' column to catalog and set centroids to zero

	Parameters
	----------

	Returns
	-------
	"""

	if VERBOSE_ING: print 'Making catalog compatible with `ngmixer.gaussap.get_gauss_aper_flux_cat()`...\n'

	__overwrite = False
	if __overwrite: raw_input('`__overwrite=True` in `make_ngmixer_gaussap_compatible_catalog()`. Press enter to proceed, ctrl+c to stop.')

	__fn_ngmix_cat = outputs.get_ngmix_compatible_catalog_filename(balrog_run=balrog_run, match_type=match_type, output_directory=output_directory, realization=realization, tile=tile, fn_unmatched_cat=fn_unmatched_cat)

	# Create new catalog #
	if os.path.isfile(__fn_ngmix_cat) is False or __overwrite:

		if VERBOSE_ING: print ' Adding column to', cat_type, 'catalog...'

		table = Table.read(fn_unmatched_cat)
		cmParams = get_cm_parameters(fn_unmatched_cat=fn_unmatched_cat, cat_type=cat_type) 
		# If the catalog already has a 'cm_pars' header, delete it and create a new column with cen1=cen2=0 #

		# Cannot overwrite table columns, must delete #
		try:
			del table['cm_pars']
			__cm_params_col = Column(cmParams, name='cm_pars')
		except:
			__cm_params_col = Column(cmParams, name='cm_pars')

		table.add_column(__cm_params_col, index=0)
		table.write(__fn_ngmix_cat, overwrite=__overwrite)

	return __fn_ngmix_cat	




def get_cm_parameters(fn_unmatched_cat, cat_type):
	"""Get parameters necessary for `ngmixer.gaussap.get_gauss_aper_flux_cat()`.

	Parameters
	----------

	Returns
	-------
	""" 

	if 'y3_gold' in cat_type:
		cmParameters = get_cm_parameters_from_y3_gold_catalog(fn_unmatched_cat=fn_unmatched_cat) 

	if cat_type in ('gal_truth', 'mof', 'sof'):
		cmParameters = get_cm_parameters_from_mof_catalog(fn_unmatched_cat=fn_unmatched_cat)

	return cmParameters 




def get_cm_parameters_from_y3_gold_catalog(fn_unmatched_cat):
	"""Get parameters necessary for `ngmixer.gaussap.get_gauss_aper_flux_cat()`.

	Parameters
	----------
	fn (str)
		Must be a FITS file.

	Returns
	-------
	"""

	# `hdu=1` is flux TODO #
	__catalog = fitsio.read(fn_unmatched_cat, hdu=1)

	# These are unmatched Y3 Gold catalogs with known headers #
	# Get data #
	__cm_g1 = __catalog['{}_cm_g_1'.format(Y3_FIT).upper()]
	__cm_g2 = __catalog['{}_cm_g_2'.format(Y3_FIT).upper()]
	__cm_t = __catalog['{}_cm_T'.format(Y3_FIT).upper()]
	__cm_flux_g = __catalog['{}_cm_flux_g'.format(Y3_FIT).upper()]
	__cm_flux_r = __catalog['{}_cm_flux_r'.format(Y3_FIT).upper()]
	__cm_flux_i = __catalog['{}_cm_flux_i'.format(Y3_FIT).upper()]
	__cm_flux_z = __catalog['{}_cm_flux_z'.format(Y3_FIT).upper()]

	# There are 9 parameters in cm_params #
	__cm_params = np.empty((len(__cm_g1), 9))

	for i in np.arange(0, len(__cm_g1)):
		# cen1, cen2 #
		__cm_params[i][0], __cm_params[i][1] = 0.0, 0.0
		# g1, g2 #
		__cm_params[i][2], __cm_params[i][3]  = __cm_g1[i], __cm_g2[i] 
		# cm_T #
		__cm_params[i][4] = __cm_t[i]
		# g-, r-, i-, z-band fluxes #
		__cm_params[i][5], __cm_params[i][6], __cm_params[i][7], __cm_params[i][8] = __cm_flux_g[i], __cm_flux_r[i], __cm_flux_i[i], __cm_flux_z[i] 

	return __cm_params




def get_cm_parameters_from_mof_catalog(fn_unmatched_cat):
	"""Get parameters necessary for `ngmixer.gaussap.get_gauss_aper_flux_cat()`. For catalogs with MOF output.

	Parameters
	----------
	fn_unmatched_cat (str)
		Must be a FITS file.

	Returns
	-------
	"""

	__catalog = fitsio.read(fn_unmatched_cat, hdu=1)

	# These are unmatched ngmix catalogs with known headers so they are hardcoded, not retrieved from `HDR1` or `HDR2` constants # 
	# Get data #
	__cm_g1_g2 = __catalog['cm_g']
	__cm_t = __catalog['cm_T']
	__cm_flux_griz = __catalog['cm_flux']

	# There are 9 parameters in cm_params #
	__cm_params = np.empty((len(__cm_g1_g2), 9))

	for i in np.arange(0, len(__cm_g1_g2)):
		# cen1, cen2 #
		__cm_params[i][0], __cm_params[i][1] = 0.0, 0.0
		# g1, g2 #
		__cm_params[i][2], __cm_params[i][3]  = __cm_g1_g2[i]
		# cm_T #
		__cm_params[i][4] = __cm_t[i]
		# g-, r-, i-, z-band fluxes #
		__cm_params[i][5], __cm_params[i][6], __cm_params[i][7], __cm_params[i][8] = __cm_flux_griz[i]

	return __cm_params




def fof_matcher(realization, tile):
	"""Get catalogs to analyze. Return FOF-analysed catalogs.

	Parameters
	----------
	realization (str) 
	tile (str) 

	Returns
	-------
	__fn_ok_1and2 OR __fn_rerun_1and2 (str)
		Complete filename for catalogs of type join=1and2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned. 
	__fn_ok_1not2 OR __fn_rerun_1not2 (str)
			Complete filename for catalogs of type join=1not2. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.
	__fn_ok_2not1 OR __fn_rerun_2not1 (str)
		Complete filename for catalogs of type join=2not1. `RUN_TYPE` determines if 'ok' or 'rerun' catalog filename is returned.
	"""

	### Filenames for input catalogs used in `ms_fof_matcher` ###
	# MOF or SOF #
	#FIXME
	FOF_FIT = 'mof' or 'sof'

	smof = get_catalog_filename(cat_type=FOF_FIT, inj=False, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
	inj_smof = get_catalog_filename(cat_type=FOF_FIT, inj=True, inj_percent=inj_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

	fof = os.path.join(BASE_PATH_TO_CATS, 'y3v02', tile, 'mof', tile+'_fofslist.fits')
	inj_fof = os.path.join(BASE_PATH_TO_CATS, 'y3v02', 'balrog_images', realization, tile, 'mof', tile+'_fofslist.fits')

	fof = os.path.join(BASE_PATH_TO_CATS, 'y3v02', tile, 'sof', tile+'_fofslist.fits')
	inj_fof = os.path.join(BASE_PATH_TO_CATS, 'y3v02', 'balrog_images', realization, tile, 'sof', tile+'_fofslist.fits')

	# Coadds. Using i-band #
	coadd = get_catalog_filename(cat_type='coadd', inj=False, inj_percent=None, realization=realization, tile=tile, band='i', base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
	inj_coadd = get_catalog_filename(cat_type='coadd', inj=True, inj_percent=inj_percent, realization=realization, tile=tile, band='i', base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)


	### Filenames for outputs of fof_matcher ###
	catalogDirectory = outputs.get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)

	### Create filenames for output catalogs created by ms_fof_matcher ###
	#TODO add fn_
	fofcoadd = os.path.join(catalogDirectory, '{}_num_match_fof_coadd.csv'.format(tile))
	fofgroups = os.path.join(catalogDirectory, '{}_fofgroups.csv'.format(tile))
	inj_fofcoadd = os.path.join(catalogDirectory, '{}_{}_num_match_inj_fof_inj_coadd.csv'.format(tile, realization))
	inj_fofgroups = os.path.join(catalogDirectory, '{}_{}_inj_fofgroups.csv'.format(tile, realization))
	origfof_injfof = os.path.join(catalogDirectory, '{}_{}_inj_fofgroup_fofgroup_match1and2.csv'.format(tile, realization))
	ok = os.path.join(catalogDirectory, '{}_{}.ok'.format(tile, realization))
	rerun = os.path.join(catalogDirectory, '{}_{}.rerun'.format(tile, realization))
	ok_inj_mof = os.path.join(catalogDirectory, '{}_{}_ok_inj_mof.csv'.format(tile, realization))
	rerun_inj_mof = os.path.join(catalogDirectory, '{}_{}_rerun_inj_mof.csv'.format(tile, realization))
	ok_mof = os.path.join(catalogDirectory, '{}_{}_ok_mof.csv'.format(tile, realization))
	rerun_mof = os.path.join(catalogDirectory, '{}_{}_rerun_mof.csv'.format(tile, realization))

	__fn_ok_1and2 = os.path.join(catalogDirectory, '{}_{}_ok_inj_mof_ok_mof_match1and2.csv'.format(tile, realization))
	__fn_ok_1not2 = os.path.join(catalogDirectory, '{}_{}_ok_inj_mof_ok_mof_match1not2.csv'.format(tile, realization))
	__fn_ok_2not1 = os.path.join(catalogDirectory, '{}_{}_ok_inj_mof_ok_mof_match2not1.csv'.format(tile, realization))

	__fn_rerun_1and2 = os.path.join(catalogDirectory, '{}_{}_rerun_inj_mof_rerun_mof_match1and2.csv'.format(tile, realization))
	__fn_rerun_1not2 = os.path.join(catalogDirectory, '{}_{}_rerun_inj_mof_rerun_mof_match1not2.csv'.format(tile, realization))
	__fn_rerun_2not1 = os.path.join(catalogDirectory, '{}_{}_rerun_inj_mof_rerun_mof_match2not1.csv'.format(tile, realization))

	# Output directory for files made in par.py #
	parpy_catalogDirectory = os.path.join(catalogDirectory, '{}_{}'.format(tile, realization))


	# WARNING: May need to overwrite if matching was interupted #
	__overwrite = False 
	if __overwrite: raw_input('`__overwrite=True` in `fof_matcher()`. Press enter to procees and ctrl+c to stop.')

	### Check file existence of last file made in fof_matcher ###
	if os.path.isfile(__fn_rerun_1and2) is False or __overwrite:

		### Run fof_matcher ###
		subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_fof_matcher', fof, inj_fof, smof, inj_smof, coadd, inj_coadd, parpy_catalogDirectory, fofcoadd, fofgroups, inj_fofcoadd, inj_fofgroups, origfof_injfof, ok, rerun, ok_inj_mof, rerun_inj_mof, ok_mof, rerun_mof, __fn_ok_1and2, __fn_rerun_1and2, __fn_ok_1not2, __fn_rerun_1not2, __fn_ok_2not1, __fn_ok_2not1])


	if RUN_TYPE == 'ok':
		return __fn_ok_1and2, __fn_ok_1not2, __fn_ok_2not1
	if RUN_TYPE == 'rerun':
		return __fn_rerun_1and2, __fn_rerun_1not2, __fn_rerun_2not1

