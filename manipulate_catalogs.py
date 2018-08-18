"""
Manipulate catalogs in the some of the following ways:
Retrieve catalog filenames.
Match catalogs with a call to `stilts_matcher`.
If necessary, reformat the catalog(s) (prior to matching in the case of... and after matching in the case of...).
    Reformat e.g.:
Add a column ('cm_pars') necessary to perform Gaussian aperture measured flux (via `ngmixer`)
Add columns to coadd catalog.
    ...
Read pandas DataFrame of the matched catalogs.

Comments are ABOVE the code they refer to.
"""

import csv
import fitsio
import ngmixer
import numpy as np
import os
import pandas as pd
import subprocess
import sys
import tarfile
from astropy.io import fits
from astropy.table import Table, Column

# BalVal #
from set_constants import * 
from catalog_headers import *
import outputs



#TODO need output dir for coadds
#TODO get rid of `band`
def get_catalog_filename(cat_type, inj, inj_percent, realization, tile, band, base_path_to_catalogs, balrog_run, output_directory=None, match_type=MATCH_TYPE):
    """Get the complete filename of a catalog that will be matched by `stilts_matcher`.
    Note that coadd catalogs must be manipulated BEFORE matching, so a reformatted coadd catalog will be returned if `cat_type='coadd'`.
    Note that different Balrog runs (`balrog_run`) may have different directory structures and filename structures.
    An (incomplete) log is at https://docs.google.com/document/d/1hb3OlmU02jS4KkvCxvw80zr1fd15LbgzOK7tps1uqBM/edit?usp=sharing

    Parameters
    ----------
    cat_type
        Catalog type.
        Set by `MATCH_CAT1` or `MATCH_CAT2` (see Table of Constants in README.md for allowed values).

    inj_percent (int)

    inj (bool)

    realization (str)

    tile

    band (str)
        Only used with coadd catalogs. Ignored if `cat_type` is not 'coadd'.

    base_path_to_catalogs

    balrog_run

    output_directory
        Needed for zipped coadds. File will be unzipped in `output_directory`.

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
        if os.path.isfile(os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', '{}_i_cat.fits'.format(tile))):
            __fn_coadd_griz = []
            for b in ALL_BANDS:
                __fn_coadd_griz.append(os.path.join(base_path_to_catalogs, 'y3v02', tile, 'coadd', '{}_{}_cat.fits'.format(tile, b)))

            __fn_cat = get_coadd_catalog_for_matcher(fn_coadd_cat_g=__fn_coadd_griz[0], fn_coadd_cat_r=__fn_coadd_griz[1], fn_coadd_cat_i=__fn_coadd_griz[2], fn_coadd_cat_z=__fn_coadd_griz[3], tile=tile, realization=realization, output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)

        #FIXME
        if os.path.isfile(os.path.join(base_path_to_catalogs, 'base_cat.tgz')):
            fnCoaddG, fnCoaddR, fnCoaddI, fnCoaddZ = unzip_coadd_catalog(fn_zip=os.path.join(base_path_to_catalogs, 'base_cat.tgz'), tile=tile, realization=realization, output_directory=output_directory, balrog_run=balrog_run)
            __fn_cat = get_coadd_catalog_for_matcher(fn_coadd_cat_g=fnCoaddG, fn_coadd_cat_r=fnCoaddR, fn_coadd_cat_i=fnCoaddI, fn_coadd_cat_z=fnCoaddZ,  tile=tile, realization=realization, output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)


        if os.path.isfile(os.path.join(base_path_to_catalogs, tile, 'base_cat.tgz')):
            fnCoaddG, fnCoaddR, fnCoaddI, fnCoaddZ = unzip_coadd_catalog(fn_zip=os.path.join(base_path_to_catalogs, tile, 'base_cat.tgz'), tile=tile, realization=realization, output_directory=output_directory, balrog_run=balrog_run)
            __fn_cat = get_coadd_catalog_for_matcher(fn_coadd_cat_g=fnCoaddG, fn_coadd_cat_r=fnCoaddR, fn_coadd_cat_i=fnCoaddI, fn_coadd_cat_z=fnCoaddZ,  tile=tile, realization=realization, output_directory=output_directory, balrog_run=balrog_run, match_type=match_type)



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


    try:
        os.path.isfile(__fn_cat)
    except UnboundLocalError:
        sys.exit('UnboundLocalError. `__fn_cat` not defined. Are you sure the tile(s) and realization(s) entered at command line exist within the basepath given at command line?\nThis is one of many possible errors!')

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




def unzip_coadd_catalog(fn_zip, tile, realization, output_directory, balrog_run):
    """Handle zipped coadd catalogs.
    The coadds will be extracted in `OUTPUT_DIRECTORY`. 

    Parameters
    ----------

    Returns
    -------
    
        Complete filenames for each band.
    """

    if VERBOSE_ING: print 'Unzipping coadd catalogs...\n'
    #TODO delete extracted coadds after matching

    __overwrite = False
    if __overwrite: raw_input('`__overwrite=True` in `unzip_coadd_catalog()`. Press enter to proceed, ctrl+c to stop')


    matchedCatalogDirectory = outputs.get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)

    # Name for new catalog #
    __fn_coadd_cat_for_matcher = os.path.join(matchedCatalogDirectory, '{}_i_coadd.fits'.format(tile))

    if os.path.isfile(__fn_coadd_cat_for_matcher) and __overwrite is False:
        if VERBOSE_ING: print 'Reformatted coadd catalog already exists. Not overwriting...'

    if os.path.isfile(__fn_coadd_cat_for_matcher) is False or __overwrite:
        if VERBOSE_ING: print 'Reformatting coadd catalog...'

        # Base catalogs are zipped #
        __tar = tarfile.open(fn_zip, 'r')

    # Base catalogs are zipped #
    __fn_coadd_griz = []

    # Extract somewhere in `OUTPUT_DIRECTORY` #
    __tar.extractall(path=matchedCatalogDirectory)
    for b in ALL_BANDS:
        __fn_b = '{}_{}_cat.fits'.format(tile, b)
        __fn_coadd_griz.append(os.path.join(matchedCatalogDirectory, __fn_b))


    return __fn_coadd_griz 



def reformat_catalog(cat_type, inj, inj_percent, realization, tile, fn_match, mag_hdr, mag_err_hdr, flux_hdr, flux_err_hdr, base_path_to_catalogs, balrog_run, match_type, output_directory):
    """Reformat catalog PRIOR to matching

    Parameters
    ----------

    Returns
    -------

    """

    __overwrite = False
    if __overwrite: raw_input('To rewrite press Enter')

    #if cat_type == 'coadd':
        #__fn_cat = get_coadd_catalog_for_matcher(cat_type=cat_type, inj_percent=inj_percent,  realization=realization, tile=tile, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_err_hdr=flux_err_hdr) 


    ### Make catalog compatible for Gaussian aperture catalog ###
    if PLOT_GAUSS_APER_FLUX and PLOT_FLUX:
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
    # TODO add coadd to get_catalog_file_name(). get_catalog_filename_for_matcher add to docstring
    if PLOT_GAUSS_APER_FLUX and PLOT_FLUX:
        input_cat1 = reformat_catalog(cat_type=MATCH_CAT1, inj=INJ1, inj_percent=inj1_percent, realization=realization, tile=tile, mag_hdr=M_HDR1, mag_err_hdr=M_ERR_HDR1, flux_hdr=FLUX_HDR1, flux_err_hdr=FLUX_ERR_HDR1, fn_match=__fn_match_1and2, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory) 
        input_cat2 = reformat_catalog(cat_type=MATCH_CAT2, inj=INJ2, inj_percent=inj2_percent, realization=realization, tile=tile, mag_hdr=M_HDR2, mag_err_hdr=M_ERR_HDR2, flux_hdr=FLUX_HDR2, flux_err_hdr=FLUX_ERR_HDR2, fn_match=__fn_match_1and2, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory) 

    if PLOT_FLUX is False and PLOT_GAUSS_APER_FLUX is False:
        input_cat1 = get_catalog_filename(cat_type=MATCH_CAT1, inj=inj1, inj_percent=inj1_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, output_directory=output_directory)
        input_cat2 = get_catalog_filename(cat_type=MATCH_CAT2, inj=inj2, inj_percent=inj2_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run, output_directory=output_directory)


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

    # Check if matched catalog contains Gaussian aperture measurements #
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
        subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/stilts_matcher', input_cat1, input_cat2, __fn_match_1and2, __fn_match_1not2, __fn_match_2not1, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2])


    print ' ----->', __fn_match_1and2, '\n'


    # Delete catalogs created for `ngmixer` once their contents are in the matched catalog #
    # Catalogs created for `ngmixer` are saved in `OUTPUT_DIRECTORY` where user has permission to edit files #
    # There is no need to have `MATCH_CAT1` and `MATCH_CAT2` in two places # 
    if PLOT_GAUSS_APER_FLUX:
        if os.path.isfile(__fn_match_1and2) is False or __overwrite:
            try:
                os.remove(input_cat1); os.remove(input_cat2) 
                if VERBOSE_ING: print ' Deleting {}\n {}'.format(input_cat1, input_cat2)
            except:
                pass

    return __fn_match_1and2, __fn_match_1not2, __fn_match_2not1, input_cat1, input_cat2





def reformat_coadd_catalog_observables_and_flags(fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z): 
    """Creates a list of observables of form '({obs}_g, {obs}_r, {obs}_i, {obs}_z)' from four catalogs. Solely for use with coadd catalogs.
    Observables are: magnitude and magnitude error. If `PLOT_FLUX=True` then flux and flux error are included.
    Also flags

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
    Note that some returns might be empty numpy arrays.
    """

    if VERBOSE_ING: print 'Getting g-, r-, i-, z-band magnitudes and magnitude errors for coadd catalog...'

    # Flags #
    flag1_hdr = 'FLAGS'
    flag2_hdr = 'IMAFLAGS_ISO'

    # Files have not yet been matched, thus do not have 'hdr' --> 'hdr_1'. Headers are known # 
    mag_hdr = 'MAG_AUTO'
    mag_err_hdr = 'MAGERR_AUTO' 
    flux_hdr = 'FLUX_AUTO' 
    flux_err_hdr = 'FLUXERR_AUTO' 

    # Open FITS files and read data #
    data_g, data_r, data_i, data_z = fitsio.read(fn_coadd_cat_g, hdu=1), fitsio.read(fn_coadd_cat_r, hdu=1), fitsio.read(fn_coadd_cat_i, hdu=1), fitsio.read(fn_coadd_cat_z, hdu=1)

    flag1_g, flag1_r, flag1_i, flag1_z = data_g[flag1_hdr], data_r[flag1_hdr], data_i[flag1_hdr], data_z[flag1_hdr]
    flag2_g, flag2_r, flag2_i, flag2_z = data_g[flag2_hdr], data_r[flag2_hdr], data_i[flag2_hdr], data_z[flag2_hdr]

    # Magnitude, magnitude error #
    mag_g, mag_err_g = data_g[mag_hdr], data_g[mag_err_hdr]
    mag_r, mag_err_r = data_r[mag_hdr], data_r[mag_err_hdr]
    mag_i, mag_err_i = data_i[mag_hdr], data_i[mag_err_hdr]
    mag_z, mag_err_z = data_z[mag_hdr], data_z[mag_err_hdr]

    if PLOT_FLUX:
        # Flux, flux error #
        flux_g, flux_err_g = data_g[flux_hdr], data_g[flux_err_hdr]
        flux_r, flux_err_r = data_r[flux_hdr], data_r[flux_err_hdr]
        flux_i, flux_err_i = data_i[flux_hdr], data_i[flux_err_hdr]
        flux_z, flux_err_z = data_z[flux_hdr], data_z[flux_err_hdr]


    __coadd_mag_griz, __coadd_mag_err_griz, __coadd_flux_griz, __coadd_flux_err_griz = np.empty(len(mag_g), dtype='S100'), np.empty(len(mag_g), dtype='S100'), np.empty(len(mag_g), dtype='S100'), np.empty(len(mag_g), dtype='S100') 
    __coadd_flag1_griz, __coadd_flag2_griz = np.empty(len(mag_g), dtype='S100'), np.empty(len(mag_g), dtype='S100')


    for i in np.arange(0, len(mag_g)):
        __coadd_flag1_griz[i] = '({}, {}, {}, {})'.format(flag1_g[i], flag1_r[i], flag1_i[i], flag1_z[i])
        __coadd_flag2_griz[i] = '({}, {}, {}, {})'.format(flag2_g[i], flag2_r[i], flag2_i[i], flag2_z[i])
        __coadd_mag_griz[i] = '({}, {}, {}, {})'.format(mag_g[i], mag_r[i], mag_i[i], mag_z[i])
        __coadd_mag_err_griz[i] = '({}, {}, {}, {})'.format(mag_err_g[i], mag_err_r[i], mag_err_i[i], mag_err_z[i])

        if PLOT_FLUX:
            __coadd_flux_griz[i] = '({}, {}, {}, {})'.format(flux_g[i], flux_r[i], flux_i[i], flux_z[i])
            __coadd_flux_err_griz[i] = '({}, {}, {}, {})'.format(flux_err_g[i], flux_err_r[i], flux_err_i[i], flux_err_z[i])

    print __coadd_mag_err_griz 
    return __coadd_mag_griz, __coadd_mag_err_griz, __coadd_flux_griz, __coadd_flux_err_griz, __coadd_flag1_griz, __coadd_flag2_griz



def get_coadd_catalog_for_matcher(fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z, tile, realization, output_directory, balrog_run, match_type): 
    """Reformat coadd catalog that will be used by `stilts_matcher`.
    The i-band coadd catalog will have the g-, r-, and z-band observables (magnitude, magnitude error, and possible flux and flux error if `PLOT_FLUX=True`) and flags appended to it.
    Thus, catalogs are matched using i-band RA and Dec.

    Parameters
    ----------

    Returns
    -------
    """

    __overwrite = False
    if __overwrite: raw_input('`__overwrite=True` in `get_coadd_catalog_for_matcher()`. Press enter to proceed, ctrl+c to exit.')

    catalogDirectory = outputs.get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)

    __fn_coadd_for_matcher = os.path.join(catalogDirectory, '{}_i_cat_combo.fits'.format(tile))

    # Check if new coadd catalog has already been created #
    if os.path.isfile(__fn_coadd_for_matcher):
        print 'New coadd catalog already exists...'


    if os.path.isfile(__fn_coadd_for_matcher) is False or __overwrite:
        print 'Adding a column to i-band coadd catalog...'

        # Reformat coadd magnitude and magnitude error to be of form '(mag_g, mag_r, mag_i, mag_z)'. Recall that this is a string #
        coadd_mags, coadd_mag_errs, coadd_fluxes, coadd_flux_errs, coadd_flags1, coadd_flags2 = reformat_coadd_catalog_observables_and_flags(fn_coadd_cat_g=fn_coadd_cat_g, fn_coadd_cat_r=fn_coadd_cat_r, fn_coadd_cat_i=fn_coadd_cat_i, fn_coadd_cat_z=fn_coadd_cat_z)


       # Create columns for new table #
        coadd_flag1_col = Column(coadd_flags1, name=NEW_COADD_FLAG_GRIZ_HDR)
        coadd_flag2_col = Column(coadd_flags2, name=NEW_COADD_IMAFLAG_ISO_GRIZ_HDR)
        coadd_mag_col = Column(coadd_mags, name=NEW_COADD_MAG_GRIZ_HDR)
        coadd_mag_err_col = Column(coadd_mag_errs, name=NEW_COADD_MAG_ERR_GRIZ_HDR)
        if PLOT_FLUX:
            coadd_flux_col = Column(coadd_fluxes, name=NEW_COADD_FLUX_GRIZ_HDR)
            coadd_flux_err_col = Column(coadd_flux_errs, name=NEW_COADD_FLUX_ERR_GRIZ_HDR)

        # Add columns to i-band coadd catalog #
        table = Table.read(fn_coadd_cat_i)
        table.add_column(coadd_mag_col, index=0)
        table.add_column(coadd_mag_err_col, index=1)
        table.add_column(coadd_flag1_col, index=2)
        table.add_column(coadd_flag2_col, index=3)
        if PLOT_FLUX:
            table.add_column(coadd_flux_col, index=4)
            table.add_column(coadd_flux_err_col, index=5)


        # Save new table as FITS file #
        table.write(__fn_coadd_for_matcher, overwrite=__overwrite)

    return __fn_coadd_for_matcher






#FIXME remove
#def get_coadd_catalog_for_matcher1(cat_type, inj_percent, inj, realization, mag_hdr, mag_err_hdr, tile, flux_hdr, flux_err_hdr, output_directory, balrog_run, base_path_to_catalogs):
    """Make FITS file that includes a column of form '(mag_g, mag_r, mag_i, mag_z)' .
    Column will be added to the i-band coadd catalog. 
    New FITS file will be matched via `stilts_matcher`. Thus, catalogs are matched via i-band RA and Dec. 
    New FITS file saved in `OUTPUT_DIRECTORY`/outputs/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{realization}/catalog_compare/` 

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
    '''
    __overwrite = False 
    if __overwrite: raw_input('`__overwrite=True` in `get_coadd_catalog_for_matcher()`. Press enter to proceed, ctrl+c to exit.') 

    catalogDirectory = outputs.get_directory(tile=tile, realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)

    __fn_coadd_for_matcher = os.path.join(catalogDirectory, '{}_i_cat_combo.fits'.format(tile))

    # Check if new coadd catalog has already been created #
    if os.path.isfile(__fn_coadd_for_matcher):
        print 'New coadd catalog already exists...'


    if os.path.isfile(__fn_coadd_for_matcher) is False or __overwrite:  
        print 'Adding a column to i-band coadd catalog...'

        # There is a separate catalog for each band. Collect filenames #
        fn_coadd_griz = []
        for b in ALL_BANDS:
            # get_catalog_filename()` handles zipped coadd catalogs #
            fn_coadd_griz.append(get_catalog_filename(cat_type=cat_type, inj=inj, inj_percent=inj_percent, realization=realization, tile=tile, band=b, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run))
        fn_coadd_cat_g, fn_coadd_cat_r, fn_coadd_cat_i, fn_coadd_cat_z = fn_coadd_griz

        # Reformat coadd magnitude and magnitude error to be of form '(mag_g, mag_r, mag_i, mag_z)'. Recall that this is a string #
        coadd_mags, coadd_mag_errs, coadd_fluxes, coadd_flux_errs, coadd_flags1, coadd_flags2 = get_reformatted_coadd_catalog_observables(fn_coadd_cat_g=fn_coadd_cat_g, fn_coadd_cat_r=fn_coadd_cat_r, fn_coadd_cat_i=fn_coadd_cat_i, fn_coadd_cat_z=fn_coadd_cat_z, mag_hdr=mag_hdr, mag_err_hdr=mag_err_hdr, flux_hdr=flux_hdr, flux_err_hdr=flux_err_hdr) 
 
           # Create columns for new table #
        coadd_flag1_col = Column(__base_flag1_griz, name=NEW_COADD_FLAG_GRIZ_HDR)
        coadd_flag2_col = Column(__base_flag2_griz, name=NEW_COADD_IMAFLAG_ISO_GRIZ_HDR)
        coadd_mag_col = Column(coadd_mags, name=NEW_COADD_MAG_GRIZ_HDR)
        coadd_mag_err_col = Column(coadd_mag_errs, name=NEW_COADD_MAG_ERR_GRIZ_HDR)
        coadd_flux_col = Column(coadd_fluxes, name=NEW_COADD_FLUX_GRIZ_HDR)
        coadd_flux_err_col = Column(coadd_flux_errs, name=NEW_COADD_FLUX_ERR_GRIZ_HDR)

        # Add columns to i-band coadd catalog #
        table = Table.read(fn_coadd_cat_i)
        table.add_column(coadd_mag_col, index=0)
        table.add_column(coadd_mag_err_col, index=1)
        table.add_column(coadd_flux_col, index=2)
        table.add_column(coadd_flux_err_col, index=3)
        table.add_column(coadd_flag1_col, index=4)
        table.add_column(coadd_flag2_col, index=5)

        # Save new table as FITS file #
        table.write(__fn_coadd_for_matcher, overwrite=__overwrite)

    return __fn_coadd_for_matcher 
    '''



def add_gaussian_aperture_flux_measurements_to_catalog(fn_unmatched_cat, tile, realization, cat_type, balrog_run, match_type, output_directory): 
    """Add Gaussian aperture measured fluxes to catalog.

    Parameters
    ----------
    fn_unmatched_cat (str)
        Complete filename of unmatched catalog.
        Set by `MATCH_CAT1` or `MATCH_CAT2`. 

    Returns
    -------
    fnNgmixCat (str)
        Complete filename for catalog that is compatible with `ngmixer.gaussap`.
    """

    if VERBOSE_ING: print 'Adding Gaussian aperture measured fluxes to {}...'.format(fn_unmatched_cat)

    ### Get `ngmixer.gaussap` compatible catalog ###
    fnNgmixCat = make_ngmixer_gaussap_compatible_catalog(fn_unmatched_cat=fn_unmatched_cat, tile=tile, realization=realization, cat_type=cat_type, balrog_run=balrog_run, match_type=match_type, output_directory=output_directory)


    ### Get Gaussian aperture measurements ###
    # Units -- `pixel_scale`: arcsec/pixel, `weight_fwhm`: arcsec #
    # Note this inserts column id as first column. Ex: (511139065, 0, [1187.84212219, 3835.43660733, 6112.35110035, 8653.79708668]) #
    __gap_flags_and_flux_griz = ngmixer.gaussap.get_gauss_aper_flux_cat(cat=fitsio.read(fnNgmixCat, hdu=1), model='cm', pixel_scale=0.263, weight_fwhm=2.5, verbose=False, fit=Y3_FIT)
    #TODO/FIXME add_gauss_aper_flux_cat(cat=fitsio.read(fnNgmixCat), model='cm', pixel_scale=0.263, weight_fwhm=2.5, verbose=False)


    # Separate fluxes and flags from `__gap_flags_and_flux_griz`. `__gap_flags_and_flux_griz` returns fluxes for all griz bands #
    __gap_flux_griz, __gap_flags = np.empty((len(__gap_flags_and_flux_griz), 4)), np.empty(len(__gap_flags_and_flux_griz))

    for i in np.arange(0, len(__gap_flags_and_flux_griz)):
        __gap_flags[i] = __gap_flags_and_flux_griz[i][1]
        __gap_flux_griz[i] = __gap_flags_and_flux_griz[i][2]


    ### Add to catalog ###
    table = Table.read(fnNgmixCat)

    # Flags #
    try:
        # Cannot overwrite column of a table, msut delete if it already exists #
        del table[GAUSS_APER_FLAGS_HDR]
        __gap_flags_col = Column(__gap_flags, name=GAUSS_APER_FLAGS_HDR)
    except:
        __gap_flags_col = Column(__gap_flags, name=GAUSS_APER_FLAGS_HDR)

    # Fluxes #
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
    """Make `ngmixer.gaussap` compatible catalog by adding 'cm_pars' column to catalog.
    If input catalog (`fn_unmatched_cat`) already has a 'cm_pars' header, delete it and create a new column with centroids (`cen1`, `cen2`) set to zero.

    Parameters
    ----------

    Returns
    -------
    __fn_ngmix_cat (str)
        Complete filename to `ngmixer.gaussap` compatible catalog.
    """

    if VERBOSE_ING: print 'Making catalog compatible with `ngmixer.gaussap.get_gauss_aper_flux_cat()`...'

    __overwrite = False
    if __overwrite: raw_input('`__overwrite=True` in `make_ngmixer_gaussap_compatible_catalog()`. Press enter to proceed, ctrl+c to stop.')

    __fn_ngmix_cat = outputs.get_ngmix_compatible_catalog_filename(balrog_run=balrog_run, match_type=match_type, output_directory=output_directory, realization=realization, tile=tile, fn_unmatched_cat=fn_unmatched_cat)

    # Create new catalog #
    if os.path.isfile(__fn_ngmix_cat) is False or __overwrite:

        if VERBOSE_ING: print ' Adding column to {} catalog...'.format(cat_type)

        table = Table.read(fn_unmatched_cat)
        cmParams = get_cm_parameters(fn_unmatched_cat=fn_unmatched_cat, cat_type=cat_type) 

        # Cannot overwrite columns of table, must delete them if they exist #
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
    'cm_pars' are: (cen1, cen2, cm_g_1, cm_g_2, cm_T, flux_g, flux_r, flux_i, flux_x)

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
    """Get parameters necessary for `ngmixer.gaussap.get_gauss_aper_flux_cat()` from Y3 Gold catalogs.
    'cm_pars' are: (cen1, cen2, cm_g_1, cm_g_2, cm_T, flux_g, flux_r, flux_i, flux_x)

    Parameters
    ----------
    fn_unmatched_cat (str)
        Must be a FITS file.

    Returns
    -------
    __cm_params (nd_array)

    """

    # `hdu=1` is flux #
    __catalog = fitsio.read(fn_unmatched_cat, hdu=1)

    # These are unmatched Y3 Gold catalogs with known headers, no need to use headers established in `catalog_headers.py` #
    # Get data #
    __cm_g1 = __catalog['{}_cm_g_1'.format(Y3_FIT).upper()]
    __cm_g2 = __catalog['{}_cm_g_2'.format(Y3_FIT).upper()]
    __cm_t = __catalog['{}_cm_T'.format(Y3_FIT).upper()]
    __cm_flux_g = __catalog['{}_cm_flux_g'.format(Y3_FIT).upper()]
    __cm_flux_r = __catalog['{}_cm_flux_r'.format(Y3_FIT).upper()]
    __cm_flux_i = __catalog['{}_cm_flux_i'.format(Y3_FIT).upper()]
    __cm_flux_z = __catalog['{}_cm_flux_z'.format(Y3_FIT).upper()]

    # There are 9 parameters in 'cm_params' #
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
    """Get parameters necessary for `ngmixer.gaussap.get_gauss_aper_flux_cat()` from catalogs with MOF outputs (some Balrog truth catalogs are MOF catalogs).
    'cm_pars' are: (cen1, cen2, cm_g_1, cm_g_2, cm_T, flux_g, flux_r, flux_i, flux_x)

    Parameters
    ----------
    fn_unmatched_cat (str)
        Must be a FITS file.

    Returns
    -------
    __cm_params (ndarray)

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




def fof_matcher(realization, tile, base_path_to_catalogs, balrog_run):
    """Get filenames of catalogs to analyze.
    Analyse their FOF groups.

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

    ### Filenames for input catalogs used in `fof_stilts_matcher` ###

    smof = get_catalog_filename(cat_type=FOF_FIT, inj=False, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)
    inj_smof = get_catalog_filename(cat_type=FOF_FIT, inj=True, inj_percent=inj_percent, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

    #TODO update how FOF lists are stored
    fof = os.path.join(base_path_to_catalogs, 'y3v02', tile, FOF_FIT.lower(), '{}_fofslist.fits'.format(tile))
    inj_fof = os.path.join(base_path_to_catalogs, 'y3v02', 'balrog_images', realization, tile, FOF_FIT.lower(), '{}_fofslist.fits'.format(tile))

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

        ### Run `fof_stilts_matcher` ###
        subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/fof_stilts_matcher', fof, inj_fof, smof, inj_smof, coadd, inj_coadd, parpy_catalogDirectory, fofcoadd, fofgroups, inj_fofcoadd, inj_fofgroups, origfof_injfof, ok, rerun, ok_inj_mof, rerun_inj_mof, ok_mof, rerun_mof, __fn_ok_1and2, __fn_rerun_1and2, __fn_ok_1not2, __fn_rerun_1not2, __fn_ok_2not1, __fn_ok_2not1])


    if RUN_TYPE == 'ok':
        return __fn_ok_1and2, __fn_ok_1not2, __fn_ok_2not1
    if RUN_TYPE == 'rerun':
        return __fn_rerun_1and2, __fn_rerun_1not2, __fn_rerun_2not1






#def get_zipped_coadd_magnitudes(fn_base_g, fn_base_r, fn_base_i, fn_base_z):
    """Get magnitudes, magnitude errors, and flags from zipped coadd catalogs.

    Parameters
    ----------
    fn*
        Filenames

    Returns
    -------
    (ndarray)
    """
    '''
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


    __base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz = np.empty(len(m_g), dtype='S100'), np.empty(len(m_g), dtype='S100'), np.empty(len(m_g), dtype='S100'), np.empty(len(m_g), dtype='S100')

    for i in np.arange(0, len(m_g)):
        
        __base_mag_griz[i] = '({}, {}, {}, {})'.format(m_g[i], m_r[i], m_i[i], m_z[i])
        __base_mag_err_griz[i] = '({}, {}, {}, {})'.format(err_g[i], err_r[i], err_i[i], err_z[i])
        __base_flag1_griz[i] = '({}, {}, {}, {})'.format(flag1_g[i], flag1_r[i], flag1_i[i], flag1_z[i])
        __base_flag2_griz[i] = '({}, {}, {}, {})'.format(flag2_g[i], flag2_r[i], flag2_i[i], flag2_z[i])

    return __base_mag_griz, __base_mag_err_griz, __base_flag1_griz, __base_flag2_griz
    '''



def stack_tiles(realization, output_directory, balrog_run, base_path_to_catalogs, all_tiles):
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
    stack_dir = outputs.get_directory(tile='stack', realization=realization, low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE) 

    # Filename for stacked catalogs #
    __fn_stack_tiles_1and2 = os.path.join(stack_dir, 'stacked_{}_{}_match1and2.csv'.format(realization, MATCH_TYPE))
    __fn_stack_tiles_1not2 = os.path.join(stack_dir, 'stacked_{}_{}_match1not2.csv'.format(realization, MATCH_TYPE))
    __fn_stack_tiles_2not1 = os.path.join(stack_dir, 'stacked_{}_{}_match2not1.csv'.format(realization, MATCH_TYPE))

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

        for t in all_tiles:

            if RUN_TYPE is None:
                fn_1and2, fn_1not2, fn_2not1 = match_catalogs(realization=realization, tile=t, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT, output_directory=output_directory, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)[:3]

            if RUN_TYPE is not None:
                fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=realization, tile=t, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs) #TODO update with inj?

            all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

        if VERBOSE_ING: print 'Stacking tiles. Stacking {} catalogs...'.format(len(all_fn_1and2))
        df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
        df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
        df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
        if VERBOSE_ING: print 'Stacking complete ... \n'


        # Save stacked catalogs as DataFrame #
        df1and2.to_csv(__fn_stack_tiles_1and2, sep=','); df1not2.to_csv(__fn_stack_tiles_1not2, sep=','); df2not1.to_csv(__fn_stack_tiles_2not1, sep=',')
        print 'Saving stacked tile catalogs as....'
        print '----->', __fn_stack_tiles_1and2 
        print '----->', __fn_stack_tiles_1not2 
        print '----->', __fn_stack_tiles_2not1 

    __number_of_stacked_tiles = len(all_tiles)

    return __fn_stack_tiles_1and2, __fn_stack_tiles_1not2, __fn_stack_tiles_2not1, __number_of_stacked_tiles



def stack_realizations(tile, balrog_run, output_directory, base_path_to_catalogs, all_realizations):
    """Concatenate catalogs with multiple realizations and fixed tile.

    Parameters
    ----------
    tile (str)
        One stacked realization catalog created per tile.

    Returns
    -------
    fn_stack_1and2 (str)
        Complete filename for stacked catalog of join=1and2.
    fn_stack_1not2 (str)
        Complete filename for stacked catalog of type join=1not2.
    fn_stack_2not1 (str)
        Complete filename for stacked catalog of type join=2not1. 
    len(all_realizations) (int)
        Number of catalogs stacked.
    """

    # Directory for stacked catalog #
    stack_dir = outputs.get_directory(tile=tile, realization='stack', low_level_dir='catalog_compare', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE) 

    # Filename for stacked catalogs #
    __fn_stack_reals_1and2 = os.path.join(stack_dir, '{}_stacked_{}_match1and2.csv'.format(tile, MATCH_TYPE))
    __fn_stack_reals_1not2 = os.path.join(stack_dir, '{}_stacked_{}_match1not2.csv'.format(tile, MATCH_TYPE)) 
    __fn_stack_reals_2not1 = os.path.join(stack_dir, '{}_stacked_{}_match2not1.csv'.format(tile, MATCH_TYPE)) 

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

        for r in all_realizations:

            if RUN_TYPE is None:
                fn_1and2, fn_1not2, fn_2not1 = match_catalogs(realization=r, tile=tile, inj1=INJ1, inj2=INJ2, inj1_percent=INJ1_PERCENT, inj2_percent=INJ2_PERCENT, output_directory=output_directory, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)[:3]

            if RUN_TYPE is not None:
                fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=r, tile=tile, balrog_run=balrog_run, output_directory=output_directory, base_path_to_catalogs=base_path_to_catalogs) #FIXME

            all_fn_1and2.append(fn_1and2); all_fn_1not2.append(fn_1not2); all_fn_2not1.append(fn_2not1)

        if VERBOSE_ING: print 'Stacking realizations. {} files...'.format(len(all_fn_1and2))
        df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1and2))
        df1not2 = pd.concat((pd.read_csv(fn) for fn in all_fn_1not2))
        df2not1 = pd.concat((pd.read_csv(fn) for fn in all_fn_2not1))
        if VERBOSE_ING: print 'Stacking complete...'


        # Save stacked catalog as DataFrame #
        df1and2.to_csv(__fn_stack_reals_1and2, sep=','); df1not2.to_csv(__fn_stack_reals_1not2, sep=','); df2not1.to_csv(__fn_stack_reals_2not1, sep=',')
        print 'Saving stacked realization catalogs as...'
        print '----->', __fn_stack_reals_1and2 
        print '----->', __fn_stack_reals_1not2 
        print '----->', __fn_stack_reals_2not1 

    __number_of_stacked_realizations = len(all_realizations)

    return __fn_stack_reals_1and2, __fn_stack_reals_1not2, __fn_stack_reals_2not1, __number_of_stacked_realizations




def get_dataframe_and_headers(realization, tile, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, inj1, inj2, inj1_percent, inj2_percent, base_path_to_catalogs, output_directory, balrog_run, all_realizations, all_tiles):
    """Get pandas DataFrame and read it.
    Rewrite headers if necessary (that is, if catalogs have been reformatted).

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
        for r in all_realizations:
            if STACK_TILES and STACK_REALIZATIONS is False: 
                fn_1and2, fn_1not2, fn_2not1, __number_of_stacked_tiles = stack_tiles(realization=r, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs, all_tiles=all_tiles)
                __number_of_stacked_realizations = None


        ### Stack realizations ###
        for t in all_tiles:
            if STACK_REALIZATIONS and STACK_TILES is False:
                fn_1and2, fn_1not2, fn_2not1, __number_of_stacked_realizations = stack_realizations(tile=t, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs, all_realizations=all_realizations)
                __number_of_stacked_tiles = None


    
    if (STACK_REALIZATIONS is False and STACK_TILES is False) or (PLOT_COMPLETENESS and STACK_TILES):
        __number_of_stacked_tiles, __number_of_stacked_realizations = None, None
        # Filenames for catalogs #
        if RUN_TYPE is None:
            #fn_1and2, fn_1not2, fn_2not1 = matcher(realization=realization, tile=tile, inj1=inj1, inj2=inj2, inj1_percent=inj1_percent, inj2_percent=inj2_percent)[:3]
            fn_1and2, fn_1not2, fn_2not1 = match_catalogs(realization=realization, tile=tile, inj1=inj1, inj2=inj2, inj1_percent=inj1_percent, inj2_percent=inj2_percent, output_directory=output_directory, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs)[:3]
        if RUN_TYPE is not None:
            fn_1and2, fn_1not2, fn_2not1 = fof_matcher(realization=realization, tile=tile, output_directory=output_directory, balrog_run=balrog_run, base_path_to_catalogs=base_path_to_catalogs) #FIXME


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

        star_mag = analysis.get_star_truth_catalog_magnitude(df_1and2=__df1and2, suf=suf)
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
        y3_gold_mag = analysis.get_y3_gold_catalog_magnitude(df_1and2=__df1and2, mag_hdr=mhdr)
        # Add new column to df #
        __df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_MAG_GRIZ_HDR, y3_gold_mag)
        if PLOT_FLUX:
            y3_gold_flux_griz, y3_gold_flux_cov_griz = analysis.get_y3_gold_catalog_flux_griz(df_1and2=__df1and2, flux_hdr=fhdr, flux_cov_hdr=fchdr)
            __df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_FLUX_GRIZ_HDR, y3_gold_flux_griz)
            __df1and2.insert(len(__df1and2.columns), NEW_Y3_GOLD_CM_FLUX_COV_HDR, y3_gold_flux_cov_griz)
        


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
    logFileDirectory = outputs.get_directory(tile=tile, realization=realization, low_level_dir='log_files', output_directory=output_directory, balrog_run=balrog_run, match_type=MATCH_TYPE)
    fn_match_log = os.path.join(logFileDirectory, '{}_{}_{}_matched_catalogs.log'.format(tile, realization, MATCH_TYPE)) 

    print 'Saving match log as:\n----->', fn_match_log

    with open(fn_match_log, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['TILE', 'REALIZATION', 'IN1', 'IN2', 'TOTAL_MATCHES_1AND2', 'TOTAL_MATCHES_1NOT2', 'TOTAL_MATCHES_2NOT1'])
        writer.writerow([tile, realization, MATCH_CAT1, MATCH_CAT2, len(__df1and2.index), len(__df1not2.index), len(__df2not1.index)])
    csvfile.close()


    return __df1and2, __df1not2, __df2not1, __mag_hdr1, __mag_hdr2, __flux_hdr1, __flux_hdr2, __mag_err_hdr1, __mag_err_hdr2, __flux_err_hdr1, __flux_err_hdr2, __number_of_stacked_tiles, __number_of_stacked_realizations
