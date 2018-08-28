"""
Get injection percent.
Note that star truth catalogs and galaxy truth catalogs in the same `base_path_to_catalogs` may not have the same number of injected objects.

Comments are ABOVE the code they refer to.
"""

import fitsio
import os

# From BalVal #
import manipulate_catalogs




def get_injection_percent(balrog_run, base_path_to_catalogs, cat_types, realization, tile, truth_cat_type='gal_truth'):
    """Get the injection percent of a truth catalog in `base_path_to_catalogs`.
    If both galaxy truth catalogs and star truth catalogs exist in `base_path_to_catalogs`, the catalog type must be set by ``truth_cat_type`.
    If one truth catalog exists in `base_path_to_catalogs`, this catalog is chosen by default.

    Parameters
    ----------
    balrog_run (str)
        Names the Balrog test. This is a substring of `base_path_to_catalogs`.
        Ex: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.
    base_path_to_catalogs (str)
        Complete base path to catalogs set by `MATCH_CAT1` and `MATCH_CAT2`.
    cat_types (list of str)
        Types of catalogs to be analysed.
        Set by `MATCH_CAT1` and `MATCH_CAT2`.
    realization (str)
        A realization refers to a Balrog injection. Realizations are numbered starting with 0.
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).
    tile (str)
        'A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181
    truth_cat_type (str)
        Allowed values: 'gal_truth' 'star_truth'
        Considered if `base_path_to_catalogs` contains both a galaxy truth catalog and a star truth catalog.

    Returns
    -------
    __inj_percent (float)
        Injection density percentage for a truth catalog in `base_path_to_catalogs`. 
        Note that a 10% injection density is 5,000 injected objects.
    """

    __objs_per_tile = 50000.0

    ### Get truth catalog ###
    if 'star_truth' in cat_types:
        # `inj_percent` is only used for TAMU Balrog runs, in which case `INJ1_PERCENT` is force set and this function will not be called. `band` is only used for coadd catalogs. #
        __fn_cat = get_catalog.get_catalog_filename(cat_type='star_truth', inj=True, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

    if 'gal_truth' in cat_types:
        __fn_cat = get_catalog.get_catalog_filename(cat_type='gal_truth', inj=True, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

    if 'star_truth' in cat_types and 'gal_truth' in cat_types:
        __fn_cat = get_catalog.get_catalog_filename(cat_type=truth_cat_type, inj=True, inj_percent=None, realization=realization, tile=tile, band=None, base_path_to_catalogs=base_path_to_catalogs, balrog_run=balrog_run)

    ### Read truth catalog ###
    data = fitsio.read(__fn_cat, hud=1)
    __num_rows = len(data)

    # Approximate objects per tile: 50,000 #
    __inj_percent = 100*(1.0*__num_rows/__objs_per_tile)

    __inj_percent = int(__inj_percent)

    
    return __inj_percent
