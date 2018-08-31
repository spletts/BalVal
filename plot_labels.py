"""
Get plot labels.
Plot labels will be overwritten if `OVERWRITE_AXLABELS=True` in `set_constants.py`.

Comments are ABOVE the code they refer to.
"""

from difflib import SequenceMatcher

from set_constants import *
from catalog_headers import TITLE_PIECE1, TITLE_PIECE2




def get_magnitude_axlabel(inj, mag_hdr, meas_or_true_cat, match_cat, band, inj_percent):
    """Create label for the axis of the magnitude plot that describes the model used to measure the observable (e.g. 'cm'),
    the injection percent (if applicable),
    and the catalog type if the catalog (given by `match_cat`) is a Y3 Gold catalog. 
    Note that `'{}'.format()` does not render TeX reliably so `%` is used at times.

    Parameters
    ----------
    inj (bool)
        If `True` the catalog given by `match_cat` is Balrog-injected.
        If `False` the catalog is not Balrog-injected (also called 'base'). 
        Set by `INJ1` or `INJ2`.

    inj_percent (int)
        Injection density percentage for the catalog given by `match_cat`.
        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
        Note that a 10% injection density corresponds to 5,000 injected objects.

    mag_hdr (str)
        Magnitude header.
        Headers refer to the catalog matched via `join=1and2` in `stilts_matcher`.
        This may set by `M_HDR1`, `M_HDR2`, or changed if catalogs are reconstructed in the case of `star_truth`, `y3_gold`, `coadd`.

    meas_or_true_cat (str)
        Describes plot axes and catalog type.
        Allowed values: 'meas', 'true'.
        Set by `AXLABEL1` or `AXLABEL2`.

    match_cat (str)
        One of the two catalogs that was matched in `stilts_matcher`.
        Set by `MATCH_CAT1` or `MATCH_CAT2`.

    band (str)
        Allowed values: 'g', 'r', 'i', 'z'.
        Set by an element in `set_constants.ALL_BANDS`.

    Returns
    -------
    __mag_axlabel (str)
        Label for the horizontal axis of the magnitude plot.
        Example: '10%_inj_cm_mag_g_true'
        Includes LaTeX \bf{} formatting so that that band is bolded.
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
    """Create label for the axis of the color plot that describes the model used to measure the observable (e.g. 'cm'),
    the injection percent (if applicable),
    and the catalog type if the catalog (given by `match_cat`) is a Y3 Gold catalog. 
    Note that `'{}'.format()` does not render TeX reliably so `%` is used at times.

    Parameters
    ----------
    inj (bool)
        If `True` the catalog given by `match_cat` is Balrog-injected.
        If `False` the catalog is not Balrog-injected (also called 'base'). 
        Set by `INJ1` or `INJ2`.

    inj_percent (int)
        Injection density percentage for the catalog given by `match_cat`.
        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
        Note that a 10% injection density corresponds to 5,000 injected objects.

    meas_or_true_cat (str)
        Describes plot axes and catalog type.
        Allowed values: 'meas', 'true'.
        Set by `AXLABEL1` or `AXLABEL2`.

    match_cat (str)
        One of the two catalogs that was matched in `stilts_matcher`.
        Set by `MATCH_CAT1` or `MATCH_CAT2`.

    band (str)
        Allowed values: 'g', 'r', 'i', 'z'.
        Set by an element in `set_constants.ALL_BANDS`.
 
    Returns
    -------
    __color_axlabel (str)
        Label for the horizontal axis of the magnitude plot.
        Example: '10%_inj_(g-r)_true'
        Includes LaTeX \bf{} formatting so that that color is bolded.
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
    """Create a short axis label to describe the difference (`axlabel_a` minus `axlabel_b`) between two axis labels.
    The short label is created by finding the longest common substring between `axlabel1` and `axlabel`;
    therefore, if the longest common substring is short, user is recommended to set `OVERWRITE_AXLABELS=True` in set_constants.py`.
    If the order of subtraction is not `axlabel_a` minus `axlabel_b`, reverse the order of `axlabel` inputs;
    that is, `axlabel_a=label1, axlabel_b=label2` --> `axlabel_a=label2, axlabel_b=label1`.

    Parameters
    ----------
    axlabel_a, axlabel_b (str)
        Axis label for an observable in `MATCH_CAT1` or `MATCH_CAT2`.
        `axlabel_a` does not need to correspond to `MATCH_CAT1`. `axlabel_b` does not need to correspond to `MATCH_CAT2`.

    band (str)
        Allowed values: 'g', 'r', 'i', 'z'.
        Set by an element in `set_constants.ALL_BANDS`.

    Returns
    -------
    __short_axlabel (str)
        Short label that describes the difference (`axlabel_a` minus `axlabel_b`) between two observables.
        Example: '10%_inj_cm_mag_g true-meas' from `axlabel_a='10%_inj_cm_mag_g_true', axlabel_b='10%_inj_cm_mag_g_meas'`.
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
        __spare_label_begin_b = '{}0%'.format(__spare_label_begin_b)

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
        

    if len(__short_axlabel) < 5: 
        print 'User may wish to set `OVERWRITE_AXLABEL=True` in `set_constants.py`.'

    # Default labels #
    if OVERWRITE_AXLABELS:
        if PLOT_COLOR: __short_axlabel = '$\Delta(\\bf{%s})$' % WRITE_COLORS[band]
        if PLOT_FLUX: __short_axlabel = '$\Delta$flux_$\\bf{%s}/$meas_flux_err' % band
        if PLOT_MAG: __short_axlabel = '$\Delta$mag_$\\bf{%s}$' % band 


    return __short_axlabel




def get_plot_suptitle(realization, tile, number_of_stacked_realizations, number_of_stacked_tiles):
    """Create plot title that describes the two catalogs matched (in `stilts_matcher`) and their Balrog-injection percents (or lack of -- called 'base').
    The plot title reflects the order in which the catalogs were matched in `stilts_matcher`; the first catalog mentioned in the plot title was used as STILTS parameter `in1`. 
    The return will be used in `matplotlib.pyplot.suptitle` for subplots, and `matplotlib.pyplot.title()` for plots.

    Parameters
    ----------
    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    number_of_stacked_realizations (int)
        Number of catalogs in the stacked realization catalog.
        This is the same as the number of realizations entered at the command line.
        Can be `None`.

    number_of_stacked_tiles (int)
        Number of catalogs in stacked tile catalog.
        The is the same as the number of tiles entered at the command line.
        Can be `None`.

    Returns
    -------
    __plot_title (str)
        Plot title that describes various properties about the two catalogs matched in `stilts_matcher`.
        Example: '10% Inj MOF Cat & 10% Inj Truth Cat' 
        In this example, STILTS parameter `in1` was a Balrog-injected MOF catalog.
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

    if NORMALIZE_MAG:
        __plot_title = 'Normalized. ' + __plot_title 

    return __plot_title 




def get_colorbar_for_magnitude_plot_axlabel(df_1and2, cm_t_hdr, cm_t_err_hdr, idx_good, meas_or_true_cat, inj_percent, inj):
    """Create the label for the colorbar of the magnitude plot.
    This label describes the property assigned to the colorbar (either size squared ('cm_T') or size squared error).

    Parameters
    ----------
    df_1and2 (pandas DataFrame)
        DataFrame that contains objects in both `MATCH_CAT1` and `MATCH_CAT2`.
        This catalog is created in `stilts_matcher` using `join=1and2`.

    cm_t_hdr (str)
        Header for the size squared (cm_T) of objects. 
        Headers refer to the catalog matched via `join=1and2` in `stilts_matcher`.
        Set by `CM_T_HDR1` or `CM_T_HDR2`.

    cm_t_err_hdr (str)
        Header for the error in 'cm_T' (size squared of objects).
        Headers refer to the catalog matched via `join=1and2` in `stilts_matcher`.
        Set by `CM_T_ERR_HDR1` or `CM_T_ERR_HDR2`.

    idx_good (1D array_like, dtype=str)
        Indices of objects without flags* (for a description of flags* see https://github.com/spletts/BalVal/blob/master/README.md#flagged-objects).
        Indices correspond to the catalog matched via `join=1and2` in `stilts_matcher`.

    meas_or_true_cat (str)
        Describes plot axes and catalog type with allowed values: 'meas', 'true'.
        Set by `AXLABEL1` or `AXLABEL2`.

    inj (bool)
        If `True` the catalog given by `match_cat` is Balrog-injected.
        If `False` the catalog is not Balrog-injected (also called 'base'). 
        Set by `INJ1` or `INJ2`.

    inj_percent (int)
        Injection density percentage for the catalog given by `match_cat`.
        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
        Note that a 10% injection density corresponds to 5,000 injected objects.

    Returns
    -------
    __cbar_label (str)
        Label for colorbar.
        Includes LaTeX \bf{} formatting.
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
