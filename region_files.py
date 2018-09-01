"""
Make DS9 region files.
Three region files are created:
(1) a region file for the objects in both catalogs specified in `set_constants.py`,
(2) a region file for the objects uniquely in the catalog set by `MATCH_CAT1` (`INJ1`, `INJ1_PERCENT`, etc) and not `MATCH_CAT2`
(3) a region file for the objects uniquely in the catalog set by `MATCH_CAT2` and not `MATCH_CAT1`.

Comments are ABOVE the code they refer to.
"""

import numpy as np
import os

# BalVal #
import outputs
from catalog_headers import MATCH_TYPE, RA_HDR1, RA_HDR2, DEC_HDR1, DEC_HDR2, MAJOR_AX_HDR1, MAJOR_AX_HDR2, MINOR_AX_HDR1, MINOR_AX_HDR2, ANGLE_HDR1, ANGLE_HDR2, CM_T_HDR1, CM_T_HDR2
from set_constants import RUN_TYPE, MATCH_CAT1, MATCH_CAT2, VERBOSE_ING


REGION_COLORS = {'1and2':'green', '1not2':'blue', '2not1':'green'}




def write_region_file_line(ra, dec, fn, color, shape, overwrite, a=None, b=None, radius=None, angle=None):
    """Write to a single DS9 region file.
    J2000 coordinate system is assumed.

    Parameters
    ----------
    ra (1D array_like, dtype=float)
        Right ascension of each object. Units: degrees.

    dec (1D array_like, dtype=float)
        Declination of each object. Units: degrees.

    fn (str)
        Complete filename for the region file.

    color (str)
        Color for regions.
        Allowed values: see http://ds9.si.edu/doc/ref/region.html#RegionProperties
            
    shape (str)
        Shape of each region.
        Allowed values: 'cirlce', 'ellipse'.

    overwrite (bool)
        If `True` the region file (given by `fn`) is overwritten if it already exists.
        If `False` the region file (given by `fn`) is not overwritten if it already exists.  

    a (1D array_like, dtype=float)
        Semi-major axis of each object. Units: arcseconds.
        Note that `A_IMAGE` and `B_IMAGE` are in units of pixels, but the conversion from DECam pixels to arcseconds is 1:1.
        Only considered if `shape=ellipse`.

    b (1D array_like, dtype=float)
        Semi-minor axis of each object. Units: arcseconds.
        Note that `A_IMAGE` and `B_IMAGE` are in units of pixels, but the conversion from DECam pixels to arcseconds is 1:1.
        Only considered if `shape=ellipse`.

    angle (1D array_like, dtype=float)
        Rotation angle of each object. Units: degrees. 
        Note that `THETA_IMAGE = THETA_J2000-90`
        Only considered if `shape=ellipse`.

    radius (1D array_like, dtype=float)
        Radius of each object. Units: arcseconds.
        Only considered if `shape=circle`.
        If `shape=cirlce` and `radius=None` a default radius of 2 arcseconds is used.

    Returns
    -------
    fn (str)
        Complete filename of the region file that has been written or rewritten to.
    """

    if VERBOSE_ING: print 'Writing to region file...'

    if overwrite: raw_input('`__overwrite=True` in `write_to_region_file()`. Press enter to procees and ctrl+c to stop.')

    if os.path.isfile(fn) is False or overwrite:

        # File descriptor #
        fd = open(fn, 'w')
        
        # Coordinate system #
        fd.write('J2000\n')

        # Write #
        if shape == 'ellipse':
            for i in np.arange(0, len(ra)):
                fd.write('ellipse {} {} {}" {}" {} #color={} width=3\n'.format(ra[i], dec[i], a[i], b[i], 90+angle[i], color))

        if shape == 'circle' and radius is not None:
            for i in np.arange(0, len(ra)):
                fd.write('circle {} {} {}" #color={} width=3\n'.format(ra[i], dec[i], radius[i], color))

        if shape == 'circle' and radius is None:
            fd.write('circle {} {} 2" #color={} width=3\n'.format(ra[i], dec[i], color))

        fd.close()

    return fn




def write_region_files(df_1and2, df_1not2, df_2not1, realization, tile, balrog_run, output_directory):
    """Write to three DS9 region files:
    (1) region file that contains objects in both `MATCH_CAT1` and `MATCH_CAT2`
    (2) region file that contains objects uniquely in `MATCH_CAT1` and not in `MATCH_CAT2`
    and (3) region file that contains objects uniquely in `MATCH_CAT2` and not in `MATCH_CAT1`.

    Parameters
    ----------
    df_1and2 (pandas DataFrame)
        DataFrame that contains objects in both `MATCH_CAT1` and `MATCH_CAT2`.
        This catalog is created in `stilts_matcher` using `join=1and2`.

    df_1not2 (pandas DataFrame)
        DataFrame that contains objects uniquely in `MATCH_CAT1` and not `MATCH_CAT2`.
        This catalog is created in `stilts_matcher` using `join=1not2`.

    df_2not1 (pandas DataFrame)
        DataFrame that contains objects uniquely in `MATCH_CAT2` and not `MATCH_CAT1`.
        This catalog is created in `stilts_matcher` using `join=2not1`.

    realization (str)
        A realization refers to a Balrog injection.
        Realizations are numbered starting with 0. Realization can be 'None' (str).
        The realization will be applied to `MATCH_CAT1` if `INJ1=True`, and `MATCH_CAT2` if `INJ2=True`.
        If both `INJ1=False` and `INJ2=False` this should be set to `None` at the command line.
        Note that although `realization` is an int, it is of type str because it is set using `sys.argv[]` (whose contents are strings).

    tile (str)
        A sky area unit used by DESDM [Dark Energy Survey Data Management] to parcel the DES footprint and organize the coadd outputs. Each tile is 0.7306 degrees on a side'. 
        See Appendix A of https://arxiv.org/abs/1801.03181

    balrog_run (str)
        Names the Balrog test. This is a substring of `BASE_PATH_TO_CATS` entered at the command line.
        Example: `base_path_to_catalogs=/data/des71.a/data/kuropat/des2247-4414_sof/` --> `balrog_run=des2247-4414_sof`.

    output_directory (str)
        All outputs of BalVal will be saved in `output_directory`.
        For specific directory structure, see https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure

    Returns
    -------
    fnRegion1and2 (str)
        Complete filename for the region file that contains objects in both `MATCH_CAT1` and `MATCH_CAT2`.

    fnRegion1not2 (str)
        Complete filename for the region file that contains objects uniquely in `MATCH_CAT1` and not `MATCH_CAT2`.

    fnRegion2not1 (str)
        Complete filename for the region file that contains objects uniquely in `MATCH_CAT2` and not `MATCH_CAT1`.
    """

    __overwrite = False

    fnRegion1and2, fnRegion1not2, fnRegion2not1 = outputs.get_region_filenames(balrog_run=balrog_run, match_type=MATCH_TYPE, output_directory=output_directory, realization=realization, tile=tile)


    ### Get RA and Dec ###
    # Handle matched catalog #
    if RUN_TYPE is None:
        ra_hdr1, dec_hdr1 = '{}_1'.format(RA_HDR1), '{}_1'.format(DEC_HDR1)
        ra_hdr2, dec_hdr2 = '{}_2'.format(RA_HDR2), '{}_2'.format(DEC_HDR2)
    if RUN_TYPE is not None:
    # MOF or SOF catalogs #
        ra1 = 'ra'; dec1 = 'dec'
        ra2 = 'ra_2'; dec2 = 'dec_2' 

    ### Get position. Arbitrarily using MATCH_CAT1 for RA and DEC for `1and2` (are similar) ###
    ra_1and2, dec_1and2 = df_1and2[ra_hdr1], df_1and2[dec_hdr1] 
    ra_1not2, dec_1not2 = df_1not2[ra_hdr1], df_1not2[dec_hdr1]
    ra_2not1, dec_2not1 = df_2not1[ra_hdr2], df_2not1[dec_hdr2]


    ### Write to region file for matched catalog. Units are arcseconds. ###

    if MATCH_CAT1 != 'coadd' and 'y3_gold' not in MATCH_CAT1:
        size_sq_1and2 = np.array(df_1and2[CM_T_HDR1])
        size_sq_1and2 = size_sq_1and2[size_sq_1and2 > 0]
        size_1and2 = np.sqrt(np.array(size_sq_1and2))

        size_sq_1not2 = np.array(df_1not2[CM_T_HDR1])
        size_sq_1not2 = size_sq_1not2[size_sq_1not2 > 0]
        size_1not2 = np.sqrt(np.array(size_sq_1not2))

        fnRegion1and2 = write_region_file_line(ra=ra_1and2, dec=dec_1and2, fn=fnRegion1and2, color='green', shape='cirlce', radius=size_1and2, overwrite=__overwrite)
        fnRegion1not2 = write_region_file_line(ra=ra_1not2, dec=dec_1not2, fn=fnRegion1not2, color='blue', shape='circle', radius=size_1not2, overwrite=__overwrite)


    if MATCH_CAT2 != 'coadd' and 'y3_gold' not in MATCH_CAT2:
        size_sq_1and2 = np.array(df_1and2[CM_T_HDR1])
        size_sq_1and2 = size_sq_1and2[size_sq_1and2 > 0]
        size_1and2 = np.sqrt(np.array(size_sq_1and2))

        size_sq_2not1 = np.array(df_2not1[CM_T_HDR2])
        size_sq_2not1 = size_sq_2not1[size_sq_2not1 > 0]
        size_2not1 = np.sqrt(np.array(size_sq_2not1))

        fnRegion1and2 = write_region_file_line(ra=ra_1and2, dec=dec_1and2, fn=fnRegion1and2, color='green', shape='cirlce', radius=size_1and2, overwrite=__overwrite)
        fnRegion2not1 = write_region_file_line(ra=ra_2not1, dec=dec_2not1, fn=fnRegion2not1, color='yellow', shape='circle', radius=size_2not1, overwrite=__overwrite)


    # Rewrite if catalogs can make more descriptive region files #
    # Coadds allow for elliptical regions #
    if MATCH_CAT1 == 'coadd' or 'y3_gold' in MATCH_CAT1:
        a_1and2, b_1and2= df_1and2[MAJOR_AX_HDR1], df_1and2[MINOR_AX_HDR1]
        a_1not2, b_1not2 = df_1not2[MAJOR_AX_HDR1], df_1not2[MINOR_AX_HDR1]
        orientation_1and2, orientation_1not2 = df_1and2[ANGLE_HDR1], df_1not2[ANGLE_HDR1]

        fnRegion1and2 = write_region_file_line(ra=ra_1and2, dec=dec_1and2, fn=fnRegion1and2, color='green', shape='ellipse', a=a_1and2, b=b_1and2, radius=None, angle=orientation_1and2, overwrite=__overwrite)
        fnRegion1not2 = write_region_file_line(ra=ra_1not2, dec=dec_1not2, fn=fnRegion1not2, color='blue', shape='ellipse', a=a_1not2, b=b_1not2, radius=None, angle=orientation_1not2, overwrite=__overwrite)


    if MATCH_CAT2 == 'coadd' or 'y3_gold' in MATCH_CAT2:
        a_1and2, b_1and2= df_1and2[MAJOR_AX_HDR2], df_1and2[MINOR_AX_HDR2]
        a_2not1, b_2not1 = df_2not1[MAJOR_AX_HDR2], df_2not1[MINOR_AX_HDR2]
        orientation_1and2, orientation_2not1 = df_1and2[ANGLE_HDR2], df_2not1[ANGLE_HDR2]

        fnRegion1and2 = write_region_file_line(ra=ra_1and2, dec=dec_1and2, fn=fnRegion1and2, color='green', shape='ellipse', a=a_1and2, b=b_1and2, radius=None, angle=orientation_1and2, overwrite=__overwrite)
        fnRegion2not1 = write_region_file_line(ra=ra_2not1, dec=dec_2not1, fn=fnRegion2not1, color='yellow', shape='ellipse', a=a_2not1, b=b_2not1, radius=None, angle=orientation_2not1, overwrite=__overwrite)


    return fnRegion1and2, fnRegion1not2, fnRegion2not1 
