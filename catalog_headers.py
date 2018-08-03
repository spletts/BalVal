"""
Get headers for the catalogs.
Get other catalog properties dependent on `class` (e.g. `title_piece` and `MATCH_TYPE`) 
Set new header names that will be used if catalogs are reformatted.
"""

# From BalVal #
import catalog_classes
from set_constants import RUN_TYPE, MATCH_CAT1, MATCH_CAT2, INJ1, INJ2, Y3_FIT, Y3_MODEL 


#FIXME
INJ1_PERCENT, INJ2_PERCENT = 20, 20




def get_class(cat_type, cat_type_pair, inj, inj_percent, suf):
	"""Get the appropriate class for the catalog type.

	Parameters
	----------
	cat_type (str)
		Catalog type. This is set by `MATCH_CAT1` or `MATCH_CAT2`. 

	inj_percent (int)

	inj (bool)

	suf (str)
		Refers to the order in which `cat_type` was matched (via join=1and2) in ms_matcher (order set by `in1` and `in2` in STILTS script). Allowed values: '_1' '_2'.

	Returns
	-------
	__class (class)
		Points to the appropriate class and class constants.
	"""

	if cat_type == 'gal_truth':
		__cat_class = catalog_classes.GalaxyTruthCatalog(inj_percent=inj_percent, inj=inj, suf=suf)

	if cat_type == 'mof':
		__cat_class = catalog_classes.MOFCatalog(inj_percent=inj_percent, inj=inj, suf=suf)

	if cat_type == 'deep_sn_mof':
		__cat_class = catalog_classes.DeepSNMOFCatalog(inj_percent=inj_percent, inj=inj, suf=suf)

	if cat_type == 'star_truth':
		__cat_class = catalog_classes.StarTruthCatalog(inj_percent=inj_percent, inj=inj, suf=suf)

	if cat_type == 'sof':
		__cat_class = catalog_classes.SOFCatalog(inj_percent=inj_percent, inj=inj, suf=suf, cat_type_pair=cat_type_pair)

	if cat_type == 'deep_sn_sof':
		__cat_class = catalog_classes.DeepSNSOFCatalog(inj_percent=inj_percent, inj=inj, suf=suf, cat_type_pair=cat_type_pair)

	if cat_type == 'coadd':
		__cat_class = catalog_classes.CoaddCatalog(inj_percent=inj_percent, inj=inj, suf=suf)

	if 'y3_gold' in cat_type:
		__cat_class = catalog_classes.Y3GoldCatalog(inj_percent=inj_percent, inj=inj, suf=suf, y3_fit=Y3_FIT, y3_model=Y3_MODEL, cat_type_pair=cat_type_pair, cat_type=cat_type)


	return __cat_class 




def get_match_type(title_piece1, title_piece2):
	"""Transform two strings of form '10% Inj MOF Cat' and '10% Inj Truth Cat' to '10%_inj_mof_cat_10%_inj_truth_cat'.

	Parameters
	----------
	title_piece1, title_piece2 (str)
		Set by `self.title_piece` in appropriate class. Class set by `MATCH_CAT*` `INJ*PERCENT`.
		Example: 10% Inj MOF Cat.

	Returns
	-------
	match_type (str)
		Reflects the order in which catalogs were matched (via join=1and2). 
		Example: 10%_inj_mof_cat_10%_inj_truth_cat, where the injected MOF catalog was STILTS parameter `in1`.
	"""

	__match_type = '{}_{}'.format(title_piece1.lower().replace(' ', '_'), title_piece2.lower().replace(' ', '_'))

	return __match_type




#def get_injection_percent(cat_types, tile, realization, base_path_to_catalogs, balrog_run):




### Get class for `MATCH_CAT1` and `MATCH_CAT2` ###
# `MATCH_CAT1` #
if RUN_TYPE is not None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, cat_type_pair=MATCH_CAT2, inj=INJ1, inj_percent=INJ1_PERCENT, suf='')
if RUN_TYPE is None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, cat_type_pair=MATCH_CAT2, inj=INJ1, inj_percent=INJ1_PERCENT, suf='_1')
# `MATCH_CAT2` #
CLASS2 = get_class(cat_type=MATCH_CAT2, cat_type_pair=MATCH_CAT1, inj=INJ2, inj_percent=INJ2_PERCENT, suf='_2')


# For plot title and filenames #
TITLE_PIECE1, TITLE_PIECE2 = CLASS1.title_piece, CLASS2.title_piece

MATCH_TYPE = get_match_type(title_piece1=TITLE_PIECE1, title_piece2=TITLE_PIECE2)

### Headers ###

# RA and Dec will be passed to `ms_matcher`. `ms_matcher` accepts 'ra' and input, not 'ra_1' (for example). #
# 'ra_1' --> 'ra' #
RA_HDR1, RA_HDR2 = CLASS1.ra_hdr[:-2], CLASS2.ra_hdr[:-2]
DEC_HDR1, DEC_HDR2 = CLASS1.dec_hdr[:-2], CLASS2.dec_hdr[:-2]

AXLABEL1, AXLABEL2 = CLASS1.axlabel, CLASS2.axlabel

# Magnitudes #
M_HDR1, M_HDR2 = CLASS1.mag_hdr, CLASS2.mag_hdr

# For magnitude error calculation #
M_ERR_HDR1, M_ERR_HDR2 = CLASS1.mag_err_hdr, CLASS2.mag_err_hdr
CM_FLUX_HDR1, CM_FLUX_HDR2 = CLASS1.cm_flux_hdr, CLASS2.cm_flux_hdr
CM_FLUX_COV_HDR1, CM_FLUX_COV_HDR2 = CLASS1.cm_flux_cov_hdr, CLASS2.cm_flux_cov_hdr

FLUX_HDR1, FLUX_HDR2 = CLASS1.flux_hdr, CLASS2.flux_hdr

# For signal to noise calculation #
CM_T_HDR1, CM_T_HDR2 = CLASS1.cm_t_hdr, CLASS2.cm_t_hdr
CM_T_ERR_HDR1, CM_T_ERR_HDR2 = CLASS1.cm_t_err_hdr, CLASS2.cm_t_err_hdr

# Flags #
FLAGS_HDR1, FLAGS_HDR2 = CLASS1.flags_hdr, CLASS2.flags_hdr
OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2 = CLASS1.obj_flags_hdr, CLASS2.obj_flags_hdr
# psf_flags is a string of form '(0,0,0,0)'; must pass through get_floats_from_string() if used. #
PSF_FLAGS_HDR1, PSF_FLAGS_HDR2 = CLASS1.psf_flags_hdr, CLASS2.psf_flags_hdr
CM_FLAGS_HDR1, CM_FLAGS_HDR2 = CLASS1.cm_flags_hdr, CLASS2.cm_flags_hdr
CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2 = CLASS1.cm_max_flags_hdr, CLASS2.cm_max_flags_hdr
CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 = CLASS1.cm_flags_r_hdr, CLASS2.cm_flags_r_hdr
CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2 = CLASS1.cm_mof_flags_hdr, CLASS2.cm_mof_flags_hdr




### For reformatting catalogs ###
# For base catalogs from COSMOS,... #
NEW_BASE_FLAG_GRIZ_HDR, NEW_BASE_IMAFLAG_ISO_GRIZ_HDR = 'FLAGS_GRIZ', 'IMAFLAGS_ISO_GRIZ'
NEW_BASE_MAG_GRIZ_HDR, NEW_BASE_MAG_ERR_GRIZ_HDR = 'MAG_AUTO_GRIZ', 'MAGERR_AUTO_GRIZ'

# For coadds #
NEW_COADD_MAG_GRIZ_HDR, NEW_COADD_MAG_ERR_GRIZ_HDR = 'MAG_AUTO_GRIZ', 'MAGERR_AUTO_GRIZ' #?
NEW_COADD_FLUX_GRIZ_HDR, NEW_COADD_FLUX_ERR_GRIZ_HDR = 'FLUX_AUTO_GRIZ', 'FLUXERR_AUTO_GRIZ'

# For star truth catalogs #
NEW_STAR_TRUTH_MAG_GRIZ_HDR = 'MAG_FROM_CORR_GRIZ'

# For Guassian aperture flux measurements #
GAUSS_APER_FLAGS_HDR = 'CM_GAP_FLAGS'
GAUSS_APER_FLUX_GRIZ_HDR = 'CM_GAP_FLUX_GRIZ'


if 'y3_gold' in MATCH_CAT1 or 'y3_gold' in MATCH_CAT2:
	NEW_Y3_GOLD_CM_FLUX_COV_HDR = '{}_CM_FLUX_COV'.format(Y3_FIT)

# For Y3 Gold #
if 'y3_gold' in MATCH_CAT1:
	NEW_Y3_GOLD_MAG_GRIZ_HDR = M_HDR1[:-2]+'_GRIZ'+'_1'
	#TODO make clear this is a new header
	NEW_Y3_GOLD_FLUX_GRIZ_HDR = '{}_CM_FLUX_GRIZ_1'.format(Y3_FIT)
if 'y3_gold' in MATCH_CAT2:
	NEW_Y3_GOLD_MAG_GRIZ_HDR = M_HDR2[:-2]+'_GRIZ'+'_2'
	NEW_Y3_GOLD_FLUX_GRIZ_HDR = '{}_CM_FLUX_GRIZ_2'.format(Y3_FIT)




# For region file #
MAJOR_AX_HDR1, MAJOR_AX_HDR2 = CLASS1.a_hdr, CLASS2.a_hdr 
MINOR_AX_HDR1, MINOR_AX_HDR2 = CLASS1.b_hdr, CLASS2.b_hdr
ANGLE1, ANGLE2 = CLASS1.angle, CLASS2.angle

FLAG_HDR_LIST = [ FLAGS_HDR1, FLAGS_HDR2, CM_FLAGS_HDR1, CM_FLAGS_HDR2, CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2, OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2, PSF_FLAGS_HDR1, PSF_FLAGS_HDR2, CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2, CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 ]
