"""
Store catalog headers and properties (e.g. injection percent, measure or true).
"""




class CoaddCatalog():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{band}_cat.fits. There is a separate catalog for each band."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj (bool)
		inj_percent (int)
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# For plot title #
		if inj:
			self.title_piece = '{}% Inj Coadd Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base Coadd Cat'

		self.axlabel = 'meas'
		# Magnitude, is one number #
		self.mag_hdr = 'MAG_AUTO{}'.format(suf)
		self.mag_axlabel = 'MAG_AUTO_meas'
		# For error calculation #
		self.mag_err_hdr = 'MAGERR_AUTO{}'.format(suf)
		self.cm_flux_hdr = 'FLUX_AUTO{}'.format(suf) #FIXME not CM flux 
		self.cm_flux_cov_hdr = None #'FLUXERR_AUTO{}'.format(suf) 
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = 'FLAGS{}'.format(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'ALPHAWIN_J2000{}'.format(suf)
		self.dec_hdr = 'DELTAWIN_J2000{}'.format(suf)
		self.angle = 'THETA_J2000{}'.format(suf)
		# Units: pixels #
		self.a_hdr = 'A_IMAGE{}'.format(suf)
		self.b_hdr = 'B_IMAGE{}'.format(suf)
		# Identifying indices corresponding to objects in matched catalog #
		self.identifier = 'NUMBER{}'.format(suf)




class GalaxyTruthCatalog():
	"""Declare headers for galaxy truth catalogs. Note that as of May 2018 galaxy truth catalogs are created using MOF (thus have the same headers).""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj (bool)
		inj_percent (int)
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# For plot title #
		# `inj` forced True for truth catalogs #
		if inj:
			self.title_piece = '{}% Inj Gal Truth Cat'.format(inj_percent)
		self.axlabel = 'true'
		# Headers are the same as `MOFCatalog` class. Reproduced below for clarity in ms_plotter.py #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag{}'.format(suf)
		self.mag_axlabel = 'cm_mag_true'
		# For error calculation #
		self.mag_err_hdr = None
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)
		# Size. cm_T units: arcseconds squared. #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_true'
		# Flags #
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags{}'.format(suf)
		# For region file #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)
		self.a_hdr, self.b_hdr = None, None
		self.angle = None




class SOFCatalog():
	"""Declare headers for SOF catalogs.""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(cat_type_pair, self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj (bool)
		inj_percent (int)
		cat_type_pair (str)
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		if inj:
			self.title_piece = '{}% Inj SOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base SOF Cat'

		self.axlabel = 'meas'
		# Headers are the same as MOFCatalog class with the exception of cm_mof_flags_hdr. Reproduced below for clarity #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag{}'.format(suf)
		self.mag_axlabel = 'cm_mag_meas'
		if cat_type_pair == 'star_truth':
			self.mag_hdr = 'psf_mag{}'.format(suf)
			self.mag_axlabel = 'psf_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None
		#TODO there is psf_flux 
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)
		# Size #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = None
		# For region file #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)
		self.a_hdr, self.b_hdr = None, None
		self.angle = None




class MOFCatalog():
	"""Declare headers for MOF catalogs."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj_percent (int)
		inj (bool)
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""

		# For plot title #
		if inj:
			self.title_piece = '{}% Inj MOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base MOF Cat'
		self.axlabel = 'meas'
		# Magnitude, is string of form (mag_g, mag_r, mag_i, mag_z)  #
		self.mag_hdr = 'cm_mag{}'.format(suf)
		self.mag_axlabel = 'cm_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None 
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)
		# Size #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags{}'.format(suf)
		# For region file #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)
		self.a_hdr, self.b_hdr = None, None
		self.angle = None




class StarTruthCatalog(): 
	"""Declare headers for star truth catalogs."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj_percent (int)
		inj (bool)
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""	

		# `inj` forced True for truth catalogs #	
		if inj:
			self.title_piece = '{}% Star Truth Cat'.format(inj_percent)
		self.axlabel = 'true'
		# Magnitude #
		self.mag_hdr = 'g_Corr{}'.format(suf) 
		self.mag_axlabel = 'mag_true' 
		# For error calculation #
		# Is of form 'PSF_MAG_ERR_{band}{}'.format(suf) #
		self.mag_err_hdr = 'PSF_MAG_ERR{}'.format(suf)
		self.cm_flux_hdr = None
		self.cm_flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = None
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'RA_new{}'.format(suf)
		self.dec_hdr = 'DEC_new{}'.format(suf)
		self.a_hdr, self.b_hdr = None, None # Use cm_T
		self.angle = None




class Y3GoldCatalog():
	"""Declare headers for Y3 Gold catalogs (https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Full_bins_of_Y3_GOLD_2_0_Columns).""" 

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf, y3_fit, y3_model, cat_type_pair, cat_type):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.

		Parameters
		----------
		inj_percent (int)

		inj (bool)

		y3_fit (str)
			e.g. `mof` or `sof`.
		y3_model (str)
			e.g. `cm` or `psf`
		suf (str)
			Refers to the order in which catalog was matched (via join=1and2) in ms_matcher (order set by STILTS parameters `in1` and `in2`). Allowed values: '_1' '_2'.
		"""


		if cat_type_pair == 'star_truth':
			y3_model = 'PSF'
		if cat_type_pair != 'star_truth':
			y3_model = 'CM'

		print 'Using {} and {} for Y3 Gold catalog...'.format(y3_fit, y3_model)

		if cat_type == 'y3_gold_2_2':
			self.title_piece = 'Y3 Gold 2-2 Cat' 
		if cat_type == 'y3_gold_2_0':
			self.title_piece = 'Y3 Gold 2-0 Cat'

		self.axlabel = 'meas'

		# Note: band dependent headers -- MAG, MAG_ERR, FLUX, FLUX_COV, PSF_FLAGS #

		# Magnitude # 
		self.mag_hdr = '{}_{}_MAG{}'.format(y3_fit, y3_model, suf).upper()
		self.mag_axlabel = 'mag_meas'
		# For error calculation #
		self.mag_hdr = '{}_{}_MAG_ERR{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_flux_hdr = '{}_{}_FLUX{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_flux_cov_hdr = '{}_{}_FLUX_COV{}'.format(y3_fit, y3_model, suf).upper()
		# Size #
		self.cm_t_hdr = '{}_{}_T{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_t_err_hdr = '{}_{}_T_ERR{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'FLAGS_GOLD{}'.format(suf).upper()
		self.obj_flags_hdr = '{}_OBJ_FLAGS{}'.format(y3_fit, suf).upper()
		self.psf_flags_hdr = '{}_PSF_FLAGS{}'.format(y3_fit, suf).upper()
		self.cm_flags_hdr = '{}_{}_FLAGS{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_max_flags_hdr = None
		# MOF_CM_MOF_FLAGS has no SOF equivalent #
		# Note: is duplicate of `cm_flags_hdr` #
		self.cm_mof_flags_hdr = '{}_{}_FLAGS{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_flags_r_hdr = '{}_{}_FLAGS_R{}'.format(y3_fit, y3_model, suf).upper()
		# For region file #
		# Note: there is also an ALPHAWIN_J2000 and DELTAWIN_J2000 #
		self.ra_hdr = 'ra{}'.format(suf).upper()
		self.dec_hdr = 'dec{}'.format(suf).upper()
		self.a_hdr = 'A_IMAGE{}'.format(suf).upper()
		self.b_hdr = 'B_IMAGE{}'.format(suf).upper()
		self.angle = 'THETA_J2000{}'.format(suf).upper()
