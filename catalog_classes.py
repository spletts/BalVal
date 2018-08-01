"""
Store catalog headers and properties.
"""

class CoaddCat():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{band}_cat.fits. There is a separate catalog for each band."""

        # Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj_percent, inj, suf):
		"""Declare headers for matched (via join=1and2) catalog. Declare descriptive constants.
		Parameters
		----------
		inj (bool)
		inj_percent (int)
			If True refers to 20% Balrog-injected catalog. If False refers to base (non-Balrog-injected) catalog. 
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
		self.mag_hdr = 'MAG_AUTO' + str(suf)
		self.mag_axlabel = 'MAG_AUTO_meas'
		# For error calculation #
		self.mag_err_hdr = 'MAGERR_AUTO' + str(suf)
		self.cm_flux_hdr = 'FLUX_AUTO' + str(suf) #FIXME not CM flux 
		self.cm_flux_cov_hdr = None #'FLUXERR_AUTO' + str(suf) 
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = 'FLAGS' + str(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'ALPHAWIN_J2000' + str(suf)
		self.dec_hdr = 'DELTAWIN_J2000' + str(suf)
		self.angle = 'THETA_J2000' + str(suf)
		# Units: pixels #
		self.a_hdr = 'A_IMAGE' + str(suf)
		self.b_hdr = 'B_IMAGE' + str(suf)
		# Identifying indices corresponding to objects in matched catalog #
		self.identifier = 'NUMBER' + str(suf)




class GalTruthCat():
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
		# Headers are the same as MOFCat class. Reproduced below for clarity in ms_plotter.py #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
	    	self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_true'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size. cm_T units: arcseconds squared. #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_true'
                # Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
		self.angle = None
		# #
		self.identifier = 'number' + str(suf)



	
class SOFCat():
        """Declare headers for SOF catalogs.""" 

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

		if inj:
			self.title_piece = '{}% Inj SOF Cat'.format(inj_percent)
		if inj is False:
                        self.title_piece = 'Base SOF Cat'

		self.axlabel = 'meas'
                # Headers are the same as MOFCat class with the exception of cm_mof_flags_hdr. Reproduced below for clarity #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
		#FIXME test?
		if 'star_truth' in (MATCH_CAT1, MATCH_CAT2):
			self.mag_hdr = 'psf_mag' + str(suf)
			self.mag_axlabel = 'psf_mag_meas'
                # For error calculation #
                self.mag_err_hdr = None
		#TODO there is psf_flux 
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
                self.cm_mof_flags_hdr = None
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
		self.angle = None




class MOFCat():
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
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None 
		self.cm_flux_hdr = 'cm_flux' + str(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
		# Size #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
		# For region file #
		self.ra_hdr = 'ra' + str(suf)
		self.dec_hdr = 'dec' + str(suf)
		self.a_hdr, self.b_hdr = None, None
		self.angle = None




class StarTruthCat(): #are there sep headers for MOFStarTruthCat and SOFStarTruthCat?
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
		self.mag_hdr = 'g_Corr' + str(suf) 
		self.mag_axlabel = 'mag_true' 
		# For error calculation #
                # Is of form 'PSF_MAG_ERR_{band}' + str(suf) #
		self.mag_err_hdr = 'PSF_MAG_ERR' + str(suf)
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
		self.ra_hdr = 'RA_new' + str(suf)
		self.dec_hdr = 'DEC_new' + str(suf)
		self.a_hdr, self.b_hdr = None, None # Use cm_T
		self.angle = None




class Y3Gold():
	"""Declare headers for Y3 Gold catalogs (https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Full_bins_of_Y3_GOLD_2_0_Columns).""" 

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

		if inj: sys.exit('error in class Y3Gold')

		if BALROG_RUN in ('grid_bkg', 'grid_bkg_noise', 'seed_test_base', 'seed_test_same_seeds' 'seed_test_same_seeds1', 'seed_test_new_seeds', 'des2247-4414_sof', 'grid_test_noise_sof') or Y3_MOF is False: 
			pref1 = 'SOF_'

		if Y3_MOF is True:
			pref1 = 'MOF_'

		if 'star' in MATCH_CAT1 or 'star' in MATCH_CAT2:
			pref2 = 'PSF_'
		if 'star' not in MATCH_CAT1 and 'star' not in MATCH_CAT2:
			pref2 = 'CM_'

		print 'Using ', pref1, '&', pref2, ' for Y3 Gold catalog ... \n'
	
		if 'y3_gold_2_2' in (MATCH_CAT1, MATCH_CAT2): 
			self.title_piece = 'Y3 Gold 2-2 Cat' 
		if 'y3_gold_2_0' in (MATCH_CAT1, MATCH_CAT2):
                        self.title_piece = 'Y3 Gold 2-0 Cat'
                self.axlabel = 'meas'

		# Note: band dependent headers -- MAG, MAG_ERR, FLUX, FLUX_COV, PSF_FLAGS #

                # Magnitude # 
                self.mag_hdr = pref1 + pref2 + 'MAG' + str(suf) 
                self.mag_axlabel = 'mag_meas'
                # For error calculation #
                self.mag_err_hdr = pref1 + pref2 + 'MAG_ERR' + str(suf) 
                self.cm_flux_hdr = pref1 + pref2 + 'FLUX' + str(suf) 
                self.cm_flux_cov_hdr = pref1 + pref2 + 'FLUX_COV' + str(suf) 
                # Size #
                self.cm_t_hdr = pref1 + pref2 + 'T' + str(suf) 
                self.cm_t_err_hdr = pref1 + pref2 + 'T_ERR' + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
                # Flags #
                self.flags_hdr = 'FLAGS_GOLD' + str(suf)
                self.obj_flags_hdr = pref1 + 'OBJ_FLAGS' + str(suf)
                self.psf_flags_hdr = pref1 + 'PSF_FLAGS_filt' + str(suf) 
                self.cm_flags_hdr = pref1 + pref2 + 'FLAGS' + str(suf) 
                self.cm_max_flags_hdr = None
		# MOF_CM_MOF_FLAGS has no SOF equivalent #
		# Note: is duplicate #
                self.cm_mof_flags_hdr = pref1 + pref2 + 'FLAGS' + str(suf)
                self.cm_flags_r_hdr = pref1 + pref2 + 'FLAGS_R' + str(suf)
                # For region file #
		# Note: there is also an ALPHAWIN_J2000 and DELTAWIN_J2000 #
		self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr = 'A_IMAGE'+ str(suf)
		self.b_hdr = 'B_IMAGE' + str(suf) 
                self.angle = 'THETA_J2000' + str(suf)
