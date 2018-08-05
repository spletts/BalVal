"""
Store catalog headers and properties.
e.g. is the catalog: Balrog-injected (if so get its injection percent), Balrog-base, measured, truth, etc.


Docstrings are repeated often, but included in each function for clarity.

Comments are ABOVE the code they refer to.
"""

import sys




class CoaddCatalog():
	"""Declare headers for coadd catalogs.
	There is a separate catalog for each band.
	Headers must be capitalized.
	Data for the headers .... are single numbers (not arrays)."""

	def __init__(self, inj, inj_percent, suf):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
		inj (bool)
			If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base'). 
			Set by `INJ1` or `INJ2`.

		inj_percent (int)
			Considered if `inj=True`, otherwise ignored.
			Injection density percentage for the catalog.
			Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
			Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
			Allowed values: '_1', '_2'
			Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
			'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
			The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj Coadd Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base Coadd Cat'

		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

		# Magnitude of each object is a single number (not an array) #
		self.mag_hdr = 'MAG_AUTO{}'.format(suf)

		# Magnitude error of each object is a single number (not an array) #
		self.mag_err_hdr = 'MAGERR_AUTO{}'.format(suf)

		# Flux of each object is a single number (not an array) #
		self.cm_flux_hdr = 'FLUX_AUTO{}'.format(suf) #FIXME not CM flux

		self.flux_hdr = 'FLUX_AUTO{}'.format(suf)

                # Flux covariance matrix #
		self.cm_flux_cov_hdr = None

		#TODO
		self.flux_err_hdr = 'FLUXERR_AUTO{}'.format(suf)

		# Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None

		#TODO can use kron radius

		### Flags ###
		self.flags_hdr = 'FLAGS{}'.format(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None

		# Units for both RA and Dec: degrees #
		self.ra_hdr = 'ALPHAWIN_J2000{}'.format(suf)
		self.dec_hdr = 'DELTAWIN_J2000{}'.format(suf)

		# Semi-major axis. Units: pixels. There is a 1:1 conversion between pixels and arcseconds #
		self.a_hdr = 'A_IMAGE{}'.format(suf)
		# Semi-minor axis #
		self.b_hdr = 'B_IMAGE{}'.format(suf)

                # `angle` required by DS9 region files for ellipses (http://ds9.si.edu/doc/ref/region.html) is given by 90+THETA_J2000 #
                # Units: degrees #
                self.angle = 'THETA_J2000{}'.format(suf)




class DeepSNMOFCatalog():
	"""Declare headers for MOF catalogs.
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z)."""

	def __init__(self, inj, inj_percent, suf, model='cm'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
		inj (bool)
			If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
			Set by `INJ1` or `INJ2`.

		inj_percent (int)
			Considered if `inj=True`, otherwise ignored.
			Injection density percentage for the catalog.
			Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
			Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
			Allowed values: '_1', '_2'
			Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
			'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
			The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj Deep SN MOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Deep SN MOF Cat'
		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

		# Magnitude of each object is an array (mag_g, mag_r, mag_i, mag_z) #
		self.mag_hdr = 'cm_mag{}'.format(suf)

		self.mag_err_hdr = None

		self.cm_flux_hdr = 'cm_flux{}'.format(suf)


		# Flux covariance matrix #
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)

		# Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)

		#TODO
		self.flux_hdr = '{}_flux{}'.format(model, suf)
		self.mag_hdr = '{}_mag{}'.format(model, suf)
		self.flux_err_hdr = '{}_flux_err{}'.format(model, suf)

		### Flags ###
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags{}'.format(suf)

		### For region file ###
		# Units for both RA and Dec: degrees #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)

		# Semi-major axis #
		self.a_hdr = None
		# Semi-minor axis #
		self.b_hdr = None

		# For DS9 region file #
		self.angle = None






class DeepSNSOFCatalog():
	"""Declare headers for SOF catalogs.
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z)."""

	def __init__(self, cat_type_pair, inj, inj_percent, suf, model='cm'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
		cat_type_pair (str)
			Refers to the catalog type that will be matched with the SOF catalog in `ms_matcher`.
			Set by `MATCH_CAT1` or `MATCH_CAT2`.

		inj (bool)
			If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
			Set by `INJ1` or `INJ2`.

		inj_percent (int)
			Considered if `inj=True`, otherwise ignored.
			Injection density percentage for the catalog.
			Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
			Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
			Allowed values: '_1', '_2'
			Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
			'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
			The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).

		model (str)
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj Deep SN SOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Deep SN SOF Cat'

		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

		# Magnitude of each object is an array (mag_g, mag_r, mag_i, mag_z) #
		self.mag_hdr = 'cm_mag{}'.format(suf)

		if cat_type_pair == 'star_truth':
			self.mag_hdr = 'psf_mag{}'.format(suf)

		self.mag_err_hdr = None

		# Flux of each object is a single number (not an array) #
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)


		#TODO
		self.flux_hdr = '{}_flux{}'.format(model, suf)
		self.mag_hdr = '{}_mag{}'.format(model, suf)
		self.flux_err_hdr = '{}_flux_err{}'.format(model, suf)

		# Flux covariance matrix #
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)

		# Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)

		### Flags ###
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = None

                # Units for both RA and Dec: degrees #
                self.ra_hdr = 'ra{}'.format(suf)
                self.dec_hdr = 'dec{}'.format(suf)

                # Semi-major axis #
                self.a_hdr = None
                # Semi-minor axis #
                self.b_hdr = None

                # For DS9 region file #
                self.angle = None




class GalaxyTruthCatalog():
	"""Declare headers for galaxy truth catalogs.
	Galaxy truth catalogs are created from MOF catalogs or COSMOS catalogs (as of August 2018).
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z).""" 

	def __init__(self, inj, inj_percent, suf, model='cm'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
                inj (bool)
			`inj=True` for truth catalogs.
                        If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
                        Set by `INJ1` or `INJ2`.

                inj_percent (int)
                        Considered if `inj=True`, otherwise ignored.
                        Injection density percentage for the catalog.
                        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
                        Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
                        Allowed values: '_1', '_2'
                        Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
                        'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
                        The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).

		model (str)
			Allowed values: 'cm', 'psf'
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj Gal Truth Cat'.format(inj_percent)

		# Axes labels will include {observable}_true #
		self.axlabel = 'true'

		# Magnitude of each object is an array (mag_g, mag_r, mag_i, mag_z) #
		self.mag_hdr = 'cm_mag{}'.format(suf)

		self.mag_err_hdr = None

		# Flux of each object is an array (flux_g, flux_r, flux_i, flux_z) #
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)

		#TODO
		self.flux_hdr = '{}_flux{}'.format(model, suf)
		self.mag_hdr = '{}_mag{}'.format(model, suf)
		self.flux_err_hdr = '{}_flux_err{}'.format(model, suf)

		# Flux covariance matrix #
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)

		# Size. cm_T units: arcseconds squared. #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)

		### Flags ###
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags{}'.format(suf)

                # Units for both RA and Dec: degrees #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)

		# Semi-major axis #
		self.a_hdr = None
		# Semi-minor axis #
		self.b_hdr = None

                # For DS9 region file #
		self.angle = None




class SOFCatalog():
	"""Declare headers for SOF catalogs.
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z).""" 

	def __init__(self, cat_type_pair, inj, inj_percent, suf, model='cm'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
		cat_type_pair (str)
			Refers to the catalog type that will be matched with the SOF catalog in `ms_matcher`.
			Set by `MATCH_CAT1` or `MATCH_CAT2`.

                inj (bool)
                        If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
                        Set by `INJ1` or `INJ2`.

                inj_percent (int)
                        Considered if `inj=True`, otherwise ignored.
                        Injection density percentage for the catalog.
                        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
                        Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
                        Allowed values: '_1', '_2'
                        Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
                        'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
                        The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj SOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base SOF Cat'

		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

                # Magnitude of each object is an array (mag_g, mag_r, mag_i, mag_z) #
		self.mag_hdr = 'cm_mag{}'.format(suf)

		if cat_type_pair == 'star_truth':
			self.mag_hdr = 'psf_mag{}'.format(suf)

		self.mag_err_hdr = None

                # Flux of each object is a single number (not an array) #
		self.cm_flux_hdr = 'cm_flux{}'.format(suf)

		#TODO
		self.flux_hdr = '{}_flux{}'.format(model, suf)
                self.mag_hdr = '{}_mag{}'.format(model, suf)
                self.flux_err_hdr = '{}_flux_err{}'.format(model, suf)

                # Flux covariance matrix #
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)

                # Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)

                ### Flags ###
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = None

                # Units for both RA and Dec: degrees #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)

                # Semi-major axis #
                self.a_hdr = None
                # Semi-minor axis #
                self.b_hdr = None

                # For DS9 region file #
		self.angle = None




class MOFCatalog():
	"""Declare headers for MOF catalogs.
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z)."""

	def __init__(self, inj, inj_percent, suf, model='cm'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
                inj (bool)
                        If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
                        Set by `INJ1` or `INJ2`.

                inj_percent (int)
                        Considered if `inj=True`, otherwise ignored.
                        Injection density percentage for the catalog.
                        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
                        Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
                        Allowed values: '_1', '_2'
                        Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
                        'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
                        The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).

		model (str)
			Allowed values: 'cm',
		"""

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Inj MOF Cat'.format(inj_percent)
		if inj is False:
			self.title_piece = 'Base MOF Cat'

		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

                # Magnitude of each object is an array (mag_g, mag_r, mag_i, mag_z) #
		self.mag_hdr = 'cm_mag{}'.format(suf)

		self.mag_err_hdr = None 

		self.cm_flux_hdr = 'cm_flux{}'.format(suf)

                # Flux covariance matrix #
		self.cm_flux_cov_hdr = 'cm_flux_cov{}'.format(suf)

                # Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)

		#TODO
                self.flux_hdr = '{}_flux{}'.format(model, suf)
                self.mag_hdr = '{}_mag{}'.format(model, suf)
                self.flux_err_hdr = '{}_flux_err{}'.format(model, suf)

                ### Flags ###
		self.flags_hdr = 'flags{}'.format(suf)
		self.obj_flags_hdr = 'obj_flags{}'.format(suf)
		self.psf_flags_hdr = 'psf_flags{}'.format(suf)
		self.cm_flags_hdr = 'cm_flags{}'.format(suf)
		self.cm_max_flags_hdr = 'cm_max_flags{}'.format(suf)
		self.cm_flags_r_hdr = 'cm_flags_r{}'.format(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags{}'.format(suf)

                ### For region file ###
                # Units for both RA and Dec: degrees #
		self.ra_hdr = 'ra{}'.format(suf)
		self.dec_hdr = 'dec{}'.format(suf)

                # Semi-major axis #
                self.a_hdr = None
                # Semi-minor axis #
                self.b_hdr = None

		# For DS9 region file #
		self.angle = None




class StarTruthCatalog(): 
	"""Declare headers for star truth catalogs.
	Data for the headers .... are arrays of form ({data}_g, {data}_r, {data}_i, {data}_z)."""

	def __init__(self, inj, inj_percent, suf, model='psf'):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
                inj (bool)
                        If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
                        Set by `INJ1` or `INJ2`.

                inj_percent (int)
                        Considered if `inj=True`, otherwise ignored.
                        Injection density percentage for the catalog.
                        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
                        Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
                        Allowed values: '_1', '_2'
                        Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
                        'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
                        The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).

		model (str)
			Allowed values: 'psf'
		"""	

		if model == 'cm':
			sys.exit("Star truth catalogs do not have CM magnitudes. Require `model='psf'`")

		# For plot title and filenames #
		if inj:
			self.title_piece = '{}% Star Truth Cat'.format(inj_percent)

		# Axes labels will include {observable}_true #
		self.axlabel = 'true'

		# For magnitude calculation #
		self.mag_hdr = 'g_Corr{}'.format(suf) 

                # Magnitude error of each object is a single number (not an array) #
		#TODO there is a band dependence
		self.mag_err_hdr = 'PSF_MAG_ERR{}'.format(suf)

		self.cm_flux_hdr = None

                # Flux covariance matrix #
		self.cm_flux_cov_hdr = None

		#TODO
		band = 'r'
                self.flux_hdr = '{}_flux_{}{}'.format(model, band, suf)
                self.mag_hdr = '{}_mag_{}{}'.format(model, band, suf)
                self.flux_err_hdr = '{}_flux_err_{}{}'.format(model, band, suf)
		#self.mag_err_hdr = '{}_mag_err_{}{}'.format(model, band, suf)


                # Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None

                ### Flags ###
		self.flags_hdr = None
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None

                ### For region file ###
                # Units for both RA and Dec: degrees #
		self.ra_hdr = 'RA_new{}'.format(suf)
		self.dec_hdr = 'DEC_new{}'.format(suf)

                # Semi-major axis #
                self.a_hdr = None
                # Semi-minor axis #
                self.b_hdr = None

                # For DS9 region file #
		self.angle = None




class Y3GoldCatalog():
	"""Declare headers for Y3 Gold catalogs.
	Headers must be in capitalized.
	Data for the headers .... are single numbers (not arrays) 
	The full list of headers is at https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Full_bins_of_Y3_GOLD_2_0_Columns""" 

	def __init__(self, cat_type, cat_type_pair, inj, inj_percent, suf, y3_fit, y3_model):
		"""Declare headers for the catalog matched via `join=1and2`. Declare descriptive constants.

		Parameters
		----------
		cat_type (str)
			Allowed values: 'y3_gold_2_0', 'y3_gold_2_2'

		cat_type_pair (str)
                        Refers to the catalog type that will be matched with the Y3 Gold catalog in `ms_matcher`.
                        Set by `MATCH_CAT1` or `MATCH_CAT2`.

                inj (bool)
                        If `True` the catalog is Balrog-injected. If `False` the catalog is not Balrog-injected ('base').
                        Set by `INJ1` or `INJ2`.

                inj_percent (int)
                        Considered if `inj=True`, otherwise ignored.
                        Injection density percentage for the catalog.
                        Set by `INJ1_PERCENT` or `INJ2_PERCENT`.
                        Note that a 10% injection density corresponds to 5,000 injected objects.

		suf (str)
                        Allowed values: '_1', '_2'
                        Once the catalogs `MATCH_CAT1` and `MATCH_CAT2` are matched via `join=1and2` with STILTS, each header changes.
                        'hdr' --> 'hdr_1' or 'hdr' --> 'hdr_2'
                        The suffix `suf` refers to whether the catalog served as `in1` or `in2` in `ms_matcher` (where `in1`, `in2`, and `join` are STILTS parameters).

		y3_fit (str)
                        Allowed values: `mof`, `sof`.
			Set by `Y3_FIT`.

                y3_model (str)
                        Allowed values: `cm` or `psf`
			Set by `Y3_MODEL`.
		"""

		y3_model = y3_model.upper()
		y3_fit = y3_fit.upper()

		if cat_type_pair == 'star_truth':
			y3_model = 'PSF'
		if cat_type_pair != 'star_truth':
			y3_model = 'CM'

		print 'Using {} and {} for Y3 Gold catalog...'.format(y3_fit, y3_model)

		# For plot title and filenames #
		if cat_type == 'y3_gold_2_2':
			self.title_piece = 'Y3 Gold 2-2 Cat' 
		if cat_type == 'y3_gold_2_0':
			self.title_piece = 'Y3 Gold 2-0 Cat'

		# Axes labels will include {observable}_meas #
		self.axlabel = 'meas'

		# Note: band dependent headers -- MAG, MAG_ERR, FLUX, FLUX_COV, PSF_FLAGS #

                # Magnitude of each object is a single number (not an array) #
		self.mag_hdr = '{}_{}_MAG{}'.format(y3_fit, y3_model, suf).upper()

		# Magnitude error of each object is a single number (not an array) #
		self.mag_err_hdr = '{}_{}_MAG_ERR{}'.format(y3_fit, y3_model, suf).upper()

		self.cm_flux_hdr = '{}_{}_FLUX{}'.format(y3_fit, y3_model, suf).upper()

                # Flux covariance matrix #
		self.cm_flux_cov_hdr = '{}_{}_FLUX_COV{}'.format(y3_fit, y3_model, suf).upper()

		#TODO find and replace `bandplacehold` with actual band
                band = 'r'
		self.flux_hdr = '{}_flux_bandplacehold{}'.format(y3_model, suf)
                #self.mag_hdr = '{}_mag_{}{}'.format(y3_model, band, suf)
                self.flux_err_hdr = '{}_flux_err_{}{}'.format(y3_model, band, suf)
                #self.mag_err_hdr = '{}_mag_err_{}{}'.format(model, band, suf)

                # Size squared of object. Units for cm_T: arcseconds squared #
		self.cm_t_hdr = '{}_{}_T{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_t_err_hdr = '{}_{}_T_ERR{}'.format(y3_fit, y3_model, suf).upper()

                ### Flags ###
		self.flags_hdr = 'FLAGS_GOLD{}'.format(suf).upper()
		self.obj_flags_hdr = '{}_OBJ_FLAGS{}'.format(y3_fit, suf).upper()
		self.psf_flags_hdr = '{}_PSF_FLAGS{}'.format(y3_fit, suf).upper()
		self.cm_flags_hdr = '{}_{}_FLAGS{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_max_flags_hdr = None
		# MOF_CM_MOF_FLAGS has no SOF equivalent #
		# Note: is duplicate of `cm_flags_hdr` #
		self.cm_mof_flags_hdr = '{}_{}_FLAGS{}'.format(y3_fit, y3_model, suf).upper()
		self.cm_flags_r_hdr = '{}_{}_FLAGS_R{}'.format(y3_fit, y3_model, suf).upper()

                ### For region file ###
                # Units for both RA and Dec: degrees #
		# Note: there is also an ALPHAWIN_J2000 and DELTAWIN_J2000 #
		self.ra_hdr = 'ra{}'.format(suf).upper()
		self.dec_hdr = 'dec{}'.format(suf).upper()

                # Semi-major axis. Units: pixels. There is a 1:1 conversion between pixels and arcseconds #
		self.a_hdr = 'A_IMAGE{}'.format(suf).upper()
                # Semi-minor axis #
		self.b_hdr = 'B_IMAGE{}'.format(suf).upper()

                # `angle` required by DS9 region files for ellipses (http://ds9.si.edu/doc/ref/region.html) is given by 90+THETA_J2000 #
		# Units: degrees #
		self.angle = 'THETA_J2000{}'.format(suf).upper()
