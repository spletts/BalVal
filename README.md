# BalVal

Conducts some [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`ms_plotter.py` calls `ms_matcher`, analyzes the matched catalog, and produces various magnitude versus Delta-magnitude plots.

___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realization tile`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/ /BalVal/ 0 DES2247-4414`

`all` is an allowed value for both `realization` and `tile`.

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment `NOTICE` to remove prompt. User should set plot attributes and catalog attributes within `ms_plotter.py`.

**Plot attributes**

Examples of plot attributes, which are explained in `ms_plotter.py` are:

`NORMALIZE` `HEXBIN` `CM_T_S2N_COLORBAR` `CM_T_COLORBAR`  `PLOT_1SIG` (Booleans).


**Catalog attributes**

Catalogs are specified using:

`MATCH_CAT1` `MATCH_CAT2` with allowed values: `mof` `sof` `star_truth` `gal_truth` `coadd`.

`INJ1` `INJ2` (Booleans) which describe if the catalog is Balrog-injected `INJ=True` or non-injected/original `INJ=False`.

**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`.

---

# Directory structure

User can change directory structure in `ms_plotter.py`; search for `User may wish to edit directory structure`. 

Directory structures depend on `outdir` `BALROG_RUN` `MATCH_TYPE`, which are defined below.

`outdir` is supplied at the command line.

`BALROG_RUN` `MATCH_TYPE` are defined in `ms_plotter.py` as follows:

If `/data/des71.a/data/kuropat/des2247-4414_sof/y3v02/` is issued at the command line, `BALROG_RUN=des2247-4414_sof`.

If `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2 = 'gal_truth', 'sof', True, True` then `MATCH_TYPE=inj_gal_truth_cat_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.

**Matched catalogs**

Matched catalogs are saved in: `outdir/catalog_compare/BALROG_RUN/MATCH_TYPE/`


**Log Files**

Log files are saved in: `outdir/log_files/BALROG_RUN/MATCH_TYPE/`


**Plots**

Plots are saved in: `outdir/plots/BALROG_RUN/MATCH_TYPE/{tile}/{plot_type}/{realization}/`

Allowed values for `{plot_type}`: `normalized` `scatter`.

Allowed values for `{realization}`: `0` `1` ... `stacked`.


**Defaults**

...

___

**Stacking multiple realizations**

Matching is performed first, then catalogs are stacked.

___


