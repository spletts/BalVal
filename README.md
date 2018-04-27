# BalVal

Conducts some [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`fof_matcher` and `par.py` analyse FOF groups. `fof_matcher` uses STILTS.

`ms_plotter.py` calls `ms_matcher` or `fof_matcher` (which calls `par.py`), analyses the matched catalog, and produces various magnitude versus Delta-magnitude plots.

___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realizations tiles`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/ /BalVal/ 0,1 DES2247-4414`

`all` is an allowed value for both `realizations` and `tiles`. If `all` is used ensure that `ALL_TILES` and `ALL_REALIZATIONS` are set correctly in `ms_plotter.py`. Alternatively, one can list realizations at command line with commas and no spaces separating the realizations (similarly for tiles).

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment `NOTICE` to remove prompt. 

User sets plot attributes and catalog attributes within `ms_plotter.py`.

**Plot attributes**

Examples of plot attributes, which are explained in `ms_plotter.py` are:

`NORMALIZE` `HEXBIN` `CM_T_S2N_COLORBAR` `CM_T_COLORBAR`  `PLOT_1SIG` (Booleans).


**Catalog attributes**

Catalogs are specified using:

`MATCH_CAT1` `MATCH_CAT2` with allowed values: `mof` `sof` `star_truth` `gal_truth` `coadd`.

`INJ1` `INJ2` (Booleans) which describe if the catalog is Balrog-injected `INJ=True` or non-injected/original `INJ=False`.


**FOF analysis**

FOF groups changed after Balrog-injection are plotted if `RUN_TYPE='rerun'`. 

FOF groups *un*changed after Balrog-injection are plotted examined if `RUN_TYPE='ok'`

If `RUN_TYPE='rerun'` or `RUN_TYPE=ok` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.

To ignore FOF analysis set `RUN_TYPE=None`. 

**Region files**

Region files are created if `MAKE_REG=True`. 

**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`. Similarly for `fof_matcher`.

---

# Directory structure

If `NO_DIR_MAKE=True` the directories below will be created if they do not already exist.

If `NO_DIR_EXIT=True` the script will exit when a nonexistent directory is encountered to remind the user to change the directory structure to their preference. User can change directory structure in `ms_plotter.py`; search for `User may wish to edit directory structure`. 

Directory stucture depends on `outdir` `BALROG_RUN` `MATCH_TYPE`, which are defined below.

`outdir` is supplied at the command line.

`BALROG_RUN` `MATCH_TYPE` are defined in `ms_plotter.py` as follows:

If `/data/des71.a/data/kuropat/des2247-4414_sof/y3v02/` is issued at the command line, `BALROG_RUN=des2247-4414_sof`.

If `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2 = 'gal_truth', 'sof', True, True` then `MATCH_TYPE=inj_gal_truth_cat_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.

**Matched catalogs**

Matched catalogs are saved in: `/outdir/catalog_compare/BALROG_RUN/MATCH_TYPE/`

Matched catalogs used for FOF analysis are saved in: `/outdir/'fof_analysis_catalog_compare'/`


**Log files**

Log files are saved in: `/outdir/log_files/BALROG_RUN/MATCH_TYPE/`

Log files for `ok` and `rerun` FOF groups are saved in: `/outdir/log_files/BALROG_RUN/MATCH_TYPE/'fof_analysis'`


**Plots**

Plots are saved in: `/outdir/plots/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/{plot_type}/`

Plots for `ok` and `rerun` FOF groups are saved in: `/outdir/plots/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/'fof_analysis'/{plot_type}/`

Allowed values for `{realization}`: `0` `1` ... `stacked`.

Allowed values for `{plot_type}`: `normalized` `scatter`.


**Region files**

Region files are saved in: `/outdir/'region_files'/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/`

Region files for `ok` and `rerun` FOF groups are saved in: `/outdir/'region_files'/BALROG_RUN/MATCH_TYPE/{tile}/{realization}/'fof_analysis'/`


**Defaults**

...

___

**Stacking multiple realizations**

Matching is performed first, then catalogs are stacked.

___


