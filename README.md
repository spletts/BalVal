# BalVal

Conducts some [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`ms_fof_matcher` and `ms_par.py` analyse FOF groups. `ms_fof_matcher` uses STILTS.

`ms_plotter.py` calls `ms_matcher` or `ms_fof_matcher` (which calls `ms_par.py`), analyses the matched catalog, and produces various plots.

___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realizations tiles`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /BalVal/ 0,1 DES2247-4414`

`None` is an allowed value for `realizations`. `all` is an allowed value for both `realizations` and `tiles`. If `all` is used ensure that `ALL_TILES` and `ALL_REALIZATIONS` are set correctly in `ms_plotter.py`. Alternatively, one can list realizations at the command line with commas and _no_ spaces separating the realizations (similarly for tiles).

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment `NOTICE` to remove prompt. 

User sets plot attributes and catalog attributes within `ms_plotter.py`. A table of user-set attributes is below.

**Plot and catalog attributes**

Parameter(s) | Type | Allowed values (if Type not bool) | Description
--- | --- | --- | ---
|`MATCH_CAT1` `MATCH_CAT2` | str | `mof` `sof` `star_truth` `gal_truth` `y3_gold` `coadd`  | Type of catalogs to analyse
|`INJ1` `INJ2` | bool | | Are `MATCH_CAT1` `MATCH_CAT2` Balrog-injected?  If `realizations=None` then the following is forced: `INJ1, INJ2 = False, False`
|`INJ1_20PERCENT` `INJ2_20PERCENT` | bool | 
| `PLOT_COLOR` | bool | | If `True` colors g-r, r-i, and i-z are plotted. If `False` magnitudes are plotted. `PLOT_COLOR` creates a 2x2 subplot with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24)). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for g-r, etc). 
| `RUN_TYPE` | str | `None` `'ok'` `'rerun'` | `'ok'`: FOF groups *un*changed after Balrog-injection. `'rerun'`: FOF groups changed after Balrog-injection. `None`: FOF analysis not conducted. If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
| `NORMALIZE` | bool | | Normalize plot to 1-sigma magnitude error? Error calculation uses measured catalogs only.
| `HIST_2D` | bool | | Plot a 2D histogram?
| `CORNER_HIST_2D` | bool | | If `True` `corner.hist2d` plots are created using [corner.py](https://github.com/dfm/corner.py).
| `PLOT_DELTA_VAX` | bool | | If `True` a difference is plotted on the vertical axis to produce a plot of form `x` versus `x-y`. If `False` a plot of form `x` versus `y` is produced.
| `SCATTER` | bool | | If `True` a scatter plot is produced.
|`HEXBIN` | bool | | If `True` a density plot via `hexbin()` is produced.
|`CM_T_S2N_COLORBAR` | bool | | If `True` a colorbar that displays cm_T signal-to-noise is added to the scatter plot. `SCATTER` must be `True`...
|`CM_T_COLORBAR` | bool | | Plot a colorbar according to cm_T?
| `CM_T_ERR_COLORBAR` | bool | | Plot a colorbar according to cm_T error?
| `BIN_CM_T_S2N` | bool | | Bin cm_T signal-to-noise? Default bins are `[0, 1, 9, 20, max(cm_t_s2n)]`
| `PLOT_1SIG` | bool | | Plot the 1-sigma curve? Errors refer to magnitude errors
| `YLOW` `YHIGH` | int or float | `None` and any real number | Limits for the vertical axis of plot. `None` results in default scaling
| `STACK_REALIZATIONS` | bool | | If `True` catalogs are matched then stacked. Plotting resumes with stacked catalog. Must be used with `realizations=all` at command line...
| `PLOT_68P` | bool | | Only considered if `NORMALIZE=True`. Plot the 68th percentile of the vertical axis centered about zero? Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `PLOT_34P_SPLIT` | bool | | Only considered if `NORMALIZE=True`. Plot the 34th percentile of the positive and negative vertical axis? Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `SUBPLOT` | bool | | If `True` four subplots (one for each griz filter) are created in a 2x2 grid. If `False` plots are made individually
| `MOF` | bool | | Only used if `RUN_TYPE` is not `None`. Does `BASEPATH` entered at command line contain MOF (`MOF=True` or SOF `MOF=False` catalogs?
| `MAKE_REG`| bool | | If `True`, three DS9 region files created: 1) objects in both catalogs, 2) objects in the first not second catalog, 3) objects in second not first  
| `NO_DIR_MAKE` | bool| | If `True` nonexistent dirs will be created. If `False`, `sys.exit()` invoked when nonexistent dirs encountered.
| `SWAP_HAX` | bool | | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis.
| `SAVE_PLOT` | bool | | Save plot?
| `SHOW_PLOT` | bool | | Show plot?
| `EH_CUTS` | bool | | Apply quality cuts introduced by Eric Huff?
| `overwrite` | bool | | Exists within the following: `ms_plotter.matcher()` `ms_plotter.fof_matcher()` `ms_plotter.make_region_files()`. If `True` existing matched catalogs or region files are overwritten.


**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`. Similarly for `fof_matcher`.

---

# Directory structure

If `NO_DIR_MAKE=True` the directories below will be created if they do not already exist.

If `NO_DIR_MAKE=False` the script will exit when a nonexistent directory is encountered to remind the user to change the directory structure to their preference. User can change directory structure in `ms_plotter.py`; search for `User may wish to edit directory structure`. 

Directory stucture depends on `OUTDIR` `BALROG_RUN` `MATCH_TYPE`, which are defined below.

`OUTDIR` is supplied at the command line.

`BALROG_RUN` `MATCH_TYPE` are defined in `ms_plotter.py` as follows:

If `/data/des71.a/data/kuropat/des2247-4414_sof/y3v02/` is issued at the command line, `BALROG_RUN=des2247-4414_sof`.

If `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2 = 'gal_truth', 'sof', True, True` then `MATCH_TYPE=inj_gal_truth_cat_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.

**Matched catalogs**

Matched catalogs are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/catalog_compare/`

Matched catalogs used for FOF analysis are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/fof_analysis_catalog_compare/`


**Log files**

Log files are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/log_files/`

Log files for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/log_files/fof_analysis/`


**Plots**

Plots are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/plots/{plot_type}/`

Plots for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/plots/fof_analysis/{plot_type}/`

Allowed values for `{realization}`: `0` `1` ... `stacked`.

Allowed values for `{plot_type}`: `normalized` `scatter`.


**Region files**

Region files are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/`

Region files for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/fof_analysis/`

___

**Stacking multiple realizations**

Matching is performed first, then catalogs are stacked.

___


# `Flags*`
Docstrings in `ms_plotter.py` frequently make reference to `flags*` which are described here...
