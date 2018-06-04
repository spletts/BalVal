# BalVal

Conducts various [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`ms_fof_matcher` and `ms_par.py` analyse FOF groups. `ms_fof_matcher` uses STILTS.

`ms_plotter.py` calls `ms_matcher` or `ms_fof_matcher` (which calls `ms_par.py`), analyses the matched catalog, and produces various comparison plots.

___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realizations tiles`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /BalVal/ 0,1 DES2247-4414`

`None` is an allowed value for `realizations`. `all` is an allowed value for both `realizations` and `tiles`. If `all` is used ensure that `ALL_TILES` and `ALL_REALIZATIONS` are set correctly in `ms_plotter.py`. Alternatively, one can list realizations at the command line with commas and _no_ spaces separating the realizations (similarly for tiles).

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment `NOTICE` to remove prompt. 

User sets plot attributes and catalog attributes within `ms_plotter.py`. A table of user-set attributes is below.


# Table of Parameters

These parameters are set within `ms_plotter.py`.

Parameter(s) | Type | Description & _allowed values_ (if Type not bool)
:---: | :---: | ---
|`MATCH_CAT1` `MATCH_CAT2` | str | Type of catalogs to analyse. _Allowed values_:  `mof` `sof` `star_truth` `gal_truth` `y3_gold` `coadd`
|`INJ1_10PERCENT` `INJ2_10PERCENT` | bool | If `True` then `MATCH_CAT1` `MATCH_CAT2` refer to 10% Balrog-injected catalog.
|`INJ1_20PERCENT` `INJ2_20PERCENT` | bool | If `True` then `MATCH_CAT1` `MATCH_CAT2` refer to 20% Balrog-injected catalog. If all the following are `False` then `MATCH_CAT1` `MATCH_CAT2` refer to base (non-Balrog-injected) catalog: `INJ1_10PERCENT` `INJ2_10PERCENT` `INJ1_20PERCENT` `INJ2_20PERCENT`. If `realizations=None` at the command line the following are force set to False: `INJ1_10PERCENT` `INJ2_10PERCENT` `INJ1_20PERCENT` `INJ2_20PERCENT`.
| `PLOT_COLOR` | bool | If `True` colors g-r, r-i, and i-z are plotted. If `False` magnitudes are plotted. `PLOT_COLOR` creates a 2x2 subplot with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24)). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for g-r, etc). 
| `RUN_TYPE` | str | _Allowed values_: `None` `'ok'` `'rerun'`. `'ok'`: FOF groups *un*changed after Balrog-injection. `'rerun'`: FOF groups changed after Balrog-injection. `None`: FOF analysis not conducted. If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
| `NORMALIZE` | bool | If `True` the magnitude plot is normalized according to the *measured* 1sigma magnitude error.
| `HIST_2D` | bool | If `True` a `matplotlib.pyplot` 2D histogram is plotted.
| `CORNER_HIST_2D` | bool | If `True` `corner.hist2d` plots are created using [corner.py](https://github.com/dfm/corner.py).
| `PLOT_DELTA_VAX` | bool | If `True` a difference is plotted on the vertical axis to produce a plot of form `x` versus `x-y`. If `False` a plot of form `x` versus `y` is produced.
| `PLOT_COMPLETENESS` | bool | If `True` a 1x2 plot grid is produced with 10% injection {magnitude/color} completeness plot and 20% {magnitude/color} completeness plot, respectively. `PLOT_COLOR` determines if the completeness plot displays magnitude or color completeness.
| `PLOT_FLUX_HIST` | bool | If `True` a 1D histogram of DeltaFlux/SigmaFlux is plotted along with a standard Gaussian (mean=0, standard_deviation=1) and a fit to the 1D histogram. DeltaFlux is computed as truth_flux minus {SOF/MOF/coadd}_flux.
| `SCATTER` | bool | If `True` a scatter plot is produced.
|`HEXBIN` | bool | If `True` a density plot via `hexbin()` is produced.
|`CM_T_S2N_COLORBAR` | bool | If `True` a colorbar that displays the *measured* cm_T signal-to-noise is added to the magnitude versus delta magnitude plot. `NORMALIZE` must be False.
|`CM_T_COLORBAR` | bool | If `True` a colorbar is added to the magnitude versus delta magnitude plot according to the *measured* cm_T. `NORMALIZE` must be False.
| `CM_T_ERR_COLORBAR` | bool | If `True` a colorbar is added to the magnitude versus delta magnitude plot according to the *measured* cm_T error. `NORMALIZE` must be False.
| `BIN_CM_T_S2N` | bool | If `True` the *measured* cm_T signal-to-noise is binned using `[0, 1, 9, 20, max(cm_t_s2n)]`
| `PLOT_1SIG` | bool | If `True` the 1sigma magnitude error curve is plotted. `NORMALIZE` must be False.
| `YLOW` `YHIGH` | int, float or `None` | Limits for the vertical axis of plot. `None` results in default scaling.
| `STACK_REALIZATIONS` | bool | If `True` catalogs are matched then stacked. Plotting resumes with stacked catalog. Must be used with `realizations=all` at command line...
| `CENTER_ERR_ABT_ZERO`| bool | If `True` the plot of the magnitude error is centered about zero. This (minorly) affects the number of objects within 1sigma_mag. If `False` the plot of the magnitude error is centered about the median of the vertical axis data each bin.
| `PLOT_68P` | bool | Only considered if `NORMALIZE=True`. If `True` the 68th percentile of the vertical axis data in each bin are plotted. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `PLOT_34P_SPLIT` | bool | Only considered if `NORMALIZE=True`. If `True` the 34th percentile of the positive and negative vertical axis data in each bin are plotted separately. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `SUBPLOT` | bool | If `True` four subplots are created in a 2x2 grid. If `False` plots are created individually.
| `MOF` | bool | Only used if `RUN_TYPE` is not `None`. Does `BASEPATH` entered at command line contain MOF (`MOF=True` or SOF `MOF=False` catalogs?
| `MAKE_REG`| bool | If `True`, three DS9 region files created containing 1) objects in both catalogs, 2) objects in the first not second catalog, 3) objects in second not first.
| `NO_DIR_MAKE` | bool| If `True` nonexistent directories will be created. If `False`, `sys.exit()` is invoked when nonexistent directories are encountered.
| `SWAP_HAX` | bool | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis.
| `SAVE_PLOT` | bool | If `True` plot is saved using name assigned by `ms_plotter.get_plot_save_name()`.
| `SHOW_PLOT` | bool | If `True` plot is displayed after it is created.
| `EH_CUTS` | bool | If `True` quality cuts introduced by Eric Huff are applied to clean the data.
| `overwrite` | bool | If `True` existing matched catalogs or region files are overwritten. Exists within `ms_plotter.matcher()`, `ms_plotter.fof_matcher()`, and `ms_plotter.make_region_files()`.


**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`. Similarly for `ms_fof_matcher`.

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

Allowed values for `{realization}`: `0` `1` ... `stack`.

Allowed values for `{tile}`: ... `stack`.

Allowed values for `{plot_type}`: `normalized` `scatter`.


**Region files**

Region files are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/`

Region files for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/fof_analysis/`

___

**Stacking multiple realizations or multiple tiles**

Matching is performed first, then catalogs are stacked.

___


# `Flags*`
Docstrings in `ms_plotter.py` frequently make reference to `flags*` which refers to specific flag cuts that are employed in `ms_plotter.get_good_index_using_primary_flags()` but described below for clarity.

In general, the following catalog headers are considered: 'flags' and 'cm_flags'. Exceptions are listed below. 

Star truth catalogs (`star_truth`) do not contain any flags.

Coadd catalogs (`coadd`) have a 'FLAGS' header but do not have a 'cm_flags' header (or equivalent), so only 'FLAGS' is used.

Y3 catalogs (`y3_gold`) are checked for 'FLAGS_GOLD' (replaced 'flags') and '{sof/mof}_cm_flags' replaces 'cm_flags'. In addition, the additional flags are examined: 'SEXTRACTOR_FLAGS_{GRIZ}', 'IMAFLAGS_ISO_{GRIZ}', and, if a Y3 catalog is compared to a MOF catalog, 'MOF_CM_MOF_FLAGS'.
