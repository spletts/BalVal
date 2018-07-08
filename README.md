# BalVal

Conducts various [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`ms_fof_matcher` and `ms_par.py` analyse FOF groups. `ms_fof_matcher` uses STILTS.

`ms_plotter.py` calls `ms_matcher` or `ms_fof_matcher` (which calls `ms_par.py`), analyses the matched catalog, and produces various comparison plots.

Note that `ms_matcher` and `ms_fof_matcher` return CSV file formats. Further, arrays in FITS files of form `(1 2 3 4)` are converted to strings of form `'(1, 2, 3, 4)'` in matched catalogs. Similarly for matrices.
___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realizations tiles`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /BalVal/ 0,1 DES2247-4414`

`None` is an allowed value for `cmd_line_realizations`. One can list realizations at the command line with commas and _no_ spaces separating the realizations (similarly for tiles).

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment `NOTICE` to remove prompt. 

User sets plot attributes and catalog attributes within `ms_plotter.py`. A table of user-set attributes is below.


# Table of Constants

These constants are set within `ms_plotter.py`. Docstrings describe parameters but a [Table of Function Parameters](https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing) is also available.

Parameter(s) | Type | Description & _allowed values_ (if Type not bool)
:---: | :---: | ---
|`MATCH_CAT1` `MATCH_CAT2` | str | Type of catalogs to analyse. _Allowed values_:  `mof` `sof` `star_truth` `gal_truth` `y3_gold` `coadd`
| `INJ1` `INJ2` | bool | If `True` then `MATCH_CAT1` `MATCH_CAT2` are Balrog-injected. If `False` then `MATCH_CAT1` `MATCH_CAT2` are base catalogs.
|`INJ1_PERCENT` `INJ2_PERCENT` | int | Injection percent for `MATCH_CAT1` `MATCH_CAT2`. It is important to specify these parameters even in the case of `INJ1=False` `INJ2=False` because these parameters impact directory structure.
| `PLOT_MAG` | bool | If `True` plots of g-, r-, i-, and z-band magnitude are created.
| `PLOT_COLOR` | bool | If `True` colors g-r, r-i, and i-z are plotted. Creates a 2x2 grid with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24)). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for g-r, etc).
| `PLOT_FLUX` | bool | If `True` a 1D histogram of Delta_Flux/Sigma_Flux is plotted along with a standard Gaussian (mean=0, standard_deviation=1) and a fit to the 1D histogram. Delta_Flux is computed as truth_flux minus {SOF/MOF/coadd}\_flux. Sigma_Flux is computed using only the measured catalog (SOF/MOF/coadd).
| `GAUSS_APER` | bool | If `True` and `PLOT_FLUX=True` a Gaussian aperture method is used to measure the flux.
| `TRIM_FLUX` | bool | If `True` and `PLOT_FLUX=True` histograms of Delta_Flux/Sigma_Flux (even those created using `GAUSS_APER`) include the 2<sup>nd</sup>-98<sup>th</sup> percentiles.
| `SAVE_PLOT` | bool | If `True` plot is saved using name assigned by `ms_plotter.get_plot_save_name()`.
| `SHOW_PLOT` | bool | If `True` plot is displayed after it is created.
| `NORMALIZE` | bool | If `True` the magnitude plot is normalized according to the *measured* 1sigma magnitude error.
| `HIST_2D` | bool | If `True` a `matplotlib.pyplot` 2D histogram is plotted.
| `CORNER_HIST_2D` | bool | If `True` `corner.hist2d` plots are created using [corner.py](https://github.com/dfm/corner.py).
| `PLOT_DELTA_VAX` | bool | If `True` a difference is plotted on the vertical axis to produce a plot of form `x` versus `x-y`. If `False` a plot of form `x` versus `y` is produced.
| `PLOT_COMPLETENESS` | bool | If `True` a 1x2 plot grid is produced with 10% injection {magnitude/color} completeness plot and 20% {magnitude/color} completeness plot, respectively (provided both injections are available for the particular `BALROG_RUN`). `PLOT_MAG` and `PLOT_COLOR` determines if the completeness plot displays magnitude or color completeness.
| `SCATTER` | bool | If `True` a `matplotlib.pyplot.scatter()` plot is produced.
|`HEXBIN` | bool | If `True` a density plot via `matplotlib.pyplot.hexbin()` is produced.
|`CM_T_COLORBAR` | bool | If `True` and `PLOT_MAG` a colorbar is added to the plot according to the *measured* cm_T. `NORMALIZE` must be False.
| `CM_T_ERR_COLORBAR` | bool | If `True` and `PLOT_MAG` a colorbar is added to the plot according to the *measured* cm_T error. `NORMALIZE` must be False.
| `BIN_CM_T_S2N` | bool | If `True` and `PLOT_MAG` the *measured* cm_T signal-to-noise is binned using `[0, 1, 9, 20, max(cm_t_s2n)]`
| `PLOT_MAG_ERR` | bool | If `True` and and `PLOT_MAG` the 1sigma magnitude error curve is plotted. `NORMALIZE` must be False.
| `MAG_YLOW`, `MAG_YHIGH` | int, float or `None` | Limits for the vertical axis of the magnitude plot. `None` results in default scaling.
| `STACK_REALIZATIONS` | bool | If `True` catalogs are matched then stacked. One stacked realization catalog is produced per tile. Plotting resumes with stacked catalog. Must be used with, for example, `cmd_line_realizations=0,1,2` at the command line.
| `STACK_TILES` | bool | If `True` catalogs are matched then stacked. One stacked tile catalog is produced per realization. Must be used with, for example, `tiles=DES0220-0207,DES0222+0043` at the command line.
| `PLOT_DIFF_ON_VAX` | bool | If `True` a difference is plotted on the vertical axis of sorts. If `False` 
| `CENTER_ERR_ABT_ZERO`| bool | If `True` the plot of the magnitude error is centered about zero. This (minorly) affects the number of objects within 1sigma_mag. If `False` the plot of the magnitude error is centered about the median of the vertical axis data each bin.
| `PLOT_68P` | bool | Considered if `NORMALIZE=True`. If `True` the 68th percentile of the vertical axis data in each bin are plotted. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `PLOT_34P_SPLIT` | bool | Considered if `NORMALIZE=True`. If `True` the 34th percentile of the positive and negative vertical axis data in each bin are plotted separately. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `SUBPLOT` | bool | If `True` four subplots are created in a 2x2 grid. If `False` plots are created individually.
| `MOF` | bool | Only used if `RUN_TYPE` is not `None`. Does `BASEPATH` entered at command line contain MOF (`MOF=True` or SOF `MOF=False` catalogs?
| `MAKE_REG`| bool | If `True`, three DS9 region files created containing 1) objects in both catalogs, 2) objects in the first not the second catalog, 3) objects in the second not the first catalog.
| `NO_DIR_MAKE` | bool| If `True` nonexistent directories will be created. If `False`, `sys.exit()` is invoked when nonexistent directories are encountered.
| `SWAP_HAX` | bool | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis.
| `RUN_TYPE` | str | _Allowed values_: `None` `'ok'` `'rerun'`. `'ok'`: FOF groups *un*changed after Balrog-injection. `'rerun'`: FOF groups changed after Balrog-injection. `None`: FOF analysis not conducted. If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
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

If `/data/des71.a/data/kuropat/des2247-4414_sof/` is issued at the command line, `BALROG_RUN=des2247-4414_sof`.

If `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2, INJ1_PERCENT, INJ2_PERCENT = 'gal_truth', 'sof', True, True`, `10`, `10` then `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.

**Matched catalogs**

Matched catalogs are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/catalog_compare/`

Matched catalogs used for FOF analysis are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/fof_analysis_catalog_compare/`


**Log files**

Log files are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/log_files/`

Log files for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/log_files/fof_analysis/`


**Plots**

Plots are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{match_type}/{tile}/{realization}/plots/{plot_obs}/` where `plot_obs` is `color` `flux` or `magnitude`. Plot names are prefaced with the plot type (for example, 'scatter_' if `SCATTER=True`).

Plots for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/plots/fof_analysis/{plot_obs}/`

Allowed values for `{realization}`: `0` `1` ... `stack`.

Allowed values for `{tile}`: ... `stack`.

`match_type` is `MATCH_TYPE` unless `PLOT_COMPLETENESS=True` in which case both 10% and 20% Balrog-injected matched catalogs are plotted in the same window, so `match_type` removes the percent injected from `MATCH_TYPE`. That is, `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat` results in `match_type=inj_gal_truth_cat_inj_sof_cat`.

**Region files**

Region files are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/`

Region files for `ok` and `rerun` FOF groups are saved in: `/{OUTDIR}/outputs/{BALROG_RUN}/{MATCH_TYPE}/{tile}/{realization}/region_files/fof_analysis/`

___

**Stacking multiple realizations or multiple tiles**

Matching is performed first, then catalogs are stacked.

___


# `flagged* objects`
Docstrings in `ms_plotter.py` frequently make reference to `flagged* objects` which refer to specific flag cuts employed in `ms_plotter.get_good_index_using_primary_flags()` but described below for clarity.

In general, the following catalog headers are considered: 'flags' and 'cm_flags'. Exceptions are listed below. 

Star truth catalogs (`star_truth`) do not contain any flags.

Coadd catalogs (`coadd`) have a 'FLAGS' header but do not have a 'cm_flags' header (or equivalent), so only 'FLAGS' is used.

Y3 catalogs (`y3_gold`) are checked for 'FLAGS_GOLD' (replaced 'flags') and '{sof/mof}\_cm_flags' replaces 'cm_flags'. In addition, the additional flags are examined: 'SEXTRACTOR_FLAGS_{GRIZ}', 'IMAFLAGS_ISO_{GRIZ}', and, if a Y3 catalog is compared to a MOF catalog, 'MOF_CM_MOF_FLAGS'.
