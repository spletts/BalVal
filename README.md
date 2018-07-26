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

**Dependencies**

`ms_plotter.py` requires [`ngmixer`](https://github.com/esheldon/ngmixer) to measure fluxes using a Gaussian aperture. 

If user has access to DES machines at FNAL, user must ssh to `@des70` or `@des71` (for `matplotlib v2.2.2`) and run the following setup script: `$source /home/s1/mspletts/setup_ngmixer_gaussap.sh`.

Note that `/home/s1/mspletts/setup_ngmixer_gaussap.sh` points to `ngmixer` as installed in `/home/s1/mspletts`, where minor changes have been made to `ngmixer.gaussap.get_gauss_aper_flux_cat()`, as follows:

~~`output['id']=cat['id']`~~
```
try:
    output['id']=cat['id']
except:
    output['id']=cat['coadd_objects_id']
```



# Table of Constants

These constants are set within `ms_plotter.py`. Docstrings describe parameters but a [Table of Function Parameters](https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing) is also available.

Parameter(s) | Type | Description & _allowed values_ (if Type not bool)
:---: | :---: | ---
|`MATCH_CAT1` `MATCH_CAT2` | str | Type of catalogs to analyse. _Allowed values_:  `coadd` `gal_truth` `mof` `sof` `star_truth` `y3_gold_2_0` `y3_gold_2_2`
| `INJ1` `INJ2` | bool | If `True` then `MATCH_CAT1` `MATCH_CAT2` are Balrog-injected. If `False` then `MATCH_CAT1` `MATCH_CAT2` are base catalogs.
|`INJ1_PERCENT` `INJ2_PERCENT` | int | Injection percent for `MATCH_CAT1` `MATCH_CAT2`. It is important to specify these parameters even in the case of `INJ1=False` `INJ2=False` because these parameters impact directory structure.
| `PLOT_MAG` | bool | If `True` plots of g-, r-, i-, and z-band magnitude are created.
| `PLOT_COLOR` | bool | If `True` colors g-r, r-i, and i-z are plotted. Creates a 2x2 grid with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24)). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for g-r, etc).
| `PLOT_FLUX` | bool | If `True` a 1D histogram of Delta_Flux/Sigma_Flux is plotted along with a standard Gaussian (mean=0, standard_deviation=1) and a fit to the 1D histogram. Delta_Flux is computed as truth_flux minus {SOF/MOF/coadd}\_flux. Sigma_Flux is computed using only the measured catalog (SOF/MOF/coadd).
| `GAUSS_APER` | bool | If `True` and `PLOT_FLUX=True` a Gaussian aperture method is used to measure the flux.
| `SIGMA_CLIP_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the normalized (to measured flux error) flux differences are sigma clipped to `N`-sigma.
| `N` | float | Used if `PLOT_FLUX=True` and `SIGMA_CLIP_NORM_FLUX_DIFF=True` to perform `N`-sigma clip.
| `NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY` | bool | Used if `PLOT_FLUX=True`. If `True` histograms are normalized (meaning area under the curve is 1).
| `RAW_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the normalized (to measured flux error) flux differences are plotted within percentile trimming nor sigma clipping. These normalized (to measured flux error) flux differences can be from Gaussian aperture measurements if `GAUSS_APER=True`.
| `TRIM_NORM_FLUX_DIFF` | bool | If `True` and `PLOT_FLUX=True` histograms of Delta_Flux/Sigma_Flux (even those created using `GAUSS_APER`) include the 2<sup>nd</sup>-98<sup>th</sup> percentiles.
| `NORMALIZE` | bool | Used if `PLOT_MAG=True`. If `True` the magnitude plot is normalized according to the *measured* 1sigma magnitude error.
| `HIST_2D` | bool | Used if `PLOT_MAG=True`. If `True` a `matplotlib.pyplot` 2D histogram is plotted.
| `CORNER_HIST_2D` | bool | Used if `PLOT_MAG=True` or Used if `PLOT_COLOR=True`. If `True` `corner.hist2d` plots are created using [corner.py](https://github.com/dfm/corner.py).
| `PLOT_DELTA_VAX` | bool | If `True` a difference is plotted on the vertical axis to produce a plot of form `x` versus `x-y`. If `False` a plot of form `x` versus `y` is produced.
| `PLOT_COMPLETENESS` | bool | Used if `PLOT_MAG=True`. If `True` a 1x2 plot grid is produced with 10% injection {magnitude/color} completeness plot and 20% {magnitude/color} completeness plot, respectively (provided both injections are available for the particular `BALROG_RUN`). `PLOT_MAG` and `PLOT_COLOR` determines if the completeness plot displays magnitude or color completeness.
| `SCATTER` | bool | Used if `PLOT_MAG=True`. If `True` a `matplotlib.pyplot.scatter()` plot is produced.
|`HEXBIN` | bool | Used if `PLOT_MAG=True`. If `True` a density plot via `matplotlib.pyplot.hexbin()` is produced.
|`CM_T_CBAR` | bool | Used if `PLOT_MAG=True`. If `True` a colorbar is added to the plot according to the *measured* cm_T. `NORMALIZE` must be False.
| `CM_T_ERR_CBAR` | bool | Used if `PLOT_MAG=True`. If `True` a colorbar is added to the plot according to the *measured* cm_T error. `NORMALIZE` must be False.
| `PLOT_MAG_ERR` | bool | If `True` and and `PLOT_MAG` the 1sigma magnitude error curve is plotted. `NORMALIZE` must be False.
| `MAG_YLOW`, `MAG_YHIGH` | int, float or `None` | Limits for the vertical axis of the magnitude plot. `None` results in default scaling.
| `SAVE_PLOT` | bool | If `True` plot is saved using name assigned by `ms_plotter.get_plot_save_name()`.
| `SHOW_PLOT` | bool | If `True` plot is displayed after it is created.
| `STACK_REALIZATIONS` | bool | If `True` catalogs are matched then stacked. One stacked realization catalog is produced per tile. Plotting resumes with stacked catalog. Must be used with, for example, `cmd_line_realizations=0,1,2` at the command line.
| `STACK_TILES` | bool | If `True` catalogs are matched then stacked. One stacked tile catalog is produced per realization. Must be used with, for example, `tiles=DES0220-0207,DES0222+0043` at the command line.
| `PLOT_DIFF_ON_VAX` | bool | If `True` a difference is plotted on the vertical axis of sorts. If `False` 
| `CENTER_ERR_ABT_ZERO`| bool | If `True` the plot of the magnitude error is centered about zero. This (minorly) affects the number of objects within 1sigma_mag. If `False` the plot of the magnitude error is centered about the median of the vertical axis data each bin.
| `PLOT_68P` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE=True`. If `True` the 68th percentile of the vertical axis data in each bin are plotted. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `PLOT_34P_SPLIT` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE=True`. If `True` the 34th percentile of the positive and negative vertical axis data in each bin are plotted separately. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `ms_plotter.normalized_delta_magnitude_plotter()`.
| `SUBPLOT` | bool | If `True` four subplots are created in a 2x2 grid. If `False` plots are created individually.
| `MOF` | bool | Only used if `RUN_TYPE` is not `None`. Does `BASEPATH` entered at command line contain MOF (`MOF=True` or SOF `MOF=False` catalogs?
| `MAKE_REG`| bool | If `True`, three DS9 region files created containing 1) objects in both catalogs, 2) objects in the first not the second catalog, 3) objects in the second not the first catalog.
| `NO_DIR_MAKE` | bool| If `True` nonexistent directories will be created. If `False`, `sys.exit()` is invoked when nonexistent directories are encountered.
| `SWAP_HAX` | bool | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis.
| `RUN_TYPE` | str | _Allowed values_: `None` `'ok'` `'rerun'`. `'ok'`: FOF groups *un*changed after Balrog-injection. `'rerun'`: FOF groups changed after Balrog-injection. `None`: FOF analysis not conducted. If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
| `VERBOSE_ING` | bool | If `True` status and progress of script will print to screen. For example, 'Plotting...', 'calculating...', 'overwriting...', etc.
| `VERBOSE_ED` | bool | If `True` results will print to screen. For example, 'got...', 'flagged...', 'checked...', etc. Note that many of these printouts are also saved in a log file.

The following functions contain `__overwrite`: `get_coadd_catalog_for_matcher()`, `make_ngmixer_gaussap_compatible_catalog()`, `get_and_reformat_base_catalog()`, `matcher()`, `fof_matcher()`, `stack_tiles()`, `stack_realizations()`, `get_coadd_catalog_for_matcher()`, `write_to_region_files()`. 
If `__overwrite=True` a `raw_input()` prompt will appear to ensure that files are to be overwritten. In these instances, press 'Enter' to continue and 'Control+c' to stop the process.

Log files, matched catalogs, and plot names are printed with a proceeding `----->` for ease of finding and opening these files.


**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`. Similarly for `ms_fof_matcher`.

---

# Condensed directory structure

Below are solely the directories and not the files that will be created within them.

```
{OUTDIR}
|    
+---- outputs
    |  
    +---- {BALROG_RUN}
        |  
        +---- {match_type}
            |
            +---- {tile}
                |
                +---- {realization}
                    |
                    +---- catalog_compare
                    |   |
                    |   +---- fof_analysis
                    |
                    +---- log_files
                    |   |
                    |   +---- fof_analysis
                    |
                    +---- plots
                    |   |
                    |   +---- color
                    |   |
                    |   +---- flux
                    |   |
                    |   +---- magnitude
                    |   |   |
                    |   |   +---- normalized
                    |   |
                    |   +---- fof_analysis
                    |       |
                    |       +---- color
                    |       |
                    |       +---- flux
                    |       |
                    |       +---- magnitude
                    |           |
                    |           +---- normalized
                    |
                    +---- region_files
                        |
                        +---- fof_analysis

```

If `NO_DIR_MAKE=True` the directories below will be created if they do not already exist.

If `NO_DIR_MAKE=False` the script will exit when a nonexistent directory is encountered to remind the user to change the directory structure to their preference. User can change directory structure in `ms_plotter.py`; search for `User may wish to edit directory structure`. 

Directory stucture depends on `OUTDIR` `BALROG_RUN` `MATCH_TYPE`, which are defined below.

`OUTDIR` is supplied at the command line.

`BALROG_RUN` `MATCH_TYPE` are defined in `ms_plotter.py` as follows:

If `/data/des71.a/data/kuropat/des2247-4414_sof/` is issued at the command line, `BALROG_RUN=des2247-4414_sof`.

If `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2, INJ1_PERCENT, INJ2_PERCENT = 'gal_truth', 'sof', True, True`, `10`, `10` then `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.


# Full directory structure 

```
{OUTDIR}
|    
+---- outputs
    |  
    +---- {BALROG_RUN}
        |  
        +---- {match_type}
            |
            +---- {tile}
                |
                +---- {realization}
                    |
                    +---- catalog_compare
                    |   |   {tile}_{realization}_{match_type}_match1and2.csv
                    |   |   {tile}_{realization}_{match_type}_match1not2.csv
                    |   |   {tile}_{realization}_{match_type}_match2not1.csv
                    |   |
                    |   +---- fof_analysis
                    |           {tile}_num_match_fof_coadd.csv
                    |           {tile}_fofgroups.csv
                    |           {tile}_{realization}_num_match_inj_fof_inj_coadd.csv
                    |           {tile}_{realization}_inj_fofgroups.csv
                    |           {tile}_{realization}_inj_fofgroup_fofgroup_match1and2.csv
                    |           {tile}_{realization}.ok
                    |           {tile}_{realization}.rerun
                    |           {tile}_{realization}_ok_inj_mof.csv
                    |           {tile}_{realization}_rerun_inj_mof.csv
                    |           {tile}_{realization}_ok_mof.csv
                    |           {tile}_{realization}_rerun_mof.csv
                    |           {tile}_{realization}_ok_inj_mof_ok_mof_match1and2.csv
                    |           {tile}_{realization}_ok_inj_mof_ok_mof_match1not2.csv
                    |           {tile}_{realization}_ok_inj_mof_ok_mof_match2not1.csv
                    |           {tile}_{realization}_rerun_inj_mof_rerun_mof_match1and2.csv
                    |           {tile}_{realization}_rerun_inj_mof_rerun_mof_match1not2.csv
                    |           {tile}_{realization}_rerun_inj_mof_rerun_mof_match2not1.csv
                    |
                    +---- log_files
                    |   |   {tile}_{realization}_{MATCH_TYPE}_color_from_mag.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_flags.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_flux_from_gaussian_aperture.log
                    |   |   {tile}_{realization}_{modified_match_type}_mag_completeness.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_mag_error_computation.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_matched_catalogs.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_main.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_outlier_mag_diffs.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_sigma_clip_flux.log
                    |   |
                    |   +---- fof_analysis
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_color_from_mag.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_flags.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_flux_from_gaussian_aperture.log
                    |           {tile}_{realization}_{mod_match_type}_{RUN_TYPE}_mag_completeness.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_mag_error_computation.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_matched_catalogs.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_main.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_outlier_mag_diffs.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_sigma_clip_flux.log
                    |
                    +---- plots
                    |   |
                    |   +---- color
                    |   |       {tile}_{realization}_{match_type}_{color}_{plot_type}_{ylim}.png
                    |   |       {tile}_{realization}_{match_type}_{color}_corner2dhist_.png
                    |   |
                    |   +---- flux
                    |   |       {tile}_{realization}_{match_type}_griz_{plot_type}_.png
                    |   |       {tile}_{realization}_{match_type}_griz_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{match_type}_griz_gauss_aper_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{match_type}_griz_{N}sigma_clip_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{match_type}_griz_{N}sigma_clip_gauss_aper_norm_flux_diff_histogram_.png
                    |   |
                    |   +---- magnitude
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_{plot_type}_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_cbar_cm_t_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_cbar_cm_t_err_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_completeness_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_cornerhist2d_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_hexbin_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_hist2d_{ylim}.png
                    |   |   |   {tile}_{realization}_{match_type}_{band(s)}_scatter_{ylim}.png
                    |   |   |
                    |   |   +---- normalized
                    |   |           norm_"
                    |   |
                    |   +---- fof_analysis
                    |       |
                    |       +---- color
                    |       |
                    |       +---- flux
                    |       |
                    |       +---- magnitude
                    |           |
                    |           +---- normalized
                    |
                    +---- region_files
                        |   {tile}_{realization}_{MATCH_TYPE}_match1and2.reg
                        |   {tile}_{realization}_{MATCH_TYPE}_match1not2.reg
                        |   {tile}_{realization}_{MATCH_TYPE}_match2not1.reg
                        |
                        +---- fof_analysis
                              {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_match1and2.reg
                              {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_match1not2.reg
                              {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_match2not1.reg
                
```
`{realization}` and `{tile}` can be `stack`.


Log files are CSVs. Not all log files will be written to; for example, if `PLOT_MAG=False`, 'mag_completeness.log' will not be created and an empty 'dummy.log' will replace it.

Plot save names ending with '\_.png' have default axes.

`match_type` is `MATCH_TYPE` unless `PLOT_COMPLETENESS=True` in which case both 10% and 20% Balrog-injected matched catalogs are plotted in the same window, so `match_type` removes the percent injected from `MATCH_TYPE`. That is, `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat` results in `match_type=inj_gal_truth_cat_inj_sof_cat`.

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

`ngmixer.gaussap.get_gauss_aper_flux_cat()` gives flags for the Gaussian aperture flux measurement. If `PLOT_FLUX=True` and `GAUSS_APER=True`, objects with these Gaussian aperture flags are removed.
