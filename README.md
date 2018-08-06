 # BalVal

Conducts various [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
### Table of Contents

[Running BalVal](http://github.com)

[Table of Constants](https://github.com/spletts/BalVal/blob/master/README.md#table-of-constants)

[Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#directory-structure)

  - [Condensed Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#condensed-directory-structure)

  - [Full Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#full-directory-structure)



___
# Running BalVal

General: `$python plotter.py {BASE_PATH_TO_CATS} {OUTPUT_DIRECTORY} {realization(s)} {tile(s)}`

Ex: `$python plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /BalVal/ 0,1 DES2247-4414`

`None` is an allowed value for `{realization(s)}`. User can list realizations at the command line with commas and _no_ spaces separating the realizations (similarly for tiles). `{tile(s)}` can also be specified with a `.dat` file.

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can set `NOTICE=False` in `set_constants.py` to remove this prompt. 

User sets plot attributes and catalog attributes within `plotter.py`. A table of user-set attributes is in 'Table of Constants'.

##### Dependencies

Requires [`ngmixer`](https://github.com/esheldon/ngmixer) to measure fluxes using a Gaussian aperture. 

If user has access to the DES machines at FNAL: 
1. Use `matplotlib v2.2.2`, via 
  - `$ssh {user}@des70.fnal.gov` OR
  - `$ssh {user}@des71.fnal.gov`
2. `$source /home/s1/mspletts/setup_ngmixer_gaussap.sh`

Note that `/home/s1/mspletts/setup_ngmixer_gaussap.sh` points to `ngmixer` as cloned in `/home/s1/mspletts/`, where minor changes have been made to `ngmixer.gaussap.get_gauss_aper_flux_cat()`. These changes are indicated by
```
### MS ### 
{modified code} 
### - ###
```

___
##### Contents of Repository

- `catalog_classes.py`: Each catalog type (see `MATCH_CAT1` in 'Table of Constants') and their relevant and allowed headers.

- `catalog_headers.py`: Contains headers for `MATCH_CAT1` and `MATCH_CAT2` (established in `set_constants.py`).

- `set_constants.py`: Set catalog and plot attributes. User interacts with this. 

- `stilts_matcher`: matches two catalogs on RA and Dec using [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/). Returns a CSV.

- `fof_stilts_matcher`: creates a number of catalogs needed to analyse FOF groups. Catalogs are matched using STILTS. Returns CSVs.

- `ms_par.py`: analyses FOF groups changed and unchanged after Balrog-injection. This is called by `fof_stilts_matcher`. 

- `outputs.py`: Creates and writes to various outputs.

- `analysis.py`: Conducts analysis for plots. E.g. remove flags, error analysis, ...

- `plotter.py`: Produces various comparison plots (see 'Table of Constants').

- `region_files.py`: Writes to region files.

Note that the CSVs created after matching two catalogs with STILTS have converted arrays in FITS files. For example, if a cell contains an array of form `(1 2 3 4)` it is converted to a string of form `'(1, 2, 3, 4)'`. Similarly for matrices, etc.



# Table of Constants

These constants are set within `set_constants.py`. Docstrings describe parameters but a [Table of Function Parameters](https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing) is also available.


Parameter(s) | Type | Description <br> `allowed values` (if Type not bool)
:---: | :---: | ---
|`MATCH_CAT1` `MATCH_CAT2` | str | Type of catalogs to analyse. <br>`coadd` `gal_truth` `mof` `sof` `star_truth` `y3_gold_2_0` `y3_gold_2_2`
| `INJ1` `INJ2` | bool | If `True` then `MATCH_CAT1` `MATCH_CAT2` are Balrog-injected. If `False` then `MATCH_CAT1` `MATCH_CAT2` are base catalogs.
|`INJ1_PERCENT` `INJ2_PERCENT` | int | These are set by `...get_injection_percent()`. User should hardcode these if `{BASE_PATH_TO_CATS}` contains catalogs with different injection percents (for example, the Balrog run for TAMU had both 10% and 20% injections).
| `PLOT_MAG` | bool | If `True` plots of g-, r-, i-, and z-band magnitude are created.
| `PLOT_COLOR` | bool | If `True` colors g-r, r-i, and i-z are plotted. Creates a 2x2 grid with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24)). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for g-r, etc). By default (`SWAP_HAX=False`) the color difference is calculated via color_from_match_cat1 - color_from_match_cat2. This order is reversed if `SWAP_HAX=True`.
| `PLOT_FLUX` | bool | If `True` a 1D histogram of Delta_Flux/Sigma_Flux is plotted along with a standard Gaussian (mean=0, standard_deviation=1) and a fit to the 1D histogram. Delta_Flux is computed as truth_flux minus {SOF/MOF/coadd}\_flux. Sigma_Flux is computed using only the measured catalog (SOF/MOF/coadd).
| `PLOT_CM_FLUX` | bool | If `True` the cm_flux is plotted. Note that this and `PLOT_GAUSS_APER_FLUX` can be `True` and both will be plotted on the same plot.
| `PLOT_GAUSSIAN_FIT` | bool | If `True` and `PLOT_FLUX=True` a fit to the distribution(s) is plotted. Distrubtions plural if `PLOT_CM_FLUX=True` and `PLOT_GAUSS_APER_FLUX=True`.
| `PLOT_PEAKS` | bool | If `True` and `PLOT_FLUX=True` a dashed vertical line representing the median of the distribution(s) is plotted. 
| `PLOT_GAUSS_APER_FLUX` | bool | If `True` and `PLOT_FLUX=True` a Gaussian aperture method is used to measure the flux.
| `SIGMA_CLIP_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the normalized (to measured flux error) flux differences are sigma clipped to `N`-sigma.
| `N` | float | Used if `PLOT_FLUX=True` and `SIGMA_CLIP_NORM_FLUX_DIFF=True` to perform `N`-sigma clip.
| `NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY` | bool | Used if `PLOT_FLUX=True`. If `True` histograms are normalized (meaning area under the curve is 1).
| `NUM_ITERS` | float OR `None` | Considered if `SIGMA_CLIP_NORM_FLUX_DIFF=True`. Number of iterations to use in `N`-sigma clipping. See [`iters`](http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html) documentation.
| `RAW_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the normalized (to measured flux error) flux differences are plotted within percentile trimming nor sigma clipping. These normalized (to measured flux error) flux differences can be from Gaussian aperture measurements if `GAUSS_APER=True`.
| `TRIM_NORM_FLUX_DIFF` | bool | If `True` and `PLOT_FLUX=True` histograms of Delta_Flux/Sigma_Flux (even those created using `GAUSS_APER`) include the 2<sup>nd</sup>-98<sup>th</sup> percentiles.
| `FLUX_XLOW` `FLUX_XHIGH` | float OR `None` | Limits for the horizontal axis of the plot created if `PLOT_FLUX=True`. If these values are `None` default `matplotlib` scaling is used. Must both be float or both be `None`.
| `COLOR_XLOW` `COLOR_XHIGH` | float OR `None` |
| `Y3_MODEL` | str | <br> `'CM'`, `'PSF'`
| `Y3_FIT` | str | Considered if `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog<br> `'MOF'`, `'SOF'`
| `NORMALIZE` | bool | Used if `PLOT_MAG=True`. If `True` the magnitude plot is normalized according to the *measured* 1sigma magnitude error.
| `HIST_2D` | bool | Used if `PLOT_MAG=True`. If `True` a `matplotlib.pyplot` 2D histogram is plotted.
| `CORNER_HIST_2D` | bool | Used if `PLOT_MAG=True` or Used if `PLOT_COLOR=True`. If `True` `corner.hist2d` plots are created using [`corner.py`](https://github.com/dfm/corner.py).
| `PLOT_COMPLETENESS` | bool | Used if `PLOT_MAG=True`. If `True` a 1x2 plot grid is produced with 10% injection {magnitude/color} completeness plot and 20% {magnitude/color} completeness plot, respectively (provided both injections are available for the particular `BALROG_RUN`). `PLOT_MAG` and `PLOT_COLOR` determines if the completeness plot displays magnitude or color completeness.
| `SCATTER` | bool | Used if `PLOT_MAG=True`. If `True` a `matplotlib.pyplot.scatter()` plot is produced.
|`HEXBIN` | bool | Used if `PLOT_MAG=True`. If `True` a density plot via `matplotlib.pyplot.hexbin()` is produced.
|`CM_T_CBAR` | bool | Used if `PLOT_MAG=True`. If `True` a colorbar is added to the plot according to the *measured* cm_T. `NORMALIZE` must be False.
| `CM_T_ERR_CBAR` | bool | Used if `PLOT_MAG=True`. If `True` a colorbar is added to the plot according to the *measured* cm_T error. `NORMALIZE` must be False.
| `PLOT_MAG_ERR` | bool | If `True` and and `PLOT_MAG` the 1sigma magnitude error curve is plotted. `NORMALIZE` must be False.
| `MAG_YLOW`, `MAG_YHIGH` | int, float or `None` | Limits for the vertical axis of the magnitude plot. `None` results in default scaling. Must both be float or both be `None`.
| `COLOR_YLOW`, `COLOR_YHIGH` | float OR `None` | Limits for the vertical axis of the plot created if `PLOT_COLOR=True`. If `False`, semi-default `matplotlib` scaling is used, with symmetric limits forced. Must both be float or both be `None`.
| `SAVE_PLOT` | bool | If `True` plot is saved using name assigned by `outputs.get_plot_filename()`.
| `SHOW_PLOT` | bool | If `True` plot is displayed after it is created.
| `STACK_REALIZATIONS` | bool | If `True` catalogs are matched then stacked. One stacked realization catalog is produced per tile. Plotting resumes with stacked catalog. Must be used with, for example, `0,1,2` at the command line.
| `STACK_TILES` | bool | If `True` catalogs are matched then stacked. One stacked tile catalog is produced per realization. Must be used with, for example, `DES0220-0207,DES0222+0043` at the command line.
| `CENTER_ERR_ABT_ZERO`| bool | If `True` the plot of the magnitude error is centered about zero. This (minorly) affects the number of objects within 1sigma_mag. If `False` the plot of the magnitude error is centered about the median of the vertical axis data each bin.
| `PLOT_68P` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE=True`. If `True` the 68th percentile of the vertical axis data in each bin are plotted. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `analysis.normalized_delta_magnitude_plotter()`.
| `PLOT_34P_SPLIT` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE=True`. If `True` the 34th percentile of the positive and negative vertical axis data in each bin are plotted separately. Bins refer to the magnitude bins used in the magnitude error calculation. Exists in `analysis.normalized_delta_magnitude_plotter()`.
| `SUBPLOT` | bool | If `True` four subplots are created in a 2x2 grid. If `False` plots are created individually.
| `MOF` | bool | Only used if `RUN_TYPE` is not `None`. Does `BASEPATH` entered at command line contain MOF (`MOF=True` or SOF `MOF=False` catalogs?
| `MAKE_REGION_FILES`| bool | If `True`, three DS9 region files created containing 1) objects in both `MATCH_CAT1` and `MATCH_CAT2`, 2) objects uniquely in `MATCH_CAT1` and not `MATCH_CAT2`, 3) objects uniquely in `MATCH_CAT2` and not `MATCH_CAT1`.
| `NO_DIR_MAKE` | bool| If `True` nonexistent directories will be created. If `False`, `sys.exit()` will be invoked when nonexistent directories are encountered.
| `SWAP_HAX` | bool | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis. Considered if `PLOT_COLOR` or `PLOT_MAG`.
| `SWAP_ORDER_OF_SUBTRACTION` | bool | If `False`, plots the observable from `MATCH_CAT1` minus the observable from `MATCH_CAT2`. IF `True`, plots the observable from `MATCH_CAT2` minus the observable from `MATCH_CAT1`. Only considered if `PLOT_COLOR` or `PLOT_MAG`. 
| `OVERRIDE_AXLABELS` | bool | If `True` the axes labels created in `plot_labels.get_short_difference_axlabel()` will be overwritten with a generic `$\Delta$ flux_{band}/meas_flux_err` for the horizontal axis if `PLOT_FLUX`, `$\Delta$ ({color})` for the vertical axis if `PLOT_COLOR`, and `$\Delta$ mag_{band}` for the vertical axis if `PLOT_MAG`. <br>`plot_labels.get_short_difference_axlabel()` attempts to produce an axis label that describes `$\Delta$ mag` by including the order of subtraction (of the observables from `MATCH_CAT1` and `MATCH_CAT2`), whether the observables are from Balrog-base or Balrog-injected catalogs, the injection percent (if applicable), and whether the observables are from truth or measured catalogs. 
| `RUN_TYPE` | str or `None` | FOF groups (if any) to analyse. <br> `None` `'ok'` `'rerun'` <br>`'ok'`: FOF groups *un*changed after Balrog-injection. <br>`'rerun'`: FOF groups changed after Balrog-injection. <br>`None`: FOF analysis not conducted. <br>If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
| `FOF_FIT` | str or `None` | Considered if `RUN_TYPE != None`. <br> `'mof'`, `'sof'`. 
| `VERBOSE_ING` | bool | If `True` status and progress of script will print to screen. For example, 'Plotting...', 'calculating...', 'overwriting...', etc.
| `VERBOSE_ED` | bool | If `True` results will print to screen. For example, 'got...', 'flagged...', 'checked...', etc. Note that many of these printouts are also saved in a log file.
| `NOTICE` | bool | If `True` after issuing `$python plotter.py ...` at the command line a summary of the catalog and prompt attributes is printed to the screen, followed by `raw_input()` so that user can review before running the script.

If `CORNER_2D_HIST`, the following are relevant:
`LVLS`
`CLRS_FOR_LVLS`

For color-coding:
`CMAPS`
`PT_COLORS`
`GAUSS_FIT_COLORS`
...


The following functions contain `__overwrite`:

`...get_and_reformat_base_catalog()`

`...get_coadd_catalog_for_matcher()`

`...fof_matcher()`

`...make_ngmixer_gaussap_compatible_catalog()`, 

`...matcher()`

`...stack_realizations()`

`...stack_tiles()`

`...write_to_region_files()`. 

If `__overwrite=True` a `raw_input()` prompt will appear to ensure that files are to be overwritten. In these instances, press 'Enter' to continue and 'Control+c' to stop the process.

Log files, matched catalogs, and plot names are printed with a proceeding `----->` for ease of finding and opening these files.


##### Warnings

If user does not have access to DES machines at FNAL:

Replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/stilts_matcher` (in `manipulate_catalogs.match_catalogs()`) with the correct path to `stilts_matcher`. Similarly for `fof_stilts_matcher`.

User should edit the paths to catalogs in `manipulate_catalogs.get_catalog_filename()`. 

...


##### Note about Y3 Gold catalogs

The following is relevant if `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog.

By default, `manipulate_catalogs.get_catalog_filename()` tries to access Y3 Gold catalogs saved in `/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/`. 

If the user has access to FNAL DES machines, and necessary Y3 Gold catalogs are not already in the directory above, user will need to download the Y3 Gold catalogs with the headers found in: `/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/y3_gold_general_column_query.sql` and edit `manipulate_catalogs.get_catalog_filename()` with the correct path to the Y3 Gold catalogs.

In general, user can edit path to catalogs in `manipulate_catalogs.get_catalog_filename()`.


---

# Output Directory Structure

If user wishes to change the directory structure imposed below, make changes to `outputs.get_directory()`

Directory names are determined by the constants that follow.

`OUTPUT_DIRECTORY`, which is supplied at the command line.

`BALROG_RUN`, which describes the name of the Balrog test. It is a substring of `BASE_PATH_TO_CATS`. For example, if `BASE_PATH_TO_CATS=/data/des71.a/data/kuropat/des2247-4414_sof/` then `BALROG_RUN=des2247-4414_sof`.

`MATCH_TYPE`, which describes the types of catalogs examined. Thus, it is a combination of `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2, INJ1_PERCENT, INJ2_PERCENT`. For example, if 
`MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'sof'`, `INJ1, INJ2 = True, True`, and `INJ1_PERCENT, INJ2_PERCENT = 10, 10` 
then `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `ms_matcher`.


### Condensed Output Directory Structure

Below are solely the directories and not the files that will be created within them.

```
{OUTPUT_DIRECTORY}
|    
+---- outputs
    |  
    +---- {BALROG_RUN}
        |  
        +---- {MATCH_TYPE}
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


### Full Output Directory Structure 

```
{OUTPUT_DIRECTORY}
|    
+---- outputs
    |  
    +---- {BALROG_RUN}
        |  
        +---- {MATCH_TYPE}
            |
            +---- {tile}
                |
                +---- {realization}
                    |
                    +---- catalog_compare
                    |   |   {tile}_{realization}_{MATCH_TYPE}_match1and2.csv
                    |   |   {tile}_{realization}_{MATCH_TYPE}_match1not2.csv
                    |   |   {tile}_{realization}_{MATCH_TYPE}_match2not1.csv
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
                    |   |       {tile}_{realization}_{MATCH_TYPE}_{color}_{plot_type}_{ylim}.png
                    |   |       {tile}_{realization}_{MATCH_TYPE}_{color}_corner2dhist_.png
                    |   |
                    |   +---- flux
                    |   |       {tile}_{realization}_{MATCH_TYPE}_griz_{plot_type}_.png
                    |   |       {tile}_{realization}_{MATCH_TYPE}_griz_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{MATCH_TYPE}_griz_gauss_aper_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{MATCH_TYPE}_griz_{N}sigma_clip_norm_flux_diff_histogram_.png
                    |   |       {tile}_{realization}_{MATCH_TYPE}_griz_{N}sigma_clip_gauss_aper_norm_flux_diff_histogram_.png
                    |   |
                    |   +---- magnitude
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_{plot_type}_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_cbar_cm_t_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_cbar_cm_t_err_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_completeness_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_cornerhist2d_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_hexbin_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_hist2d_{ylim}.png
                    |   |   |   {tile}_{realization}_{MATCH_TYPE}_{band(s)}_scatter_{ylim}.png
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

`{realization}` and `{tile}` will be `stack` if `STACK_REALIZATIONS=True` and `STACK_TILES=True`, respectively.

Log files are CSVs although they have `.log` extensions. Not all log files will be written to; for example, if `PLOT_MAG=False`, 'mag_completeness.log' will not be created and an empty 'dummy.log' will replace it.

Plot save names ending with '\_.png' have default axes.

`modified_match_type` is `MATCH_TYPE` unless `PLOT_COMPLETENESS=True` in which case both 10% and 20% Balrog-injected matched catalogs are plotted in the same window (if applicable). `modified_match_type` removes the percent injected from `MATCH_TYPE`. That is, `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat` results in `match_type=inj_gal_truth_cat_inj_sof_cat`.

___
##### Stacking multiple realizations or multiple tiles

Matching is performed first, then catalogs are stacked.

___
# `flagged* objects`
Docstrings frequently make reference to `flagged* objects` which refer to specific flag cuts employed in `analysis.get_good_index_using_primary_flags()` but described below for clarity.

In general, the following catalog headers are considered: 'flags' and 'cm_flags'. Exceptions are listed below. 

Star truth catalogs (`star_truth`) do not contain any flags.

Coadd catalogs (`coadd`) have a 'FLAGS' header but do not have a 'cm_flags' header (or equivalent), so only 'FLAGS' is used.

Y3 catalogs (`y3_gold`) are checked for 'FLAGS_GOLD' (replaced 'flags') and '{sof/mof}\_cm_flags' replaces 'cm_flags'. In addition, the additional flags are examined: 'SEXTRACTOR_FLAGS_{GRIZ}', 'IMAFLAGS_ISO_{GRIZ}', and, if a Y3 catalog is compared to a MOF catalog, 'MOF_CM_MOF_FLAGS'.

`ngmixer.gaussap.get_gauss_aper_flux_cat()` gives flags for the Gaussian aperture flux measurement. If `PLOT_FLUX=True` and `GAUSS_APER=True`, objects with these Gaussian aperture flags are removed.
