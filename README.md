# BalVal: [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) Validation

Conducts various [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

Validation tests include comparing the color, magnitude, and flux from input (truth) catalogs to output (Balrog-injected) catalogs. See below for more capabilities.

___
### Table of Contents

[Running BalVal](https://github.com/spletts/BalVal/blob/master/README.md#running-balval)
 
  - [Dependencies](https://github.com/spletts/BalVal/blob/master/README.md#dependencies)
  
  - [Contents of Repository](https://github.com/spletts/BalVal/blob/master/README.md#contents-of-repository)

[Table of Constants](https://github.com/spletts/BalVal/blob/master/README.md#table-of-constants)

  - [Plot Labels](https://github.com/spletts/BalVal/blob/master/README.md#plot-labels)

  - [Plot Colors](https://github.com/spletts/BalVal/blob/master/README.md#plot-colors)

[Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#directory-structure)

  - [Expected Catalog Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#expected-catalog-directory-structure)

    - [Y3 Gold Catalogs and Deep SN Coadds](https://github.com/spletts/BalVal/blob/master/README.md#y3-gold-catalogs-and-deep-sn-coadds)
    
  - [Condensed Output Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#condensed-output-directory-structure)

  - [Full Output Directory Structure](https://github.com/spletts/BalVal/blob/master/README.md#full-output-directory-structure)
  
    - [Output Notes](https://github.com/spletts/BalVal/blob/master/README.md#output-notes)
  
  
[Flagged* Objects](https://github.com/spletts/BalVal/blob/master/README.md#flagged-objects)



___
# Running BalVal

##### Setup
Prior to running, user with access to the DES machines at FNAL must:

1. Use `matplotlib v2.2.2` via `$ssh {user}@des70.fnal.gov` or `$ssh {user}@des71.fnal.gov`, and
2. `$source /home/s1/mspletts/setup_ngmixer_gaussap.sh`.

##### Use
General: `$python plotter.py {BASE_PATH_TO_CATS} {OUTPUT_DIRECTORY} {realization(s)} {tile(s)}`

Ex: `$python plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/ /BalVal/ 0,1 DES2247-4414`

`None` is an allowed value for `{realization(s)}`. User can list realizations at the command line with commas and _no_ spaces separating the realizations (similarly for tiles). `{tile(s)}` can also be specified with a `.dat` file.

User sets plot attributes and catalog attributes within `set_constants.py`. User-set attributes are listed in a [Table of Constants](https://github.com/spletts/BalVal/blob/master/README.md#table-of-constants).


### Dependencies

Non-standard dependencies:

* Requires [`ngmixer`](https://github.com/esheldon/ngmixer) to measure fluxes using a Gaussian aperture. 

  * Note that `/home/s1/mspletts/setup_ngmixer_gaussap.sh` points to `ngmixer` as cloned in `/home/s1/mspletts/`, where minor changes have been made to `ngmixer.gaussap.get_gauss_aper_flux_cat()`. These changes are indicated by
   ```
   ### MS ### 
   {modified code} 
   ### - ###
   ```

* Requires [`corner.py`](https://github.com/dfm/corner.py) to make corner plots.


### Contents of Repository

See docstring and the beginning of each script for more detail.

- `set_constants.py`: Set catalog and plot attributes. User interacts with this.

- `plotter.py`: Produces plots. User runs this script.

- `catalog_classes.py`: Contains the relevant and permitted headers for each catalog type (see `MATCH_CAT1`, `MATCH_CAT2` in [Table of Constants](https://github.com/spletts/BalVal/blob/master/README.md#table-of-constants) for allowed catalog types).

- `catalog_headers.py`: Contains headers for `MATCH_CAT1` and `MATCH_CAT2` (`MATCH_CAT1` and `MATCH_CAT2` are established in `set_constants.py`). 

- `stilts_matcher`: Matches two catalogs on RA and Dec using [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/). Returns a CSV.

- `fof_stilts_matcher`: Creates a number of catalogs needed to analyse FOF groups. Catalogs are matched using [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/). Returns CSVs.

- `par.py`: Finds FOF groups changed and unchanged after Balrog-injection. This is called by `fof_stilts_matcher`. Written by Brian Yanny.

- `outputs.py`: Creates and writes to various output files.

- `analysis.py`: Conducts analysis needed to make plots. E.g. remove flags, error analysis, etc.

- `region_files.py`: Writes to region files.

Note that after matching two catalogs with STILTS (in `stilts_matcher` and `fof_stilts_matcher`), the output CSV contains *strings*, not arrays. For example, if in an input FITS file a cell contains an array of form `(1 2 3 4)`, the array is converted to a string of form `'(1, 2, 3, 4)'` in the matched CSVs. Similarly for matrices, etc.



___
# Table of Constants

These constants are set within `set_constants.py`. Alphabetized below.

Docstrings describe parameters but a [Table of Function Parameters and Returns](https://docs.google.com/spreadsheets/d/1utJdA9SpigrbDTsmtHcqcECP9sARHeFNUXdp47nyyro/edit?usp=sharing) is also available.


Parameter(s) | Type | Description <br> `allowed values` (if Type not bool)
:---: | :---: | ---
| `CENTER_MAG_ERR_ABT_ZERO`| bool | If `True` and `PLOT_MAG=True` and `NORMALIZE=True` the magnitude error is centered about zero in the plot. If `False` the magnitude error is centered about the median of the magnitude difference in each bin (which is close to zero). The distinction (minorly) affects the number of objects within one-sigma of the magnitude error. 
|`CM_T_CBAR` | bool | Used if `PLOT_MAG=True` and `NORMALIZE_MAG=False`. If `True` a colorbar of the measured cm_T (size squared) is added to the plot.
| `CM_T_ERR_CBAR` | bool | Used if `PLOT_MAG=True` and `NORMALIZE_MAG=False`. If `True` a colorbar of the measured cm_T error (error of size squared) is added to the plot.
| `COLOR_YLOW` `COLOR_YHIGH` | float OR `None` | Limits for the vertical axis of the plot created if `PLOT_COLOR=True`. If `None` symmetric limits are forced based on the default plot scaling of `matplotlib`. Must both be float or both be `None`.
| `CORNER_HIST_2D` | bool | Used if `PLOT_MAG=True` or `PLOT_COLOR=True`. If `True` `corner.hist2d` plots are created using [`corner.py`](https://github.com/dfm/corner.py).
| `CORNER_HIST_2D_LVLS` | 1D array | If `CORNER_2D_HIST=True` this is used as `levels` in [`corner.hist2d()`](https://github.com/dfm/corner.py/blob/master/corner/corner.py).
| `CORNER_HIST_2D_LVLS_CLRS` | list of str | If `CORNER_HIST_2D=True` this is used as the color(s) for the `levels` above.
| `CORNER_HIST_2D_LVLS_CLRS_LABEL` | list of str | If `CORNER_HIST_2D=True` this is used to label `levels`.
| `FLUX_XLOW` `FLUX_XHIGH` | float OR `None` | Limits for the horizontal axis of the plot created if `PLOT_FLUX=True`. If these values are `None` default `matplotlib` scaling is used. Must both be float or both be `None`.
| `FOF_FIT` | str | Used if `RUN_TYPE` is not `None` to specify the catalogs to be used in FOF analysis. User should set this if `RUN_TYPE` is not `None`. <br> `'mof'` `'sof'`.
| `HEXBIN` | bool | Used if `PLOT_MAG=True`. If `True` a density plot is produced via `matplotlib.pyplot.hexbin()`.
| `HIST_2D` | bool | Used if `PLOT_MAG=True`. If `True` a 2D histogram is plotted via `matplotlib.pyplot.hist2d`.
| `INJ1` `INJ2` | bool | If `INJ1=True` then `MATCH_CAT1` is Balrog-injected. If `False` then `MATCH_CAT1` is a base (non-injected) catalog. Similarly for `INJ2` and `MATCH_CAT2`.
|`INJ1_PERCENT` `INJ2_PERCENT` | int | These are set by `calculate_injection_percent.get_injection_percent()`. User can hardcode these, and should do so if `{BASE_PATH_TO_CATS}` contains catalogs with multiple injection densities (for example, the Balrog run for TAMU had both 10% and 20% injections).
| `MAG_YLOW` `MAG_YHIGH` | float OR `None` | Limits for the vertical axis of the plot created if `PLOT_MAG=True`. `None` results in default `matplotlib` scaling. Must both be float or both be `None`.
| `MAKE_REGION_FILES`| bool | If `True`, three DS9 region files are created, containing 1) objects in both `MATCH_CAT1` and `MATCH_CAT2`, 2) objects uniquely in `MATCH_CAT1` and not `MATCH_CAT2`, 3) objects uniquely in `MATCH_CAT2` and not `MATCH_CAT1`.
|`MATCH_CAT1` `MATCH_CAT2` | str | Type of catalogs to analyse. <br>`'coadd'` `'gal_truth'` `'mof'` `'sof'` `'star_truth'` `'y3_gold_2_0'` `'y3_gold_2_2'` `'deep_sn_mof'` `'deep_sn_sof'`
| `N` | float | Used if `PLOT_FLUX=True` and `SIGMA_CLIP_NORM_FLUX_DIFF=True` to perform `N`-sigma clip.
| `NO_DIR_MAKE` | bool | If `True` nonexistent directories will be created. If `False`, `sys.exit()` will be invoked when nonexistent directories are encountered to remind user to edit directory structure in `outputs.get_directory()`.
| `NORMALIZE_MAG` | bool | Used if `PLOT_MAG=True`. If `True` the magnitude plot is normalized according to the measured one-sigma magnitude error. Red dashed lines representing +/- one-sigma are also plotted.
| `NORMALIZE_NORM_FLUX_DIFF_VIA_DENSITY` | bool | Used if `PLOT_FLUX=True`. If `True` histograms are normalized (meaning, in this instance, that the area under the histogram is 1).
| `NOTICE` | bool | If `True` after issuing `$python plotter.py ...` at the command line, a summary of the catalog and plot attributes is printed to the screen, followed by `raw_input()` so that user can review before running the script.
| `NUM_ITERS` | float OR `None` | Considered if `SIGMA_CLIP_NORM_FLUX_DIFF=True` and `PLOT_FLUX=True`. Number of iterations to use in `N`-sigma clip. See [`iters`](http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html) documentation.
| `__overwrite` | bool | If `True` a `raw_input()` prompt will appear to ensure that files are to be overwritten. In these instances, press 'Enter' to continue and 'Control+c' to stop the process. This constant exisits within the following functions: `...get_and_reformat_base_catalog()`, `...get_coadd_catalog_for_matcher()`, `...fof_matcher()`, `...make_ngmixer_gaussap_compatible_catalog()`, `manipulate_catalogs.match_catalogs()`, `...stack_realizations()`, `...stack_tiles()`, `...write_to_region_files()`. 
| `OVERWRITE_AXLABELS` | bool | If `True` the axes labels created in `plot_labels.get_short_difference_axlabel()` will be overwritten with a generic `$\Delta$ flux_{band}/meas_flux_err` for the horizontal axis if `PLOT_FLUX=True`, `$\Delta$ ({color})` for the vertical axis if `PLOT_COLOR=True`, and `$\Delta$ mag_{band}` for the vertical axis if `PLOT_MAG=True`.
| `PLOT_34P_SPLIT` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE_MAG=True`. If `True` the 34<sup>th</sup> percentiles of the binned positive magnitude differences are plotted. One percentile is computed per bin. Similarly for the negative magnitude differences. Bins refer to the magnitude bins used in the magnitude error calculation.
| `PLOT_68P` | bool | Considered if `PLOT_MAG=True` and `NORMALIZE_MAG=True`. If `True` the 68<sup>th</sup> percentiles of the binned magnitude differences are plotted. One percentile is computed per bin. Bins refer to the magnitude bins used in the magnitude error calculation.
| `PLOT_CM_FLUX` | bool | If `True` and `PLOT_FLUX=True` the CM flux is plotted (as a 1D histogram). Note that both this and `PLOT_GAUSS_APER_FLUX` can be `True`, and both will be plotted on the same plot.
| `PLOT_COLOR` | bool | If `True` colors (g-r), (r-i), and (i-z) (horizontal axis) versus color difference (vertical axis) are plotted. Creates a 2x2 grid for each color, with subplots corresponding to different magnitude bins (currently \[20,21), \[21,22), \[22,23), and \[23,24), but user can change this via `__mag_bins_for_color` in `analysis.get_color_from_binned_magnitude()`). Magnitudes are binned according to values in `MATCH_CAT1` for the leading filter (g for (g-r), etc). By default (`SWAP_ORDER_OF_SUBTRACTION=False`) the color difference is calculated via the color from `MATCH_CAT1` minus the color from `MATCH_CAT2`. This order is reversed if `SWAP_ORDER_OF_SUBTRACTION=True`.
| `PLOT_COMPLETENESS` | bool | Used if `PLOT_MAG=True`. If `True` one completeness plot per band is produced. The completeness of objects with a magnitude difference less than 2 (set by `__cutoff_mag_diff` in `analysis.get_magnitude_completeness()`) is considered. This must be used with either `MATCH_CAT1` or `MATCH_CAT2` set to `gal_truth` or `star_truth` (because the completeness is determined by the recovery of magnitudes in a truth catalog). Magnitude bins are set via `COMPLETENESS_PLOT_MAG_BINS` in `set_constants.py`.
| `PLOT_FLUX` | bool | If `True` a 1D histogram of the g-, r-, i-, and z-band flux difference normalized to the measured flux error are plotted in a 2x2 plot grid. This must be used with either `PLOT_CM_FLUX=True` or `PLOT_GAUSS_APER_FLUX=True` or both. The flux difference is computed as the flux from the truth catalog minus the flux from the measured catalog (if `MATCH_CAT1` or `MATCH_CAT2` is a truth catalog). Otherwise, if `SWAP_ORDER_OF_SUBTRACTION=False` the flux difference is calculated via the flux from `MATCH_CAT1` minus the flux from `MATCH_CAT2`. A standard Gaussian distribution (mean=0, standard_deviation=1) is also plotted.
| `PLOT_GAUSS_APER_FLUX` | bool | If `True` and `PLOT_FLUX=True` a Gaussian aperture method is used to measure the fluxes that are to be plotted. Note that both this and `PLOT_CM_FLUX` can be `True`, and both will be plotted on the same plot.
| `PLOT_GAUSSIAN_FIT` | bool | If `True` and `PLOT_FLUX=True` a Gaussian fit to the distribution(s) is plotted. If `PLOT_CM_FLUX=True` and `PLOT_GAUSS_APER_FLUX=True` the Gaussian fits to both will be plotted.
| `PLOT_MAG` | bool | If `True` plots of g-, r-, i-, and z-band magnitudes (horizontal axis) versus magnitude differences (vertical axis) are created. If `SWAP_ORDER_OF_SUBTRACTION=False` the magnitude difference is computed via the magnitude from `MATCH_CAT1` minus the magnitude from `MATCH_CAT2`.
| `PLOT_MAG_ERR` | bool | If `True` and `PLOT_MAG=True` the one-sigma measured magnitude error curve is plotted, centered according to `CENTER_MAG_ERR_ABT_ZERO`. `NORMALIZE_MAG` must be False.
| `PLOT_PEAKS` | bool | If `True` and `PLOT_FLUX=True` a dashed vertical line representing the median of the normalized (to measured flux error) flux difference is plotted. If `PLOT_CM_FLUX=True` and `PLOT_GAUSS_APER_FLUX=True` the medians of both are plotted.
| `RAW_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the normalized (to measured flux error) flux differences (determined by `PLOT_CM_FLUX` and `PLOT_GAUSS_APER_FLUX`) are plotted without percentile trimming and without sigma clipping.
| `RUN_TYPE` | str or `None` | FOF groups (if any) to analyse. <br> `None` `'ok'` `'rerun'` <br>`None`: FOF analysis not conducted. <br>`'ok'`: FOF groups *un*changed after Balrog-injection. <br>`'rerun'`: FOF groups changed after Balrog-injection.<br>If `RUN_TYPE='rerun'` or `RUN_TYPE='ok'` then `MATCH_CAT1` `MATCH_CAT2` `INJ1` and `INJ2` will be overwritten.
| `SAVE_BAD_IDX` | bool | If `True` the indices of objects with [flags*](https://github.com/spletts/BalVal/blob/master/README.md#flagged-objects) are stored.
| `SAVE_PLOT` | bool | If `True` the plot is saved using then filename assigned by `outputs.get_plot_filename()`.
| `SHOW_PLOT` | bool | If `True` the plot is displayed.
| `SCATTER` | bool | Used if `PLOT_MAG=True` or `PLOT_COLOR=True`. If `True` a `matplotlib.pyplot.scatter()` plot is produced.
| `SIGMA_CLIP_NORM_FLUX_DIFF` | bool | Used if `PLOT_FLUX=True`. If `True` the distribution of normalized (to measured flux error) flux differences is sigma clipped using `N`-sigma. This distribution can be one or both of `PLOT_CM_FLUX` and `PLOT_GAUSS_APER_FLUX`.
| `STACK_REALIZATIONS` | bool | If `True` catalogs are matched then stacked. One stacked realization catalog is produced per tile. Plotting resumes with stacked catalog. Must be used with, for example, `0,1,2` at the command line.
| `STACK_TILES` | bool | If `True` catalogs are matched then stacked. One stacked tile catalog is produced per realization. Must be used with, for example, `DES0220-0207,DES0222+0043` at the command line.
| `SUBPLOT` | bool | If `True` and `PLOT_MAG=True` four subplots are created in a 2x2 grid. If `False` plots are created individually.
| `SWAP_HAX` | bool | If `False` (default) `MATCH_CAT1` values are plotted on the horizontal axis. If `True` `MATCH_CAT2` values are plotted on the horizontal axis. Considered if `PLOT_COLOR=True` or `PLOT_MAG=True`.
| `SWAP_ORDER_OF_SUBTRACTION` | bool | If `False` differences are computed via the observable from `MATCH_CAT1` minus the observable from `MATCH_CAT2`. If `True` differences are computed via the observable from `MATCH_CAT2` minus the observable from `MATCH_CAT1`. Considered if `PLOT_COLOR` or `PLOT_MAG`.
| `TRIM_NORM_FLUX_DIFF` | bool | If `True` and `PLOT_FLUX=True` histograms of the distribution(s) of normalized (to measured flux error) flux difference are trimmed to include the 2<sup>nd</sup>-98<sup>th</sup> percentiles. This distribution can be one or both of `PLOT_CM_FLUX` and `PLOT_GAUSS_APER_FLUX`.
| `VERBOSE_ED` | bool | If `True` results will print to screen. For example, 'Computed...', 'Flagged...', 'Checked...', etc. Note that many of these printouts are also saved in a log file.
| `VERBOSE_ING` | bool | If `True` status and progress of script will print to screen. For example, 'Plotting...', 'calculating...', 'overwriting...', etc.
| `Y3_MODEL` | str | Considered if `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog.<br> `'CM'` `'PSF'`
| `Y3_FIT` | str | Considered if `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog. `set_constants.py` attempts to set this using `MATCH_CAT1` or `MATCH_CAT2`, but user can hardcode this. <br> `'MOF'` `'SOF'`

##### Plot Labels

`plot_labels.get_short_difference_axlabel()` attempts to produce an axis label that describes: the order of subtraction (of the observables from `MATCH_CAT1` and `MATCH_CAT2`), which of the observables are from Balrog-base and Balrog-injected catalogs, the injection percent (if applicable), and whether the observables are from truth or measured catalogs. But the return of `plot_labels.get_short_difference_axlabel()` doesn't always look nice (particularly if its inputs have a short common substring). User can set `OVERWRITE_AXLABEL=True` in `set_constants.py` if axes labels are difficult to read.

##### Plot Colors

Colormaps and colors for plots are set by the following in `set_constants.py`: `CMAPS`, `PT_COLORS`, `GAP_FLUX_HIST`, `FIT_TO_GAP_FLUX`, `FLUX_HIST`, and `FIT_TO_FLUX`, where `GAP` represents the Gaussian aperture measured flux. These variables are explicitly stated here so user can more easily alter the plot colors.


___
# Directory Structure

### Expected Catalog Directory Structure

In general, user can edit path to and filename of catalogs in `manipulate_catalogs.get_catalog_filename()`.

##### Y3 Gold Catalogs and Deep SN Coadds

The following is relevant if either `MATCH_CAT1` or `MATCH_CAT2` is a Y3 Gold catalog or a deep SN coadd chip.

By default, `manipulate_catalogs.get_catalog_filename()` tries to access Y3 Gold catalogs and deep SN coadds in `/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/` and `/data/des71.a/data/mspletts/balrog_validation_tests/deep_sn_catalogs/`, respectively.

If user has access to FNAL DES machines, and necessary Y3 Gold catalogs are not already in the directory above, user will need to download the Y3 Gold catalogs with the headers found in: `/data/des71.a/data/mspletts/balrog_validation_tests/y3_gold_catalogs/y3_gold_general_column_query.sql` and edit `manipulate_catalogs.get_catalog_filename()` with the correct path to the Y3 Gold catalogs.


### Condensed Output Directory Structure

If user wishes to change the output directory structure imposed below, make changes to `outputs.get_directory()`.

Directory names are determined by the constants that follow.

`OUTPUT_DIRECTORY`, which is supplied at the command line.

`BALROG_RUN`, which describes the name of the Balrog test. This is a substring of `BASE_PATH_TO_CATS`. For example, if `BASE_PATH_TO_CATS=/data/des71.a/data/kuropat/des2247-4414_sof/` then `BALROG_RUN=des2247-4414_sof`.

`MATCH_TYPE`, which describes the types of catalogs examined. Thus, it is a combination of `MATCH_CAT1, MATCH_CAT2, INJ1, INJ2, INJ1_PERCENT, INJ2_PERCENT`. For example, if 
`MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'sof'`, `INJ1, INJ2 = True, True`, and `INJ1_PERCENT, INJ2_PERCENT = 10, 10` 
then `MATCH_TYPE=10%_inj_gal_truth_cat_10%_inj_sof_cat`. Note that `MATCH_TYPE` reflects the order in which the catalogs were matched in `stilts_matcher` (that is, the catalog that served as `in1` is listed first in `MATCH_TYPE`).


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
                    |   |   {tile}_{realization}_{MATCH_TYPE}_num_objs_in_1sig_mag.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_flags.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_mag_completeness.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_mag_error_computation.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_matched_catalogs.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_outlier_mag_diffs.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_sigma_clip_flux.log
                    |   |   {tile}_{realization}_{MATCH_TYPE}_flux_from_gaussian_aperture.log
                    |   |
                    |   +---- fof_analysis
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_color_from_mag.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_num_objs_in_1sig_mag.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_flags.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_mag_completeness.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_mag_error_computation.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_matched_catalogs.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_outlier_mag_diffs.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_sigma_clip_flux.log
                    |           {tile}_{realization}_{MATCH_TYPE}_{RUN_TYPE}_flux_from_gaussian_aperture.log
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

##### Output Notes

  - `{realization}` and `{tile}` will be `'stack'` if `STACK_REALIZATIONS=True` and `STACK_TILES=True`, respectively.
  
  - When stacking multiple realizations or multiple tiles, matching is performed first, then catalogs are stacked.

  - Log files are CSVs although they have `.log` extensions. Not all log files will be written to; for example, if `PLOT_MAG=False`, 'mag_completeness.log' will not be created and an empty 'dummy.log' will replace it.
  
  - Log filenames, matched catalog filenames, and plot filenames are printed with a preceding `----->` for ease of finding and opening these files.




___
# Flagged* Objects
Docstrings frequently make reference to 'flagged* objects' which refer to specific flag cuts employed in `analysis.get_good_index_using_primary_flags()` and described below for clarity.

In general, the following catalog headers are considered: 'flags' and 'cm_flags'. Exceptions are listed below. 

Star truth catalogs (`star_truth`) do not contain any flags.

Coadd catalogs (`coadd`) have a 'FLAGS' header but do not have a 'cm_flags' header (or equivalent), so only 'FLAGS' is used.

Y3 Gold catalogs (`y3_gold_2_0` `y3_gold_2_2`) are checked for 'FLAGS_GOLD' (which replaces 'flags' mentioned above), and 'SOF_CM_FLAGS' or 'MOF_CM_FLAGS' replaces 'cm_flags' (depending on `Y3_FIT`). In addition, the following flags are examined: 'SEXTRACTOR_FLAGS_{GRIZ}', 'IMAFLAGS_ISO_{GRIZ}', and, if a Y3 Gold catalog is compared to a MOF catalog, 'MOF_CM_MOF_FLAGS'.

`ngmixer.gaussap.get_gauss_aper_flux_cat()` gives flags for the Gaussian aperture flux measurement. If `PLOT_FLUX=True` and `PLOT_GAUSS_APER_FLUX=True`, objects with these Gaussian aperture flags are removed.


___
If user does not have access to DES machines at FNAL:

Replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/stilts_matcher` (in `manipulate_catalogs.match_catalogs()`) with the correct path to `stilts_matcher`. Similarly for `fof_stilts_matcher`.

There's more for the non-DES machine user!...
