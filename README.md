# BalVal

Conducts some [Balrog-GalSim](https://github.com/sweverett/Balrog-GalSim) validation tests.

___
**Contents of the repository**

`ms_matcher` matches two catalogs using STILTS.

`ms_plotter.py` calls `ms_matcher`, analyzes the matched catalog, and produces various magnitude versus Delta-magnitude plots.

___

# Running BalVal

General: `$python ms_plotter.py base_path_to_catalogs output_directory realization tile`

Ex: `$python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/ /plots/ all DES2247-4414`

After the above command is issued, a prompt will appear so that the user can confirm plot attributes. This is to prevent plots from being overwritten when testing new additions to the script. User can comment out `NOTICE` to remove prompt.

Examples of plot attributes, which are explained in `ms_plotter.py` are:

`NORMALIZE` `HEXBIN` `CM_T_S2N_COLORBAR` `CM_T_COLORBAR`  `PLOT_1SIG`

**Warnings**

User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/ms_matcher` (in `ms_plotter.matcher()`) with the correct path to `ms_matcher`.

___

**Stacking multiple realizations**

Matching is performed first, then catalogs are stacked.

___

**Directory structure for saving plots**

Directory structure: `outdir/{tile}/{plot_type}/{realization}/`.

Allowed values for `{plot_type}` are: `normalized` `scatter`

Allowed values for `{realization}` are: `0` `1` ... `stacked`

User can change directory structure in `ms_plotter.get_plot_save_name()`
