# superstamp_subtraction
Utilities for performing image subtraction photometry on the Kepler superstamps, with a minimal list of cluster members provided for testing.

To run this code, you will need a local copy of the Kepler superstamp FITS files for NGC 6791 and NGC 6819. It is suggested that you store these in a subdirectory of the code directory called `superstamps`, however there is an option when running the code to specfiy any file path.

Before running the code using the ensemble centroid (centroid option `e`), you will need to run `globalcentroid.py`, which takes as its options the location of the superstamp directory (`-s`/`--sfilepath`) and the cluster's NGC number (`-n`/`--ngc`). Data files containing a list of the bright stars used to calculate the ensemble centroid have been included. These measurements will be saved in a subdirectory called `centroids`.

To perform image subtraction on a given target, run `subtract.py`. This code takes the following options:

* `-k`, `--kic`: target KIC ID, integer, required
* `-s`, `--stampfilepath`: file location of superstamps, string, default `./superstamps/`
* `-p`, `--centroidfilepath`:  file location of ensemble centroid data, string, default `./centroids/`
* `-r`, `--resamplefactor`: resampling factor, integer, default 2
* `-c`, `--centroid`: choosing between ensemble or target centroid, string, required with options `e` or `t`
* `-d`, `--dims`: dimensions of cutout around target from central pixel to edge, integer, default `3`, with options from `2` to `19`
* `-f`, `--figures`: toggle to make summary plots, bool, default True
* `-i`, `--interactive`: toggle to show interactive plots when code has run, bool, default False

For example, a typical run of the code might look like:

`>> python subtract.py -k 2437171 -c e -r 21 -p ../otherfolder/ensemblecentroids/`

Note that the code determines which cluster the star is in based on its KIC ID, so that does not need to be specified here.
