# superstamp_subtraction
Utilities for performing image subtraction photometry on the Kepler superstamps

Full description of parameters to come. Minimal list of input members provided for testing

To run this code, you will need a local copy of the Kepler superstamp FITS files for NGC 6791 and NGC 6819. It is suggested that you store these in a subdirectory of the code directory called `superstamps`, however there is an option when running the code to specfiy any file path.

Before running the code using the ensemble centroid (centroid option `e`), you will need to run globalcentroid.py, which takes as its options the location of the superstamp directory and the cluster's NGC number. Data files containing a list of the bright stars used to calculate the ensemble centroid have been included. These measurements will be saved in a subdirectory called `centroids`.
