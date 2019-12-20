import numpy as np
import matplotlib.pyplot as plt
from functools import partial

def regrid(flux, factor=2):
    """
    Takes a 2D array and expands it

    The flux array is expanded into a series of factor by factor submatricies. These
    submatricies are arranged into the same ordering as the original array.

    Arguments:
    :: flux (2D ndarray) :: The 2D array to be regridded

    Keyword arguments:
    :: factor (int) :: The scaling factor for the regridding. 
        Default is 2

    Returns:
    :: 2D ndarray :: Containing the regridded values

    Function courtesy of Alan Robertson: http://github.com/alan-robertson
    """

    # Increasing function; don't try to multiply lists as it copies the reference 
    inc = lambda factor, value: np.zeros((factor, factor)) + (value/(factor**2))
    inc = partial(inc, factor)
    
    # Some useful formats
    block_shape = tuple(list(flux.shape) + [factor, factor])
    final_shape = tuple(np.array(flux.shape) * factor) 
    transpose_ordering = [0,3,1,2]

    # The Map operation 
    re_flux = np.array(list(map(inc, flux.reshape(flux.size))))

    # Reshaping to the correct view
    re_flux = re_flux.reshape(block_shape).transpose(transpose_ordering)
    re_flux = re_flux.reshape(final_shape)
    
    return re_flux

def regrid_slow(arr, factor):
   new_arr = np.zeros((arr.shape[0]*factor, arr.shape[1]*factor))

   for index, val in np.ndenumerate(new_arr):
      new_arr[index] = arr[int(index[0]/factor),int(index[1]/factor)]

   return new_arr
