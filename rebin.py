import numpy as np

# IDL rebin in python
# Thanks to Andrea Zonca
# https://gist.github.com/zonca/1348792

def rebin(a_, n, m):
    """
    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions

    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array

    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated

    See Also
    --------
    resize : Return a new array with the specified shape.

    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])

    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])

#    """
    a = a_.copy()
    M, N = a.shape
#    m, n = new_shape
    if M % m or N % n:
      print('REBIN: Result dimensions must be integer factor of original dimensions.')
      return
    if m<M:
        return a.reshape((m, int(M/m),n, int(N/n))).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, int(m/M), axis=0), int(n/N), axis=1)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
