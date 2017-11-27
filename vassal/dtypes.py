""" Routines to validate data consistency

See Also
--------



"""
import numpy as np
import pandas as pd
from pandas.core.dtypes.common import is_list_like

def is_1darray_like(data):
    """
    Parameters
    ----------
    data : object
        the object to be tested

    Returns
    -------
    bool
        True if data is 1d list-like object
    
    References
    ----------
    .. _Source:
        https://github.com/pandas-dev/pandas/blob/master/pandas/core/dtypes/common.py
        
    Examples
    --------
    
    Scalar
    >>> is_1darray_like(1)
    False
    
    List
    >>> is_1darray_like([1])
    True
    
    np.array 1d
    >>> is_1darray_like(np.arange(10))
    True
    
    np.array nd
    >>> is_1darray_like(np.zeros(shape=(4,3)))
    False
    
    pd.Series
    >>> s = pd.Series(np.arange(10))
    >>> is_1darray_like(s)
    True
    
    pd.DataFrame 
    >>> df = pd.DataFrame(s)
    >>> is_1darray_like(df)
    False
    
    dictionary 1d
    >>> d1 = dict(zip(s,s))
    >>> is_1darray_like(d1)
    True
    
    dictionary nd
    >>> d2 = dict(zip(s, [[1,2]]*len(s)))
    >>> is_1darray_like(d2)
    False
    
    generator
    >>> gen = (i for i in range(10))
    >>> is_1darray_like(gen)
    False
    
        
    """
    test = True

    # Test fails if data has no attribute __getitem__

    if not hasattr(data, '__getitem__'):
        test = False

    else:
        item = data[0]

        # Test fails if data item has attribute __len__ (ie multi dimensional)

        if hasattr(item, '__len__'):
            test = False

    return test


if __name__ == '__main__':
    import doctest
    doctest.testmod()
