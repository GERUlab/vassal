""" Routines to validate data consistency and convert types

"""
import numpy as np
import pandas as pd

# -------------------------------------------------------------------------------
# Type checker

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

        # Test fails if 1 data item has attribute __len__ (ie multi dimensional)

        itemkeys = [i for i, __ in enumerate(iter(data))]

        itemhaslen = [hasattr(data.__getitem__(i), '__len__') for i in itemkeys]

        if any(itemhaslen):
            test = False

    return test


def is_valid_group_dict(grpdict):
    """Test if grpdict has str keys and int or list of int values
    
    Examples
    --------
    
    >>> d = {'foo': 1, 'bar': [1, 2]}
    >>> is_valid_group_dict(d)
    True
    
    >>> d = {0:1, 1:2}
    >>> is_valid_group_dict(d)
    False
    
    >>> d = {'foo':1, 'bar':[1,'2']}
    >>> is_valid_group_dict(d)
    False
    
    """

    keys = grpdict.keys()
    values = grpdict.values()

    test_keys = all(isinstance(k, str) for k in keys)

    test_values = all(is_int_or_list_of_int(v) for v in values)

    return test_keys and test_values


def is_int_or_list_of_int(item):
    """Test if an item is int or list of int
    
    Examples
    --------
    
    >>> is_int_or_list_of_int(10)
    True
    
    >>> is_int_or_list_of_int('10')
    False
    
    >>> is_int_or_list_of_int([1,2])
    True
    
    >>> is_int_or_list_of_int([1, '2'])
    False
    
    """

    if isinstance(item, int):

        test = True

    elif isinstance(item, list):

        test = all(isinstance(i, int) for i in item)

    else:

        test = False

    return test


# -------------------------------------------------------------------------------
# Type converters


def arraylike_to_nparray(arraylike):
    """Return np array from 1d array like
    
    Examples
    --------
    
    >>> d = {0:3, 1:4}
    >>> arraylike_to_nparray(d)
    array([3, 4])
     
     >>> s = pd.Series([1,2])
     >>> arraylike_to_nparray(s)
     array([1, 2], dtype=int64)
     
     >>> lst = [1,2]
     >>> arraylike_to_nparray(lst)
     array([1, 2])
    
    """

    # I use assert because this function is always called after is_1darray_like

    assert (is_1darray_like(arraylike))

    # The only problem comes with dictionnaries. Dict values sould be passed to
    # np.array.

    if isinstance(arraylike, dict):

        nparr = np.array(arraylike.values())

    else:

        nparr = np.array(arraylike, dtype=None)

    return nparr

def nested2d_to_flatlist(nestedlist):
    """Returns a flat list from a nested 2d list regardless item is iterable
    
    Examples
    --------
    
    >>> a = [1, [2,3], 4]
    >>> nested2d_to_flatlist(a)
    [1, 2, 3, 4]
    
    """

    flatlist = []

    for item in nestedlist:
        if hasattr(item, '__len__'):
            [flatlist.append(subitem) for subitem in item]
        else:
            flatlist.append(item)

    return flatlist



if __name__ == '__main__':
    import doctest
    doctest.testmod()
