"""Singular Spectrum Analysis with numpy

"""
import warnings

import numpy as np
import matplotlib.pyplot as plt

from vassal.dtypes import (
    is_1darray_like,
    is_valid_group_dict)

from vassal.base import BaseSSA
from vassal.plot import PlotSSA

try:
    import pandas as pd

    __TS_DEFAULT_TYPE__ = 'pdseries'
except:
    warnings.warn('pandas module missing: __TS__DEFAULT_TYPE__ set to np.array')
    __TS_DEFAULT_TYPE__ = 'nparray'


# TODO: test behavior with np.nan in series

class BasicSSA(BaseSSA, PlotSSA):
    """A class for basic Singular Spectrum Analysis 
    
    Singular Spectrum Analysis (SSA) is a non-parametric method
    to decompose and recompose a signal into specific components:
    trend, seasonality, noise, ... 
    
    The method consists in embedding the time series into a
    trajectory matrix. The matrix is then decomposed into 
    eigentriples using singular value decomposition. Eigentriples 
    can be grouped by user. Each group can be recomposed into a 
    time series component.
    
    Grouping can be user defined based on the interpretation of
    the singular value plot or the w-corr plot.
    
    Parameters
    ----------
    #TODOC
    
    Examples
    --------
    
    Loading time series
    
    >>> co2 = pd.read_csv("test/co2.csv", index_col=0, header=None)
    >>> co2 = co2[co2.columns[0]]
    >>> co2.describe()
    count    468.000000
    mean     337.053526
    std       14.966220
    min      313.180000
    25%      323.530000
    50%      335.170000
    75%      350.255000
    max      366.840000
    Name: 1, dtype: float64
    
    Decomposition
    
    >>> co2_ssa = BasicSSA(co2)
    
    Reconstruction
    
    >>> groups = { 'Trend': [0, 3], 'Season': [1,2,4,5] }
    >>> co2_ssa.reconstruct(groups)
    >>> print co2_ssa.groups
    ['Original', 'Trend', 'Season', 'Residuals']

    >>> co2_ssa['Season'].describe()
    count    468.000000
    mean      -0.002320
    std        2.055533
    min       -3.461722
    25%       -1.883337
    50%        0.291613
    75%        1.757231
    max        3.301371
    dtype: float64
    
    Weighted correlation of components
    
    >>> co2_ssa.wcorr(components=3)
    array([[ 1.        , -0.0124797 ,  0.00814555],
           [-0.0124797 ,  1.        ,  0.97536049],
           [ 0.00814555,  0.97536049,  1.        ]])
           

        
    References
    ----------
    
    [1] Singular Spectrum Analysis for Time Series | Nina Golyandina | Springer. 
    Accessed November 19, 2017. //www.springer.com/gp/book/9783642349126.
        
    """

    def __init__(self, ts=None, window=None, svdmethod='nplapack',
                 usetype=__TS_DEFAULT_TYPE__):

        # we pass the svd method to the base class init in order to wrap the
        # proper decompose method

        super(BasicSSA, self).__init__(ts=ts, svdmethod=svdmethod,
                                       usetype= usetype)

        # define window length if none

        if window is None:
            window = self._n_ts // 2

        self.window = window

        # define number of trajectory vectors

        self._k = self._n_ts - self.window + 1

    # --------------------------------------------------------
    # Properties


    # --------------------------------------------------------
    # Public methods


    def wcorr(self, components=None):
        """Compute the weighted correlation matrix

        See equation in ref [1], paragraph separability

        Returns
        -------

        References
        ----------

        [1] Hassani, Hossein. "Singular Spectrum Analysis: Methodology and Comparison."
        MPRA Paper, April 1, 2007. https://mpra.ub.uni-muenchen.de/4991/.


        """
        # TODO: optimize for loop

        # check for components type

        if isinstance(components, int):

            comp_idx = range(components)

        elif is_1darray_like(components):

            comp_idx = components

        elif components is None:

            comp_idx = self._xi.keys()

        else:

            raise TypeError('components should be either None, int or '
                            'array-like.')

        # check if components exists

        if not set(comp_idx).issubset(self._xi.keys()):
            raise IndexError('Components are out of range.')

        k = self._k  # number of lagged vectors
        n = self._n_ts  # series length
        w = self.window  # window parameter
        cn = len(comp_idx)  # number of components

        w_k = min((k, w, n - k)) * np.ones(n)

        # intialize w-corr upper right matrix

        wcorr_ur = np.zeros(shape=(cn, cn))

        # reconstruction of selected components

        tsn = np.array([self._hankelmatrix_to_ts(x) for x in self._xi.values()[:cn]])

        # diag offsets

        offsets = range(cn)

        # we compute wcorr using an indices based lagged convolution
        # doing so we compute the upper right correlation matrix
        # then we perform the symmetry

        for offset in offsets:
            lagged_tsn = np.array([tsn[offset:, :], tsn[:cn - offset, :]])

            # weighted sum for components i j lagged by offset
            # see reference for equation

            wsum_ij = np.sum(w_k * np.product(lagged_tsn, axis=0), axis=1)
            sqrtwsum_i = np.sqrt(np.sum(w_k * lagged_tsn[0] ** 2, axis=1))
            sqrtwsum_j = np.sqrt(np.sum(w_k * lagged_tsn[1] ** 2, axis=1))

            rho = wsum_ij / (sqrtwsum_i * sqrtwsum_j)

            # wcorr offset diagonal axis 0 coordinate

            drng = np.array(offsets[:len(rho)])
            wcorr_ur[drng, drng + offset] = rho

        # get lower left symmetry

        wcorr_ll = wcorr_ur.T.copy()

        # remove diagonal values of lower left triangular matrix

        np.fill_diagonal(wcorr_ll, 0.)

        # restore symmetry by addition of upperright and lower left

        wcorr = wcorr_ur + wcorr_ll

        return wcorr

    # --------------------------------------------------------
    # Private methods

    def _embedseries(self):
        """Embed a time series into a L-trajectory matrix
        
        Returns
        -------
        x : np.matrix
            the trajectory matrix of size (window, k)
        
        """

        ts = self.ts
        w = self.window
        n = self._n_ts
        k = self._k

        x = np.zeros(shape=(w, k))

        for i in range(k):
            x[:, i] = ts[i:i + w]

        return np.matrix(x)


    def _reconstruct_group(self, idx):

        u, s, v = self.svd

        x = self._embedseries()

        m, n = x.shape

        x_grp = np.matrix(np.zeros(shape=(m,n)))

        if isinstance(idx, int):
            idx = [idx]

        for i in idx:
            s_i = np.sqrt(s[i]) # TODO: check if I really need to do the square root !!
            u_i = u[:, i]  # eigenvector i corresponding to si
            v_i = x.T * u_i / s_i
            x_i = s_i * u_i * v_i.T
            x_grp += x_i

        # anti diagonal averaging

        ts = self._hankelmatrix_to_ts(x_grp)

        return ts

    @staticmethod
    def _hankelmatrix_to_ts(x):
        """Average the antidiagonal of Hankel matrix to return 1d time series
        
        Parameters
        ----------
        matrix : np.matrix

        Returns
        -------
        
        timeseries: np.array

        """

        ts = [np.mean(x[::-1, :].diagonal(i)) for i in range(-x.shape[0] + 1, x.shape[1])]

        return np.array(ts)


class ToeplitzSSA(BaseSSA):
    pass



if __name__ == '__main__':
    import doctest

    doctest.testmod()
