"""Singular Spectrum Analysis Class

"""
import warnings

import numpy as np

from vassal.base import BaseSSA
from vassal.plot import PlotSSA

try:
    import pandas as pd

    __TS_DEFAULT_TYPE__ = 'pdseries'
except ImportError:
    pd = None
    warnings.warn('pandas module missing: __TS__DEFAULT_TYPE__ set to np.array')
    __TS_DEFAULT_TYPE__ = 'nparray'

class BasicSSA(BaseSSA, PlotSSA):
    """A class for basic Singular Spectrum Analysis 
    """

    def __init__(self, ts=None, window=None, svdmethod='nplapack',
                 usetype=__TS_DEFAULT_TYPE__):
        """Basic Singular Spectrum Analysis
        
        Basic Singular Spectrum Analysis embed the times series into an Hankel
        matrix (i.e. anti-diagonal elements are all equal). For a time series
        :math:`Y =\{y_1, y_2, ..., y_N\}` of length :math:`N`, 
        :math:`Y` is embedded with respect to the window parameter :math:`L` 
        into a matrix :math:`X` of shape :math:`(L,K)`:
        
        .. math::
        
            X= \\begin{bmatrix}
                    y_1     & y_2     & y_3     & \dots   & y_K      \\\\
                    y_2     & y_3     & y_4     & \dots   & y_{K+1}  \\\\
                    y_3     & y_4     & y_5     & \dots   & y_{K+2}  \\\\
                    \\vdots & \\vdots & \\vdots & \\ddots & \\vdots   \\\\
                    y_L     & y_{L+1} & y_{L+2} & \dots   & y_N
                \end{bmatrix}
        
        Parameters
        ----------
        ts : arraylike
            One dimensional array-like object (np.array, dict, list, pd.Series) 
            holding the time series values. If ts is as pandas.Series object and
            if svdmethod is set to 'pdseries' index is kept and pass the the 
            results of SSA.
        window : int, optionnal
            The window parameter of basic SSA.
        svdmethod : str, optionnal
            Lorem Ipsum.
        usetype : str, optionnal
            Lorem Ipsums


        Examples
        --------
        
        * Decomposition
        
        >>> import numpy as np
        >>> import pandas as pd
        >>> pd.options.display.float_format = '{:,.2f}'.format
        >>> np.random.seed(0)
        >>> timeseries = np.random.randint(low=0, high=10, size=100)
        >>> myssa = BasicSSA(ts=timeseries)
        >>> u, s, v = myssa.decompose()
        >>> # Print the 3 first singular values
        >>> print s[:3]
        [ 220.04832942   42.43314868   42.11375246]
        
        * Grouping and reconstruction
        
        >>> groups = {'trend':0, 'signal':range(1,10)}
        >>> myssa.reconstruct(groups)
        >>> myssa['trend'].describe()
        count   100.00
        mean      4.35
        std       0.09
        min       4.09
        25%       4.32
        50%       4.36
        75%       4.40
        max       4.51
        dtype: float64
        
        >>> myssa.to_frame().describe()
               trend  signal  ssa_original  ssa_reconstruction  ssa_residuals
        count 100.00  100.00        100.00              100.00         100.00
        mean    4.35   -0.08          4.33                4.33           0.05
        std     0.09    2.03          2.82                2.82           1.83
        min     4.09   -4.80          0.00                0.00          -5.00
        25%     4.32   -1.49          2.00                2.00          -1.16
        50%     4.36   -0.10          4.00                4.00           0.10
        75%     4.40    1.47          7.00                7.00           1.17
        max     4.51    4.64          9.00                9.00           4.31

        """

        # we pass the svd method to the base class init in order to wrap the
        # proper decompose method

        super(BasicSSA, self).__init__(ts=ts, svdmethod=svdmethod,
                                       usetype=usetype)

        # define window length if none

        if window is None:
            window = self._n_ts // 2

        self.window = window

        # define number of trajectory vectors

        self._k = self._n_ts - self.window + 1

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
        k = self._k

        x = np.zeros(shape=(w, k))

        for i in range(k):
            x[:, i] = ts[i:i + w]

        return np.matrix(x)

    def _reconstruct_group(self, idx):

        u, s, v = self.svd

        x = self._embedseries()

        m, n = x.shape

        x_grp = np.matrix(np.zeros(shape=(m, n)))

        if isinstance(idx, int):
            idx = [idx]

        for i in idx:
            s_i = np.sqrt(
                s[i])  # TODO: check if I really need to do the square root !!
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
        x : np.matrix

        Returns
        -------
        
        timeseries: np.array

        """

        ts = [np.mean(x[::-1, :].diagonal(i)) for i in
              range(-x.shape[0] + 1, x.shape[1])]

        return np.array(ts)


class ToeplitzSSA(BaseSSA):
    def _reconstruct_group(self, grpidx):
        pass

    def _embedseries(self):
        pass


if __name__ == '__main__':
    import doctest

    doctest.testmod()
