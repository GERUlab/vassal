"""Singular Spectrum Analysis with numpy

"""
import warnings

import matplotlib.pyplot as plt
import numpy as np

from vassal.devutil.performance import mytimer
from vassal.dtypes import is_1darray_like
from vassal.base import BaseSSA

try:
    import pandas as pd

    __TS_DEFAULT_TYPE__ = pd.Series
except:
    warnings.warn('pandas module missing: __TS__DEFAULT_TYPE__ set to np.array')
    __TS_DEFAULT_TYPE__ = np.array


# TODO: test behavior with np.nan in series
# TODO: add option full_matrices



class BasicSSA(BaseSSA):
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

    def __init__(self, ts=None, window=None, svdmethod='nplapack', tstype=__TS_DEFAULT_TYPE__):



        super(BasicSSA, self).__init__(svdmethod)

        self.ts = np.array(ts)
        self.tstype = tstype
        self._n = len(ts)

        self._matrixgrp = dict()

        # define window length if none

        if window is None:
            window = self._n // 2

        self.window = window

        # define number of trajectory vectors

        self._k = self._n - self.window + 1

        # compute trajectory matrix

        self._x = self._embedseries()

        # set rank of trajectory matrix

        self._xrank = np.linalg.matrix_rank(self._x)

        # add original matrix to matrix group

        self._matrixgrp['Original'] = self._x

        # run decomposition

        #self.decompose()

    def __getitem__(self, item):
        # TODO : error handling
        ts = self._getseries(item)
        return self.tstype(ts)

    # --------------------------------------------------------
    # Properties

    @property
    def groups(self):
        """List of reconstructed group names
        
        Group names are ordered so that the first one corresponds
        the the 'Original' time series group and the last one to
        the 'Residuals'.
        
        """

        groups = self._matrixgrp.keys()

        if 'Original' in groups:

            firstgroup = ['Original']
            others = [n for n in groups if n != 'Original']
            sortedgroups = firstgroup + others

        return sortedgroups

    # --------------------------------------------------------
    # Public methods

    def reconstruct(self, groups=None):

        # TODO: DOC
        print self._xi

        # Define a list of group indexes

        if groups is None:
            idx_list = [range(len(self._xi))]
            names = ['reconstruction']
        else:
            idx_list = [i for i in groups.values()]
            names = [name for name in groups.keys()]

        for name, idx_grp in zip(names, idx_list):
            x_group = [self._xi[key] for key in idx_grp]
            x_sum = np.sum(x_group, axis=0)

            self._matrixgrp[name] = x_sum

        all_grp_idx = [ix for sublist in idx_list for ix in sublist]

        residual_idx = [ix for ix in range(len(self._xi)) if ix not in all_grp_idx]

        x_res = [self._xi[ix] for ix in residual_idx]
        x_res_sum = np.sum(x_res, axis=0)

        if bool(x_res_sum.any()):
            self._matrixgrp['Residuals'] = x_res_sum

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

            raise TypeError('components should be either None, int or array-like.')

        # check if components exists

        if not set(comp_idx).issubset(self._xi.keys()):
            raise IndexError('Components are out of range.')

        k = self._k  # number of lagged vectors
        n = self._n  # series length
        w = self.window  # window parameter
        cn = len(comp_idx)  # number of components

        w_k = min((k, w, n - k)) * np.ones(n)

        # intialize w-corr upper right matrix

        wcorr_ur = np.zeros(shape=(cn, cn))

        # reconstruction of selected components

        tsn = np.array([self._antidiagmean(x) for x in self._xi.values()[:cn]])

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
        n = self._n
        k = self._k

        x = np.zeros(shape=(w, k))

        for i in range(k):
            x[:, i] = ts[i:i + w]

        return np.matrix(x)

    """
    def decompose(self):

        # rank of the trajectory matrix x
        d = self._xrank
        x = self._x
        assert (d == min(self._x.shape))

        # decomposition of the trajectory matrix x
        # u and v are unitary and s is a 1-d array of d singular values.

        u, s, v = np.linalg.svd(self._x, full_matrices=False)

        # note: types are all np.matrix

        self._xi = dict()

        for i in range(d):
            si = np.sqrt(s[i])  # square root of eigenvalue i
            ui = u[:, i]  # eigenvector i corresponding to si
            vi = x.T * ui / si

            self._xi[i] = si * ui * vi.T

        self.svd = [u, s, v]
    """
    @property
    def _xi(self):
        d = self._xrank
        x = self._x
        xi = dict()
        u, s, v = self.svd
        for i in range(d):
            si = np.sqrt(s[i])  # square root of eigenvalue i
            ui = u[:, i]  # eigenvector i corresponding to si
            vi = x.T * ui / si
            self._xi[i] = si * ui * vi.T
        return xi

    def _getseries(self, name):

        x = self._matrixgrp[name]

        # anti diagonal averaging

        ts = self._antidiagmean(x)

        return ts

    @staticmethod
    def _antidiagmean(x):
        """Average the antidiagonal of matrix
        
        Parameters
        ----------
        matrix : np.matrix

        Returns
        -------
        
        timeseries: np.array

        """

        ts = [np.mean(x[::-1, :].diagonal(i)) for i in range(-x.shape[0] + 1, x.shape[1])]

        return np.array(ts)

    # --------------------------------------------------------
    # Plotting methods

    def plot(self, pltname='values', show=True, **pltkw):

        if pltname not in self._plotnames:
            names = ','.join(self._plotnames)
            raise AttributeError('Unknown plot name \'{}\'. Name should be on of {}.'.format(pltname, names))

        elif pltname == 'values':
            fig, ax = self._value_plot(**pltkw)

        elif pltname == 'reconstruction':
            fig, ax = self._reconstruction_plot(**pltkw)

        elif pltname == 'wcorr':
            fig, ax = self._wcorr_plot(**pltkw)

        elif pltname == 'vectors':
            fig, ax = self._vectors_plot(**pltkw)

        elif pltname == 'paired':
            fig, ax = self._paired_plot(**pltkw)

        else:
            raise NotImplementedError('{} not implemented'.format(pltname))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if show is True:
            plt.show()

        return fig, ax

    @property
    def _plotnames(self):
        names = [
            'values',
            'reconstruction',
            'wcorr',
            'vectors',
            'paired'
        ]
        return names

    def _value_plot(self, n=50, **pltkw):

        # eigenvalues
        eigenvalues = self.svd[1] # TODO: check if needed to raise power

        #
        fig = plt.figure()
        ax = fig.gca()
        ax.semilogy(eigenvalues[:n], '-ok', markersize=4., alpha=0.5)
        ax.set_ylabel('Component Norms')
        ax.set_xlabel('Index')

        return fig, ax

    def _reconstruction_plot(self, **pltkw):

        groups = self.groups

        if len(groups) > 1:

            fig, axarr = plt.subplots(len(groups), 1, sharex=True)

            for i, g in enumerate(groups):
                ts = self._getseries(g)
                axarr[i].plot(ts, **pltkw)
                axarr[i].set_title(g)
                if i == len(groups) - 1:
                    axarr[i].set_xlabel('Time')
        else:
            pass

        return fig, axarr

    def _wcorr_plot(self, n=20, *args, **kwargs):

        wcorr = self.wcorr(components=n)

        fig = plt.figure()
        ax = fig.gca()
        im = ax.pcolor(wcorr, vmin=-1, vmax=1, cmap='PiYG')
        ax.set_aspect('equal')

        # set ticks

        ticks = np.arange(wcorr.shape[0])
        ax.set_xticks(ticks + 0.5, minor=False)
        ax.set_yticks(ticks + 0.5, minor=False)

        ax.set_xticklabels(int(x) for x in ticks)
        ax.set_yticklabels(int(x) for x in ticks)

        ax.set_title('w-correlation matrix')

        fig.colorbar(im)

        return fig, ax

    def _vectors_plot(self, n = 10, **pltkw):
        """
        The rows of v are the eigenvectors of a.H a. The columns of u are the 
        eigenvectors of a a.H. For row i in v and column i in u, the 
        corresponding eigenvalue is s[i]**2.
        
        Parameters
        ----------
        n

        Returns
        -------

        """
        # TODO: type error
        u = self.svd[0]
        s = self.svd[1]**2 # TODO: check if power is needed

        fig = plt.figure()

        # grid size

        m = int(np.ceil(np.sqrt(n)))

        ax = None

        for i in range(n):
            ax = plt.subplot(m, m, i+1, sharey=ax)
            ax.plot(u[:,i], **pltkw)
            ax.set_xticks([])
            ax.set_yticks([])

            contribution =  s[i]/np.sum(s)*100

            title = 'EV{0} ({1:.0f} %)'.format(i+1, contribution)

            ax.set_title(title, {'fontsize':10.})

        fig.suptitle('Eigenvectors plot')

        return fig, fig.get_axes()

    def _paired_plot(self, pairs=zip(range(0,9),range(1,10)), **pltkw):

        # TODO: check type pairs list of tuple of size 2

        print pairs
        u = self.svd[0]
        s = self.svd[1]**2 # TODO: check if power is needed

        fig = plt.figure()

        m = int(np.ceil(np.sqrt(len(pairs))))

        ax = None

        for i, j in pairs:
            ax = plt.subplot(m, m, i+1, sharey=ax)
            ax.plot(u[:,j], u[:,i], **pltkw)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect('equal')

            ssum = np.sum(s)

            contribution1 =  s[i]/ssum*100
            contribution2 =  s[j]/ssum*100

            title = 'EV{0} ({1:.0f}%) vs EV{2} ({3:.0f}%)'.format(
                i+1,
                contribution1,
                i+2,
                contribution2
            )

            ax.set_title(title, {'fontsize':10.})

            fig.suptitle('Pairs of eigenvectors')

        return fig, fig.get_axes()



if __name__ == '__main__':
    import doctest

    doctest.testmod()
