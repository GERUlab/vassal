"""Base class for SSA objects

"""

# Get svd algorithm from numpy scipy and sklearn

import abc
import numpy as np

# Get svd algorithm from numpy scipy and sklearn
from numpy.linalg import svd as nplapack
from scipy.linalg import svd as splapack
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import svds as sparpack
from sklearn.utils.extmath import randomized_svd
# utils
from sklearn.utils.extmath import svd_flip


class BaseSSA(object):
    """Base class of SSA object
    
    BaseSSA is a base class for SSA objects. When initiated, BaseSSA map the svd 
    method selected by the user to the decomposition method.
    
    Existing SVD algorithm in python are wrapped to this base class to ensure
    consisency in the outputs generated by the algoritm: Singular vectors are 
    assign to the object as np.matrix type while singular values are 1d np.array.
    
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, svdmethod):
        """Base class init method
        """

        # reference singular value decomposition results
        # 0: eigenvectors
        # 1: eignevalues
        # 2: factorvectors #TODO: check for more precise references

        self.svd = [None, None, None]

        self.decompose = self._SVD_METHODS_MAP[svdmethod]

    @property
    def _SVD_METHODS_MAP(self):
        """Map user-selected svd method to the proper wrapper
        """
        svdmap = {
            'nplapack': self._nplapack_wrapper,
            'splapack': self._splapack_wrapper,
            'sparpack': self._sparpack_wrapper,
            'skrandom': self._skrandom_wrapper
        }
        return svdmap

    @abc.abstractmethod
    def _embedseries(self):
        """Time Series Embedding abstract method"""
        pass

    def _nplapack_wrapper(self, full_matrices=True):
        """Wrapper for numpy.linalg.svd
               
        Apply SVD to the embedding matrix of shape (`M`, `N`) using the 
        `numpy.linalg.svd`_ algorithm based on LAPACK implementation of SVD. 
        
        Parameters
        ----------
        
        full_matrices : bool, optional
            If True (default), `u` and `v` have the shapes (`M`, `M`) and
            (`N`, `N`), respectively.  Otherwise, the shapes are (`M`, `K`)
            and (`K`, `N`), respectively, where `K` = min(`M`, `N`).
            
        See Also
        --------
        
        .. _`numpy.linalg.svd`:
           https://docs.scipy.org/doc/numpy-1.11.0/reference/generated/numpy.linalg.svd.html
            
        """

        # Matrix to be decomposed

        x = self._embedseries()

        # Apply decomposition

        u, s, v = nplapack(x, full_matrices=full_matrices, compute_uv=True)

        # saving svd

        self.svd = [np.matrix(u), s, np.matrix(v)]

        return self.svd

    def _splapack_wrapper(self, full_matrices=True, check_finite=False,
                          lapack_driver='gesdd'):
        """Wrapper for scipy.linalg.svd
               
        Apply SVD to the embedding matrix of shape (`M`, `N`) using the 
        `scipy.linalg.svd`_ algorithm based on LAPACK implementation of SVD. 
        
        Parameters
        ----------
        full_matrices : bool, optional
            If True (default), `U` and `Vh` are of shape ``(M, M)``, ``(N, N)``.
            If False, the shapes are ``(M, K)`` and ``(K, N)``, where
            ``K = min(M, N)``.
        check_finite : bool, optional
            Whether to check that the input matrix contains only finite numbers.
            Disabling may give a performance gain, but may result in problems
            (crashes, non-termination) if the inputs do contain infinities or NaNs.
        lapack_driver : {'gesdd', 'gesvd'}, optional
            Whether to use the more efficient divide-and-conquer approach
            (``'gesdd'``) or general rectangular approach (``'gesvd'``)
            to compute the SVD. MATLAB and Octave use the ``'gesvd'`` approach.
            Default is ``'gesdd'``.

        See Also
        --------
        
        .. _`scipy.linalg.svd`:
           https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html
            

        """
        # Matrix to be decomposed

        x = self._embedseries()

        # Decomposition

        u, s, v = splapack(x, full_matrices=full_matrices,
                           compute_uv=True,
                           overwrite_a=True,
                           check_finite=check_finite,
                           lapack_driver=lapack_driver)

        self.svd = [np.matrix(u), s, np.matrix(v)]

        return self.svd

    def _sparpack_wrapper(self, k=None, ncv=None, tol=0, v0=None, maxiter=None):
        """Wrapper for scipy.sparse.linalg.svds

        Apply Singular Value Decomposition to the embedding matrix of shape 
        (`M`, `N`) using the `scipy.sparse.linalg.svds`_ algorithm. 

        Parameters
        ----------


        See Also
        --------
        
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html

        """
        # Matrix to be decomposed

        x = csc_matrix(self._embedseries())

        # Default k value is full svd

        if k is None:
            k = min(x.shape) - 1

        u, s, v = sparpack(x, k=k, ncv=ncv, tol=tol, which='LM', v0=v0,
                           maxiter=maxiter, return_singular_vectors=True)

        # with this implementation vectors needs to be flipped to match lapack
        # other format and sign ambiguities have to be solved to force
        # deterministic output with svd_flip

        u, v = svd_flip(u[:, ::-1], v[::-1, :])

        self.svd = [np.matrix(u), s[::-1], np.matrix(v)]

        return self.svd

    def _skrandom_wrapper(self, k=None, n_oversamples=10, n_iter='auto',
                          power_iteration_normalizer='auto', random_state=None):
        """Wrapper to sklearn.utils.extmath.randomized_svd
        
        Apply Singular Value Decomposition to the embedding matrix of shape 
        (`M`, `N`) using the `sklearn.utils.extmath.randomized_svd`_ algorithm. 
            
        Parameters
        ----------
        
        k : int
            Number of singular values and vectors to extract.
        n_oversamples : int (default is 10)
            Additional number of random vectors to sample the range of M so as
            to ensure proper conditioning. The total number of random vectors
            used to find the range of M is n_components + n_oversamples. Smaller
            number can improve speed but can negatively impact the quality of
            approximation of singular vectors and singular values.
        n_iter : int or 'auto' (default is 'auto')
            Number of power iterations. It can be used to deal with very noisy
            problems. When 'auto', it is set to 4, unless `n_components` is small
            (< .1 * min(X.shape)) `n_iter` in which case is set to 7.
            This improves precision with few components.
            .. versionchanged:: 0.18
        power_iteration_normalizer : 'auto' (default), 'QR', 'LU', 'none'
            Whether the power iterations are normalized with step-by-step
            QR factorization (the slowest but most accurate), 'none'
            (the fastest but numerically unstable when `n_iter` is large, e.g.
            typically 5 or larger), or 'LU' factorization (numerically stable
            but can lose slightly in accuracy). The 'auto' mode applies no
            normalization if `n_iter`<=2 and switches to LU otherwise.
            .. versionadded:: 0.18
        random_state : int, RandomState instance or None, optional (default=None)
            The seed of the pseudo random number generator to use when shuffling
            the data.  If int, random_state is the seed used by the random number
            generator; If RandomState instance, random_state is the random number
            generator; If None, the random number generator is the RandomState
            instance used by `np.random`.
            
        See Also
        -------
        
        *  .. _`sklearn.utils.extmath.randomized_svd`:
            https://github.com/scikit-learn/scikit-learn/blob/439bf1ac896aacafe90177bc28bf036cc8e8e4d9/sklearn/utils/extmath.py#L228
            
        * .. _`sklearn.decomposition.TruncatedSVD`:
            http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html
            

        """

        # Matrix to be decomposed

        x = csc_matrix(self._embedseries())

        # if k is None get the maximum

        if k is None:
            k = min(x.shape) - 1

        if random_state is None:
            random_state = np.random.RandomState()

        # Sklearn randomized svd decomposition

        u, s, v = randomized_svd(x, n_components=k, n_oversamples=n_oversamples,
                                 n_iter=n_iter,
                                 power_iteration_normalizer=power_iteration_normalizer,
                                 transpose=False, flip_sign=True,
                                 random_state=random_state)

        # store output

        self.svd = [np.matrix(u), s, np.matrix(v)]

        return self.svd


if __name__ == '__main__':
    from vassal.ssa import BasicSSA
    import numpy as np

    np.random.seed(0)
    s = np.random.randint(low=0, high=10, size=200)

    npssa = BasicSSA(s, svdmethod='nplapack')
    spssa = BasicSSA(s, svdmethod='splapack')
    sp2ssa = BasicSSA(s, svdmethod='sparpack')
    skssa = BasicSSA(s, svdmethod='skrandom')
    u1, s1, v1 = npssa.decompose()
    u2, s2, v2 = spssa.decompose()
    u3, s3, v3 = sp2ssa.decompose()
    u4, s4, v4 = skssa.decompose()
    print np.allclose(s1, s2)
    print np.allclose(s2[:-1], s3)
    print np.allclose(s3, s4)
    skssa['Original'].plot()

    groups = {
        'trend': 0,
        'season': [1, 2, 3]
    }

    skssa.reconstruct(groups)
    skssa.plot('reconstruction')