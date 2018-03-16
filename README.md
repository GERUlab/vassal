# VASSAL: VArious Singular Spectrum AnaLysis with python

UNDER DEVELOPMENT

Author: Damien Delforge

Email: damien.delforge@uclouvain.be

## Todo

* Compatibility with python 3
* Computing the weighted correlation matrix is too long for large time series
* Reconstructing series is too long for long time series
* Make 'paired' plot able to plot any pair, not just consecutive one

## Introduction

### What is Singular Spectrum Analysis?

Singular Spectrum Analysis (SSA) is a non-parametric and model free method for time series decomposition, reconstruction (and foracasting). The general walktrhough of SSA consists in (1) **embedding** the time series into a trajectory matrix of lagged vectors, (2) **decomposing** the trajectory matrix using singular value decomposition (SVD), (3) **grouping** the resulting components based on similarities between their singular values or eigenvectors to reconstruct interpretable components of the original time series. The later is usually supervised.

The main hypothesis behind SSA is separability of the components.

### Potential application of SSA

* Smoothing, filtering, noise reduction
* Structured components extraction (ie. trend or seasonality)
* Filling missing values
* Forecasting
* Nonlinear time series analysis

### Variants of SSA

Different variants of SSA could be declined based either on the embedding method, the decomposition method or the dimensional nature of the dataset of interest.

* **Basic SSA**: The basic 1d SSA algorithm also known as the Broomhead-King variant of SSA (or BK-SSA). The time series is embedded in a Hankel matrix.
* **Toeplitz SSA**: The Toeplitz variant of SSA also known as Vautard-Ghil variant of SSA (or VG-SSA). The time series is embedded in a Toeplitz matrix. Toeplitz SSA should be used when time series is known to be stationary.
* **SSA-ICA**: ICA refers to Independent Component Analysis (ICA) and replace SVD. This variants helps to separate components in case of weak separability. As a less stable procedure than SVD, SSA-ICA is best used in a two stages procedure, a first separation is done using a basic SVD method, then remaining mixed-up components are decomposed using SSA-ICA.

Some other variants are not 1-dimensional such as: **Multichannel SSA** (multiple time series), **2d-SSA** (arrays and images), **nd-SSA** (nd arrays).

## SVD algorithms
TODOC

see: https://stats.stackexchange.com/a/159602/87558

### Comparison with rSSA package: 
see: https://cran.r-project.org/web/packages/Rssa/Rssa.pdf

rSSA package uses either 'nutrlan', 'propack', 'full svd' with the dgesdd routine
and 'eigen' as full SVD via eigendecompsition of the cross-product matrix

see: https://code.lbl.gov/pipermail/trlan-users/2009-May/000007.html

>Here is a little longer answer to your question on comparing ARPACK 
with TRLan.  TRLan (and nuTRLan) implements a restarted version of 
Lanczos algorithm, just like ARPACK implements a restarted version of 
Arnoldi algorithm.  In this regard, the user has control over the 
memory usage by controlling the maximum basis size.  Another 
similarity is that both can keep an arbitrary number of basis vectors 
when restarting -- this is the key advantage of these methods over 
earlier restarted versions.
On symmetric (or Hermitian) problems, when the basis vectors 
corresponding to the same Ritz values are saved during restarting, 
TRLan and ARPACK are theoretically equivalent.  One difference is that 
TRLan uses Ritz vectors while ARPACK uses the vectors produced by the 
implicit QR procedure.  This makes TRLan a little easier to understand 
and implement.  This difference is mainly useful for software 
implementors -- it is of no consequence to the end users.
What do have some consequence are the following.  TRLan can take 
advantage of the symmetry in the original problem as Ichi has pointed 
out.  TRLan and especially nuTRLan use more advanced strategies to 
decide what Ritz values to save during restarting.  These strategies 
have been demonstrated to be very effective.  In general, the 
restarted version of Lanczos would need more matrix-vector 
multiplications than the un-restarted version.  In cases where the 
un-restarted Lanczos can be used, TRLan was shown to use nearly the 
same number of matrix-vector multiplications.  On more difficult 
eigenvalue problems, TRLan usually performed better because of the new 
restarting strategies.

### References for sklearn TruncatedSVD:   

* Finding structure with randomness: Stochastic algorithms for constructing
  approximate matrix decompositions
  Halko, et al., 2009 http://arxiv.org/abs/arXiv:0909.4061
* A randomized algorithm for the decomposition of matrices
  Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert
* An implementation of a randomized algorithm for principal component
  analysis
  A. Szlam et al. 2014

## SVD algorithms in python scientific librairies

* [`numpy.linalg.svd`](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.linalg.svd.html)

The decomposition is performed using LAPACK, with option ´full_matrices´, 

* [`scipy.linalg.svd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html)

The decompotion is performed usin LAPACK, with option ´full_matrices´, additionnaly the user can choose a lapack_driver.

* [`scipy.sparse.linalg.svds`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html)

Compute the largest k singular values/vectors for a sparse matrix. This is a naive implementation using ARPACK as an eigensolver on A.H * A or A * A.H, depending on which one is more efficient.

* [`sklearn.decomposition.TruncatedSVD`](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html#sklearn.decomposition.TruncatedSVD)

This estimator supports two algorithms: a fast randomized SVD solver, and a “naive” algorithm that uses ARPACK as an eigensolver on (X * X.T) or (X.T * X), whichever is more efficient.

## SVD algorithms selection

Both `nplapack` and `splapack` use the LAPACK algorithm for full svd decomposition but the scipy implementation allows more flexibility. Hence it is kept. The same goes with `sparpack` and `skarpack`, as `skarpack` is just a wrapper to `sparpack` with fewer arguments allowed. To my knowledge, `skrandom` has no equivalent. 

Are kept: `splapack`,`sparpack` and `skrandom`.








## Acknowledgment

TODOC

## References

[1] Singular Spectrum Analysis for Time Series | Nina Golyandina | Springer. Accessed November 19, 2017. //www.springer.com/gp/book/9783642349126.
