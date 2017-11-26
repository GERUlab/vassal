# A Python implementation of various Singular Spectrum Analysis Algorithms

Important Note:

_This package is under development and follows a documentation driven development strategy. This is a DRAFT README and information related to the package corresponds to the end-stage of development._

## Introduction

### What is Singular Spectrum Analysis?

Singular Spectrum Analysis (SSA) is a non-parametric and model free method for time series decomposition, reconstruction (and foracasting). The general walktrhough of SSA consists in (1) **embedding** the time series into a trajectory matrix of lagged vectors, (2) **decomposing** the trajectory matrix using singular value decomposition (SVD), (3) **grouping** the resulting components based on similarities between their singular values or eigenvectors to reconstruct interpretable components of the original time series. The later is usually supervised.

The main hypothesis behind SSA is separability of the components.

### Potential application of SSA?

* Smoothing, filtering, noise reduction
* Structured components extraction (ie. trend or seasonality)
* Filling missing values
* Forecasting
* Nonlinear time series analysis

### Variants of SSA

Different variants of SSA could be declined based either on the embedding method or the decomposition method.

* **Basic SSA**: The basic 1d SSA algorithm also known as the Broomhead-King variant of SSA (or BK-SSA). The time series is embedded in a Hankel matrix.
* **Toeplitz SSA**: The Toeplitz variant of SSA also known as Vautard-Ghil variant of SSA (or VG-SSA). The time series is embedded in a Toeplitz matrix. Toeplitz SSA should be used when time series is known to be stationary.
* **SSA-ICA**: ICA refers to Independent Component Analysis (ICA) and replace SVD. This variants helps to separate components in case of weak separability. As a less stable procedure than SVD, SSA-ICA is best used in a two stages procedure, a first separation is done using a basic SVD method, then remaining mixed-up components are decomposed using SSA-ICA.


## Acknowledgment

TODOC

## References

[1] Singular Spectrum Analysis for Time Series | Nina Golyandina | Springer. Accessed November 19, 2017. //www.springer.com/gp/book/9783642349126.


