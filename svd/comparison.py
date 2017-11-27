"""Comparison of existing svd algorithm

"""

import time
import numpy as np
import sklearn.decomposition
import seaborn
import matplotlib.pyplot as plt
import pandas as pd

RANK = 50
N_COLS = 1000


def evaluate_svd(svd_fn, reconstruct_fn, min_rows=100, max_rows=5000,
                 n_samples=100, n_cols=N_COLS, rank=RANK, random_seed=0):
    np.random.seed(random_seed)
    elapsed_times = []
    errors = []
    n_rows_array = np.random.randint(low=min_rows, high=max_rows + 1,
                                     size=n_samples)
    n_rows_array.sort()

    for n_rows in n_rows_array:
        # construct a low-rank matrix
        left = np.random.randn(n_rows, rank)
        right = np.random.randn(rank, n_cols)
        full = np.dot(left, right)

        # how long does it take to perform the SVD?
        start_t = time.time()
        svd_outputs = svd_fn(full)
        end_t = time.time()
        elapsed_t = end_t - start_t
        elapsed_times.append(elapsed_t)

        # compute mean absolte error of reconstruction
        reconstructed = reconstruct_fn(svd_outputs)
        diff = full - reconstructed
        mae = np.mean(np.abs(diff))
        errors.append(mae)
        print(
        "n_rows=%d ==> time = %0.4f, MAE = %0.8f" % (n_rows, elapsed_t, mae))
    max_error = np.max(errors)
    print("Max Error=%f" % max_error)
    assert max_error < 0.0000001
    return n_rows_array, elapsed_times, errors


# Full SVD with NumPy

def np_svd(X):
    return np.linalg.svd(X, full_matrices=False, compute_uv=True)


def np_inv_svd(svd_outputs):
    U, s, V = svd_outputs
    return np.dot(U, np.dot(np.diag(s), V))


# Truncated SVD with scikit-learn

def sklearn_svd(X, rank=RANK):
    tsvd = sklearn.decomposition.TruncatedSVD(rank)
    X_reduced = tsvd.fit_transform(X)
    return (tsvd, X_reduced)


def sklearn_inv_svd(svd_outputs):
    tsvd, X_reduced = svd_outputs
    return tsvd.inverse_transform(X_reduced)


# Perform timings

n_rows, np_times, np_errors = evaluate_svd(np_svd, np_inv_svd)


def sklearn_randomized_svd(X, rank=RANK):
    tsvd = sklearn.decomposition.TruncatedSVD(rank, algorithm="randomized",
                                              n_iter=1)
    X_reduced = tsvd.fit_transform(X)
    return (tsvd, X_reduced)


def sklearn_arpack_svd(X, rank=RANK):
    tsvd = sklearn.decomposition.TruncatedSVD(rank, algorithm="arpack")
    X_reduced = tsvd.fit_transform(X)
    return (tsvd, X_reduced)


def sklearn_inv_svd(svd_outputs):
    tsvd, X_reduced = svd_outputs
    return tsvd.inverse_transform(X_reduced)


n_rows, sklearn_randomized_times, sklearn_randomized_errors = evaluate_svd(
    sklearn_randomized_svd, sklearn_inv_svd)
n_rows, sklearn_arpack_times, sklearn_arpack_errors = evaluate_svd(
    sklearn_arpack_svd, sklearn_inv_svd)

figure = seaborn.plt.figure(figsize=(10, 10))
seaborn.plt.xlim(0, 5000)
seaborn.plt.ylim(0, 3)
seaborn.regplot(x=pd.Series(n_rows, name="n_rows"),
                y=pd.Series(np_times, name="elapsed time (s)"))
seaborn.regplot(x=pd.Series(n_rows, name="n_rows"),
                y=pd.Series(sklearn_randomized_times, name="elapsed time (s)"))
seaborn.regplot(x=pd.Series(n_rows, name="n_rows"),
                y=pd.Series(sklearn_arpack_times, name="elapsed time (s)"))

seaborn.plt.legend(
    ("numpy.linalg.svd", "TruncatedSVD (randomized)", "TruncatedSVD (arpack)"))
seaborn.plt.title(
    "Time to perform SVD (on matrices of rank 50 with 1000 columns)")


plt.show()