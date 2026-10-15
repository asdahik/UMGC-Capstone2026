# this is a helper script that is kept from the main ml script for readability and modularity for other projects

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class PCAImputer(BaseEstimator, TransformerMixin):
    """Iterative low-rank imputation using singular value decomposition (SVD), imitating NIPALS-PCA imputation from reference paper.
        fit() learns the loadings from training data while transform() fills NaNs in new rows using those fixed loadings"""

    # initialize with established parameters
    def __init__(self, k=5, iters=200, tol=1e-6):
        print("initializing imputation")
        self.k, self.iters, self.tol = k, iters, tol

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        miss = np.isnan(X)
        Z = np.where(miss, np.nanmean(X, axis=0), X)
        for _ in range(self.iters):
            m = Z.mean(0)
            U, S, Vt = np.linalg.svd(Z-m, full_matrices=False)
            R = (U[:, :self.k] * S[:self.k] @ Vt[:self.k] + m)
            Zn = np.where(miss, R, X)
            done = np.abs(Zn - Z).max() < self.tol
            Z = Zn
            if done:
                break
        self.col_mean_ = np.nanmean(X, axis=0)
        self.mean_ = Z.mean(0)
        self.V_ = Vt[:self.k]
        self.train_filled_ = Z
        return self
    
    def transform(self, X):
        print("running transform")
        X = np.asarray(X, float)
        miss = np.isnan(X)
        Z = np.where(miss, self.col_mean_, X)
        for _ in range(self.iters):
            R = ((Z - self.mean_) @ self.V_.T) @ self.V_ + self.mean_
            Zn = np.where(miss, R, X)
            done = np.abs(Zn - Z).max() < self.tol
            Z = Zn
            if done:
                break
        return Z