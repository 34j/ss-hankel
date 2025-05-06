from typing import Any

import numpy as np
from numpy.typing import NDArray


def score(
    feigval: NDArray[Any], eigvec: NDArray[Any], /, ord: Any = None
) -> NDArray[Any]:
    """
    Compute the accuracy of the solution of the eigenvalue problem.

    |F(λ)v|/|F(λ)||v|

    Parameters
    ----------
    feigval : NDArray[Any]
        The function evaluated at the eigenvalues (..., n_eig, n, n).
    eigvec : NDArray[Any]
        The eigenvectors of shape (..., n, n_eig).
    ord : Any, optional
        Order of the norm to be passed to np.linalg.norm, by default None

    Returns
    -------
    NDArray[Any]
        The score of the solution.

    """
    if (
        feigval.shape[-2] == feigval.shape[-1] == eigvec.shape[-2]
        and feigval.shape[-3] == eigvec.shape[-1]
        and feigval.shape[:-3] == eigvec.shape[:-2]
    ):
        pass
    else:
        raise ValueError(
            "feigval should be of shape (..., n_eig, n, n)"
            " and eigvec should be of shape (..., n, n_eig), "
            f"but got {feigval.shape=} and {eigvec.shape=}"
        )

    eigvec = np.swapaxes(eigvec, -1, -2)
    return (
        np.linalg.norm(np.matmul(feigval, eigvec[..., None])[..., 0], axis=-1)
        / np.linalg.norm(eigvec, axis=-1)
        / np.linalg.norm(feigval, axis=(-2, -1))
    )
