import warnings
from typing import Any

import attrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import NDArray


@attrs.frozen(kw_only=True)
class Multiplicity:
    value: complex
    """The eigenvalue."""
    algebraic_multiplicity: int
    """The algebraic multiplicity of the eigenvalue.
    The number of times the eigenvalue appears
    as a root of the characteristic polynomial."""
    geometric_multiplicity: int | None
    """The geometric multiplicity of the eigenvalue.
    The dimension of the eigenspace of the eigenvalue.
    Less than or equal to the algebraic multiplicity."""


def get_multiplicity(
    eigval: "NDArray[Any]",
    eigvec: "NDArray[Any] | None" = None,
    atol: float | None = None,
) -> list[Multiplicity]:
    """
    Get the multiplicity of the eigenvalue.

    Does not support batched eigenvalues.

    Parameters
    ----------
    eigval : Array | NativeArray
        The eigenvalues.
    eigvec : Array | NativeArray | None, optional
        The eigenvectors, by default None.
    atol : float | None, optional
        The threshold to treat eigenvalues as the same.

    Returns
    -------
    int
        The multiplicity of the eigenvalue.

    """
    atol = atol or 0
    if eigval.ndim != 1:
        raise ValueError("eigval should be 1D array.")
    if eigvec is not None:
        if eigvec.ndim != 2:
            raise ValueError("eigvec should be 2D array.")
        if eigval.shape[0] != eigvec.shape[1]:
            raise ValueError(
                f"{eigval.shape[0]=} should be equal to {eigvec.shape[1]=}."
            )
    eigval_dists = np.abs(eigval[:, None] - eigval[None, :])
    eigval_dists_close = eigval_dists < atol
    eigval_index = list(range(eigval.shape[0]))
    result = []
    while eigval_index:
        i = eigval_index.pop(0)
        close_index = eigval_dists_close[i, :].nonzero()[0]
        for j in close_index:
            if i == j:
                continue
            elif j in eigval_index:
                eigval_index.remove(j)
            else:
                warnings.warn(
                    "atol is too large or too small.", RuntimeWarning, stacklevel=2
                )
        algebraic_multiplicity = len(close_index)
        if eigvec is None:
            geometric_multiplicity = None
        else:
            geometric_multiplicity = np.linalg.matrix_rank(
                eigvec[:, close_index], tol=atol
            )
        result.append(
            Multiplicity(
                value=eigval[i],
                algebraic_multiplicity=algebraic_multiplicity,
                geometric_multiplicity=geometric_multiplicity,
            )
        )
    return result


def plot_eigval(
    eigval: "NDArray[Any]",
    eigvec: "NDArray[Any] | None" = None,
    atol: float | None = None,
    ax: Axes | None = None,
) -> None:
    """
    Plot eigenvalues with annotation of the multiplicity.

    Does not support batched eigenvalues.

    Parameters
    ----------
    eigval : Array | NativeArray
        The eigenvalues.
    eigvec : Array | NativeArray | None, optional
        The eigenvectors, by default None.
    atol: float | None, optional
        The threshold to treat eigenvalues as the same.
    ax : plt.Axes | None, optional
        The axes to plot, by default None.

    """
    multiplicity = get_multiplicity(eigval, eigvec, atol)
    ax_ = ax or plt.gca()
    ax_.plot(eigval.real, eigval.imag, "x")
    for m in multiplicity:
        ax_.text(
            m.value.real,
            m.value.imag,
            f"{m.algebraic_multiplicity}"
            + (
                f",{m.geometric_multiplicity}"
                if m.geometric_multiplicity is not None
                else ""
            ),
        )
    ax_.set_xlabel("Re")
    ax_.set_ylabel("Im")
