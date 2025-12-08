from collections.abc import Callable
from typing import Any

import attrs
import numpy as np
from numpy.typing import NDArray

from ._main import CircleResult

ResultType = tuple[NDArray[Any], NDArray[Any]]


@attrs.frozen(kw_only=True)
class GridResult[TResult: ResultType](CircleResult):
    result: TResult
    """The result returned by `solve_nep` for each grid circle."""


def solve_nep_grid[TResult: ResultType](
    f: Callable[["NDArray[Any]"], "NDArray[Any]"],
    solve_nep: Callable[
        [
            Callable[["NDArray[Any]"], "NDArray[Any]"],
            NDArray[np.complexfloating[Any]],
            NDArray[np.floating[Any]],
        ],
        TResult,
    ],
    search_range: tuple[complex, complex],
    window: complex,
    /,
) -> GridResult[TResult]:
    """
    Sakurai-Sugiura method for the circle.

    Parameters
    ----------
    f : Callable[[np.ndarray], np.ndarray]
        An analytic function (F(z)).
        Array of shape [circle_n_points] will be passed
        and should return [circle_n_points, ..., n, n] array.
    solve_nep : Callable[Callable[[np.ndarray], np.ndarray],
        np.ndarray, np.ndarray, TResult]
        Given a matrix-valued function, center, and radius,
        returns eigenvalues and eigenvectors
        of the function inside the circle defined by the center and radius.
    search_range : tuple[complex, complex]
        The range of the grid.
    window : complex
        The window of the grid.
    num_vectors : int, optional
        Number of linearly independent vectors (L), by default None.

    Returns
    -------
    GridResult[TResult]
        The eigenvalues and eigenvectors.

    """
    if len(search_range) != 2:
        raise ValueError("search_range should be tuple of length 2.")
    if window.real <= 0 or window.imag <= 0:
        raise ValueError(
            "The real and imaginary part of window should be both positive."
        )
    grid_count_real = np.ceil(
        abs((search_range[1] - search_range[0]).real) / window.real
    )
    grid_count_imag = np.ceil(
        abs((search_range[1] - search_range[0]).imag) / window.imag
    )
    window_real = abs((search_range[1] - search_range[0]).real) / grid_count_real
    window_imag = abs((search_range[1] - search_range[0]).imag) / grid_count_imag
    circle_radius = np.sqrt(window_real**2 + window_imag**2) / 2
    circle_center_real = (
        min(search_range[0].real, search_range[1].real)
        + (0.5 + np.arange(grid_count_real)[:, None]) * window_real
    )
    circle_center_imag = (
        min(search_range[0].imag, search_range[1].imag)
        + (0.5 + np.arange(grid_count_imag)[None, :]) * window_imag
    )
    circle_center = (circle_center_real + 1j * circle_center_imag).flatten()
    # array of arrays
    n = None

    def f_(*args: Any, **kwargs: Any) -> NDArray[Any]:
        nonlocal n
        res = f(*args, **kwargs)
        n = res.shape[-1]
        return res

    res = solve_nep(
        f_,
        circle_center,
        circle_radius[None],
    )

    if n is None:
        raise AssertionError("n should be set in f_.")

    eigvals, eigvecs = res
    eigvals_merged = []
    eigvecs_merged = []
    for i in range(eigvals.shape[-1]):
        eigval = eigvals[..., i].item()
        eigvec = eigvecs[..., i].item()
        if eigval is None:
            eigval = np.empty(0)
            eigvec = np.empty((n, 0))
        mask = (
            (eigval.real < circle_center[..., i].real + window_real / 2)
            & (eigval.real > circle_center[..., i].real - window_real / 2)
            & (eigval.imag < circle_center[..., i].imag + window_imag / 2)
            & (eigval.imag > circle_center[..., i].imag - window_imag / 2)
        )
        eigval = eigval[mask]
        eigvec = eigvec[:, mask]
        eigvals_merged.append(eigval)
        eigvecs_merged.append(eigvec)
    return GridResult(
        eigval=np.concatenate(eigvals_merged, axis=-1)
        if eigvals_merged
        else np.empty(0),
        eigvec=np.concatenate(eigvecs_merged, axis=-1)
        if eigvecs_merged
        else np.empty((n, 0)),
        circle_center=circle_center,
        circle_radius=circle_radius,
        n=n,
        result=res,
    )
