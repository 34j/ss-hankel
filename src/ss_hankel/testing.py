from typing import Any

import numpy as np
from numpy.typing import NDArray


def asakura_example_1(z: NDArray[Any]) -> NDArray[Any]:
    """
    Example 1 from Asakura 2009.

    Parameters
    ----------
    z : NDArray[Any]
        The input array of shape (...).

    Returns
    -------
    NDArray[Any]
        Matrix of shape (..., 3, 3).

    """
    return np.stack(
        [
            np.stack(
                [
                    2 * np.exp(z) + np.cos(z) - 14,
                    (z**2 - 1) * np.sin(z) + (2 * np.exp(z) + 14) * np.cos(z),
                    2 * np.exp(z) - 14,
                ],
                axis=-1,
            ),
            np.stack(
                [
                    (z + 3) * (np.exp(z) - 7),
                    np.sin(z) + (z + 3) * (np.exp(z) - 7) * np.cos(z),
                    (z + 3) * (np.exp(z) - 7),
                ],
                axis=-1,
            ),
            np.stack(
                [np.exp(z) - 7, (np.exp(z) - 7) * np.cos(z), np.exp(z) - 7], axis=-1
            ),
        ],
        axis=-2,
    )


def asakura_example_1_eigvals() -> NDArray[Any]:
    """
    The eigenvalues of example 1 from Asakura 2009.

    Returns
    -------
    NDArray[Any]
        Eigenvalues of shape (6,).

    """
    return np.asarray(
        [
            -np.pi,
            -np.pi / 2,
            0,
            np.pi / 2,
            np.log(7),
            np.pi,
        ],
    )
