from typing import Any, Literal

import cyclopts
import numpy as np
import sympy as sp
from numpy.typing import NDArray
from rich import print
from sympy.parsing.mathematica import parse_mathematica

from ._main import ss_h_circle
from ._score import score

app = cyclopts.App(name="ss-hankel")


@app.default()
def main(
    expr: str,
    /,
    *,
    num_vectors: int | None = None,
    max_order: int | None = 8,
    circle_n_points: int = 256,
    circle_center: complex = 0,
    circle_radius: float = 1,
    rtol: float | Literal["auto"] = "auto",
    atol: float = 1e-6,
    seed: int | None = None,
) -> None:
    """
    Compute the eigenvalues and eigenvectors of a matrix-valued function.

    Some default parameters (circle_n_points, max_order) are set larger
    than the default in the Python API for convenience.

    Parameters
    ----------
    expr : str
        The matrix-valued function in Mathematica format.
    num_vectors : int, optional
        Number of linearly independent vectors (L), by default None.
    max_order : int, optional
        Maximum order of the moments μ_k and s_k, by default 8.
        The size of hankel matrix is num_vectors * max_order.
    circle_n_points : int, optional
        Number of integration points on the circle, by default 256
    circle_center : complex, optional
        The center of the circle, by default 0
    circle_radius : float, optional
        The radius of the circle, by default 1
    rtol : float, optional
        The relative threshold to treat eigenvalues as zero, by default "auto"
        If "auto", the threshold is determined by
        searching the largest gap of singular values
    atol : float, optional
        The absolute threshold to treat eigenvalues as zero, by default 1e-6
    seed : int, optional
        The seed for the random number generator, by default None.
        If None, the random number generator is not seeded.

    """
    mat = parse_mathematica(expr)
    mat_lambda = sp.lambdify("x", mat, modules=["numpy"])

    def f(x: NDArray[Any]) -> NDArray[Any]:
        return np.stack(
            [
                np.stack([np.broadcast_to(item, x.shape) for item in row], axis=-1)
                for row in mat_lambda(x)
            ],
            axis=-2,
        )

    eig = ss_h_circle(
        f,
        num_vectors=num_vectors,
        max_order=max_order,
        circle_n_points=circle_n_points,
        circle_radius=circle_radius,
        circle_center=circle_center,
        rtol=rtol,
        atol=atol,
        rng=np.random.default_rng(seed),
    )
    result = {
        "eigenvalues": eig.eigval,
        "eigenvectors (columns)": eig.eigvec,
        "|F(λ)v|/|F(λ)||v|": score(f(eig.eigval), eig.eigvec),
        "singular_values": eig.s,
    }
    for k, v in result.items():
        print(f"{k}:")
        print(v)
