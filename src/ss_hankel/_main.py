import warnings
from collections.abc import Callable
from typing import Any, Literal, TypedDict

import attrs
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eig
from typing_extensions import NotRequired


def _get_random_matrix(size: tuple[int, ...], rng: np.random.Generator) -> np.ndarray:
    """
    Get a column fullrank matrix.

    Parameters
    ----------
    size : tuple[int, ...]
        The shape of the matrix.
    rng : np.random.Generator
        The random number generator.

    Returns
    -------
    np.ndarray
        A random regular matrix.

    """
    while True:
        np.random.random()
        m = rng.random(size) + 1j * rng.random(size)
        s = np.linalg.svd(m, compute_uv=False)
        if s[-1] / s[0] > 1e-4:
            return m


def broadcast_without_repeating(*arrays: NDArray[Any]) -> tuple[NDArray[Any], ...]:
    """
    Broadcast arrays without repeating arrays.

    Returns
    -------
    tuple[NDArray[Any], ...]
        The broadcasted arrays.

    """
    np.broadcast_shapes(*[a.shape for a in arrays])
    max_dim = max(a.ndim for a in arrays)
    return tuple(array[(None,) * (max_dim - array.ndim) + (...,)] for array in arrays)


def hankel_matrix(a: NDArray[Any], /) -> NDArray[Any]:
    """
    Create a Hankel matrix from a vector.

    Parameters
    ----------
    a : NDArray[Any]
        The input vector of shape [..., n]

    Returns
    -------
    NDArray[Any]
        The Hankel matrix of shape [..., n, n]

    """
    i = np.arange(a.shape[-1])
    return a[..., i[:, None] + i[None, :]]


@attrs.frozen(kw_only=True)
class SSHCircleResult:
    eigval: NDArray[Any]
    """The eigenvalues, an array of shape [...] of
        (array of shape [neigval] if max_neigval >= neigval, else None)"""
    eigvec: NDArray[Any]
    """The eigenvectors, an array of shape [...] of
        (array of shape [n, neigval] if max_neigval >= neigval, else None)"""
    s: NDArray[Any]
    """Array of shape [..., K * L]"""
    s_valid: NDArray[np.bool]
    """Boolearn array of shape [..., K * L]"""

    # input parameters
    n: int
    """The size of the matrix."""
    circle_center: NDArray[Any]
    """The center of the circle of shape [...]."""
    circle_radius: NDArray[Any]
    """The radius of the circle of shape [...]."""

    def __iter__(self) -> Any:
        yield self.eigval
        yield self.eigvec


class SSHKwargs(TypedDict):
    num_vectors: NotRequired[int | None]
    """Number of linearly independent vectors (L), by default None."""
    max_order: NotRequired[int | None]
    """Maximum order of the moments μ_k and s_k, by default None."""
    circle_n_points: NotRequired[int]
    """Number of integration points on the circle, by default 16"""
    circle_center: NotRequired[ArrayLike]
    """The center of the circle of shape [...], by default 0"""
    circle_radius: NotRequired[ArrayLike]
    """The radius of the circle of shape [...], by default 1"""
    rtol: NotRequired[float | Literal["auto"]]
    """The relative threshold to treat eigenvalues as zero, by default "auto"
    If "auto", the threshold is determined by
    searching the largest gap of singular values"""
    atol: NotRequired[float]
    """The absolute threshold to treat eigenvalues as zero, by default 1e-6"""
    max_neigval: NotRequired[int | None]
    """The maximum number of eigenvalues to proceed calculation, by default None"""
    rng: NotRequired[np.random.Generator | None]
    """The random number generator, by default None"""


class MaxOrderTooSmallWarning(RuntimeWarning):
    pass


class NEigvalExceedMaxWarning(RuntimeWarning):
    pass


class EigvalsOutsidePathWarning(RuntimeWarning):
    pass


def ss_h_circle(
    f: Callable[["NDArray[Any]"], "NDArray[Any]"],
    /,
    *,
    num_vectors: int | None = None,
    max_order: int | None = None,
    circle_n_points: int = 16,
    circle_center: ArrayLike = 0,
    circle_radius: ArrayLike = 1,
    rtol: float | Literal["auto"] = "auto",
    atol: float = 1e-6,
    max_neigval: int | None = None,
    rng: np.random.Generator | None = None,
) -> SSHCircleResult:
    """
    Sakurai-Sugiura method for the circle.

    Parameters
    ----------
    f : Callable[[np.ndarray], np.ndarray]
        An analytic function (F(z)).
        Array of shape [circle_n_points] will be passed
        and should return [circle_n_points, ..., n, n] array.
    num_vectors : int, optional
        Number of linearly independent vectors (L), by default None.
    max_order : int, optional
        Maximum order of the moments μ_k and s_k, by default None.
        The size of hankel matrix is num_vectors * max_order.
    circle_n_points : int, optional
        Number of integration points on the circle, by default 16
    circle_center : complex, optional
        The center of the circle of shape [...], by default 0
    circle_radius : float, optional
        The radius of the circle of shape [...], by default 1
    rtol : float, optional
        The relative threshold to treat eigenvalues as zero, by default "auto"
        If "auto", the threshold is determined by
        searching the largest gap of singular values
    atol : float, optional
        The absolute threshold to treat eigenvalues as zero, by default 1e-6
    max_neigval : int | None, optional
        The maximum number of eigenvalues to proceed calculation, by default None
    rng : np.random.Generator | None, optional
        The random number generator, by default None

    Returns
    -------
    SakuraiSugiuraCircleResult
        The eigenvalues and eigenvectors.

    Warnings
    --------
    MaxOrderTooSmallWarning
        The maximum order is too small against the number of eigenvalues.
    EigvalsOutsidePathWarning
        Some eigenvalues are outside the path.
    NEigvalExceedMaxWarning
        The number of eigenvalues is larger than ``max_neigval``.

    References
    ----------
    Asakura, J., Sakurai, T., Tadano, H., Ikegami, T., & Kimura, K. (2009).
    A numerical method for nonlinear eigenvalue problems using contour integrals.
    JSIAM Letters, 1, 52–55.
    https://doi.org/10.14495/jsiaml.1.52

    Kravanja, P., & Van Barel, M. (1999).
    A Derivative-Free Algorithm for Computing Zeros of Analytic Functions.
    Computing (Vienna/New York), 63, 69–91.
    https://doi.org/10.1007/s006070050051

    Xiao, J., Meng, S., Zhang, C., & Zheng, C. (2016).
    Resolvent sampling based Rayleigh-Ritz method
    for large-scale nonlinear eigenvalue problems.
    Computer Methods in Applied Mechanics and Engineering, 310, 33–57.
    https://doi.org/10.1016/j.cma.2016.06.018

    """
    num_vectors = num_vectors or 1
    if num_vectors < 1:
        raise ValueError(f"{num_vectors=} should be greater than 0")
    max_order = max_order or 2
    if max_order < 1:
        raise ValueError(f"{max_order=} should be greater than 0")
    if rtol == "auto" and max_order * num_vectors <= 1:
        raise ValueError(
            f"To set rtol to 'auto', {max_order=} * {num_vectors=} "
            "should be greater than 1"
        )
    rng_ = np.random.default_rng() if rng is None else rng
    circle_center = np.asarray(circle_center)
    circle_radius = np.asarray(circle_radius)
    if not np.all(np.isfinite(circle_center)):
        raise ValueError(f"circle_center should be finite, but got {circle_center}")
    if not np.all(np.isfinite(circle_radius)):
        raise ValueError(f"circle_radius should be finite, but got {circle_radius}")
    circle_center, circle_radius = broadcast_without_repeating(
        circle_center, circle_radius
    )
    additional_dims = circle_center.ndim

    # use scaled and shifted weights (8)
    # [circle_n_points, ...]
    w = np.exp(2 * np.pi * 1j * np.arange(circle_n_points) / circle_n_points)[
        (...,) + (None,) * additional_dims
    ]
    x = circle_center + circle_radius * w
    # [circle_n_points, ..., n, n]
    fvals = f(x)

    # check fvals.shape
    if fvals.shape[0] != circle_n_points:
        raise ValueError(
            f"f should return array of shape [{circle_n_points=}, ..., n, n], "
            f"but got {fvals.shape}"
        )

    # check all additional dimensions are broadcastable
    if fvals.shape[-1] != fvals.shape[-2]:
        raise ValueError(
            "f should return array which last two dimensions are the same "
            f"(batched matrix), but got {fvals.shape}"
        )
    try:
        additional_shape = np.broadcast_shapes(
            circle_center.shape, circle_radius.shape, fvals.shape[1:-2]
        )
        if not circle_center.ndim == circle_radius.ndim == fvals.ndim - 3:
            raise ValueError
    except ValueError as e:
        raise ValueError(
            f"{circle_center.shape=}, {circle_radius.shape=}, {fvals.shape[1:-2]=} "
            "should be the same length and broadcastable"
        ) from e
    circle_center = np.broadcast_to(circle_center, additional_shape)
    circle_radius = np.broadcast_to(circle_radius, additional_shape)

    # check num_vectors
    n: int = fvals.shape[-1]
    if num_vectors > n:
        raise ValueError(
            f"num_vectors {num_vectors} should be less than or equal to"
            f" the matrix width and height {fvals.shape[-2:]=}"
        )

    # check fvals are finite
    if not np.all(np.isfinite(fvals)):
        raise ValueError(
            "f should return finite values, but got "
            f"{(~np.isfinite(fvals)).sum()} non-finite values "
            f"out of {fvals.size}"
        )

    # [n, L]
    V = _get_random_matrix((n, num_vectors), rng_)
    U = _get_random_matrix((n, num_vectors), rng_)
    # [circle_n_points, ..., n, L]
    # from [circle_n_points, ..., n, n] @ [circle_n_points, n, L]
    # -> [circle_n_points, ..., n, L]
    fvalsinvV = np.linalg.solve(fvals, V[(None,) * (additional_dims + 1) + (...,)])
    # fvalsinvV = np.linalg.inv(fvals) @ V[(None,) * (additional_dims + 1)]
    k = np.arange(2 * max_order)
    # [k, circle_n_points, ..., n, L] -> [k, ..., n, L]
    S = np.mean(
        (w[None, ...] ** (k[(slice(None),) + (None,) * (additional_dims + 1)] + 1))[
            (...,) + (None,) * (2)
        ]
        * fvalsinvV[None, ...],
        axis=1,
    )
    # [L, n] @ [k, ..., n, L] -> [k, ..., L, L]
    M = U.conj().T @ S
    k = np.arange(max_order)
    # [K, K, ..., L, L]
    Hs = M[k[:, None] + k[None, :] + 1, ...]
    H = M[k[:, None] + k[None, :], ...]
    # [..., K, L, K, L]
    Hs = np.moveaxis(Hs, (0, 1), (-4, -2))
    H = np.moveaxis(H, (0, 1), (-4, -2))
    # [..., K * L, K * L]
    Hs = Hs.reshape((*Hs.shape[:-4], Hs.shape[-3] * Hs.shape[-2], -1))
    H = H.reshape((*H.shape[:-4], H.shape[-3] * H.shape[-2], -1))
    # [..., K * L] (step 6)
    s = np.linalg.svd(H, compute_uv=False)
    # Omit small singular value components (step 7)
    # https://arxiv.org/pdf/1510.07522 p.11 (28)
    if rtol == "auto":
        s_valid_count = np.diff(np.log(s), axis=-1).argmin(axis=-1, keepdims=True) + 1
        s_valid = np.arange(s.shape[-1]) < s_valid_count
    else:
        s_valid = s > rtol * s[..., [0]]
    s_valid = s_valid & (s > atol)
    neigvals = np.sum(s_valid, axis=-1)
    eigvals = np.empty(neigvals.shape, dtype=object)
    eigvecs = np.empty(neigvals.shape, dtype=object)
    if (neigvals >= max_order).any():
        warnings.warn(
            f"Max order {max_order} is too small against"
            f" number of eigenvalues {neigvals}",
            MaxOrderTooSmallWarning,
            stacklevel=2,
        )
    it = np.nditer(neigvals, flags=["multi_index"])
    for neigval in it:
        neigval: int  # type: ignore
        if neigval == 0:
            eigval = np.empty((0,), dtype=H.dtype)
            eigvec = np.empty((n, 0), dtype=H.dtype)
        elif max_neigval is not None and neigval > max_neigval:
            warnings.warn(
                f"Number of eigenvalues {neigval} is larger than"
                f" the maximum number of eigenvalues {max_neigval}",
                NEigvalExceedMaxWarning,
                stacklevel=2,
            )
            eigval = None
            eigvec = None
        else:
            # slice hankel matrices
            # [..., neigval, neigval]
            Hscut = Hs[it.multi_index + (slice(neigval),) * 2]
            Hcut = H[it.multi_index + (slice(neigval),) * 2]

            # eigenvalues
            # [..., neigval], [..., neigval(vector), neigval(eigenvalue)]
            eigvalh, eigvech = eig(Hscut, b=Hcut)
            eigvalh_out = (
                np.abs(eigvalh - circle_center[it.multi_index])
                > circle_radius[it.multi_index] + atol
            )
            if np.any(eigvalh_out):
                warnings.warn(
                    "Some eigenvalues are outside the path.",
                    EigvalsOutsidePathWarning,
                    stacklevel=2,
                )
            eigval = (
                circle_radius[it.multi_index] * eigvalh + circle_center[it.multi_index]
            )

            # eigenvectors
            # [k, ..., n, L] -> [neigval, ..., n, L] -> [..., n, L * neigval]
            # -> [..., n, neigval]
            Scut = np.concatenate(S[:neigval], axis=-1)[it.multi_index][:, :neigval]
            eigvec = Scut @ eigvech
            eigvec /= np.linalg.vector_norm(eigvec, axis=-2, keepdims=True)

        # assign
        eigvals[it.multi_index] = eigval
        eigvecs[it.multi_index] = eigvec
    if neigvals.shape == ():
        eigvals = eigvals.item()
        eigvecs = eigvecs.item()
    return SSHCircleResult(
        eigval=eigvals,
        eigvec=eigvecs,
        s=s,
        s_valid=s_valid,
        n=n,
        circle_center=circle_center,
        circle_radius=circle_radius,
    )
