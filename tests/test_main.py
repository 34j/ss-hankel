import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_less

from ss_hankel import (
    EigvalsOutsidePathWarning,
    MaxOrderTooSmallWarning,
    SSHKwargs,
    score,
    ss_h_circle,
)
from ss_hankel.testing import asakura_example_1, asakura_example_1_eigvals


def test_main() -> None:
    kwargs = SSHKwargs(
        num_vectors=2,
        max_order=8,
        circle_n_points=256,
        circle_radius=[3.3, 0.01],
        circle_center=[0, 2.8],
        rng=np.random.default_rng(0),
    )
    eig = ss_h_circle(
        asakura_example_1,
        **kwargs,
    )
    eigval, eigvec = eig

    # global array shape
    assert eigval.shape == (2,)
    assert eigvec.shape == (2,)
    assert eig.s.shape == (2, 2 * 8)
    assert eig.s_valid.shape == (2, 2 * 8)

    # first circle
    assert eigval[0].shape == (6,)
    assert eigvec[0].shape == (3, 6)

    # check if eigval is correct
    assert_allclose(
        eigval[0][np.argsort(np.real(eigval[0]))],
        asakura_example_1_eigvals(),
        atol=1e-10,
    )
    # check if score is small
    assert_array_less(
        score(asakura_example_1(eigval[0]), eigvec[0]), 1e-8, strict=False
    )

    # second circle
    # no eigenvalues
    assert eigval[1].shape == (0,)
    assert eigvec[1].shape == (3, 0)


def test_tol_too_small():
    with (
        pytest.warns(EigvalsOutsidePathWarning),
        pytest.warns(MaxOrderTooSmallWarning),
        pytest.warns(RuntimeWarning, match="invalid value encountered in multiply"),
    ):
        eigval, _ = ss_h_circle(
            asakura_example_1,
            num_vectors=2,
            max_order=8,
            circle_n_points=256,
            circle_radius=3.3,
            circle_center=0,
            rng=np.random.default_rng(0),
            rtol=1e-20,
            atol=1e-20,
        )
    eigval = eigval[np.isfinite(eigval)]
    assert_allclose(
        eigval[np.argsort(np.real(eigval))],
        asakura_example_1_eigvals(),
        atol=1e-10,
    )
