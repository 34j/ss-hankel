# SS Hankel

<p align="center">
  <a href="https://github.com/34j/ss-hankel/actions/workflows/ci.yml?query=branch%3Amain">
    <img src="https://img.shields.io/github/actions/workflow/status/34j/ss-hankel/ci.yml?branch=main&label=CI&logo=github&style=flat-square" alt="CI Status" >
  </a>
  <a href="https://ss-hankel.readthedocs.io">
    <img src="https://img.shields.io/readthedocs/ss-hankel.svg?logo=read-the-docs&logoColor=fff&style=flat-square" alt="Documentation Status">
  </a>
  <a href="https://codecov.io/gh/34j/ss-hankel">
    <img src="https://img.shields.io/codecov/c/github/34j/ss-hankel.svg?logo=codecov&logoColor=fff&style=flat-square" alt="Test coverage percentage">
  </a>
</p>
<p align="center">
  <a href="https://github.com/astral-sh/uv">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json" alt="uv">
  </a>
  <a href="https://github.com/astral-sh/ruff">
    <img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Ruff">
  </a>
  <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white&style=flat-square" alt="pre-commit">
  </a>
</p>
<p align="center">
  <a href="https://pypi.org/project/ss-hankel/">
    <img src="https://img.shields.io/pypi/v/ss-hankel.svg?logo=python&logoColor=fff&style=flat-square" alt="PyPI Version">
  </a>
  <img src="https://img.shields.io/pypi/pyversions/ss-hankel.svg?style=flat-square&logo=python&amp;logoColor=fff" alt="Supported Python versions">
  <img src="https://img.shields.io/pypi/l/ss-hankel.svg?style=flat-square" alt="License">
</p>

---

**Documentation**: <a href="https://ss-hankel.readthedocs.io" target="_blank">https://ss-hankel.readthedocs.io </a>

**Source Code**: <a href="https://github.com/34j/ss-hankel" target="_blank">https://github.com/34j/ss-hankel </a>

---

**Derivative-free** method to find zeros of analytic (holomorphic) functions / solve nonlinear (polynomial / generalized) eigenvalue problems using contour integration. (Block SS-Hankel method, Block Sakurai Sugiura method)

## Installation

Install this via pip (or your favourite package manager):

```shell
pip install ss-hankel
```

## Usage

Below is a simple example of solving the following nonlinear eigenvalue problem from [NEP-PACK Tutorial](https://nep-pack.github.io/NonlinearEigenproblems.jl/dev/tutorial_python_call/#Tutorial:-Using-NEP-PACK-from-python).

$$
f(x) = \begin{pmatrix}
3+e^{0.5x} & 2+2x+e^{0.5x} \\
3+e^{0.5x} & -1+x+e^{0.5x}
\end{pmatrix}
, \quad f(x) v = 0
$$

```python
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ss_hankel import score, ss_h_circle


def f(x: NDArray[Any]) -> NDArray[Any]: # deriative is not needed!
    return np.stack(
        [
            np.stack([3 + np.exp(0.5 * x), 2 + 2 * x + np.exp(0.5 * x)], axis=-1),
            np.stack([3 + np.exp(0.5 * x), -1 + x + np.exp(0.5 * x)], axis=-1),
        ],
        axis=-2,
    )


eig = ss_h_circle(
    f,
    num_vectors=2,
    max_order=8,
    circle_n_points=128, # number of integration points
    circle_radius=3, # radius of the contour
    circle_center=0, # center of the contour
)
print(f"Eigenvalues: {eig.eigval}") # eigenvalues inside the contour
print(f"Eigenvectors: {eig.eigvec}") # corresponding eigenvectors
print(f"|f(λ)v|/|f(λ)||v|: {score(f(eig.eigval), eig.eigvec)}")
```

```text
Eigenvalues: [-3.+3.19507247e-15j]
Eigenvectors: [[-0.45042759-0.61296714j]
 [-0.38438888-0.52309795j]]
|f(λ)v|/|f(λ)||v|: [1.37836544e-15]
```

- **Batch calculation** (function and/or contour) is supported.
  - Steps until SVD are batched. Function evaluations are batched (called only once).
  - Only the final step (solving small generalized eigenvalue problem) is not batched because the size of the eigenvalue problem (the number of eigenvalues in the contour) might be different and moreover `scipy.linalg.eig` does not support batch calculation.
- Since random matrices `U,V` are used in the algorithm, the results may vary slightly on each run. `np.random.Generator` can be passed to control the randomness.
- To get **zeros of an analytic function**, set `lambda x: f(x)[..., None, None]` as an argument. The SS-Hankel method for 1x1 matrix is completely equivalent to the Kravanja (1999)'s derivative-free root-finding method.
- The default parameters are set to be impractically small. Consider increasing `circle_n_points` and `max_order` based on the problem and `num_vectors` based on the matrix size.
- The number of eigenvalues (zeros) inside the contour is estimated by evaluating the numerical rank of the Hankel matrix. By default the singular values below the largest gap between singular values are considered meaningless, as propsed in Xiao (2016), but the behaviour can be controlled by manually setting `rtol`. `atol` (default: `1e-6`) is useful in the case when no eigenvalues are inside the contour.

## CLI Usage

```shell
> ss-hankel "{{3+Exp[x/2],2+2x+Exp[x/2]},{3+Exp[x/2],-1+x+Exp[x/2]}}" --circle-radius 4
eigenvalues:
[-3.-2.29788612e-15j]
eigenvectors (columns):
[[0.35283836-0.67388339j]
 [0.30110753-0.57508306j]]
|F(λ)v|/|F(λ)||v|:
[9.82824873e-16]
singular_values:
[1.36659229e-01 5.51578001e-17 3.11252713e-17 2.25070948e-17
 1.05446714e-17 9.42202841e-18 6.28427578e-18 2.84988862e-18]
```

## References

- [Asakura, J., Sakurai, T., Tadano, H., Ikegami, T., & Kimura, K. (2009). A numerical method for nonlinear eigenvalue problems using contour integrals. JSIAM Letters, 1, 52–55.](https://doi.org/10.1007/s006070050051)
- [Kravanja, P., & Van Barel, M. (1999). A Derivative-Free Algorithm for Computing Zeros of Analytic Functions. Computing (Vienna/New York), 63, 69–91.](https://doi.org/10.14495/jsiaml.1.52)
- [Xiao, J., Meng, S., Zhang, C., & Zheng, C. (2016). Resolvent sampling based Rayleigh-Ritz method for large-scale nonlinear eigenvalue problems. Computer Methods in Applied Mechanics and Engineering, 310, 33–57.](https://doi.org/10.1016/j.cma.2016.06.018)

## Alternatives

- [nep-pack/NonlinearEigenproblems.jl: Nonlinear eigenvalue problems in Julia: Iterative methods and benchmarks](https://github.com/nep-pack/NonlinearEigenproblems.jl)

### Zeros of analytic functions

- [rparini/cxroots: Find all the roots (zeros) of a complex analytic function within a given contour in the complex plane.](https://github.com/rparini/cxroots)
- [nennigb/polze: A python package to locate poles and zeros of a meromorphic function with their multiplicities](https://github.com/nennigb/polze)
- [Spectral-Analysis-UPB/PyZEAL: Project dealing with the numerical calculation of zeros, poles and residues of holomorphic and meromorphic functions. It aspires to be a PYthon ZEAL (ZEros of AnaLytic functions, a Fortran90 package) successor.](https://github.com/Spectral-Analysis-UPB/PyZEAL)

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- prettier-ignore-start -->
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- markdownlint-disable -->
<!-- markdownlint-enable -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- prettier-ignore-end -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

## Credits

[![Copier](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-grayscale-inverted-border-orange.json)](https://github.com/copier-org/copier)

This package was created with
[Copier](https://copier.readthedocs.io/) and the
[browniebroke/pypackage-template](https://github.com/browniebroke/pypackage-template)
project template.
