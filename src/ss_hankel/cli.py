import numpy as np
import typer
from rich import print

from ._main import ss_h_circle
from ._score import score
from .testing import asakura_example_1, asakura_example_1_eigvals

app = typer.Typer()


@app.command()
def main() -> None:
    """Compute the eigenvalues and eigenvectors of the Example 1 from Asakura 2009."""
    eigval, eigvec = ss_h_circle(
        asakura_example_1,
        num_vectors=2,
        max_order=8,
        circle_n_points=256,
        circle_radius=3.3,
        circle_center=0,
        rng=np.random.default_rng(0),
    )
    print(f"diff = {np.abs(eigval - asakura_example_1_eigvals())}")
    print(f"|F(λ)v|/|F(λ)||v| = {score(asakura_example_1(eigval), eigvec)}")
