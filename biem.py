import seaborn as sns
from array_api._2024_12 import Array
from array_api_compat import array_namespace
from array_api_compat import numpy as np
from axis_batch import AxisBatch
from biem_helmholtz_sphere import BIEMKwargs, biem
from matplotlib import pyplot as plt
from tqdm import tqdm
from ultrasphere import create_standard

from ss_hankel import SSHKwargs
from ss_hankel._main import ss_h_hankel_matrix

xp = np
c = create_standard(1)
sakurai_kwargs = SSHKwargs(
    num_vectors=12,
    max_order=12,
    circle_n_points=32768,
    circle_center=2 - 2j,
    circle_radius=2,
)
biem_kwargs = BIEMKwargs(
    n_end=20,
    centers=xp.asarray([[-2, 0], [2, 0]]),
    radii=xp.asarray([1, 1]),
    k=xp.asarray(xp.nan),
    eta=xp.asarray(0),
    kind="outer",
    alpha=xp.asarray(0),
    beta=xp.asarray(1),
    force_matrix=True,
)
for i in ["centers", "radii", "eta"]:
    if i in biem_kwargs:
        biem_kwargs[i] = biem_kwargs[i][None, ...]
xp = array_namespace(biem_kwargs["centers"], biem_kwargs["radii"])


def f(z: Array) -> Array:
    # [z, ...]
    z = xp.asarray(z)
    b = AxisBatch(z, axis=0, size=1000000)
    for zc in tqdm(b, disable=len(b) == 1, desc="Finding resonance"):
        biem_kwargs["k"] = zc
        b.send(
            biem(
                c,
                **biem_kwargs,
            ).matrix
        )
    val = b.value
    val = xp.reshape(
        val,
        (
            *val.shape[:-4],
            val.shape[-4] * val.shape[-3],
            val.shape[-2] * val.shape[-1],
        ),
    )
    if xp.any(xp.isnan(val)):
        raise ValueError("NaN in the matrix", val, biem_kwargs)
    return np.asarray(val)


H, Hs = ss_h_hankel_matrix(
    f,
    **sakurai_kwargs,
)
svdvals = []
for i in range(H.shape[0]):
    svdvals.append(np.linalg.svdvals(H[:i, :i]))
sns.set_theme()
fig, ax = plt.subplots()
data = [(i + 1, j) for i, sv in enumerate(svdvals) for j in sv]
ax.scatter(*zip(*data, strict=False), s=1, marker="x")
ax.set_yscale("log")
ax.set_xlabel("Matrix size")
ax.set_ylabel("Singular values")
ax.set_title("Singular values of the biem matrix vs matrix size")
fig.savefig("biem_singular_values.png")
