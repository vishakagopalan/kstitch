import logging
import os

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy
import shapely
import time
from itertools import repeat
from multiprocess import Pool

# Adds a wall-clock timestamp to every log message automatically, so no
# per-call changes are needed inside frechet() or kendall_tpca().
# Guarded so re-importing this module (e.g. via reticulate::import_from_path
# across multiple R sessions, or repeated devtools::load_all()) never
# attaches duplicate handlers, which would otherwise cause each log line to
# print multiple times.
logger = logging.getLogger("kendall_tpca")
if not logger.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logger.addHandler(_handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False


# ── Shape geometry helpers ─────────────────────────────────────────────────────

def cell_features(x):
    """Compute scalar shape descriptors for a single cell contour.

    Parameters
    ----------
    x : np.ndarray, shape (2, L)
        Raw (x, y) landmark coordinates.

    Returns
    -------
    list
        [scale, major_axis, minor_axis, eccentricity, area_in_px,
         perimeter, roundness, convexity_score]
    """
    x_copy = x.copy()
    mu = x_copy.mean(axis=1)
    for i in range(x_copy.shape[1]):
        x_copy[:, i] -= mu
    scale = np.linalg.norm(x_copy, ord='fro')
    _, sigma, _ = np.linalg.svd(x_copy)
    major_axis, minor_axis = sigma[0], sigma[1]
    eccentricity = np.sqrt(1 - (minor_axis ** 2 / major_axis ** 2))

    poly = shapely.Polygon(x.T)
    area_in_px = shapely.area(poly)
    perimeter = shapely.length(poly)
    hull = shapely.convex_hull(shapely.MultiPoint(x_copy.T))
    convexity_score = hull.length / perimeter
    roundness = perimeter ** 2 / (4 * np.pi * area_in_px)

    return [scale, major_axis, minor_axis, eccentricity,
            area_in_px, perimeter, roundness, convexity_score]


def preprocess(x):
    """Centre and scale a shape to the pre-shape sphere.

    Parameters
    ----------
    x : np.ndarray, shape (2, L)

    Returns
    -------
    np.ndarray, shape (2, L)
    """
    mu = x.mean(axis=1)
    for i in range(x.shape[1]):
        x[:, i] -= mu
    return x / np.linalg.norm(x, ord='fro')


def interp_shape_and_preprocess(x_original, num_vertices):
    """Interpolate a contour to a fixed number of landmarks and preprocess.

    Parameters
    ----------
    x_original : np.ndarray, shape (2, L)
    num_vertices : int

    Returns
    -------
    np.ndarray, shape (2, num_vertices)
    """
    x = preprocess(x_original)
    x_u, unique_mask = np.unique(x, axis=1, return_index=True)
    keep = sorted(set(np.arange(x.shape[1])) & set(unique_mask))
    x_u = x[:, keep]

    tck, _ = scipy.interpolate.splprep(
        [x_u[0, :], x_u[1, :]], s=0, k=1, per=True
    )
    xi, yi = scipy.interpolate.splev(np.linspace(0, 1, num_vertices), tck)
    return preprocess(np.array([xi, yi]))


# ── Riemannian geometry ────────────────────────────────────────────────────────

def OPA(A, B):
    """Ordinary Procrustes alignment: rotationally align A to B.

    Parameters
    ----------
    A, B : np.ndarray, shape (2, L)

    Returns
    -------
    np.ndarray, shape (2, L)
        Rotated A.
    """
    u, _, v_t = np.linalg.svd(B @ A.T)
    return (u @ v_t) @ A


def reparam_OPA(A, B):
    """Reparameterisation + OPA: find the cyclic shift of A closest to B.

    Parameters
    ----------
    A, B : np.ndarray, shape (2, L)

    Returns
    -------
    np.ndarray, shape (2, L)
        Best-aligned version of A.
    """
    A_best = A
    best_error = np.linalg.norm(A - B)
    for shift in range(A.shape[1]):
        gamma = np.roll(np.arange(A.shape[1]), shift)
        A_rot = OPA(A[:, gamma], B)
        error = np.linalg.norm(A_rot - B)
        if error < best_error:
            best_error = error
            A_best = A_rot
    return A_best


def exp(p, v, theta=None):
    """Exponential map on the pre-shape sphere.

    Parameters
    ----------
    p : np.ndarray, shape (2, L)
        Base point (mean shape).
    v : np.ndarray, shape (2, L)
        Tangent vector at p.
    theta : float, optional
        Geodesic step size. Defaults to ``||v||``.

    Returns
    -------
    np.ndarray, shape (2, L)
    """
    if theta is None:
        theta = np.linalg.norm(v)
    p2 = np.cos(theta) * p + np.sin(theta) * v / np.linalg.norm(v)
    return p2 / np.linalg.norm(p2, ord='fro')


def exp_safe(p, v, theta=None, eps=1e-12):
    """Exponential map with a zero-norm guard.

    Parameters
    ----------
    p : np.ndarray, shape (2, L)
    v : np.ndarray, shape (2, L)
    theta : float, optional
    eps : float
        Minimum norm below which v is treated as zero.

    Returns
    -------
    np.ndarray, shape (2, L)
    """
    v_norm = np.linalg.norm(v)
    if v_norm < eps:
        return p.copy()
    if theta is None:
        theta = v_norm
    p2 = np.cos(theta) * p + np.sin(theta) * v / v_norm
    return p2 / np.linalg.norm(p2, ord='fro')


def log(p1, p2, reparam=True):
    """Logarithmic map: lift p2 into the tangent space at p1.

    Parameters
    ----------
    p1 : np.ndarray, shape (2, L)
        Base point.
    p2 : np.ndarray, shape (2, L)
        Target point on the manifold.
    reparam : bool
        If True, apply reparameterisation + OPA before mapping; otherwise
        OPA only.

    Returns
    -------
    np.ndarray, shape (2, L)
        Tangent vector at p1 pointing towards p2.
    """
    s = p2.shape
    p2 = reparam_OPA(p2, p1) if reparam else OPA(p2, p1)
    p2_v, p1_v = p2.flatten(), p1.flatten()
    dot = np.clip(np.dot(p1_v, p2_v), -1, 1)
    theta = np.arccos(dot)
    frac = 1 if theta == 0 else theta / np.sin(theta)
    return ((p2_v - np.dot(p2_v, p1_v) * p1_v) * frac).reshape(s)


def theta_value(p1, p2):
    """Geodesic distance between two pre-shapes.

    Parameters
    ----------
    p1, p2 : np.ndarray, shape (2, L)

    Returns
    -------
    float
    """
    p2 = reparam_OPA(p2, p1)
    dot = np.clip(np.dot(p1.flatten(), p2.flatten()), -1, 1)
    return np.arccos(dot)


# ── Parallel helper ────────────────────────────────────────────────────────────

def log_map_values(shapes, ref_shape, query_idxes):
    """Compute log maps from ref_shape to a subset of shapes (for parallel use).

    Parameters
    ----------
    shapes : np.ndarray, shape (N, 2, L)
    ref_shape : np.ndarray, shape (2, L)
    query_idxes : array-like of int

    Returns
    -------
    np.ndarray, shape (len(query_idxes), 2, L)
    """
    result = np.zeros((len(query_idxes), shapes.shape[1], shapes.shape[2]))
    for ctr, idx in enumerate(query_idxes):
        result[ctr] = log(ref_shape, shapes[idx])
    return result


# ── Frechet mean ───────────────────────────────────────────────────────────────

def frechet(X, eta=0.001, use_parallel=False, num_chunks=4, max_iter=100,
            tol=1e-4, store_history=False):
    """Compute the Frechet mean of a set of pre-shapes.

    Parameters
    ----------
    X : np.ndarray, shape (N, 2, L)
    eta : float
        Gradient step size.
    use_parallel : bool
    num_chunks : int
        Number of parallel workers.
    max_iter : int
    tol : float
        Convergence tolerance on average tangent vector. 
    store_history : bool
        If True, also return the mean history and convergence trace.

    Returns
    -------
    np.ndarray or list
        Frechet mean (shape (2, L)), or [mean, mean_history, history] when
        ``store_history=True``.
    """
    mu = X[0]
    history = []
    if store_history:
        mu_history = np.zeros((max_iter, X.shape[1], X.shape[2]))

    start_time = time.time()
    logger.info(
        "Frechet mean: starting on %d shapes (eta=%.4g, tol=%.4g, max_iter=%d, parallel=%s)",
        X.shape[0], eta, tol, max_iter, use_parallel
    )

    for i in range(max_iter):
        if use_parallel:
            chunk_size = int(X.shape[0] / num_chunks)
            idxes = np.arange(X.shape[0])
            chunks = [
                idxes[c * chunk_size:(c + 1) * chunk_size]
                for c in range((len(idxes) + chunk_size - 1) // chunk_size)
            ]
            with Pool(num_chunks) as pool:
                dmu_list = pool.starmap(
                    log_map_values, zip(repeat(X), repeat(mu), chunks)
                )
            dmu = sum(
                v for batch in dmu_list for v in batch
            ) / X.shape[0]
        else:
            dmu = sum(log(mu, X[j]) for j in range(X.shape[0])) / X.shape[0]

        prev_mu = mu
        mu = exp(mu, dmu * eta)

        if store_history:
            mu_history[i] = mu

        grad_norm = np.linalg.norm(dmu)

        change = theta_value(mu, prev_mu)
        logger.debug("Frechet iteration %d: grad_norm = %.6f", i, grad_norm)

        history.append(change)

        if i % 10 == 0 and i > 0:
            elapsed = time.time() - start_time
            logger.info(
                "Frechet mean: iteration %d, grad_norm = %.6f (elapsed %.1fs)",
                i, grad_norm, elapsed
            )

        if grad_norm < tol and i > 0:
            elapsed = time.time() - start_time
            logger.info(
                "Frechet mean: converged at iteration %d (grad_norm = %.6f < tol = %.4g, elapsed %.1fs)",
                i, grad_norm, tol, elapsed
            )
            break
    else:
        elapsed = time.time() - start_time
        logger.info(
            "Frechet mean: reached max_iter (%d) without converging (final grad_norm = %.6f, elapsed %.1fs)",
            max_iter, grad_norm, elapsed
        )

    if store_history:
        return [mu, mu_history[:i], history]
    return mu



# ── Tangent Principal Component Analysis ────────────────────────────────────────────────
def kendall_tpca(X, mu=None, max_frechet_mean_iter=100, use_parallel=True,
        frechet_mean_tol=1e-4, num_chunks=4, eta=1, store_history=False):
    """Run Principal Geodesic Analysis on a set of pre-shapes.

    Parameters
    ----------
    X : np.ndarray, shape (N, 2, L)
    mu : np.ndarray, optional
        Pre-computed Frechet mean. Computed from X if None.
    max_frechet_mean_iter : int
    use_parallel : bool
    frechet_mean_tol : float
    num_chunks : int
    eta : float
    store_history : bool

    Returns
    -------
    tuple
        (p, lambdas, v, mu_or_mu_info, U_flat)
        p          — embeddings, shape (N, n_pcs)
        lambdas    — PC variances
        v          — PC directions, shape (2L, n_pcs)
        mu_or_info — Frechet mean, or history list when store_history=True
        U_flat     — flattened log-map matrix, shape (N, 2L)
    """
    start_time = time.time()
    logger.info("kendall_tpca: running on %d shapes, %d landmarks", X.shape[0], X.shape[2])

    if mu is None:
        logger.info("kendall_tpca: no mu provided, computing Frechet mean")
        mu_result = frechet(
            X, eta=eta, max_iter=max_frechet_mean_iter,
            num_chunks=num_chunks, use_parallel=use_parallel,
            tol=frechet_mean_tol, store_history=store_history
        )
        mu = mu_result[0] if store_history else mu_result
    else:
        logger.info("kendall_tpca: using provided mu")

    U_flat = np.zeros((X.shape[0], 2 * X.shape[2]))
    for i in range(X.shape[0]):
        U_i = log(mu, X[i])
        U_flat[i] = U_i.flatten()

    logger.info("kendall_tpca: log-map complete, running SVD on %s matrix", U_flat.shape)

    P, sigma, R_t = np.linalg.svd(U_flat, full_matrices=False)
    lambdas = sigma ** 2 / (X.shape[0] - 1)
    embedding = P @ np.diag(sigma)

    total_var = lambdas.sum()
    top5_frac = lambdas[:5].sum() / total_var if total_var > 0 else float("nan")
    elapsed = time.time() - start_time
    logger.info(
        "kendall_tpca: done. %d PCs, top 5 explain %.1f%% of variance (elapsed %.1fs)",
        len(lambdas), top5_frac * 100, elapsed
    )

    if store_history:
        return embedding, lambdas, R_t.T, mu_result, U_flat

    return embedding, lambdas, R_t.T, mu, U_flat

# ── Reconstruction helpers ─────────────────────────────────────────────────────

def reconstruct_shapes_from_pca(v, mu, lambdas, sds_to_plot, num_pcs=8):
    """Reconstruct contour shapes along each PC at specified SD values.

    Parameters
    ----------
    v : np.ndarray, shape (2L, n_pcs)
    mu : np.ndarray, shape (2, L)
    lambdas : np.ndarray
    sds_to_plot : array-like of float
    num_pcs : int

    Returns
    -------
    dict with keys x, y, PC, sd_val (all numpy arrays)
    """
    rows = []
    for j in range(num_pcs):
        vec = v[:, j]
        for scale in sds_to_plot:
            scale_to_store = scale
            if scale == 0:
                scale = 1e-6
            theta = scale * np.sqrt(lambdas[j])
            shape = exp(mu, vec.reshape(mu.shape), theta=theta)
            df = pd.DataFrame(shape.T, columns=["x", "y"])
            df["PC"] = np.array([f"PC{j + 1}"] * len(df), dtype=object)
            df["sd_val"] = float(scale_to_store)
            rows.append(df)

    out = pd.concat(rows, axis=0, ignore_index=True)
    out["PC"] = out["PC"].to_numpy(dtype=object)
    return {
        "x":      out["x"].to_numpy(dtype=float),
        "y":      out["y"].to_numpy(dtype=float),
        "PC":     out["PC"].to_numpy(dtype=object),
        "sd_val": out["sd_val"].to_numpy(dtype=float),
    }


def reconstruct_shape_from_pca_coords(pca_coords, v, mu, lambdas=None,
                                      coords_are_sds=False, return_df=True):
    """Reconstruct a single shape from a PCA coordinate vector.

    Parameters
    ----------
    pca_coords : array-like, shape (d,)
    v : np.ndarray, shape (2L, n_pcs)
    mu : np.ndarray, shape (2, L)
    lambdas : np.ndarray, optional
        Required when ``coords_are_sds=True``.
    coords_are_sds : bool
    return_df : bool

    Returns
    -------
    pd.DataFrame or np.ndarray
    """
    pca_coords = np.asarray(pca_coords, dtype=float)
    num_pcs = len(pca_coords)
    if coords_are_sds:
        if lambdas is None:
            raise ValueError("lambdas must be provided when coords_are_sds=True.")
        pca_coords = pca_coords * np.sqrt(lambdas[:num_pcs])

    tangent_vec = (v[:, :num_pcs] @ pca_coords).reshape(mu.shape)
    shape = exp_safe(mu, tangent_vec)

    if return_df:
        df = pd.DataFrame(shape.T, columns=["x", "y"])
        df["coord_source"] = "pca_coord_vector"
        return df
    return shape


def project_onto_dataset(X, mu, v):
    """Project new shapes onto existing TPCA directions.

    Parameters
    ----------
    X : np.ndarray, shape (N, 2, L)
    mu : np.ndarray, shape (2, L)
    v : np.ndarray, shape (2L, n_pcs)

    Returns
    -------
    np.ndarray, shape (N, n_pcs)
    """
    U_flat = np.zeros((X.shape[0], 2 * X.shape[2]))
    for i in range(X.shape[0]):
        U_flat[i] = log(mu, X[i]).flatten()
    return U_flat @ v


# ── Pipeline entry points ──────────────────────────────────────────────────────

def compute_pre_shape_embedding(
        pre_shape_output_dir,
        boundary_parquet_path=None,
        df=None,
        num_vertices_to_sample=30,
        cell_ids_to_analyze=None,
        x_vertex_col="vertex_x",
        y_vertex_col="vertex_y",
        cell_id_col="cell_id"):
    """Compute pre-shape embeddings from cell contours and write to disk.

    Accepts either a parquet file path or an in-memory pandas DataFrame.
    Writes ``Pre_Shape_Space_Embedding.h5`` and ``Shape_Metadata.csv.gz``
    to ``pre_shape_output_dir``.

    Parameters
    ----------
    pre_shape_output_dir : str
    boundary_parquet_path : str, optional
    df : pd.DataFrame, optional
    num_vertices_to_sample : int
    cell_ids_to_analyze : list, optional
    x_vertex_col, y_vertex_col, cell_id_col : str
    """
    if df is not None:
        df_shapes = pa.Table.from_pandas(df)
    elif boundary_parquet_path is not None:
        df_shapes = pq.read_table(boundary_parquet_path)
    else:
        raise ValueError("One of boundary_parquet_path or df must be provided.")

    os.makedirs(pre_shape_output_dir, exist_ok=True)

    groups = df_shapes.group_by([cell_id_col]).aggregate(
        [(x_vertex_col, "list"), (y_vertex_col, "list")]
    )
    distinct_ids = groups.column(cell_id_col).to_pylist()
    if cell_ids_to_analyze is None:
        cell_ids_to_analyze = distinct_ids
    num_groups = len(cell_ids_to_analyze)

    pre_shape_mat = np.zeros((num_groups, 2, num_vertices_to_sample))
    cell_ids, idx_vals = [], []
    areas         = np.zeros(num_groups)
    roundness     = np.zeros(num_groups)
    major_axes    = np.zeros(num_groups)
    minor_axes    = np.zeros(num_groups)
    perimeters    = np.zeros(num_groups)
    eccentricity  = np.zeros(num_groups)
    convexity     = np.zeros(num_groups)

    idx = 0
    for group in groups.to_struct_array():
        cell_id = str(group[cell_id_col])
        if cell_id not in cell_ids_to_analyze:
            continue
        if idx % 1000 == 0 and idx > 0:
            logger.info("Pre-shape embedding: %d / %d", idx, num_groups)

        x_orig = np.array([
            group[f"{x_vertex_col}_list"].values.to_numpy(),
            group[f"{y_vertex_col}_list"].values.to_numpy(),
        ])
        (_, major_axes[idx], minor_axes[idx], eccentricity[idx],
         areas[idx], perimeters[idx], roundness[idx],
         convexity[idx]) = cell_features(x_orig)

        pre_shape_mat[idx] = interp_shape_and_preprocess(
            x_orig.copy(), num_vertices_to_sample
        )
        cell_ids.append(cell_id)
        idx_vals.append(idx)
        idx += 1

    meta = pd.DataFrame({
        "cell_id":    cell_ids,
        "numpy_idx":  np.array(idx_vals),
        "R_idx":      np.array(idx_vals) + 1,
        "roundness":  roundness,
        "area":       areas,
        "major_axis": major_axes,
        "minor_axis": minor_axes,
        "perimeter":  perimeters,
        "eccentricity": eccentricity,
        "convexity":  convexity,
    })
    meta.to_csv(
        os.path.join(pre_shape_output_dir, "Shape_Metadata.csv.gz"),
        index=False
    )

    with h5py.File(
        os.path.join(pre_shape_output_dir, "Pre_Shape_Space_Embedding.h5"), "w"
    ) as h5:
        h5.create_dataset("pre_shape_space_embedding", data=pre_shape_mat)


def run_kendall_tpca(pre_shape_input_dir, output_dir, cell_ids_to_analyze=None,
            cell_id_col="cell_id", max_frechet_mean_iter=1000, eta=1,
            use_parallel=False, num_threads=8, frechet_mean_tol=1e-4):
    """Run Kendall TPCA on a pre-shape embedding and write results to disk.

    Reads ``Pre_Shape_Space_Embedding.h5`` and ``Shape_Metadata.csv.gz``
    from ``pre_shape_input_dir``; writes ``TPCA_Info.h5`` to
    ``output_dir``.

    Parameters
    ----------
    pre_shape_input_dir : str
    output_dir : str
    cell_ids_to_analyze : list, optional
    cell_id_col : str
    max_frechet_mean_iter : int
    eta : float
    use_parallel : bool
    num_threads : int
    frechet_mean_tol : float
    """
    os.makedirs(output_dir, exist_ok=True)

    meta = pd.read_csv(
        os.path.join(pre_shape_input_dir, "Shape_Metadata.csv.gz")
    )
    if cell_ids_to_analyze is not None:
        mask = meta[cell_id_col].isin(cell_ids_to_analyze)
        idxes = meta.loc[mask, "numpy_idx"].values
    else:
        idxes = meta["numpy_idx"].values

    with h5py.File(
        os.path.join(pre_shape_input_dir, "Pre_Shape_Space_Embedding.h5"), "r"
    ) as h5:
        pre_shape_mat = h5["pre_shape_space_embedding"][:]

    with h5py.File(
        os.path.join(output_dir, "TPCA_Info.h5"), "w") as h5:
        p, lambdas, v, mu, U_flat = kendall_tpca(
            pre_shape_mat[idxes],
            max_frechet_mean_iter=max_frechet_mean_iter,
            eta=eta,
            use_parallel=use_parallel,
            frechet_mean_tol=frechet_mean_tol,
            num_chunks=num_threads,
            store_history=False,
        )
        h5.create_dataset("processed_idxes", data=idxes)
        h5.create_dataset("embedding",        data=p)
        h5.create_dataset("variances",        data=lambdas)
        h5.create_dataset("v_matrix",         data=v)
        h5.create_dataset("frechet_mean",     data=mu)
        h5.create_dataset("u_flattened",      data=U_flat)
        logger.info("TPCA_Info.h5 written to %s", output_dir)
