import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree, Voronoi
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics import silhouette_score
from itertools import combinations
import json
import sys

# =========================
# PARAMETERS
# =========================
PIXEL_SIZE        = 9.92   # Å/pixel (2.48 x bin4)
CUT_OFF_PIX       = 125    # pixels for RDF range (~124 nm at bin4)
BINS              = 50
N_ORIENT_CLUSTERS = 3      # number of orientation populations

# Neighbourhood radii for psi6 in Angstroms
# Set to ~1x, 1.5x, 2x, 3x your mean NND (42.6 nm = 426 Å)
PSI6_RADII = [300, 400, 500, 700]

# Distance restriction for end-on / packing-plane analysis
OUTLIER_SIGMA    = 2.0   # exclude fibrils > 2 std from cluster centre
NND_SEARCH_FACTOR = 1.5  # neighbourhood radius = 1.5 x mean NND

# Orientation clustering for end-on view
# Maximum k to try when searching for discrete orientation groups
MAX_ORIENT_K = 6
# Minimum silhouette score improvement to add another cluster
SILHOUETTE_THRESHOLD = 0.05
# Fibrils whose distance to nearest cluster centre > this multiple of
# the cluster's internal spread are kept as "isolated" (not forced into a group)
ISOLATION_SIGMA = 2.5


# =========================
# LOAD MODEL2POINT
# =========================
def load_model2point(file):
    fibrils = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = line.split()
            try:
                if len(vals) >= 4:
                    contour_id = int(vals[0])
                    xyz = [float(vals[1]), float(vals[2]), float(vals[3])]
                    fibrils.setdefault(contour_id, []).append(xyz)
            except ValueError:
                continue
    fibrils = [np.array(pts) for pts in fibrils.values() if len(pts) >= 2]
    print(f"Loaded {len(fibrils)} fibrils")
    return fibrils


# =========================
# FIBRIL PROPERTIES
# =========================
def get_fibril_properties(fibrils, pixel_size=PIXEL_SIZE):
    props = []
    for f in fibrils:
        centroid = f.mean(axis=0)           # 3D centroid in pixels
        centered = f - centroid
        _, s, vh = np.linalg.svd(centered)
        axis = vh[0]
        if axis[2] < 0:
            axis = -axis
        proj = centered @ axis
        length_pix = proj.max() - proj.min()
        length_ang = length_pix * pixel_size
        aspect_ratio = float(s[0] / s[1]) if len(s) >= 2 and s[1] > 0 else 1.0
        props.append({
            'centroid_pix': centroid,           # 3D, pixels
            'centroid_ang': centroid * pixel_size,  # 3D, Angstroms
            'axis': axis,
            'length_ang': length_ang,
            'aspect_ratio': aspect_ratio,
            'n_points': len(f)
        })
    return props


# =========================
# ORIENTATION ANALYSIS
# =========================
def orientation_analysis(props, n_clusters=N_ORIENT_CLUSTERS):
    axes    = np.array([p['axis'] for p in props])
    azimuth = np.degrees(np.arctan2(axes[:, 1], axes[:, 0]))
    polar   = np.degrees(np.arccos(np.clip(axes[:, 2], -1, 1)))

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=20)
    labels = kmeans.fit_predict(axes)

    print(f"\n--- ORIENTATION CLUSTERS (k={n_clusters}) ---")
    for k in range(n_clusters):
        members   = axes[labels == k]
        mean_axis = members.mean(axis=0)
        mean_axis /= np.linalg.norm(mean_axis)
        az  = np.degrees(np.arctan2(mean_axis[1], mean_axis[0]))
        pol = np.degrees(np.arccos(np.clip(mean_axis[2], -1, 1)))
        print(f"  Cluster {k}: {len(members)} fibrils, "
              f"mean azimuth={az:.1f}°, polar={pol:.1f}°")

    return labels, azimuth, polar


# =========================
# DISCRETE ORIENTATION CLUSTERING WITH ISOLATED FIBRILS
# =========================
def find_discrete_orientation_groups(props):
    """
    Finds discrete orientation sub-groups using orientation tensor clustering.

    Unlike orientation_analysis() which forces all fibrils into k groups,
    this function:
      1. Uses the outer-product (orientation tensor) representation so
         the clustering is direction-invariant (axis and -axis are the same)
      2. Searches k=2..MAX_ORIENT_K and picks the k with the best
         silhouette score, as long as adding k improves it by
         SILHOUETTE_THRESHOLD
      3. After clustering, fibrils that are more than ISOLATION_SIGMA
         standard deviations from their cluster centre are labelled -1
         (isolated / ungrouped) rather than forced into a group

    Returns
    -------
    group_labels : np.ndarray (n_fibrils,)
        -1 = isolated, 0..K-1 = group id
    group_axes   : list of np.ndarray
        Mean orientation axis per group (unit vector)
    k_best       : int
        Number of groups found
    """
    axes = np.array([p['axis'] for p in props])
    n    = len(axes)

    # Normalise and flip to consistent hemisphere
    norms = np.linalg.norm(axes, axis=1, keepdims=True)
    axes  = axes / np.where(norms > 0, norms, 1)
    mean_all = axes.mean(axis=0)
    for i in range(n):
        if np.dot(axes[i], mean_all) < 0:
            axes[i] = -axes[i]

    # Orientation tensor representation: outer product, flattened
    Q = np.array([np.outer(a, a).ravel() for a in axes])   # (n, 9)

    # Search for best k
    best_k    = 1
    best_sil  = -1.0
    best_lbls = np.zeros(n, dtype=int)

    print(f"\n--- DISCRETE ORIENTATION GROUP SEARCH ---")
    for k in range(2, min(MAX_ORIENT_K + 1, n)):
        km  = KMeans(n_clusters=k, random_state=42, n_init=20)
        lbl = km.fit_predict(Q)
        if len(set(lbl)) < 2:
            continue
        sil = silhouette_score(Q, lbl)
        print(f"  k={k}: silhouette={sil:.4f}")
        if sil > best_sil + SILHOUETTE_THRESHOLD:
            best_sil  = sil
            best_k    = k
            best_lbls = lbl

    print(f"  → Best k={best_k}  silhouette={best_sil:.4f}")

    # Compute per-group mean axes
    group_axes = []
    for g in range(best_k):
        mask = best_lbls == g
        m    = axes[mask].mean(axis=0)
        m   /= np.linalg.norm(m)
        group_axes.append(m)

    # Compute angle of each fibril to its cluster centre
    # Fibrils far from their centre → isolated
    group_labels = best_lbls.copy()
    for g in range(best_k):
        mask   = best_lbls == g
        ga     = group_axes[g]
        angles = np.degrees(np.arccos(np.clip(
            [abs(np.dot(axes[i], ga)) for i in np.where(mask)[0]], -1, 1)))
        thresh = angles.mean() + ISOLATION_SIGMA * angles.std()
        idxs   = np.where(mask)[0]
        for ii, ang in zip(idxs, angles):
            if ang > thresh:
                group_labels[ii] = -1   # isolated

    n_isolated = (group_labels == -1).sum()
    print(f"  Isolated fibrils (not in any group): {n_isolated}")
    for g in range(best_k):
        mask = group_labels == g
        ang_between = np.degrees(np.arccos(np.clip(
            abs(np.dot(group_axes[g],
                       group_axes[(g+1) % best_k])), -1, 1)))
        print(f"  Group {g}: n={mask.sum()}  "
              f"mean_axis={group_axes[g].round(3)}")
    if best_k >= 2:
        for ga, gb in combinations(range(best_k), 2):
            ang = np.degrees(np.arccos(np.clip(
                abs(np.dot(group_axes[ga], group_axes[gb])), -1, 1)))
            print(f"  Angle between Group {ga} and Group {gb}: {ang:.1f}°")

    return group_labels, group_axes, best_k



# =========================
# PER-GROUP CHAIN DETECTION
# =========================
def detect_chains_in_group(proj, group_name='',
                            collinearity_thresh=168.0,
                            search_factor=1.4):
    """
    Detect collinear chains within a single orientation group's
    end-on projection.  Uses the projected 2D positions of ONLY
    the fibrils belonging to that group.

    Returns chain_labels, chains, stats dict.
    """
    n = len(proj)
    if n < 3:
        return np.full(n, -1, dtype=int), [], {}

    tree = cKDTree(proj)
    nnds = [tree.query(proj[i], k=2)[0][1] for i in range(n)]
    r    = np.mean(nnds) * search_factor

    nb = {i: [j for j in tree.query_ball_point(proj[i], r=r) if j != i]
          for i in range(n)}

    row_adj = {i: set() for i in range(n)}
    for k in range(n):
        nbs = nb[k]
        if len(nbs) < 2:
            continue
        vecs = {j: proj[j] - proj[k] for j in nbs}
        nbl  = list(nbs)
        for a in range(len(nbl)):
            for b in range(a + 1, len(nbl)):
                i, j = nbl[a], nbl[b]
                vi = vecs[i] / (np.linalg.norm(vecs[i]) + 1e-9)
                vj = vecs[j] / (np.linalg.norm(vecs[j]) + 1e-9)
                ang = np.degrees(np.arccos(np.clip(np.dot(vi, vj), -1, 1)))
                if ang >= collinearity_thresh:
                    row_adj[i].add(j); row_adj[j].add(i)

    visited = np.zeros(n, dtype=bool)
    chains  = []
    for start in range(n):
        if visited[start] or not row_adj[start]:
            visited[start] = True; continue
        chain = []; queue = [start]; visited[start] = True
        while queue:
            node = queue.pop(0); chain.append(node)
            for nb_ in row_adj[node]:
                if not visited[nb_]: visited[nb_] = True; queue.append(nb_)
        if len(chain) >= 2:
            chains.append(chain)

    chain_set = {fi: ci for ci, c in enumerate(chains) for fi in c}
    for k in range(n):
        if k in chain_set: continue
        nbs = nb[k]
        if len(nbs) < 2: continue
        vecs = {j: proj[j] - proj[k] for j in nbs}
        nbl  = list(nbs)
        for a in range(len(nbl)):
            for b in range(a + 1, len(nbl)):
                i, j = nbl[a], nbl[b]
                if i not in chain_set or j not in chain_set: continue
                if chain_set[i] != chain_set[j]: continue
                vi = vecs[i] / (np.linalg.norm(vecs[i]) + 1e-9)
                vj = vecs[j] / (np.linalg.norm(vecs[j]) + 1e-9)
                ang = np.degrees(np.arccos(np.clip(np.dot(vi, vj), -1, 1)))
                if ang >= collinearity_thresh:
                    ci = chain_set[i]; chains[ci].append(k); chain_set[k] = ci; break

    labels = np.full(n, -1, dtype=int)
    chains.sort(key=len, reverse=True)
    for ci, c in enumerate(chains):
        for fi in c: labels[fi] = ci

    cl       = [len(c) for c in chains]
    n_in     = int((labels >= 0).sum())
    stats = {
        'n_chains':          len(chains),
        'n_isolated':        int((labels == -1).sum()),
        'n_in_chains':       n_in,
        'frac_in_chains':    float(n_in / n),
        'mean_chain_length': float(np.mean(cl)) if cl else 0.0,
        'max_chain_length':  int(max(cl)) if cl else 0,
        'chain_lengths':     cl,
    }
    if group_name:
        print(f"  {group_name}: chains={len(chains)}  isolated={stats['n_isolated']}  "
              f"in_chains={n_in}/{n} ({n_in/n*100:.0f}%)  "
              f"max={stats['max_chain_length']}  mean={stats['mean_chain_length']:.1f}")
    return labels, chains, stats

# =========================
# PER-GROUP END-ON PROJECTION  (THE IMOD VIEW)
# =========================
def orientation_group_endon(props, group_labels, group_axes):
    """
    For each discrete orientation group, projects the 3D centroids onto
    the plane perpendicular to that group's mean axis — i.e. exactly the
    view you get in IMOD when you rotate to look straight down the fibrils.

    Uses the TRUE 3D centroid positions (centroid_ang, all three components).

    For each group returns:
      proj_x, proj_y  — 2D coordinates in the perpendicular plane (nm)
      depth           — distance along the viewing axis (nm)
                        used to fade fibrils in the "other" group
      psi6_pp         — ψ₆ computed in that plane
      psi6_pp_per     — per-fibril ψ₆
      plane_angle     — angle of viewing axis from Z (°)
      normal          — viewing axis (unit vector)

    Also returns the "global" entry using the overall director.
    """
    centers3d = np.array([p['centroid_ang'] for p in props]) / 10.0  # nm
    n         = len(centers3d)
    n_groups  = len(group_axes)

    def _project(centers, view_axis):
        """Project 3D points onto plane perpendicular to view_axis."""
        v = view_axis / np.linalg.norm(view_axis)
        # Remove depth component
        depth    = centers.dot(v)
        in_plane = centers - depth[:, None] * v
        # Two orthonormal in-plane axes via Gram-Schmidt
        ref = np.array([1., 0., 0.])
        if abs(np.dot(ref, v)) > 0.9:
            ref = np.array([0., 1., 0.])
        e1 = ref - np.dot(ref, v) * v;  e1 /= np.linalg.norm(e1)
        e2 = np.cross(v, e1);            e2 /= np.linalg.norm(e2)
        return in_plane.dot(e1), in_plane.dot(e2), depth, e1, e2

    def _psi6(px_own, py_own, px_all, py_all):
        """
        Compute psi6 for own-group fibrils.
        Neighbour search uses ALL fibrils (real physical spacing).
        Only same-group neighbours contribute to the angle calculation.
        This avoids the artefact where using only half the fibrils
        inflates the apparent NND and distorts the search radius.
        """
        pts_own = np.column_stack([px_own, py_own])
        pts_all = np.column_stack([px_all, py_all])
        n_own   = len(pts_own)
        n_all   = len(pts_all)
        if n_own < 3:
            return 0.0, [0.0] * n_own
        # NND using ALL fibrils -> correct physical spacing
        tree_all = cKDTree(pts_all)
        nnds_all = [tree_all.query(pts_all[i], k=2)[0][1] for i in range(n_all)]
        r        = np.mean(nnds_all) * NND_SEARCH_FACTOR
        # For each own-group fibril: find ALL neighbours within r,
        # keep only those that are also own-group fibrils
        # We need a mapping: for each position in pts_own,
        # find its index in pts_all
        # pts_own is a subset of pts_all -- find matching indices
        tree_own = cKDTree(pts_own)
        vals = []
        for i in range(n_own):
            c   = pts_own[i]
            # find all fibrils within r in the full set
            idx_all = tree_all.query_ball_point(c, r=r)
            # keep only those whose position matches a pts_own position
            nbs = []
            for j_all in idx_all:
                p = pts_all[j_all]
                # is this point in pts_own?
                dist, j_own = tree_own.query(p, k=1)
                if dist < 0.1 and j_own != i:   # same point, different fibril
                    nbs.append(j_own)
            if len(nbs) < 2:
                vals.append(np.nan)
                continue
            angs = [np.arctan2(pts_own[j, 1] - c[1], pts_own[j, 0] - c[0])
                    for j in nbs]
            vals.append(float(np.abs(np.mean(np.exp(6j * np.array(angs))))))
        valid = [v for v in vals if not np.isnan(v)]
        return float(np.mean(valid)) if valid else 0.0, \
               [0.0 if np.isnan(v) else v for v in vals]

    results = {}

    print(f"\n--- ORIENTATION GROUP END-ON PROJECTIONS ---")
    print(f"  {'Group':<15} {'N':>5}  {'psi6_pp':>10}  "
          f"{'view->Z':>9}  Note")
    print(f"  {'-'*15} {'-'*5}  {'-'*10}  {'-'*9}  {'-'*25}")

    for g, gax in enumerate(group_axes):
        mask     = group_labels == g
        n_g      = mask.sum()
        if n_g < 4:
            print(f"  {'Group '+str(g):<15} {n_g:>5}  skipping (< 4)")
            continue

        c_g      = centers3d[mask]
        px, py, depth_g, e1, e2 = _project(c_g, gax)
        # All fibrils projected along this group's axis (for correct NND)
        px_g_all, py_g_all, _, _, _ = _project(centers3d, gax)
        psi6, psi6_per = _psi6(px, py, px_g_all, py_g_all)
        plane_angle    = float(np.degrees(np.arccos(
            np.clip(abs(np.dot(gax / np.linalg.norm(gax),
                               np.array([0., 0., 1.]))), 0, 1))))

        # Also project ALL fibrils with this viewing axis
        # (so we can show the "other" group faded by depth)
        px_all, py_all, depth_all, _, _ = _project(centers3d, gax)

        # Per-group chain detection on own-axis projection
        chain_lbl_g, chains_g, chain_stats_g = detect_chains_in_group(
            np.column_stack([px, py]),
            group_name=f'Group {g}')

        results[g] = {
            'label':         f'Group {g}',
            'n':             int(n_g),
            'psi6_pp':       psi6,
            'psi6_pp_per':   psi6_per,
            'plane_angle':   plane_angle,
            'normal':        gax.tolist(),
            'proj_x':        px.tolist(),
            'proj_y':        py.tolist(),
            'depth':         depth_g.tolist(),
            'proj_all_x':    px_all.tolist(),
            'proj_all_y':    py_all.tolist(),
            'depth_all':     depth_all.tolist(),
            'group_labels_all': group_labels.tolist(),
            'chain_labels':  chain_lbl_g.tolist(),
            'chain_stats':   chain_stats_g,
        }

        print(f"  {'Group '+str(g):<15} {n_g:>5}  {psi6:>10.4f}  "
              f"{plane_angle:>8.1f}°  "
              f"{'Strong' if psi6>0.6 else 'Moderate' if psi6>0.3 else 'Weak'}")

    # Global: project along overall mean director
    axes_all  = np.array([p['axis'] for p in props])
    Q_all     = np.array([np.outer(a, a) for a in axes_all]).mean(axis=0)
    evals, evecs = np.linalg.eigh(Q_all)
    director  = evecs[:, np.argmax(evals)]
    if director[2] < 0:
        director = -director

    px_g, py_g, depth_glo, _, _ = _project(centers3d, director)
    psi6_glo, psi6_per_glo       = _psi6(px_g, py_g, px_g, py_g)
    ang_glo = float(np.degrees(np.arccos(np.clip(abs(director[2]), 0, 1))))

    chain_lbl_glo, _, chain_stats_glo = detect_chains_in_group(
        np.column_stack([px_g, py_g]), group_name='Global')

    results['global'] = {
        'label':         'Global',
        'n':             n,
        'psi6_pp':       psi6_glo,
        'psi6_pp_per':   psi6_per_glo,
        'plane_angle':   ang_glo,
        'normal':        director.tolist(),
        'proj_x':        px_g.tolist(),
        'proj_y':        py_g.tolist(),
        'depth':         depth_glo.tolist(),
        'proj_all_x':    px_g.tolist(),
        'proj_all_y':    py_g.tolist(),
        'depth_all':     depth_glo.tolist(),
        'group_labels_all': group_labels.tolist(),
        'chain_labels':  chain_lbl_glo.tolist(),
        'chain_stats':   chain_stats_glo,
    }

    print(f"  {'Global':<15} {n:>5}  {psi6_glo:>10.4f}  "
          f"{ang_glo:>8.1f}°  (all fibrils, mean director)")

    return results


# =========================
# ALIGNMENT / NEMATIC ORDER
# =========================
def alignment_analysis(props):
    axes = np.array([p['axis'] for p in props])

    Q = np.zeros((3, 3))
    for a in axes:
        Q += np.outer(a, a) - np.eye(3) / 3.0
    Q /= len(axes)
    eigenvalues = np.linalg.eigvalsh(Q)
    S2 = float(1.5 * eigenvalues.max())

    eigenvectors = np.linalg.eigh(Q)[1]
    director = eigenvectors[:, np.argmax(eigenvalues)]

    print(f"\n--- ALIGNMENT / NEMATIC ORDER ---")
    print(f"  S2 nematic order parameter: {S2:.4f}")
    print(f"  {'Strong alignment' if S2 > 0.7 else 'Moderate alignment' if S2 > 0.4 else 'Weak/no alignment'}")
    print(f"  Director axis: [{director[0]:.3f}, {director[1]:.3f}, {director[2]:.3f}]")

    pairwise_angles = []
    for i in range(len(axes)):
        for j in range(i + 1, len(axes)):
            cos_a = abs(np.dot(axes[i], axes[j]))
            angle = np.degrees(np.arccos(np.clip(cos_a, 0, 1)))
            pairwise_angles.append(angle)
    pairwise_angles = np.array(pairwise_angles)

    print(f"  Mean pairwise fibril angle: {np.mean(pairwise_angles):.1f}°")
    print(f"  Std pairwise fibril angle:  {np.std(pairwise_angles):.1f}°")

    return S2, director, pairwise_angles


# =========================
# NEAREST NEIGHBOUR DISTANCES
# =========================
def nearest_neighbour_distances(props, orientation_labels):
    centers = np.array([p['centroid_ang'] for p in props])
    tree    = cKDTree(centers)
    dists, idxs = tree.query(centers, k=2)
    nnd = dists[:, 1] / 10  # nm

    same_nnd, diff_nnd = [], []
    for i in range(len(centers)):
        j = idxs[i, 1]
        (same_nnd if orientation_labels[i] == orientation_labels[j] else diff_nnd).append(nnd[i])

    print(f"\n--- NEAREST NEIGHBOUR DISTANCES ---")
    print(f"  Mean NND:   {np.mean(nnd):.1f} nm")
    print(f"  Median NND: {np.median(nnd):.1f} nm")
    print(f"  Std NND:    {np.std(nnd):.1f} nm")
    if same_nnd:
        print(f"  Same-orientation NND: {np.mean(same_nnd):.1f} nm (n={len(same_nnd)})")
    if diff_nnd:
        print(f"  Diff-orientation NND: {np.mean(diff_nnd):.1f} nm (n={len(diff_nnd)})")

    return nnd, same_nnd, diff_nnd


# =========================
# FIBRIL LENGTH DISTRIBUTION
# =========================
def length_distribution(props):
    lengths = [p['length_ang'] / 10 for p in props]
    print(f"\n--- FIBRIL LENGTHS ---")
    print(f"  Mean:   {np.mean(lengths):.1f} nm")
    print(f"  Median: {np.median(lengths):.1f} nm")
    print(f"  Std:    {np.std(lengths):.1f} nm")
    return lengths


# =========================
# AXIS-TO-AXIS DISTANCES
# =========================
def axis_to_axis_distances(fibrils, props, pixel_size=PIXEL_SIZE, sample_every=5):
    print("\n--- AXIS-TO-AXIS DISTANCES ---")
    sampled   = [f[::sample_every] * pixel_size for f in fibrils]
    min_dists = []
    for i in range(len(sampled)):
        tree      = cKDTree(sampled[i])
        local_min = []
        for j in range(len(sampled)):
            if i == j:
                continue
            dists, _ = tree.query(sampled[j])
            local_min.append(dists.min())
        if local_min:
            min_dists.append(min(local_min))
    min_dists = np.array(min_dists) / 10
    print(f"  Mean:   {np.mean(min_dists):.1f} nm")
    print(f"  Median: {np.median(min_dists):.1f} nm")
    print(f"  Std:    {np.std(min_dists):.1f} nm")
    return min_dists


# =========================
# BUNDLE DETECTION
# =========================
def bundle_analysis(props, orientation_labels):
    centers    = np.array([p['centroid_ang'] for p in props])
    clustering = DBSCAN(eps=500, min_samples=3).fit(centers)
    bundle_labels = clustering.labels_
    n_bundles  = len(set(bundle_labels)) - (1 if -1 in bundle_labels else 0)
    n_noise    = (bundle_labels == -1).sum()
    print(f"\n--- BUNDLE ANALYSIS ---")
    print(f"  Bundles found: {n_bundles}")
    print(f"  Isolated fibrils: {n_noise}")
    return bundle_labels


# =========================
# VORONOI ANALYSIS
# =========================
def voronoi_analysis(props):
    centers = np.array([p['centroid_ang'] for p in props])[:, :2]

    if len(centers) < 4:
        print(f"  Skipping Voronoi — need at least 4 points (have {len(centers)})")
        return None, np.array([])

    vor              = Voronoi(centers)
    neighbour_counts = np.zeros(len(centers), dtype=int)
    for ridge in vor.ridge_points:
        neighbour_counts[ridge[0]] += 1
        neighbour_counts[ridge[1]] += 1

    print("\n--- VORONOI COORDINATION NUMBERS ---")
    unique, counts = np.unique(neighbour_counts, return_counts=True)
    for u, c in zip(unique, counts):
        print(f"  {u} neighbours: {c} fibrils")

    frac_6 = (neighbour_counts == 6).sum() / len(neighbour_counts)
    print(f"  Fraction with exactly 6 neighbours: {frac_6:.3f} "
          f"({'consistent with hexagonal packing' if frac_6 > 0.4 else 'not hexagonally packed'})")

    return vor, neighbour_counts


# =========================
# VORONOI PER CLUSTER
# =========================
def voronoi_analysis_per_cluster(props, orientation_labels):
    centers_2d  = np.array([p['centroid_ang'][:2] for p in props])
    cluster_ids = sorted(set(orientation_labels))

    print("\n--- VORONOI COORDINATION PER CLUSTER ---")
    print(f"  {'Population':<20} {'N':>5}  {'6-fold frac':>12}  {'Mean coord':>12}  Assessment")
    print(f"  {'-'*20} {'-'*5}  {'-'*12}  {'-'*12}  {'-'*25}")

    results = {}

    if len(centers_2d) >= 4:
        vor = Voronoi(centers_2d)
        nc  = np.zeros(len(centers_2d), dtype=int)
        for ridge in vor.ridge_points:
            nc[ridge[0]] += 1
            nc[ridge[1]] += 1
        frac6 = (nc == 6).sum() / len(nc)
        print(f"  {'ALL (global)':<20} {len(centers_2d):>5}  {frac6:>12.3f}  "
              f"{nc.mean():>12.2f}  "
              f"{'consistent with hex' if frac6 > 0.4 else 'not hexagonally packed'}")
        results['global'] = {'frac_6': float(frac6), 'mean_coord': float(nc.mean()),
                             'counts': nc, 'n': len(centers_2d)}

    for k in cluster_ids:
        mask = orientation_labels == k
        c2d  = centers_2d[mask]
        if len(c2d) < 4:
            print(f"  {'Cluster '+str(k):<20} {mask.sum():>5}  — skipping (need >=4)")
            continue
        vor = Voronoi(c2d)
        nc  = np.zeros(len(c2d), dtype=int)
        for ridge in vor.ridge_points:
            nc[ridge[0]] += 1
            nc[ridge[1]] += 1
        frac6 = (nc == 6).sum() / len(nc)
        print(f"  {'Cluster '+str(k):<20} {mask.sum():>5}  {frac6:>12.3f}  "
              f"{nc.mean():>12.2f}  "
              f"{'consistent with hex' if frac6 > 0.4 else 'not hexagonally packed'}")
        results[k] = {'frac_6': float(frac6), 'mean_coord': float(nc.mean()),
                      'counts': nc, 'n': int(mask.sum())}

    return results


# =========================
# BOND ANGLE DISTRIBUTION
# =========================
def bond_angle_distribution(centers_2d, radius=500):
    tree       = cKDTree(centers_2d)
    all_angles = []
    for i, c in enumerate(centers_2d):
        idx        = tree.query_ball_point(c, r=radius)
        neighbours = [j for j in idx if j != i]
        if len(neighbours) < 2:
            continue
        vecs = [centers_2d[j] - c for j in neighbours]
        for a, b in combinations(vecs, 2):
            cos_angle = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9)
            angle     = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
            all_angles.append(angle)
    return np.array(all_angles)


# =========================
# 2D PSI6
# =========================
def psi6_2d(centers_2d, radii=None):
    if radii is None:
        radii = PSI6_RADII

    tree    = cKDTree(centers_2d)
    results = {}

    for r in radii:
        psi_vals = []
        psi_per  = []
        for i, c in enumerate(centers_2d):
            idx        = tree.query_ball_point(c, r=r)
            neighbours = [j for j in idx if j != i]
            if len(neighbours) < 2:
                psi_per.append(0.0)
                continue
            angles = [np.arctan2((centers_2d[j] - c)[1], (centers_2d[j] - c)[0])
                      for j in neighbours]
            psi = np.abs(np.mean(np.exp(6j * np.array(angles))))
            psi_vals.append(psi)
            psi_per.append(psi)
        results[r] = {
            'mean':       float(np.mean(psi_vals)) if psi_vals else 0.0,
            'per_fibril': np.array(psi_per)
        }

    compare_r = PSI6_RADII[1]
    for k, name in [(4, 'psi4'), (2, 'psi2')]:
        psi_vals = []
        for i, c in enumerate(centers_2d):
            idx        = tree.query_ball_point(c, r=compare_r)
            neighbours = [j for j in idx if j != i]
            if len(neighbours) < 2:
                continue
            angles = [np.arctan2((centers_2d[j] - c)[1], (centers_2d[j] - c)[0])
                      for j in neighbours]
            psi_vals.append(np.abs(np.mean(np.exp(k * 1j * np.array(angles)))))
        results[name] = float(np.mean(psi_vals)) if psi_vals else 0.0

    return results


# =========================
# 2D RDF
# =========================
def rdf_2d(centers_2d, cutoff_ang=None, bins=BINS):
    if cutoff_ang is None:
        cutoff_ang = CUT_OFF_PIX * PIXEL_SIZE
    if len(centers_2d) < 3:
        r = np.linspace(0, cutoff_ang, bins)
        return r, np.zeros(bins), np.zeros(bins)

    tree  = cKDTree(centers_2d)
    dists = []
    for c in centers_2d:
        idx = tree.query_ball_point(c, cutoff_ang)
        for j in idx:
            d = np.linalg.norm(c - centers_2d[j])
            if d > 0:
                dists.append(d)

    dists       = np.array(dists)
    hist, edges = np.histogram(dists, bins=bins, range=(0, cutoff_ang))
    r           = 0.5 * (edges[1:] + edges[:-1])

    shell_area = np.pi * (edges[1:]**2 - edges[:-1]**2)
    area       = float(np.prod(centers_2d.max(axis=0) - centers_2d.min(axis=0)))
    density    = len(centers_2d) / (area + 1e-9)
    g_r        = hist / (len(centers_2d) * shell_area * density + 1e-9)

    mins, maxs  = centers_2d.min(axis=0), centers_2d.max(axis=0)
    rand        = np.random.uniform(mins, maxs, size=centers_2d.shape)
    rand_tree   = cKDTree(rand)
    rand_dists  = []
    for c in rand:
        idx = rand_tree.query_ball_point(c, cutoff_ang)
        for j in idx:
            d = np.linalg.norm(c - rand[j])
            if d > 0:
                rand_dists.append(d)
    rand_hist, _ = np.histogram(rand_dists, bins=bins, range=(0, cutoff_ang))
    g_rand        = rand_hist / (len(rand) * shell_area * density + 1e-9)

    return r, g_r, g_rand


# =========================
# PER-CLUSTER PSI6 AND 2D RDF
# =========================
def per_cluster_analysis(props, orientation_labels):
    centers_2d  = np.array([p['centroid_ang'][:2] for p in props])
    cluster_ids = sorted(set(orientation_labels))

    r0 = PSI6_RADII[0]
    r1 = PSI6_RADII[2]

    print("\n--- PER-CLUSTER psi6 AND 2D RDF ---")
    print(f"  {'Population':<20} {'N':>5}  "
          f"{'psi6(r='+str(r0//10)+'nm)':>12}  "
          f"{'psi6(r='+str(r1//10)+'nm)':>12}  "
          f"{'psi4':>8}  {'1st peak nm':>12}")
    print(f"  {'-'*20} {'-'*5}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*12}")

    cluster_data = {}

    psi_g            = psi6_2d(centers_2d)
    r_g, gr_g, gr_rg = rdf_2d(centers_2d)
    peak_g           = float(r_g[np.argmax(gr_g[1:]) + 1]) / 10 if len(gr_g) > 1 else 0.0
    ba_g             = bond_angle_distribution(centers_2d)

    print(f"  {'ALL (global)':<20} {len(centers_2d):>5}  "
          f"{psi_g[r0]['mean']:>12.4f}  {psi_g[r1]['mean']:>12.4f}  "
          f"{psi_g.get('psi4', 0):>8.4f}  {peak_g:>12.1f}")

    cluster_data['global'] = {
        'label': 'Global', 'n': len(centers_2d),
        'psi6': psi_g,
        'r': r_g, 'g_r': gr_g, 'g_rand': gr_rg,
        'centers_2d': centers_2d,
        'bond_angles': ba_g
    }

    for k in cluster_ids:
        mask = orientation_labels == k
        if mask.sum() < 5:
            print(f"  Cluster {k}: only {mask.sum()} fibrils — skipping (need >=5)")
            continue
        c2d               = centers_2d[mask]
        psi_k             = psi6_2d(c2d)
        r_k, gr_k, gr_rk  = rdf_2d(c2d)
        peak_k            = float(r_k[np.argmax(gr_k[1:]) + 1]) / 10 if len(gr_k) > 1 else 0.0
        ba_k              = bond_angle_distribution(c2d)

        print(f"  {'Cluster '+str(k):<20} {mask.sum():>5}  "
              f"{psi_k[r0]['mean']:>12.4f}  {psi_k[r1]['mean']:>12.4f}  "
              f"{psi_k.get('psi4', 0):>8.4f}  {peak_k:>12.1f}")

        cluster_data[k] = {
            'label': f'Cluster {k}', 'n': int(mask.sum()),
            'psi6': psi_k,
            'r': r_k, 'g_r': gr_k, 'g_rand': gr_rk,
            'centers_2d': c2d,
            'bond_angles': ba_k
        }

    return cluster_data


# =========================
# PLOT PER-CLUSTER RESULTS
# =========================
def plot_per_cluster(cluster_data):
    COLORS = ['#00d4ff', '#ff6b35', '#7fff6b', '#c87fff', '#ffcc44']

    fig, axes = plt.subplots(1, 4, figsize=(22, 5), facecolor='#0d1117')
    fig.suptitle('Per-Orientation-Cluster Analysis', color='white',
                 fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')

    ax = axes[0]; style(ax)
    for i, (k, d) in enumerate(cluster_data.items()):
        lw = 2.5 if k == 'global' else 1.8
        ls = '--'  if k == 'global' else '-'
        ax.plot(d['r'] / 10, d['g_r'],
                color=COLORS[i % len(COLORS)], lw=lw, ls=ls,
                label=f"{d['label']} (n={d['n']})", alpha=0.9)
    ax.axhline(1, color='#444', lw=0.8, ls=':')
    ax.set_xlabel('r (nm)'); ax.set_ylabel('g(r)')
    ax.set_title('2D RDF per Cluster')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax = axes[1]; style(ax)
    keys  = list(cluster_data.keys())
    r0, r1 = PSI6_RADII[0], PSI6_RADII[2]
    psiA  = [cluster_data[k]['psi6'][r0]['mean']      for k in keys]
    psiB  = [cluster_data[k]['psi6'][r1]['mean']      for k in keys]
    psi4  = [cluster_data[k]['psi6'].get('psi4', 0)  for k in keys]
    x = np.arange(len(keys)); w = 0.25
    for offset, vals, lbl in [(-w, psiA, f'psi6 r={r0//10}nm'),
                                (0,  psiB, f'psi6 r={r1//10}nm'),
                                (w,  psi4, 'psi4')]:
        bars = ax.bar(x + offset, vals, width=w, label=lbl,
                      edgecolor='#0d1117', alpha=0.6 if 'psi4' in lbl else 1.0)
        for bar, col, val in zip(bars,
                                  [COLORS[i % len(COLORS)] for i in range(len(keys))],
                                  vals):
            bar.set_color(col)
            ax.text(bar.get_x() + bar.get_width()/2, val + 0.01,
                    f'{val:.2f}', ha='center', va='bottom', color='white', fontsize=7)
    ax.axhline(0.6, color='#7fff6b', lw=1, ls='--', alpha=0.5, label='Strong (0.6)')
    ax.axhline(0.3, color='#ff6b35', lw=1, ls='--', alpha=0.5, label='Weak (0.3)')
    ax.set_xticks(x)
    ax.set_xticklabels([cluster_data[k]['label'] for k in keys],
                       color='#8b949e', fontsize=8, rotation=15)
    ax.set_ylim(0, 1); ax.set_ylabel('Order parameter')
    ax.set_title('psi6 / psi4 per Cluster')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax = axes[2]; style(ax)
    for i, (k, d) in enumerate(cluster_data.items()):
        ba = d['bond_angles']
        if len(ba) == 0: continue
        hist, edges = np.histogram(ba, bins=36, range=(0, 180))
        centres = 0.5 * (edges[1:] + edges[:-1])
        lw = 2.5 if k == 'global' else 1.5
        ls = '--'  if k == 'global' else '-'
        ax.plot(centres, hist / (hist.max() + 1e-9),
                color=COLORS[i % len(COLORS)], lw=lw, ls=ls,
                label=d['label'], alpha=0.9)
    for deg, col, lbl in [(60, '#ffffff', '60°'), (90, '#888', '90°'), (120, '#ffffff', '120°')]:
        ax.axvline(deg, color=col, lw=0.8, ls=':', alpha=0.5)
        ax.text(deg + 1, 0.94, lbl, color=col, fontsize=7)
    ax.set_xlabel('Bond angle (°)'); ax.set_ylabel('Normalised count')
    ax.set_title('Bond Angle Distribution\n(60°=hex, 90°=square)')
    ax.set_xlim(0, 180)
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax = axes[3]; style(ax)
    for i, (k, d) in enumerate(cluster_data.items()):
        vals = [d['psi6'][r]['mean'] for r in PSI6_RADII]
        lw = 2.5 if k == 'global' else 1.8
        ls = '--'  if k == 'global' else '-'
        ax.plot([r / 10 for r in PSI6_RADII], vals, 'o' + ls,
                color=COLORS[i % len(COLORS)], lw=lw, ms=6,
                label=d['label'], alpha=0.9)
    ax.axhline(0.6, color='#7fff6b', lw=1, ls='--', alpha=0.4)
    ax.axhline(0.3, color='#ff6b35', lw=1, ls='--', alpha=0.4)
    ax.set_xlabel('Neighbourhood radius (nm)'); ax.set_ylabel('psi6')
    ax.set_ylim(0, 1); ax.set_title('psi6 vs Radius per Cluster')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig('per_cluster_analysis.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved per_cluster_analysis.png")


# =========================
# PLOT ALIGNMENT ANALYSIS
# =========================
def plot_alignment(pairwise_angles, S2, props, orientation_labels):
    COLORS = ['#00d4ff', '#ff6b35', '#7fff6b', '#c87fff', '#ffcc44']
    axes_arr = np.array([p['axis'] for p in props])

    fig, ax_list = plt.subplots(1, 3, figsize=(18, 5), facecolor='#0d1117')
    fig.suptitle('Fibril Alignment Analysis', color='white',
                 fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')

    ax = ax_list[0]; style(ax)
    ax.hist(pairwise_angles, bins=36, range=(0, 90), color='#00d4ff',
            edgecolor='#0d1117', alpha=0.8)
    ax.axvline(np.mean(pairwise_angles), color='white', lw=1.5, ls='--',
               label=f'Mean: {np.mean(pairwise_angles):.1f}°')
    ax.set_xlabel('Pairwise angle between fibril axes (°)')
    ax.set_ylabel('Count')
    ax.set_title(f'Fibril Axis Alignment\nS2 = {S2:.3f}')
    ax.legend(fontsize=8, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')
    col = '#7fff6b' if S2 > 0.7 else '#ffcc44' if S2 > 0.4 else '#ff6b35'
    ax.text(0.97, 0.95, f'S2 = {S2:.3f}', transform=ax.transAxes,
            ha='right', va='top', color=col, fontsize=12, fontweight='bold')

    ax = ax_list[1]
    ax.set_facecolor('#161b22')
    ax.remove()
    ax = fig.add_subplot(1, 3, 2, projection='polar')
    ax.set_facecolor('#161b22')
    ax.tick_params(colors='#8b949e', labelsize=7)
    ax.title.set_color('#e6edf3')

    cluster_ids = sorted(set(orientation_labels))
    for ki, k in enumerate(cluster_ids):
        mask    = orientation_labels == k
        az_k    = np.arctan2(axes_arr[mask, 1], axes_arr[mask, 0])
        doubled = (2 * az_k) % (2 * np.pi)
        hist_k, bin_edges = np.histogram(doubled, bins=36)
        theta = (bin_edges[:-1] + bin_edges[1:]) / 2
        width = bin_edges[1] - bin_edges[0]
        ax.bar(theta, hist_k, width=width, bottom=0,
               color=COLORS[ki % len(COLORS)], alpha=0.7,
               label=f'Cluster {k}')
    ax.set_title('Azimuthal Orientation\nby Cluster', pad=15)
    ax.legend(fontsize=7, loc='upper right', bbox_to_anchor=(1.3, 1.1),
              labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax = ax_list[2]; style(ax)
    polar_all = np.degrees(np.arccos(np.clip(axes_arr[:, 2], -1, 1)))
    for ki, k in enumerate(cluster_ids):
        mask    = orientation_labels == k
        polar_k = polar_all[mask]
        ax.hist(polar_k, bins=18, range=(0, 90), alpha=0.6,
                color=COLORS[ki % len(COLORS)], edgecolor='#0d1117',
                label=f'Cluster {k}')
    ax.set_xlabel('Polar angle from Z axis (°)')
    ax.set_ylabel('Count')
    ax.set_title('Fibril Tilt from Z\n(0°=vertical, 90°=horizontal)')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig('alignment_analysis.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved alignment_analysis.png")


# =========================
# PLOT ORIENTATION GROUP END-ON VIEWS
# =========================
def plot_orientation_group_endon(group_endon_results, group_labels, group_axes):
    """
    Plots the per-group end-on projections styled like IMOD:
    - Black background
    - Each fibril drawn as a circle (truly end-on = circular cross-section)
    - In-focus group = bright colour
    - Other group = grey, faded by depth (distance along viewing axis)
    - Global (all fibrils, mixed) shown first for comparison
    """
    GROUP_COLS = ['#E05C00', '#0070B8', '#15803D', '#9333EA', '#D97706']
    n_groups   = len(group_axes)
    n_panels   = n_groups + 1   # +1 for global

    fig, axes = plt.subplots(1, n_panels,
                              figsize=(6 * n_panels, 6.5),
                              facecolor='black')
    if n_panels == 1:
        axes = [axes]
    fig.suptitle('End-on projections — looking down fibril axis\n'
                 '(left = all fibrils mixed, right = per orientation group)',
                 color='white', fontsize=13, fontweight='bold')

    DIAM_NM = 8.0   # fibril diameter for circle radius

    def draw_panel(ax, px_all, py_all, depth_all, group_lbl_all,
                   focus_group, title):
        ax.set_facecolor('black')
        ax.set_aspect('equal')

        pts = np.column_stack([px_all, py_all])
        # depth-normalise for alpha
        d_min, d_max = depth_all.min(), depth_all.max()
        d_range      = d_max - d_min + 1e-9

        # Draw from back to front
        order = np.argsort(depth_all)
        for i in order:
            g = group_lbl_all[i]
            if focus_group == -99:          # global panel — all same colour
                col   = '#00FF41'
                alpha = 0.85
            elif g == focus_group:          # in-focus group
                col   = GROUP_COLS[g % len(GROUP_COLS)]
                alpha = 0.92
            else:                           # other group — depth-fade to grey
                d_norm = (depth_all[i] - d_min) / d_range
                alpha  = 0.15 + 0.30 * (1 - d_norm)
                col    = '#666666'
            circ = plt.Circle((pts[i, 0], pts[i, 1]),
                               DIAM_NM / 2,
                               color=col, alpha=alpha, zorder=3)
            ax.add_patch(circ)

        ax.set_xlim(pts[:, 0].min() - DIAM_NM * 2,
                    pts[:, 0].max() + DIAM_NM * 2)
        ax.set_ylim(pts[:, 1].min() - DIAM_NM * 2,
                    pts[:, 1].max() + DIAM_NM * 2)
        xl = ax.get_xlim(); yl = ax.get_ylim()
        # Scale bar
        sbx = xl[0] + (xl[1] - xl[0]) * 0.04
        sby = yl[0] + (yl[1] - yl[0]) * 0.04
        ax.plot([sbx, sbx + 100], [sby, sby],
                color='white', lw=2.5, zorder=9, solid_capstyle='butt')
        ax.text(sbx + 50, sby + (yl[1] - yl[0]) * 0.04,
                '100 nm', ha='center', fontsize=8,
                color='white', zorder=9, fontweight='bold')
        ax.set_xlabel('P1 (nm)', fontsize=9, color='#AAAAAA')
        ax.set_ylabel('P2 (nm)', fontsize=9, color='#AAAAAA')
        ax.tick_params(colors='#666666', labelsize=8)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(title, color='white', fontsize=9, pad=6)

    # Panel 0: global (all fibrils)
    g_res = group_endon_results.get('global', {})
    if g_res:
        draw_panel(axes[0],
                   np.array(g_res['proj_all_x']),
                   np.array(g_res['proj_all_y']),
                   np.array(g_res['depth_all']),
                   np.array(g_res['group_labels_all']),
                   focus_group=-99,
                   title=f"All fibrils mixed\n"
                         f"n={g_res['n']}  ψ₆={g_res['psi6_pp']:.3f}")

    # Per-group panels
    for g in range(n_groups):
        if g not in group_endon_results:
            continue
        res = group_endon_results[g]
        col = GROUP_COLS[g % len(GROUP_COLS)]
        # Coloured border for this group's panel
        for sp in axes[g + 1].spines.values():
            sp.set_visible(True)
            sp.set_color(col)
            sp.set_linewidth(2.5)
        draw_panel(axes[g + 1],
                   np.array(res['proj_all_x']),
                   np.array(res['proj_all_y']),
                   np.array(res['depth_all']),
                   np.array(res['group_labels_all']),
                   focus_group=g,
                   title=f"Viewing along Group {g} axis\n"
                         f"n={res['n']} in focus  ψ₆={res['psi6_pp']:.3f}"
                         f"  tilt={res['plane_angle']:.0f}°")

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig('orientation_group_endon.png', dpi=150,
                bbox_inches='tight', facecolor='black')
    print("Saved orientation_group_endon.png")


# =========================
# SUMMARY MATPLOTLIB FIGURE
# =========================
def make_summary_figure(fibrils, props, orientation_labels, bundle_labels,
                        cluster_data, nnd, lengths, vor, voronoi_counts, S2):

    psi6_results = cluster_data['global']['psi6']
    r_rdf        = cluster_data['global']['r']
    g_r          = cluster_data['global']['g_r']
    g_rand       = cluster_data['global']['g_rand']
    mid_r        = PSI6_RADII[1]
    psi6_per     = psi6_results[mid_r]['per_fibril']

    COLORS      = ['#00d4ff', '#ff6b35', '#7fff6b', '#c87fff', '#ffcc44']
    cmap_orient = cm.get_cmap('tab10')
    centers_ang = np.array([p['centroid_ang'] for p in props])

    fig = plt.figure(figsize=(20, 16), facecolor='#0d1117')
    fig.suptitle('Tau Fibril Organisation Analysis', color='white',
                 fontsize=18, fontweight='bold', y=0.98)

    def make_ax(row, col, colspan=1):
        ax = plt.subplot2grid((3, 4), (row, col), colspan=colspan, fig=fig)
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')
        return ax

    ax1 = plt.subplot2grid((3, 4), (0, 0), colspan=2, fig=fig)
    ax1.set_facecolor('#0d1117')
    for sp in ax1.spines.values(): sp.set_color('#30363d')
    ax1.tick_params(colors='#8b949e', labelsize=8)
    ax1.xaxis.label.set_color('#8b949e'); ax1.yaxis.label.set_color('#8b949e')
    ax1.title.set_color('#e6edf3')
    for i, f in enumerate(fibrils):
        c = cmap_orient(orientation_labels[i] % 10)
        ax1.plot(f[:, 0] * PIXEL_SIZE / 10, f[:, 1] * PIXEL_SIZE / 10,
                 color=c, lw=1.2, alpha=0.75)
    ax1.scatter(centers_ang[:, 0]/10, centers_ang[:, 1]/10,
                c=[cmap_orient(l % 10) for l in orientation_labels],
                s=15, zorder=5, alpha=0.9)
    ax1.set_xlabel('X (nm)'); ax1.set_ylabel('Y (nm)')
    ax1.set_title('Fibril Map — Orientation Clusters')
    ax1.set_aspect('equal')

    ax2 = plt.subplot2grid((3, 4), (0, 2), colspan=2, fig=fig)
    ax2.set_facecolor('#0d1117')
    for sp in ax2.spines.values(): sp.set_color('#30363d')
    ax2.tick_params(colors='#8b949e', labelsize=8)
    ax2.xaxis.label.set_color('#8b949e'); ax2.yaxis.label.set_color('#8b949e')
    ax2.title.set_color('#e6edf3')
    norm     = Normalize(vmin=0, vmax=1)
    cmap_psi = cm.get_cmap('plasma')
    for i, f in enumerate(fibrils):
        c = cmap_psi(norm(psi6_per[i]))
        ax2.plot(f[:, 0] * PIXEL_SIZE / 10, f[:, 1] * PIXEL_SIZE / 10,
                 color=c, lw=1.5, alpha=0.8)
    sc = ax2.scatter(centers_ang[:, 0]/10, centers_ang[:, 1]/10,
                     c=psi6_per, cmap='plasma', vmin=0, vmax=1, s=20, zorder=5)
    plt.colorbar(sc, ax=ax2, label='psi6', fraction=0.03)
    ax2.set_xlabel('X (nm)'); ax2.set_ylabel('Y (nm)')
    ax2.set_title('psi6 Hexagonal Order — Spatial Map')
    ax2.set_aspect('equal')

    ax3 = make_ax(1, 0)
    ax3.plot(r_rdf/10, g_r, color='#58a6ff', lw=2.5, label='Global', zorder=5)
    ax3.plot(r_rdf/10, g_rand, '--', color='#f85149', lw=1.5, alpha=0.6, label='Random')
    for i, (k, d) in enumerate(cluster_data.items()):
        if k == 'global': continue
        ax3.plot(d['r']/10, d['g_r'], color=COLORS[i % len(COLORS)],
                 lw=1.2, alpha=0.55, label=d['label'])
    ax3.axhline(1, color='#8b949e', lw=0.8, ls=':')
    ax3.set_xlabel('r (nm)'); ax3.set_ylabel('g(r) 2D')
    ax3.set_title('2D Radial Distribution Function')
    ax3.legend(fontsize=6, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax4 = make_ax(1, 1)
    ax4.hist(nnd, bins=25, color='#3fb950', edgecolor='#0d1117', alpha=0.8)
    ax4.axvline(np.median(nnd), color='white', lw=1.5, ls='--',
                label=f'Median: {np.median(nnd):.1f} nm')
    ax4.axvline(8, color='#f85149', lw=1.5, ls=':', label='Fibril diameter ~8 nm')
    ax4.set_xlabel('NND (nm)'); ax4.set_ylabel('Count')
    ax4.set_title('Nearest Neighbour Distances')
    ax4.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax5 = make_ax(1, 2)
    for i, (k, d) in enumerate(cluster_data.items()):
        vals = [d['psi6'][r]['mean'] for r in PSI6_RADII]
        lw = 2.5 if k == 'global' else 1.5
        ls = '--'  if k == 'global' else '-'
        ax5.plot([r/10 for r in PSI6_RADII], vals, 'o'+ls,
                 color=COLORS[i % len(COLORS)], lw=lw, ms=6, label=d['label'])
    ax5.axhline(0.6, color='#3fb950', lw=1, ls='--', alpha=0.5, label='Strong')
    ax5.axhline(0.3, color='#f85149', lw=1, ls='--', alpha=0.5, label='Weak')
    ax5.set_xlabel('Neighbourhood radius (nm)'); ax5.set_ylabel('psi6')
    ax5.set_ylim(0, 1); ax5.set_title('psi6 vs Radius per Cluster')
    ax5.legend(fontsize=6, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax6 = make_ax(1, 3)
    if voronoi_counts is not None and len(voronoi_counts) > 0:
        unique_coords, coord_counts = np.unique(voronoi_counts, return_counts=True)
        colors_bar = ['#f85149' if u == 6 else '#58a6ff' for u in unique_coords]
        ax6.bar(unique_coords, coord_counts, color=colors_bar,
                edgecolor='#0d1117', width=0.7)
        frac_6 = (voronoi_counts == 6).sum() / len(voronoi_counts)
        ax6.text(0.98, 0.95, f'6-fold: {frac_6:.1%}', transform=ax6.transAxes,
                 ha='right', va='top', color='#f85149', fontsize=9, fontweight='bold')
    ax6.set_xlabel('Voronoi neighbours'); ax6.set_ylabel('Count')
    ax6.set_title('Coordination Number\n(6 = hexagonal)')

    ax7 = make_ax(2, 0)
    ax7.hist(lengths, bins=25, color='#bc8cff', edgecolor='#0d1117', alpha=0.85)
    ax7.axvline(np.median(lengths), color='white', lw=1.5, ls='--',
                label=f'Median: {np.median(lengths):.0f} nm')
    ax7.set_xlabel('Length (nm)'); ax7.set_ylabel('Count')
    ax7.set_title('Fibril Length Distribution')
    ax7.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax8 = make_ax(2, 1)
    for i, (k, d) in enumerate(cluster_data.items()):
        ba = d['bond_angles']
        if len(ba) == 0: continue
        hist, edges = np.histogram(ba, bins=36, range=(0, 180))
        centres = 0.5 * (edges[1:] + edges[:-1])
        lw = 2.5 if k == 'global' else 1.5
        ls = '--'  if k == 'global' else '-'
        ax8.plot(centres, hist / (hist.max() + 1e-9),
                 color=COLORS[i % len(COLORS)], lw=lw, ls=ls,
                 label=d['label'], alpha=0.9)
    ax8.axvline(60,  color='white', lw=0.8, ls=':', alpha=0.6)
    ax8.axvline(90,  color='#888',  lw=0.8, ls=':', alpha=0.6)
    ax8.axvline(120, color='white', lw=0.8, ls=':', alpha=0.4)
    ax8.text(61, 0.92, '60°', color='white', fontsize=7)
    ax8.text(91, 0.92, '90°', color='#aaa',  fontsize=7)
    ax8.set_xlabel('Bond angle (°)'); ax8.set_ylabel('Norm. count')
    ax8.set_title('Bond Angle Distribution\n(60°=hex, 90°=square)')
    ax8.set_xlim(0, 180)
    ax8.legend(fontsize=6, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax9 = make_ax(2, 2)
    ax9.hist(psi6_per, bins=20, range=(0, 1), color='#da3633',
             edgecolor='#0d1117', alpha=0.85)
    ax9.axvline(float(np.mean(psi6_per)), color='white', lw=1.5, ls='--',
                label=f'Mean: {np.mean(psi6_per):.3f}')
    ax9.set_xlabel('psi6 per fibril'); ax9.set_ylabel('Count')
    ax9.set_title('psi6 Distribution')
    ax9.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax10 = make_ax(2, 3)
    ax10.axis('off')
    mean_psi  = float(psi6_results[mid_r]['mean'])
    psi4_val  = float(psi6_results.get('psi4', 0))
    frac_6c   = float((voronoi_counts == 6).sum() / len(voronoi_counts)) \
                if voronoi_counts is not None and len(voronoi_counts) > 0 else 0.0
    n_bundles = len(set(bundle_labels)) - (1 if -1 in bundle_labels else 0)
    dominant  = 'Hexagonal' if mean_psi > psi4_val else 'Square/other'
    hex_col   = '#3fb950' if mean_psi > 0.6 else '#ffa657' if mean_psi > 0.3 else '#f85149'
    s2_col    = '#3fb950' if S2 > 0.7 else '#ffa657' if S2 > 0.4 else '#f85149'

    rows = [
        ("SUMMARY",         "",                              '#e6edf3', 11, True),
        ("",                "",                              '#8b949e',  8, False),
        ("Fibrils",         f"{len(fibrils)}",               '#58a6ff',  9, False),
        ("psi6",            f"{mean_psi:.3f}",               '#ffa657',  9, False),
        ("psi4",            f"{psi4_val:.3f}",               '#bc8cff',  9, False),
        ("Dominant order",  dominant,                        hex_col,    9, False),
        ("S2 alignment",    f"{S2:.3f}",                     s2_col,     9, False),
        ("6-fold Voronoi",  f"{frac_6c:.1%}",                '#ffa657',  9, False),
        ("Mean NND",        f"{np.mean(nnd):.1f} nm",        '#3fb950',  9, False),
        ("Median length",   f"{np.median(lengths):.0f} nm",  '#bc8cff',  9, False),
        ("Bundles",         f"{n_bundles}",                  '#58a6ff',  9, False),
    ]
    y = 0.97
    for label, val, color, size, bold in rows:
        wt = 'bold' if bold else 'normal'
        if val:
            ax10.text(0.05, y, label, transform=ax10.transAxes,
                      color='#8b949e', fontsize=size-1, va='top')
            ax10.text(0.95, y, val, transform=ax10.transAxes,
                      color=color, fontsize=size, va='top', ha='right', fontweight=wt)
        else:
            ax10.text(0.5, y, label, transform=ax10.transAxes,
                      color=color, fontsize=size, va='top', ha='center', fontweight=wt)
        y -= 0.068

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig('fibril_analysis_summary.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved fibril_analysis_summary.png")


# =========================
# EXPORT JSON
# =========================
def export_json(fibrils, props, orientation_labels, bundle_labels,
                cluster_data, nnd, same_nnd, diff_nnd,
                lengths, axis_dists, voronoi_counts, voronoi_per_cluster, S2,
                psi6_perp=0.0, end_on_clusters=None, packing_plane_results=None,
                bundle_end_on=None, row_periodicity=None, row_chains=None,
                az=None, pol=None,
                group_labels=None, group_axes=None,
                group_endon_results=None):

    psi6_results = cluster_data['global']['psi6']
    mid_r        = PSI6_RADII[1]
    psi6_per     = psi6_results[mid_r]['per_fibril'].tolist()

    # Default arrays if not provided
    if az  is None: az  = [0.0] * len(fibrils)
    if pol is None: pol = [0.0] * len(fibrils)
    if group_labels is None: group_labels = orientation_labels

    cluster_export = {}
    for k, d in cluster_data.items():
        cluster_export[str(k)] = {
            'label': d['label'],
            'n': d['n'],
            'psi6_by_radius': {str(r): float(d['psi6'][r]['mean'])
                               for r in PSI6_RADII},
            'psi4':        float(d['psi6'].get('psi4', 0)),
            'psi2':        float(d['psi6'].get('psi2', 0)),
            'r_nm':        (d['r'] / 10).tolist(),
            'g_r':         d['g_r'].tolist(),
            'g_rand':      d['g_rand'].tolist(),
            'bond_angles': d['bond_angles'].tolist()
        }

    vc_list = voronoi_counts.tolist() if voronoi_counts is not None and len(voronoi_counts) > 0 else []
    vf6     = float((voronoi_counts == 6).sum() / len(voronoi_counts)) \
              if voronoi_counts is not None and len(voronoi_counts) > 0 else 0.0

    data = {
        "pixel_size":  PIXEL_SIZE,
        "n_fibrils":   len(fibrils),
        "fibrils": [
            {
                # 2D centroid — XY only (nm) — kept for backward compatibility
                "centroid_nm":         (props[i]['centroid_ang'][:2] / 10).tolist(),
                # *** NEW: full 3D centroid (nm) — needed for true end-on projection ***
                "centroid_3d_nm":      (props[i]['centroid_ang'] / 10).tolist(),
                "points_xy_nm":        (f[:, :2] * PIXEL_SIZE / 10).tolist(),
                "orientation_cluster": int(orientation_labels[i]),
                # *** NEW: discrete orientation group (-1 = isolated) ***
                "orientation_group":   int(group_labels[i]),
                "bundle":              int(bundle_labels[i]),
                "length_nm":           float(props[i]['length_ang'] / 10),
                "axis":                props[i]['axis'].tolist(),
                # *** NEW: azimuth and polar angles of fibril axis ***
                "azimuth":             float(az[i]),
                "polar":               float(pol[i]),
                "aspect_ratio":        float(props[i]['aspect_ratio']),
                "psi6":                float(psi6_per[i]) if i < len(psi6_per) else 0.0,
                "chain_label":         int(row_chains['chain_labels'][i])
                                       if row_chains and i < len(row_chains.get('chain_labels',[]))
                                       else -1
            }
            for i, f in enumerate(fibrils)
        ],
        "metrics": {
            "mean_psi6":              float(psi6_results[mid_r]['mean']),
            "psi4":                   float(psi6_results.get('psi4', 0)),
            "psi2":                   float(psi6_results.get('psi2', 0)),
            "psi6_by_radius":         {str(r): float(psi6_results[r]['mean'])
                                       for r in PSI6_RADII},
            "S2_nematic":             float(S2),
            "psi6_perp_director":     float(psi6_perp),
            "mean_nnd_nm":            float(np.mean(nnd)),
            "median_nnd_nm":          float(np.median(nnd)),
            "same_orient_nnd_nm":     float(np.mean(same_nnd)) if same_nnd else None,
            "diff_orient_nnd_nm":     float(np.mean(diff_nnd)) if diff_nnd else None,
            "mean_length_nm":         float(np.mean(lengths)),
            "median_length_nm":       float(np.median(lengths)),
            "mean_axis_axis_nm":      float(np.mean(axis_dists)),
            "n_bundles":              int(len(set(bundle_labels)) - (1 if -1 in bundle_labels else 0)),
            "n_isolated":             int((bundle_labels == -1).sum()),
            "n_orientation_clusters": N_ORIENT_CLUSTERS,
            "voronoi_frac_6":         vf6
        },
        "cluster_data":         cluster_export,
        "nnd_distribution":     {"all_nm": nnd.tolist()},
        "axis_dists_nm":        axis_dists.tolist(),
        "lengths_nm":           lengths,
        "voronoi_coordination": vc_list,
        "voronoi_per_cluster":  {
            str(k): {'frac_6': v['frac_6'], 'mean_coord': v['mean_coord'], 'n': v['n']}
            for k, v in voronoi_per_cluster.items()
        },
        "end_on_clusters": {
            str(k): {
                'label':                v['label'],
                'n':                    v['n'],
                'psi6_perp':            v['psi6_perp'],
                'psi6_perp_per_fibril': v['psi6_perp_per_fibril'],
                'mean_tilt':            v['mean_tilt'],
                'std_tilt':             v['std_tilt'],
                'azimuth_perp':         v['azimuth_perp'],
                'proj_x':               v['proj_x'],
                'proj_y':               v['proj_y'],
            }
            for k, v in (end_on_clusters or {}).items()
        },
        "end_on_global": {
            'proj_all_x':    end_on_clusters['global']['proj_all_x']   if end_on_clusters else [],
            'proj_all_y':    end_on_clusters['global']['proj_all_y']   if end_on_clusters else [],
            'orient_labels': end_on_clusters['global']['orient_labels'] if end_on_clusters else [],
        },
        "packing_plane": {
            str(k): {
                'label':       v['label'],
                'n':           v['n'],
                'n_excluded':  v.get('n_excluded', 0),
                'psi6_pp':     v['psi6_pp'],
                'psi6_pp_per': [0.0 if np.isnan(x) else x for x in v['psi6_pp_per']],
                'plane_angle': v['plane_angle'],
                'normal':      v['normal'],
                'proj_x':      v['proj_x'],
                'proj_y':      v['proj_y'],
            }
            for k, v in (packing_plane_results or {}).items()
        },
        "bundle_end_on": {
            str(k): v
            for k, v in (bundle_end_on or {}).items()
        },
        "row_periodicity": row_periodicity or {},
        "row_chains": {
            'chain_labels':  (row_chains or {}).get('chain_labels', []),
            'chain_lengths': (row_chains or {}).get('chain_lengths', []),
            'proj_x':        (row_chains or {}).get('proj_x', []),
            'proj_y':        (row_chains or {}).get('proj_y', []),
            'stats':         (row_chains or {}).get('stats', {}),
        },
        # *** NEW: discrete orientation group end-on projections ***
        # These are the true 3D projections looking along each group's axis
        # — equivalent to the IMOD rotation view
        "orientation_group_endon": {
            str(k): {
                'label':            v['label'],
                'n':                v['n'],
                'psi6_pp':          v['psi6_pp'],
                'psi6_pp_per':      v['psi6_pp_per'],
                'plane_angle':      v['plane_angle'],
                'normal':           v['normal'],
                'proj_x':           v['proj_x'],       # own-group fibrils
                'proj_y':           v['proj_y'],
                'depth':            v['depth'],
                'proj_all_x':       v['proj_all_x'],   # all fibrils in this view
                'proj_all_y':       v['proj_all_y'],
                'depth_all':        v['depth_all'],
                'group_labels_all': v['group_labels_all'],
            }
            for k, v in (group_endon_results or {}).items()
        },
        # *** NEW: orientation group metadata ***
        "orientation_groups": {
            str(g): {
                'mean_axis': ga.tolist(),
                'n':         int((np.array(group_labels) == g).sum()),
            }
            for g, ga in enumerate(group_axes or [])
        },
    }

    class NoNaNEncoder(json.JSONEncoder):
        def iterencode(self, o, _one_shot=False):
            chunks = super().iterencode(o, _one_shot)
            for chunk in chunks:
                chunk = chunk.replace('NaN', '0.0')
                chunk = chunk.replace('Infinity', '0.0')
                chunk = chunk.replace('-Infinity', '0.0')
                yield chunk

    with open("fibril_data.json", "w") as f:
        json.dump(data, f, cls=NoNaNEncoder)
    print("Exported fibril_data.json")
    return data


# =========================
# PACKING PLANE ANALYSIS
# =========================
def find_packing_plane(props, orientation_labels=None):
    centers  = np.array([p['centroid_ang'] for p in props])
    centered = centers - centers.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)

    p1     = vh[0]
    p2     = vh[1]
    normal = vh[2]

    if normal[2] < 0:
        normal = -normal

    proj = np.column_stack([centered @ p1, centered @ p2]) / 10

    z           = np.array([0, 0, 1])
    plane_angle = float(np.degrees(np.arccos(np.clip(abs(np.dot(normal, z)), 0, 1))))

    tree     = cKDTree(proj)
    nnd_vals = [tree.query(proj[i], k=2)[0][1] for i in range(len(proj))]
    search_r = float(np.mean(nnd_vals)) * NND_SEARCH_FACTOR

    psi6_vals = []
    psi6_per  = []
    for i, c in enumerate(proj):
        idx        = tree.query_ball_point(c, r=search_r)
        neighbours = [j for j in idx if j != i]
        if len(neighbours) < 3:
            psi6_per.append(np.nan)
            continue
        angles = [np.arctan2((proj[j] - c)[1], (proj[j] - c)[0]) for j in neighbours]
        psi = float(np.abs(np.mean(np.exp(6j * np.array(angles)))))
        psi6_vals.append(psi)
        psi6_per.append(psi)

    psi6_pp = float(np.nanmean(psi6_per)) if psi6_per else 0.0

    return p1, p2, normal, proj, psi6_pp, psi6_per, plane_angle


def packing_plane_per_cluster(props, orientation_labels):
    centers     = np.array([p['centroid_ang'] for p in props])
    cluster_ids = sorted(set(orientation_labels))
    results     = {}

    print(f"\n--- PACKING PLANE psi6 PER CLUSTER ---")
    print(f"  (Outlier cutoff: {OUTLIER_SIGMA}x std  |  NND factor: {NND_SEARCH_FACTOR})")
    print(f"  {'Population':<20} {'N':>5} {'Excl':>5}  {'psi6_pp':>10}  "
          f"{'Plane->Z':>9}  Assessment")
    print(f"  {'-'*20} {'-'*5} {'-'*5}  {'-'*10}  {'-'*9}  {'-'*15}")

    def run_cluster(mask):
        sub_centers = centers[mask]
        sub_props   = [props[i] for i in np.where(mask)[0]]
        if len(sub_centers) < 4:
            return None
        centre    = sub_centers.mean(axis=0)
        dists     = np.linalg.norm(sub_centers - centre, axis=1)
        std_dist  = dists.std()
        keep_mask = dists <= (dists.mean() + OUTLIER_SIGMA * std_dist)
        n_excl    = (~keep_mask).sum()
        sub_props_clean = [sub_props[i] for i in np.where(keep_mask)[0]]
        if len(sub_props_clean) < 4:
            return None
        p1, p2, normal, proj, psi6_pp, psi6_per, plane_angle = \
            find_packing_plane(sub_props_clean)
        return {
            'n':           int(keep_mask.sum()),
            'n_excluded':  int(n_excl),
            'psi6_pp':     psi6_pp,
            'psi6_pp_per': psi6_per,
            'plane_angle': plane_angle,
            'normal':      normal.tolist(),
            'p1':          p1.tolist(),
            'p2':          p2.tolist(),
            'proj_x':      proj[:, 0].tolist(),
            'proj_y':      proj[:, 1].tolist(),
        }

    all_mask = np.ones(len(props), dtype=bool)
    g = run_cluster(all_mask)
    if g:
        g['label'] = 'Global'
        results['global'] = g
        print(f"  {'ALL (global)':<20} {g['n']:>5} {g['n_excluded']:>5}  "
              f"{g['psi6_pp']:>10.4f}  {g['plane_angle']:>8.1f}  "
              f"{'Strong' if g['psi6_pp']>0.6 else 'Moderate' if g['psi6_pp']>0.3 else 'Weak'}")

    for k in cluster_ids:
        mask = orientation_labels == k
        r    = run_cluster(mask)
        if r is None:
            print(f"  {'Cluster '+str(k):<20} {mask.sum():>5}     -  too few after cleaning")
            continue
        r['label'] = f'Cluster {k}'
        results[k] = r
        print(f"  {'Cluster '+str(k):<20} {r['n']:>5} {r['n_excluded']:>5}  "
              f"{r['psi6_pp']:>10.4f}  {r['plane_angle']:>8.1f}  "
              f"{'Strong' if r['psi6_pp']>0.6 else 'Moderate' if r['psi6_pp']>0.3 else 'Weak'}")

    return results


def bundle_end_on_analysis(props, bundle_labels, director):
    centers     = np.array([p['centroid_ang'] for p in props])
    bundle_ids  = sorted([b for b in set(bundle_labels) if b != -1])

    d_unit = director / np.linalg.norm(director)
    ref    = np.array([1,0,0]) if abs(d_unit[0]) < 0.9 else np.array([0,1,0])
    e1     = np.cross(d_unit, ref); e1 /= np.linalg.norm(e1)
    e2     = np.cross(d_unit, e1);  e2 /= np.linalg.norm(e2)

    proj_all = np.column_stack([
        np.dot(centers - centers.mean(axis=0), e1),
        np.dot(centers - centers.mean(axis=0), e2)
    ]) / 10

    results = {}
    print(f"\n--- END-ON ψ6 PER BUNDLE ---")
    print(f"  {'Bundle':<12} {'N':>5}  {'ψ6_pp':>8}  {'plane->Z':>8}")
    print(f"  {'-'*12} {'-'*5}  {'-'*8}  {'-'*8}")

    for bid in bundle_ids:
        mask     = bundle_labels == bid
        n        = mask.sum()
        if n < 4:
            continue
        sub_props = [props[i] for i in np.where(mask)[0]]

        p1, p2, normal, proj_b, psi6_pp, psi6_pp_per, plane_angle = \
            find_packing_plane(sub_props)

        fidxs           = np.where(mask)[0]
        sub_proj_global = proj_all[mask]

        results[int(bid)] = {
            'label':         f'Bundle {bid}',
            'n':             int(n),
            'psi6_pp':       psi6_pp,
            'psi6_pp_per':   [0.0 if np.isnan(x) else x for x in psi6_pp_per],
            'plane_angle':   plane_angle,
            'normal':        normal.tolist(),
            'proj_x':        proj_b[:, 0].tolist(),
            'proj_y':        proj_b[:, 1].tolist(),
            'proj_global_x': sub_proj_global[:, 0].tolist(),
            'proj_global_y': sub_proj_global[:, 1].tolist(),
            'fibril_indices': fidxs.tolist(),
        }
        print(f"  {'Bundle '+str(bid):<12} {n:>5}  {psi6_pp:>8.4f}  "
              f"{plane_angle:>7.1f}°")

    results['_global_proj'] = {
        'proj_x': proj_all[:, 0].tolist(),
        'proj_y': proj_all[:, 1].tolist(),
        'bundle_labels': bundle_labels.tolist(),
    }
    return results


# =========================
# 3D HEXAGONAL COLUMNAR ANALYSIS
# =========================
def hexagonal_3d_analysis(props, director):
    centers = np.array([p['centroid_ang'] for p in props])
    tree    = cKDTree(centers)

    d = director / np.linalg.norm(director)

    nnd_vals = []
    for i in range(len(centers)):
        dists, _ = tree.query(centers[i], k=2)
        nnd_vals.append(dists[1])
    search_r = float(np.mean(nnd_vals)) * 2.0

    ref = np.array([1, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1, 0])
    e1  = np.cross(d, ref); e1 /= np.linalg.norm(e1)
    e2  = np.cross(d, e1);  e2 /= np.linalg.norm(e2)

    tilt_angles  = []
    azimuth_perp = []
    psi6_perp_per_fibril = []

    for i, c in enumerate(centers):
        idx        = tree.query_ball_point(c, r=search_r)
        neighbours = [j for j in idx if j != i]

        perp_angles = []
        for j in neighbours:
            v   = centers[j] - c
            v_n = v / (np.linalg.norm(v) + 1e-9)

            cos_tilt = abs(np.dot(v_n, d))
            tilt_angles.append(np.degrees(np.arccos(np.clip(cos_tilt, 0, 1))))

            v_perp = v - np.dot(v, d) * d
            if np.linalg.norm(v_perp) > 1e-6:
                az = np.arctan2(np.dot(v_perp, e2), np.dot(v_perp, e1))
                azimuth_perp.append(az)
                perp_angles.append(az)

        if len(perp_angles) >= 2:
            psi6_perp_per_fibril.append(
                float(np.abs(np.mean(np.exp(6j * np.array(perp_angles)))))
            )
        else:
            psi6_perp_per_fibril.append(0.0)

    tilt_angles  = np.array(tilt_angles)
    azimuth_perp = np.array(azimuth_perp)
    psi6_perp    = float(np.mean(psi6_perp_per_fibril)) if psi6_perp_per_fibril else 0.0

    print(f"\n--- 3D HEXAGONAL COLUMNAR ANALYSIS ---")
    print(f"  Search radius used: {search_r/10:.1f} nm")
    print(f"  Mean tilt of inter-fibril vectors from director: {np.mean(tilt_angles):.1f}°")
    print(f"  ψ6 in plane perpendicular to director: {psi6_perp:.4f}")

    return tilt_angles, azimuth_perp, psi6_perp, psi6_perp_per_fibril


def hexagonal_3d_per_cluster(props, orientation_labels, director):
    centers = np.array([p['centroid_ang'] for p in props])

    d   = director / np.linalg.norm(director)
    ref = np.array([1, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1, 0])
    e1  = np.cross(d, ref); e1 /= np.linalg.norm(e1)
    e2  = np.cross(d, e1);  e2 /= np.linalg.norm(e2)

    proj_all = np.column_stack([
        np.dot(centers - centers.mean(axis=0), e1),
        np.dot(centers - centers.mean(axis=0), e2)
    ]) / 10

    def run_for_subset(idx_mask):
        sub_centers = centers[idx_mask]
        if len(sub_centers) < 3:
            return None

        centre     = sub_centers.mean(axis=0)
        dists_c    = np.linalg.norm(sub_centers - centre, axis=1)
        keep       = dists_c <= (dists_c.mean() + OUTLIER_SIGMA * dists_c.std())
        n_excl     = (~keep).sum()
        sub_centers = sub_centers[keep]
        if len(sub_centers) < 3:
            return None

        tree     = cKDTree(sub_centers)
        nnd_vals = [tree.query(sub_centers[i], k=2)[0][1]
                    for i in range(len(sub_centers))]
        search_r = float(np.mean(nnd_vals)) * NND_SEARCH_FACTOR

        tilt_angles  = []
        azimuth_perp = []
        psi6_per     = []

        for i, c in enumerate(sub_centers):
            idx        = tree.query_ball_point(c, r=search_r)
            neighbours = [j for j in idx if j != i]
            perp_a     = []
            for j in neighbours:
                v   = sub_centers[j] - c
                v_n = v / (np.linalg.norm(v) + 1e-9)
                cos_tilt = abs(np.dot(v_n, d))
                tilt_angles.append(np.degrees(np.arccos(np.clip(cos_tilt, 0, 1))))
                v_perp = v - np.dot(v, d) * d
                if np.linalg.norm(v_perp) > 1e-6:
                    az = np.arctan2(np.dot(v_perp, e2), np.dot(v_perp, e1))
                    azimuth_perp.append(az)
                    perp_a.append(az)

            if len(perp_a) >= 3:
                psi6_per.append(float(np.abs(np.mean(np.exp(6j * np.array(perp_a))))))
            else:
                psi6_per.append(np.nan)

        tilt_arr = np.array(tilt_angles)
        az_arr   = np.array(azimuth_perp)
        valid    = [p for p in psi6_per if not np.isnan(p)]
        psi6v    = float(np.mean(valid)) if valid else 0.0

        sub_proj = np.column_stack([
            np.dot(sub_centers - centers.mean(axis=0), e1),
            np.dot(sub_centers - centers.mean(axis=0), e2)
        ]) / 10

        return {
            'psi6_perp':            psi6v,
            'psi6_perp_per_fibril': [0.0 if np.isnan(p) else p for p in psi6_per],
            'n_valid_psi6':         len(valid),
            'n_excluded':           int(n_excl),
            'mean_tilt':            float(np.mean(tilt_arr)) if len(tilt_arr) else 0.0,
            'std_tilt':             float(np.std(tilt_arr))  if len(tilt_arr) else 0.0,
            'azimuth_perp':         az_arr.tolist(),
            'proj_x':               sub_proj[:, 0].tolist(),
            'proj_y':               sub_proj[:, 1].tolist(),
            'n':                    int(keep.sum())
        }

    cluster_ids = sorted(set(orientation_labels))
    results     = {}

    print(f"\n--- END-ON ψ6 PER CLUSTER ---")
    print(f"  {'Population':<20} {'N':>5}  {'ψ6_perp':>10}  {'Mean tilt':>10}  Assessment")
    print(f"  {'-'*20} {'-'*5}  {'-'*10}  {'-'*10}  {'-'*20}")

    all_mask = np.ones(len(props), dtype=bool)
    g        = run_for_subset(all_mask)
    g['label']    = 'Global'
    g['proj_all_x'] = proj_all[:, 0].tolist()
    g['proj_all_y'] = proj_all[:, 1].tolist()
    g['orient_labels'] = orientation_labels.tolist()
    results['global'] = g
    print(f"  {'ALL (global)':<20} {g['n']:>5}  "
          f"{g['psi6_perp']:>10.4f}  {g['mean_tilt']:>9.1f}°  "
          f"{'Strong' if g['psi6_perp']>0.6 else 'Moderate' if g['psi6_perp']>0.3 else 'Weak'}")

    for k in cluster_ids:
        mask = orientation_labels == k
        if mask.sum() < 3:
            continue
        r = run_for_subset(mask)
        r['label'] = f'Cluster {k}'
        results[k] = r
        print(f"  {'Cluster '+str(k):<20} {r['n']:>5}  "
              f"{r['psi6_perp']:>10.4f}  {r['mean_tilt']:>9.1f}°  "
              f"{'Strong' if r['psi6_perp']>0.6 else 'Moderate' if r['psi6_perp']>0.3 else 'Weak'}")

    return results


# =========================
# PLOT 3D HEX ANALYSIS
# =========================
def plot_3d_hex(tilt_angles, azimuth_perp, psi6_perp, psi6_perp_per_fibril,
                psi6_xy, director):
    fig, axes = plt.subplots(1, 4, figsize=(22, 5), facecolor='#0d1117')
    fig.suptitle('3D Hexagonal Columnar Analysis (end-on view)',
                 color='white', fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')

    ax = axes[0]; style(ax)
    ax.hist(tilt_angles, bins=36, range=(0, 90),
            color='#00d4ff', edgecolor='#0d1117', alpha=0.8)
    ax.axvline(90, color='#7fff6b', lw=1.5, ls='--', alpha=0.7,
               label='90° = ⊥ to director')
    ax.axvline(float(np.mean(tilt_angles)), color='white', lw=1.5, ls='--',
               label=f'Mean: {np.mean(tilt_angles):.1f}°')
    ax.set_xlabel('Tilt angle from director (°)')
    ax.set_ylabel('Count')
    ax.set_xlim(0, 90)
    ax.set_title('Inter-fibril vector tilt\n(90° = columnar packing)')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    axes[1].remove()
    ax = fig.add_subplot(1, 4, 2, projection='polar')
    ax.set_facecolor('#161b22')
    ax.tick_params(colors='#8b949e', labelsize=7)
    ax.title.set_color('#e6edf3')
    hist_az, edges = np.histogram(azimuth_perp, bins=36)
    theta = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]
    ax.bar(theta, hist_az, width=width, bottom=0,
           color='#ffa657', alpha=0.8, edgecolor='#0d1117')
    for angle in np.linspace(0, 2 * np.pi, 7)[:-1]:
        ax.axvline(angle, color='#7fff6b', lw=0.8, alpha=0.4, ls='--')
    ax.set_title(f'Azimuthal distribution\naround director\nψ6_perp = {psi6_perp:.3f}', pad=15)

    ax = axes[2]; style(ax)
    labels = ['ψ6\n(XY plane)', 'ψ6\n(⊥ director)']
    vals   = [psi6_xy, psi6_perp]
    colors = ['#58a6ff', '#7fff6b']
    bars   = ax.bar(labels, vals, color=colors, edgecolor='#0d1117', width=0.5, alpha=0.85)
    ax.axhline(0.6, color='#7fff6b', lw=1, ls='--', alpha=0.4, label='Strong (0.6)')
    ax.axhline(0.3, color='#ff6b35', lw=1, ls='--', alpha=0.4, label='Weak (0.3)')
    ax.set_ylim(0, 1); ax.set_ylabel('ψ6')
    ax.set_title('ψ6 comparison:\nXY projection vs end-on view')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.02,
                f'{val:.3f}', ha='center', color='white', fontsize=11, fontweight='bold')

    ax = axes[3]; style(ax)
    ax.hist(psi6_perp_per_fibril, bins=20, range=(0, 1),
            color='#ffa657', edgecolor='#0d1117', alpha=0.85)
    ax.axvline(psi6_perp, color='white', lw=1.5, ls='--',
               label=f'Mean: {psi6_perp:.3f}')
    ax.set_xlabel('ψ6_perp per fibril'); ax.set_ylabel('Count')
    ax.set_title('Per-fibril ψ6 distribution\n(end-on view)')
    ax.legend(fontsize=7, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    plt.tight_layout(rect=[0, 0.04, 1, 0.94])
    plt.savefig('hexagonal_3d_analysis.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved hexagonal_3d_analysis.png")


# =========================
# END-ON VIEW VISUALISATION
# =========================
def plot_end_on_view(props, director, orientation_labels, psi6_perp_per_fibril,
                     pixel_size=PIXEL_SIZE):
    centers = np.array([p['centroid_ang'] for p in props])

    d   = director / np.linalg.norm(director)
    ref = np.array([1, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1, 0])
    e1  = np.cross(d, ref); e1 /= np.linalg.norm(e1)
    e2  = np.cross(d, e1);  e2 /= np.linalg.norm(e2)

    proj = np.column_stack([
        np.dot(centers - centers.mean(axis=0), e1),
        np.dot(centers - centers.mean(axis=0), e2)
    ]) / 10

    fig, axes = plt.subplots(1, 3, figsize=(21, 7), facecolor='#0d1117')
    fig.suptitle('End-on View — Looking Down the Fibril Axis',
                 color='white', fontsize=15, fontweight='bold')

    ORIENT_COLORS = ['#00d4ff', '#ff6b35', '#7fff6b', '#c87fff', '#ffcc44']
    cmap_psi      = cm.get_cmap('plasma')
    norm_psi      = Normalize(vmin=0, vmax=1)

    def style(ax, xlabel='e1 (nm)', ylabel='e2 (nm)'):
        ax.set_facecolor('#0d1117')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.set_xlabel(xlabel, color='#8b949e')
        ax.set_ylabel(ylabel, color='#8b949e')
        ax.set_aspect('equal')
        ax.grid(True, color='#1a2535', linewidth=0.5, alpha=0.5)

    ax = axes[0]; style(ax)
    ax.set_title('Orientation clusters + Voronoi', color='#e6edf3', pad=10)

    if len(proj) >= 4:
        try:
            from scipy.spatial import Voronoi as _Vor
            vor = _Vor(proj)
            for ridge_pts, ridge_verts in zip(vor.ridge_points, vor.ridge_vertices):
                if -1 not in ridge_verts:
                    vx = vor.vertices[ridge_verts, 0]
                    vy = vor.vertices[ridge_verts, 1]
                    bound = np.abs(proj).max() * 1.1
                    if np.all(np.abs(vx) < bound) and np.all(np.abs(vy) < bound):
                        ax.plot(vx, vy, '-', color='#1a2535', lw=0.8, alpha=0.7, zorder=1)
        except Exception:
            pass

    for k in sorted(set(orientation_labels)):
        mask = orientation_labels == k
        col  = ORIENT_COLORS[k % len(ORIENT_COLORS)]
        ax.scatter(proj[mask, 0], proj[mask, 1],
                   c=col, s=80, zorder=3, edgecolors='#070b12',
                   linewidths=0.8, label=f'Cluster {k} (n={mask.sum()})')

    ax.legend(fontsize=8, labelcolor='#8b949e',
              facecolor='#0d1117', edgecolor='#30363d', loc='upper right')

    ax = axes[1]; style(ax)
    ax.set_title('ψ6 end-on (looking down fibril axis)', color='#e6edf3', pad=10)

    psi6_arr = np.array(psi6_perp_per_fibril)
    sc = ax.scatter(proj[:, 0], proj[:, 1],
                    c=psi6_arr, cmap='plasma', vmin=0, vmax=1,
                    s=120, zorder=3, edgecolors='#070b12', linewidths=0.8)
    plt.colorbar(sc, ax=ax, label='ψ6_perp', fraction=0.04, pad=0.02)

    tree = cKDTree(proj)
    for i in range(len(proj)):
        dists, idxs = tree.query(proj[i], k=4)
        for dist, j in zip(dists[1:], idxs[1:]):
            if i < j:
                col = cmap_psi(norm_psi((psi6_arr[i] + psi6_arr[j]) / 2))
                ax.plot([proj[i, 0], proj[j, 0]],
                        [proj[i, 1], proj[j, 1]],
                        '-', color=col, lw=0.8, alpha=0.4, zorder=2)

    ax.text(0.02, 0.97, f'Mean ψ6_perp = {psi6_arr.mean():.3f}',
            transform=ax.transAxes, color='#ffcc44', fontsize=9,
            va='top', fontweight='bold')

    ax = axes[2]; style(ax)
    ax.set_title('NN bond network\n(colour = bond length)',
                 color='#e6edf3', pad=10)

    bond_lengths = []
    for i in range(len(proj)):
        dists, idxs = tree.query(proj[i], k=4)
        for dist, j in zip(dists[1:], idxs[1:]):
            if i < j:
                bond_lengths.append(dist)
    bond_lengths = np.array(bond_lengths)
    bmin, bmax = bond_lengths.min(), bond_lengths.max()
    cmap_bond = cm.get_cmap('RdYlGn_r')
    norm_bond = Normalize(vmin=bmin, vmax=bmax)

    for i in range(len(proj)):
        dists, idxs = tree.query(proj[i], k=4)
        for dist, j in zip(dists[1:], idxs[1:]):
            if i < j:
                col = cmap_bond(norm_bond(dist))
                ax.plot([proj[i, 0], proj[j, 0]],
                        [proj[i, 1], proj[j, 1]],
                        '-', color=col, lw=1.5, alpha=0.7, zorder=2)

    ax.scatter(proj[:, 0], proj[:, 1],
               c='white', s=60, zorder=3, edgecolors='#070b12', linewidths=0.8)

    sm = cm.ScalarMappable(cmap='RdYlGn_r', norm=norm_bond)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label('Bond length (nm)', color='#8b949e')
    cb.ax.yaxis.set_tick_params(color='#8b949e')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='#8b949e')

    ax.text(0.02, 0.97,
            f'Bond length: {bond_lengths.mean():.1f} ± {bond_lengths.std():.1f} nm\n'
            f'CV = {bond_lengths.std()/bond_lengths.mean():.3f}',
            transform=ax.transAxes, color='#7fff6b', fontsize=8,
            va='top', fontweight='bold')

    d_str = f'Director: [{d[0]:.2f}, {d[1]:.2f}, {d[2]:.2f}]'
    fig.text(0.5, 0.01, d_str, ha='center', color='#4a6080', fontsize=9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('end_on_view.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved end_on_view.png")


def measure_row_periodicity(props, packing_plane_results, orientation_labels):
    from scipy.signal import find_peaks

    pp_global = packing_plane_results.get('global', {})
    if not pp_global or not pp_global.get('proj_x'):
        print("\n--- ROW PERIODICITY: no packing plane data ---")
        return {}

    proj = np.column_stack([pp_global['proj_x'], pp_global['proj_y']])

    print(f"\n--- ROW PERIODICITY ANALYSIS ---")

    res   = 3.0
    xmin, xmax = proj[:,0].min()-10, proj[:,0].max()+10
    ymin, ymax = proj[:,1].min()-10, proj[:,1].max()+10
    nx = max(32, int((xmax-xmin)/res))
    ny = max(32, int((ymax-ymin)/res))

    density = np.zeros((ny, nx))
    for x, y in proj:
        ix = int((x-xmin)/res); iy = int((y-ymin)/res)
        ix = np.clip(ix, 0, nx-1); iy = np.clip(iy, 0, ny-1)
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                if 0<=iy+dy<ny and 0<=ix+dx<nx:
                    density[iy+dy, ix+dx] += np.exp(-(dx**2+dy**2)/2)

    fft2   = np.fft.fft2(density - density.mean())
    power  = np.abs(np.fft.fftshift(fft2))**2

    cy, cx = ny//2, nx//2
    mask   = np.ones_like(power, dtype=bool)
    for dy in range(-3, 4):
        for dx in range(-3, 4):
            if 0<=cy+dy<ny and 0<=cx+dx<nx:
                mask[cy+dy, cx+dx] = False
    power_masked        = power.copy(); power_masked[~mask] = 0
    peak_idx            = np.unravel_index(power_masked.argmax(), power.shape)
    peak_ky             = (peak_idx[0] - cy) / (ny * res)
    peak_kx             = (peak_idx[1] - cx) / (nx * res)

    k_mag    = np.sqrt(peak_kx**2 + peak_ky**2)
    if k_mag < 1e-9:
        print("  No dominant periodicity found in 2D FFT")
        return {}

    spacing_fft   = 1.0 / k_mag
    wavevec_angle = np.degrees(np.arctan2(peak_ky, peak_kx))
    row_direction = wavevec_angle + 90

    print(f"  Row spacing (FFT): {spacing_fft:.1f} nm")

    k_unit    = np.array([peak_kx, peak_ky]) / k_mag
    pos_1d    = proj @ k_unit

    bins_1d   = np.linspace(pos_1d.min()-5, pos_1d.max()+5, 500)
    dr        = bins_1d[1] - bins_1d[0]
    hist_1d, edges = np.histogram(pos_1d, bins=bins_1d)
    r_1d      = 0.5*(edges[1:]+edges[:-1])

    from scipy.ndimage import gaussian_filter1d
    hist_sm   = gaussian_filter1d(hist_1d.astype(float), sigma=2)

    h_cent    = hist_sm - hist_sm.mean()
    acf       = np.correlate(h_cent, h_cent, mode='full')
    acf       = acf[len(acf)//2:]
    lags      = dr * np.arange(len(acf))

    min_lag   = 5
    min_idx   = int(min_lag / dr)
    peaks, _  = find_peaks(acf[min_idx:], height=acf[0]*0.05, distance=int(10/dr))

    if len(peaks):
        row_spacing_acf = lags[peaks[0] + min_idx]
        print(f"  Row spacing (autocorrelation): {row_spacing_acf:.1f} nm")
    else:
        row_spacing_acf = spacing_fft
        print(f"  Autocorrelation: no clear peak, using FFT value")

    result = {
        'row_direction_deg':  float(row_direction),
        'row_spacing_fft_nm': float(spacing_fft),
        'row_spacing_acf_nm': float(row_spacing_acf),
        'n_rows':             int((pos_1d.max()-pos_1d.min())/row_spacing_acf) if len(peaks) else 0,
        'pos_1d':             pos_1d.tolist(),
        'acf':                acf[:200].tolist(),
        'acf_lags':           lags[:200].tolist(),
        'r_1d':               r_1d.tolist(),
        'hist_1d':            hist_sm.tolist(),
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), facecolor='#0d1117')
    fig.suptitle('Row Periodicity Analysis — End-On View',
                 color='white', fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e'); ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')

    ax = axes[0]; style(ax); ax.set_aspect('equal'); ax.set_facecolor('#0d1117')
    psi6_arr = np.array(pp_global.get('psi6_pp_per', [0]*len(proj)))
    sc = ax.scatter(proj[:,0], proj[:,1], c=psi6_arr, cmap='plasma',
                    vmin=0, vmax=1, s=80, edgecolors='#070b12', linewidths=0.8, zorder=3)
    plt.colorbar(sc, ax=ax, label='ψ6', fraction=0.04)
    row_vec   = np.array([np.cos(np.radians(row_direction)), np.sin(np.radians(row_direction))])
    n_rows_draw = int((pos_1d.max()-pos_1d.min()) / row_spacing_acf) + 2
    for i in range(n_rows_draw + 1):
        offset    = pos_1d.min() + i * row_spacing_acf
        centre_pt = k_unit * offset
        t         = np.linspace(-300, 300, 2)
        ax.plot(centre_pt[0] + t * row_vec[0], centre_pt[1] + t * row_vec[1],
                '-', color='#3fb950', lw=0.8, alpha=0.5, zorder=2)
    ax.set_xlabel('P1 (nm)'); ax.set_ylabel('P2 (nm)')
    ax.set_title(f'Packing plane + detected rows\nrow spacing = {row_spacing_acf:.1f} nm')

    ax = axes[1]; style(ax)
    ax.plot(r_1d, hist_sm, color='#58a6ff', lw=1.5)
    ax.fill_between(r_1d, hist_sm, alpha=0.2, color='#58a6ff')
    if len(peaks):
        ax.axvline(lags[peaks[0]+min_idx], color='#3fb950', lw=1.5, ls='--',
                   label=f'Row spacing: {row_spacing_acf:.1f} nm')
    ax.set_xlabel('Position along wavevector (nm)'); ax.set_ylabel('Fibril density')
    ax.set_title('1D fibril density profile')
    ax.legend(fontsize=8, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax = axes[2]; style(ax)
    ax.plot(lags[:200], acf[:200]/acf[0], color='#ffa657', lw=1.5)
    if len(peaks):
        ax.axvline(row_spacing_acf, color='#3fb950', lw=1.5, ls='--',
                   label=f'First peak: {row_spacing_acf:.1f} nm')
        exp_spacing = 42.6 * np.sqrt(3) / 2
        ax.axvline(exp_spacing, color='#ff6b35', lw=1, ls=':',
                   label=f'Expected hex: {exp_spacing:.1f} nm')
    ax.set_xlabel('Lag (nm)'); ax.set_ylabel('Normalised ACF')
    ax.set_title('Autocorrelation of 1D density')
    ax.legend(fontsize=8, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')
    ax.set_xlim(0, 200)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig('row_periodicity.png', dpi=150, bbox_inches='tight', facecolor='#0d1117')
    print("  Saved row_periodicity.png")

    return result


# =========================
# ROW CHAIN DETECTION
# =========================
def detect_row_chains(props, packing_plane_results,
                      collinearity_thresh=168.0,
                      search_factor=1.4):
    pp = packing_plane_results.get('global', {})
    if not pp or not pp.get('proj_x'):
        print("\n--- ROW CHAIN DETECTION: no packing plane data ---")
        return {}

    proj = np.column_stack([pp['proj_x'], pp['proj_y']])
    n    = len(proj)

    tree     = cKDTree(proj)
    nnd_vals = [tree.query(proj[i], k=2)[0][1] for i in range(n)]
    search_r = float(np.mean(nnd_vals)) * search_factor

    print(f"\n--- ROW CHAIN DETECTION ---")
    print(f"  Collinearity threshold: {collinearity_thresh}°")
    print(f"  Search radius: {search_r:.1f} nm  ({search_factor}× mean NND)")

    neighbours = {}
    for i in range(n):
        idx = tree.query_ball_point(proj[i], r=search_r)
        neighbours[i] = [j for j in idx if j != i]

    row_adj = {i: set() for i in range(n)}

    for k in range(n):
        nb = neighbours[k]
        if len(nb) < 2:
            continue
        vecs = {j: (proj[j] - proj[k]) for j in nb}
        nb_list = list(nb)
        for a in range(len(nb_list)):
            for b in range(a+1, len(nb_list)):
                i, j = nb_list[a], nb_list[b]
                vi = vecs[i] / (np.linalg.norm(vecs[i]) + 1e-9)
                vj = vecs[j] / (np.linalg.norm(vecs[j]) + 1e-9)
                cos_angle = np.dot(vi, vj)
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                if angle >= collinearity_thresh:
                    row_adj[i].add(j)
                    row_adj[j].add(i)

    visited = np.zeros(n, dtype=bool)
    chains  = []
    for start in range(n):
        if visited[start]:
            continue
        if not row_adj[start]:
            visited[start] = True
            continue
        chain   = []
        queue   = [start]
        visited[start] = True
        while queue:
            node = queue.pop(0)
            chain.append(node)
            for nb in row_adj[node]:
                if not visited[nb]:
                    visited[nb] = True
                    queue.append(nb)
        if len(chain) >= 2:
            chains.append(sorted(chain))

    chain_set = {}
    for ci, chain in enumerate(chains):
        for fi in chain:
            chain_set[fi] = ci

    for k in range(n):
        if k in chain_set:
            continue
        nb = neighbours[k]
        if len(nb) < 2:
            continue
        vecs = {j: (proj[j] - proj[k]) for j in nb}
        nb_list = list(nb)
        for a in range(len(nb_list)):
            for b in range(a+1, len(nb_list)):
                i, j = nb_list[a], nb_list[b]
                if i not in chain_set or j not in chain_set:
                    continue
                if chain_set[i] != chain_set[j]:
                    continue
                vi = vecs[i] / (np.linalg.norm(vecs[i]) + 1e-9)
                vj = vecs[j] / (np.linalg.norm(vecs[j]) + 1e-9)
                angle = np.degrees(np.arccos(np.clip(np.dot(vi, vj), -1, 1)))
                if angle >= collinearity_thresh:
                    ci = chain_set[i]
                    chains[ci].append(k)
                    chain_set[k] = ci
                    break

    chain_labels = np.full(n, -1, dtype=int)
    for ci, chain in enumerate(chains):
        for fi in chain:
            chain_labels[fi] = ci

    chains.sort(key=len, reverse=True)
    for ci, chain in enumerate(chains):
        for fi in chain:
            chain_labels[fi] = ci

    chain_lengths = [len(c) for c in chains]
    n_isolated    = int((chain_labels == -1).sum())
    n_chains      = len(chains)

    print(f"  Chains found: {n_chains}")
    print(f"  Isolated fibrils: {n_isolated}")
    if chain_lengths:
        print(f"  Chain lengths: min={min(chain_lengths)}  max={max(chain_lengths)}  "
              f"mean={np.mean(chain_lengths):.1f}")
    print(f"  Fibrils in chains: {int((chain_labels>=0).sum())}/{n} "
          f"({int((chain_labels>=0).sum())/n*100:.0f}%)")

    stats = {
        'n_chains':          n_chains,
        'n_isolated':        n_isolated,
        'mean_chain_length': float(np.mean(chain_lengths)) if chain_lengths else 0.0,
        'max_chain_length':  int(max(chain_lengths)) if chain_lengths else 0,
        'frac_in_chains':    float((chain_labels >= 0).sum() / n),
        'collinearity_thresh': collinearity_thresh,
        'search_factor':     search_factor,
    }

    CHAIN_COLORS = [
        '#00d4ff','#ff6b35','#7fff6b','#c87fff','#ffcc44',
        '#ff4d8b','#44eeff','#ffaa33','#aa44ff','#44ffaa',
        '#ff8844','#4488ff','#ffff44','#ff44ff','#44ff88',
        '#88ff44','#8844ff','#ff8888','#88ffff','#ffff88',
    ]

    fig, axes = plt.subplots(1, 3, figsize=(21, 7), facecolor='#0d1117')
    fig.suptitle('Local Row Chain Detection — Collinear Fibril Groups',
                 color='white', fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#0d1117')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e'); ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')
        ax.set_aspect('equal')

    ax = axes[0]; style(ax)
    ax.set_title(f'Row chains (n={n_chains} chains, {n_isolated} isolated)\nColour = chain identity', fontsize=9)

    for ci, chain in enumerate(chains):
        col = CHAIN_COLORS[ci % len(CHAIN_COLORS)]
        chain_proj = proj[chain]
        ctree = cKDTree(chain_proj)
        for fi_idx, fi in enumerate(chain):
            dists, idxs = ctree.query(chain_proj[fi_idx], k=min(3, len(chain)))
            for j_idx in idxs[1:]:
                if chain[fi_idx] < chain[j_idx]:
                    ax.plot([chain_proj[fi_idx, 0], chain_proj[j_idx, 0]],
                            [chain_proj[fi_idx, 1], chain_proj[j_idx, 1]],
                            '-', color=col, lw=1.5, alpha=0.5, zorder=2)

    for i in range(n):
        ci = chain_labels[i]
        col = CHAIN_COLORS[ci % len(CHAIN_COLORS)] if ci >= 0 else '#2a2a4a'
        ms  = 80 if ci >= 0 else 30
        ax.scatter(proj[i, 0], proj[i, 1], c=col, s=ms,
                   edgecolors='#070b12', linewidths=0.8, zorder=3)

    ax.set_xlabel('P1 (nm)'); ax.set_ylabel('P2 (nm)')

    ax = axes[1]; style(ax)
    ax.set_title('Chain length per fibril', fontsize=9)
    max_len = max(chain_lengths) if chain_lengths else 1
    cmap_len = cm.get_cmap('YlOrRd')
    norm_len = Normalize(vmin=0, vmax=max_len)

    for i in range(n):
        ci  = chain_labels[i]
        ln  = len(chains[ci]) if ci >= 0 else 0
        col = cmap_len(norm_len(ln))
        ax.scatter(proj[i, 0], proj[i, 1], c=[col], s=80,
                   edgecolors='#070b12', linewidths=0.8, zorder=3)

    sm = cm.ScalarMappable(cmap='YlOrRd', norm=norm_len)
    sm.set_array([])
    cb = plt.colorbar(sm, ax=axes[1], fraction=0.04, pad=0.02)
    cb.set_label('Chain length (fibrils)', color='#8b949e')
    cb.ax.yaxis.set_tick_params(color='#8b949e')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='#8b949e')
    axes[1].set_xlabel('P1 (nm)'); axes[1].set_ylabel('P2 (nm)')

    ax = axes[2]
    ax.set_facecolor('#161b22')
    for sp in ax.spines.values(): sp.set_color('#30363d')
    ax.tick_params(colors='#8b949e', labelsize=8)
    ax.xaxis.label.set_color('#8b949e'); ax.yaxis.label.set_color('#8b949e')
    ax.title.set_color('#e6edf3')

    if chain_lengths:
        bins = range(2, max(chain_lengths)+2)
        ax.hist(chain_lengths, bins=bins, color='#ffa657', edgecolor='#0d1117', alpha=0.85)
        ax.axvline(np.mean(chain_lengths), color='white', lw=1.5, ls='--',
                   label=f'Mean: {np.mean(chain_lengths):.1f}')
        ax.set_xlabel('Chain length (fibrils)'); ax.set_ylabel('Count')
        ax.legend(fontsize=8, labelcolor='#8b949e', facecolor='#0d1117', edgecolor='#30363d')

    ax.set_title(f'Chain length distribution\n'
                 f'{int((chain_labels>=0).sum())}/{n} fibrils in chains '
                 f'({int((chain_labels>=0).sum())/n*100:.0f}%)', fontsize=9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('row_chains.png', dpi=150, bbox_inches='tight', facecolor='#0d1117')
    print("  Saved row_chains.png")

    return {
        'chain_labels':  chain_labels.tolist(),
        'chains':        chains,
        'chain_lengths': chain_lengths,
        'proj_x':        proj[:, 0].tolist(),
        'proj_y':        proj[:, 1].tolist(),
        'stats':         stats,
    }


# =========================
# PLOT PACKING PLANE VIEW
# =========================
def plot_packing_plane_view(packing_plane_results, psi6_xy):
    COLORS = ['#00d4ff', '#ff6b35', '#7fff6b', '#c87fff', '#ffcc44']
    n_panels = len(packing_plane_results)
    if n_panels == 0:
        return

    fig, axes = plt.subplots(2, max(2, n_panels), figsize=(5 * n_panels, 10),
                             facecolor='#0d1117')
    fig.suptitle('Packing Plane View — Looking from True Normal to Fibril Sheet',
                 color='white', fontsize=14, fontweight='bold')

    def style(ax):
        ax.set_facecolor('#0d1117')
        for sp in ax.spines.values(): sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e', labelsize=8)
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e6edf3')
        ax.set_aspect('equal')

    cmap_psi = cm.get_cmap('plasma')
    norm_psi = Normalize(vmin=0, vmax=1)

    for i, (k, d) in enumerate(packing_plane_results.items()):
        proj_x = np.array(d['proj_x'])
        proj_y = np.array(d['proj_y'])
        psi6pp = np.array(d['psi6_pp_per'])
        psi6pp = np.array([0.0 if np.isnan(v) else v for v in psi6pp])

        col    = COLORS[i % len(COLORS)]
        label  = d['label']
        angle  = d['plane_angle']
        psi6v  = d['psi6_pp']
        n_excl = d.get('n_excluded', 0)

        ax = axes[0, i]; style(ax)
        ax.scatter(proj_x, proj_y, c=col, s=60, edgecolors='#070b12', linewidths=0.8, zorder=3)
        if len(proj_x) > 1:
            tree = cKDTree(np.column_stack([proj_x, proj_y]))
            for fi in range(len(proj_x)):
                dists, idxs = tree.query([proj_x[fi], proj_y[fi]], k=4)
                for dist, j in zip(dists[1:], idxs[1:]):
                    if fi < j:
                        ax.plot([proj_x[fi], proj_x[j]], [proj_y[fi], proj_y[j]],
                                '-', color=col, lw=0.8, alpha=0.3, zorder=2)
        ax.set_title(f'{label}\nψ6={psi6v:.3f}  normal→Z={angle:.0f}°  excl={n_excl}', fontsize=8)
        ax.set_xlabel('P1 (nm)'); ax.set_ylabel('P2 (nm)')

        ax = axes[1, i]; style(ax)
        sc = ax.scatter(proj_x, proj_y, c=psi6pp, cmap='plasma',
                        vmin=0, vmax=1, s=60, edgecolors='#070b12', linewidths=0.8, zorder=3)
        plt.colorbar(sc, ax=ax, label='ψ6', fraction=0.04)
        imp = psi6v - psi6_xy
        imp_col = '#7fff6b' if imp > 0.05 else '#ffcc44' if imp > 0 else '#f85149'
        ax.text(0.02, 0.97,
                f'ψ6_packing={psi6v:.3f}\nψ6_XY={psi6_xy:.3f}\n'
                f'{"+" if imp>=0 else ""}{imp:.3f} vs XY',
                transform=ax.transAxes, fontsize=7, va='top',
                color=imp_col, fontweight='bold')
        ax.set_xlabel('P1 (nm)'); ax.set_ylabel('P2 (nm)')

    for j in range(n_panels, axes.shape[1]):
        axes[0, j].set_visible(False)
        axes[1, j].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('packing_plane_view.png', dpi=150,
                bbox_inches='tight', facecolor='#0d1117')
    print("Saved packing_plane_view.png")


# =========================
# MAIN
# =========================
def main(file):
    fibrils = load_model2point(file)
    props   = get_fibril_properties(fibrils)

    print(f"\nTotal fibrils: {len(fibrils)}")

    orientation_labels, az, pol = orientation_analysis(props)
    lengths                     = length_distribution(props)
    nnd, same_nnd, diff_nnd     = nearest_neighbour_distances(props, orientation_labels)
    vor, voronoi_counts           = voronoi_analysis(props)
    voronoi_per_cluster           = voronoi_analysis_per_cluster(props, orientation_labels)
    axis_dists                    = axis_to_axis_distances(fibrils, props)
    bundle_labels                 = bundle_analysis(props, orientation_labels)
    S2, director, pairwise_angles = alignment_analysis(props)

    # Per-cluster psi6 + 2D RDF
    cluster_data = per_cluster_analysis(props, orientation_labels)

    # 3D hexagonal columnar analysis (end-on view) — global
    psi6_xy = float(cluster_data['global']['psi6'][PSI6_RADII[1]]['mean'])
    tilt_angles, azimuth_perp, psi6_perp, psi6_perp_per_fibril = \
        hexagonal_3d_analysis(props, director)

    # Per-cluster end-on analysis
    end_on_clusters = hexagonal_3d_per_cluster(props, orientation_labels, director)

    # Packing plane analysis
    packing_plane_results = packing_plane_per_cluster(props, orientation_labels)

    # Per-bundle end-on analysis
    bundle_end_on = bundle_end_on_analysis(props, bundle_labels, director)

    # Row periodicity
    row_periodicity = measure_row_periodicity(props, packing_plane_results,
                                              orientation_labels)

    # Row chain detection
    row_chains = detect_row_chains(props, packing_plane_results)

    # *** NEW: Discrete orientation group clustering with isolated fibrils ***
    group_labels, group_axes, n_groups = find_discrete_orientation_groups(props)

    # *** NEW: True 3D end-on projections per orientation group (the IMOD view) ***
    group_endon_results = orientation_group_endon(props, group_labels, group_axes)

    # Plots
    make_summary_figure(fibrils, props, orientation_labels, bundle_labels,
                        cluster_data, nnd, lengths, vor, voronoi_counts, S2)
    plot_per_cluster(cluster_data)
    plot_alignment(pairwise_angles, S2, props, orientation_labels)
    plot_3d_hex(tilt_angles, azimuth_perp, psi6_perp, psi6_perp_per_fibril,
                psi6_xy, director)
    plot_end_on_view(props, director, orientation_labels, psi6_perp_per_fibril)
    plot_packing_plane_view(packing_plane_results, psi6_xy)

    # *** NEW: Plot per-group end-on views (IMOD-style) ***
    plot_orientation_group_endon(group_endon_results, group_labels, group_axes)

    # JSON export — now includes centroid_3d_nm, azimuth, polar,
    # orientation_group, and orientation_group_endon projections
    export_json(fibrils, props, orientation_labels, bundle_labels,
                cluster_data, nnd, same_nnd, diff_nnd,
                lengths, axis_dists, voronoi_counts, voronoi_per_cluster, S2,
                psi6_perp=psi6_perp, end_on_clusters=end_on_clusters,
                packing_plane_results=packing_plane_results,
                bundle_end_on=bundle_end_on,
                row_periodicity=row_periodicity,
                row_chains=row_chains,
                az=az, pol=pol,
                group_labels=group_labels,
                group_axes=group_axes,
                group_endon_results=group_endon_results)

    print("\nDone. Outputs:")
    print("  fibril_analysis_summary.png")
    print("  per_cluster_analysis.png")
    print("  alignment_analysis.png")
    print("  hexagonal_3d_analysis.png")
    print("  end_on_view.png")
    print("  packing_plane_view.png")
    print("  row_periodicity.png")
    print("  row_chains.png")
    print("  orientation_group_endon.png   ← NEW: IMOD-style per-group views")
    print("  fibril_data.json")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fibrilorganisationscript_alignment_endon.py <model2point_output.txt>")
        sys.exit(1)
    main(sys.argv[1])
