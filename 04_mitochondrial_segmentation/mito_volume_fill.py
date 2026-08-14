"""
mito_volume_fill.py
===================
Fills mitochondrial volume from membrane segmentation using:
  1. Convex hull of outermost membrane voxels (interpolates across gaps)
  2. Clipped by dilated background mask (prevents over-extension)
  3. Fill interior voids
  4. Expand to include orphan granules near the boundary

Optimised for large cryo-ET volumes (400M+ voxels).

Requirements:
  pip install numpy scipy scikit-image tifffile

Usage:
  python mito_volume_fill.py --check   (identify label values first)
  python mito_volume_fill.py           (run the fill)
"""

import numpy as np
import tifffile
from scipy.spatial import ConvexHull
from scipy.ndimage import binary_fill_holes, distance_transform_edt
from skimage.morphology import ball

# =============================================================================
# CONFIGURATION
# =============================================================================

SEG_PATH          = "Pos018_ImageForProc_1.2_GranulesandMembranes_Pre-trained2DU-Net_Depth5.tiff"
OUTPUT_PATH       = "mito_filled.tif"
OUTPUT_OVERLAY    = "mito_overlay.tif"

MEMBRANE_LABEL    = 1
GRANULE_LABEL     = 2
BACKGROUND_LABEL  = 3    # set to None if no background mask

# How many voxels to expand background mask as safety margin (~15nm)
BACKGROUND_DILATION = 15  # voxels

# How far outside the filled volume an orphan granule can be and still
# be considered inside the mito (hull underestimate at sparse edges)
# 20 voxels ~ 20nm — generous enough to catch edge cases
ORPHAN_GRANULE_MARGIN = 20  # voxels

VOXEL_SIZE_NM     = 0.992

# =============================================================================


def check_labels(seg_path):
    print(f"Loading {seg_path}...")
    seg   = tifffile.imread(seg_path)
    total = seg.size
    print(f"Shape: {seg.shape}\n")
    print("Label values and voxel counts:")
    print("  (background mask = largest non-zero)")
    print("  (membrane        = medium)")
    print("  (granules        = smallest)\n")
    for val in np.unique(seg):
        count = int(np.sum(seg == val))
        print(f"  Label {val:>3}: {count:>12,} voxels  ({100*count/total:.2f}%)")


def dilate_fast(mask, radius):
    """
    Fast dilation using distance transform.
    Much faster than morphological dilation on large volumes.
    Expands mask by 'radius' voxels in all directions.
    """
    dist = distance_transform_edt(~mask)
    return dist <= radius


def fill_mito_volume(seg_path, output_path, output_overlay,
                     membrane_label, granule_label, background_label,
                     background_dilation, orphan_granule_margin,
                     voxel_size_nm):

    print(f"Loading {seg_path}...")
    seg = tifffile.imread(seg_path)
    print(f"  Shape: {seg.shape}")
    print(f"  Labels: {np.unique(seg)}")

    mem_mask  = (seg == membrane_label)
    gran_mask = (seg == granule_label)
    bg_mask   = (seg == background_label) if background_label is not None else None

    print(f"\nMembrane voxels:   {int(np.sum(mem_mask)):,}")
    if bg_mask is not None:
        print(f"Background voxels: {int(np.sum(bg_mask)):,}")

    if int(np.sum(mem_mask)) < 10:
        print(f"\nERROR: Too few membrane voxels — check MEMBRANE_LABEL={membrane_label}")
        return

    # -------------------------------------------------------------------------
    # Step 1: Convex hull of membrane voxels
    # -------------------------------------------------------------------------
    print("\nStep 1: Computing convex hull of membrane voxels...")
    coords = np.array(np.where(mem_mask)).T.astype(np.float64)

    if len(coords) > 20000:
        print(f"  Subsampling {len(coords):,} → 20,000 points...")
        idx    = np.random.choice(len(coords), 20000, replace=False)
        coords = coords[idx]

    hull   = ConvexHull(coords)
    print("  Rasterising hull (few minutes for large volumes)...")
    filled = rasterise_hull(hull, seg.shape)
    print(f"  Hull volume: {int(np.sum(filled)):,} voxels")

    # -------------------------------------------------------------------------
    # Step 2: Clip to expanded background mask
    # -------------------------------------------------------------------------
    if bg_mask is not None:
        print(f"\nStep 2: Expanding background mask by {background_dilation} voxels "
              f"(~{background_dilation * voxel_size_nm:.1f} nm)...")
        bg_expanded = dilate_fast(bg_mask, background_dilation)
        allowed     = bg_expanded | mem_mask | gran_mask
        before      = int(np.sum(filled))
        filled      = filled & allowed
        after       = int(np.sum(filled))
        clipped     = before - after
        print(f"  Clipped {clipped:,} voxels ({100*clipped/before:.1f}%) outside mask")
    else:
        print("\nStep 2: No background mask — skipping")

    # -------------------------------------------------------------------------
    # Step 3: Fill interior voids
    # -------------------------------------------------------------------------
    print("\nStep 3: Filling interior voids...")
    filled = binary_fill_holes(filled)
    print(f"  Volume after fill: {int(np.sum(filled)):,} voxels")

    # -------------------------------------------------------------------------
    # Step 4: Expand to include orphan granules near the boundary
    # Convex hull may slightly underestimate at edges where membrane is sparse,
    # leaving real granules just outside. Any granule within
    # ORPHAN_GRANULE_MARGIN voxels of the fill boundary is included.
    # -------------------------------------------------------------------------
    print(f"\nStep 4: Checking for orphan granules near boundary "
          f"(within {orphan_granule_margin} voxels / "
          f"~{orphan_granule_margin * voxel_size_nm:.0f} nm)...")

    dist_from_filled = distance_transform_edt(~filled)
    orphan_granules  = gran_mask & (dist_from_filled <= orphan_granule_margin)
    n_orphans        = int(np.sum(orphan_granules))

    if n_orphans > 0:
        print(f"  Found {n_orphans:,} orphan granule voxels — including in fill...")
        filled = filled | orphan_granules
        filled = binary_fill_holes(filled)
        print(f"  Volume after orphan inclusion: {int(np.sum(filled)):,} voxels")
    else:
        print("  No orphan granules found within margin")

    # -------------------------------------------------------------------------
    # Save filled volume
    # -------------------------------------------------------------------------
    tifffile.imwrite(output_path,
                     (filled * 255).astype(np.uint8),
                     compression="zlib")
    print(f"\nSaved filled volume: {output_path}")

    # Coloured overlay — cyan=filled, blue=membrane, red=granules
    overlay = np.zeros((*seg.shape, 3), dtype=np.uint8)
    overlay[filled]    = [0,   180, 180]   # cyan  = filled mito volume
    overlay[mem_mask]  = [50,  120, 255]   # blue  = membrane
    overlay[gran_mask] = [200,  50,  50]   # red   = granules
    tifffile.imwrite(output_overlay, overlay, compression="zlib")
    print(f"Saved overlay:       {output_overlay}")

    # -------------------------------------------------------------------------
    # Measurements
    # -------------------------------------------------------------------------
    voxel_vol_nm3 = voxel_size_nm ** 3
    filled_vox    = int(np.sum(filled))
    vol_nm3       = filled_vox * voxel_vol_nm3
    vol_um3       = vol_nm3 / 1e9
    gran_in_mito  = int(np.sum(gran_mask & filled))
    gran_vol_nm3  = gran_in_mito * voxel_vol_nm3
    gran_fraction = gran_vol_nm3 / vol_nm3 if vol_nm3 > 0 else 0

    print(f"\n{'='*48}")
    print(f"  Mito volume:       {vol_nm3:>12,.1f} nm³")
    print(f"  Mito volume:       {vol_um3:>12.6f} µm³")
    print(f"  Granule volume:    {gran_vol_nm3:>12,.1f} nm³")
    print(f"  Granule / mito:    {gran_fraction:>12.4f}")
    print(f"{'='*48}")
    print(f"\nOpen {output_overlay} in Fiji to validate visually.")
    print("Cyan should fill the full mito interior including any previously")
    print("orphaned granule regions.")


def rasterise_hull(hull, shape):
    """Rasterise convex hull into binary volume. Batched for memory efficiency."""
    zz, yy, xx = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]]
    pts         = np.column_stack([
        zz.ravel(), yy.ravel(), xx.ravel()
    ]).astype(np.float32)

    equations  = hull.equations.astype(np.float32)
    batch_size = 1_000_000
    inside     = np.zeros(len(pts), dtype=bool)
    n_batches  = int(np.ceil(len(pts) / batch_size))

    for i in range(n_batches):
        start = i * batch_size
        end   = start + batch_size
        batch = pts[start:end]
        inside[start:end] = np.all(
            batch @ equations[:, :3].T + equations[:, 3] <= 1e-6,
            axis=1
        )
        print(f"  Rasterising... {min(end, len(pts)):,}/{len(pts):,}", end="\r")

    print()
    return inside.reshape(shape)


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        check_labels(SEG_PATH)
    else:
        fill_mito_volume(
            seg_path              = SEG_PATH,
            output_path           = OUTPUT_PATH,
            output_overlay        = OUTPUT_OVERLAY,
            membrane_label        = MEMBRANE_LABEL,
            granule_label         = GRANULE_LABEL,
            background_label      = BACKGROUND_LABEL,
            background_dilation   = BACKGROUND_DILATION,
            orphan_granule_margin = ORPHAN_GRANULE_MARGIN,
            voxel_size_nm         = VOXEL_SIZE_NM,
        )
