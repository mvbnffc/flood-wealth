"""
Fix GHS-MOD so that population cells that aren't assigned a GHS-MOD value
are assigned the nearest valid one.
"""

import logging
import sys
import os

import numpy as np
import rasterio
from scipy import ndimage

if __name__ == "__main__":
    try:
        ghs_mod_path: str = snakemake.input["ghs_mod"]
        pop_path: str = snakemake.input["pop_file"]
        output_path: str = snakemake.output["fixed_ghs_mod"]
        country: str = snakemake.wildcards["ISO3"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Fixing GHS-MOD invalid values for {country} (nearest valid fill for populated cells).")

logging.info("Creating output directory (if it doesn't already exist)")
out_dir = os.path.dirname(output_path)
os.makedirs(out_dir, exist_ok=True)

VALID_CODES = np.array([11, 12, 13, 21, 22, 23, 30], dtype=np.float32)

def fill_invalid_by_nearest(arr: np.ndarray, valid_codes: np.ndarray, fill_mask: np.ndarray) -> np.ndarray:
    """
    Fill cells indicated by fill_mask with the nearest valid code in pixel space.
    Valid = value in valid_codes (and not NaN).
    """
    out = arr.astype(np.float32, copy=True)

    is_valid = np.isin(out, valid_codes) & (~np.isnan(out))

    if not np.any(fill_mask):
        logging.info("No cells require filling. Writing output unchanged.")
        return out

    if not np.any(is_valid):
        logging.warning("No valid GHS-MOD cells found anywhere. Writing output unchanged.")
        return out

    # For each pixel, get indices of the nearest valid pixel
    _, (ri, ci) = ndimage.distance_transform_edt(
        ~is_valid,
        return_distances=True,
        return_indices=True,
    )

    nearest_vals = out[ri, ci]
    out[fill_mask] = nearest_vals[fill_mask]
    return out

logging.info("Reading GHS-MOD and population rasters")
with rasterio.open(ghs_mod_path) as ghs_src, rasterio.open(pop_path) as pop_src:
    ghs = ghs_src.read(1).astype(np.float32, copy=False)
    pop = pop_src.read(1).astype(np.float32, copy=False)

    # Convert nodata to NaN for processing
    ghs_nodata = ghs_src.nodata
    pop_nodata = pop_src.nodata

    if ghs_nodata is not None:
        ghs = ghs.copy()
        ghs[ghs == ghs_nodata] = np.nan

    if pop_nodata is not None:
        pop = pop.copy()
        pop[pop == pop_nodata] = np.nan

    # Basic alignment check (fail loudly if mismatch)
    if (ghs.shape != pop.shape) or (ghs_src.transform != pop_src.transform):
        raise ValueError(
            f"GHS-MOD and population rasters are not aligned for {country} "
            "(shape and/or transform differ). Reproject/resample first."
        )

    profile = ghs_src.profile.copy()

logging.info("Building fill mask: pop > 0 AND GHS-MOD is invalid (incl. NaN)")
pop_positive = (~np.isnan(pop)) & (pop > 0)
is_valid = np.isin(ghs, VALID_CODES) & (~np.isnan(ghs))
fill_mask = pop_positive & (~is_valid)

logging.info(f"Cells to fill: {int(fill_mask.sum()):,}")

logging.info("Filling invalid GHS-MOD values by nearest valid cell")
fixed = fill_invalid_by_nearest(ghs, VALID_CODES, fill_mask)

logging.info("Preparing output raster (restoring nodata)")
# If input had no nodata, choose one that won't collide with valid codes
if ghs_nodata is None:
    ghs_nodata = -9999.0

fixed_out = fixed.copy()
fixed_out[np.isnan(fixed_out)] = np.float32(ghs_nodata)

profile.update(dtype="float32", nodata=np.float32(ghs_nodata), compress="deflate")

logging.info("Writing fixed GHS-MOD raster")
with rasterio.open(output_path, "w", **profile) as dst:
    dst.write(fixed_out.astype(np.float32), 1)

logging.info("Done.")

