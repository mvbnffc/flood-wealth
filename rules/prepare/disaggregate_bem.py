"""
Distribute GIRI BEM capital stock data across an building volume layer from GHSL
This is used to estimate the RES / NRES capital stock at risk from flooding.
"""

import argparse
import logging
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.windows import Window
from rasterio.warp import reproject
from rasterio.env import Env
from rasterio.vrt import WarpedVRT



if __name__ == "__main__":

    try:
        res_bem_path: str = snakemake.input["res_bem"]
        nres_bem_path: str = snakemake.input["nres_bem"]
        res_volume_path: str = snakemake.input["res_volume"]
        nres_volume_path: str = snakemake.input["nres_volume"]
        res_output_path: str = snakemake.output["res_capstock"]
        nres_output_path: str = snakemake.output["nres_capstock"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Disaggregating GIRI Building Exposure Model (BEM) across GHS Built Volume Layers.")

# --------------------- HELPER FUNCTIONS ---------------------

def read_f32(ds, window):
    return ds.read(1, window=window, out_dtype="float32", fill_value=0.0)

def clip_array_to_factor(a: np.ndarray, factor: int) -> np.ndarray:
    """Clip array so both dims are divisible by factor."""
    nrows = (a.shape[0] // factor) * factor
    ncols = (a.shape[1] // factor) * factor
    if nrows == 0 or ncols == 0:
        return a[:0, :0]
    return a[:nrows, :ncols]

def resample_sum_block(a: np.ndarray, factor: int) -> np.ndarray:
    """Block-sum downsample by integer factor."""
    r, c = a.shape
    rr = r // factor
    cc = c // factor
    a = a[: rr * factor, : cc * factor]
    a = a.reshape(rr, factor, cc, factor)
    return a.sum(axis=(1, 3))

def repeat_2d(a: np.ndarray, factor: int) -> np.ndarray:
    """Nearest-repeat upscale by integer factor on both axes."""
    return np.repeat(np.repeat(a, factor, axis=0), factor, axis=1)

# -----------------------------------------------------------------

logging.info("Opening rasters and preparing a snapped 3\" grid aligned to the 150\" grid.")
with rasterio.open(res_bem_path) as res_bem_src, rasterio.open(nres_bem_path) as nres_bem_src, \
        rasterio.open(res_volume_path) as res_volume_src, rasterio.open(nres_volume_path) as nres_volume_src:
    if res_bem_src.crs != res_volume_src.crs:
        raise ValueError("GIRI and GHS rasters must have the same CRS.")
    
    SCALE = 50 # 150" / 3"
    # Build an exact snapped 3" grid that nests perfectly in the 150" grid
    dx3 = res_bem_src.transform.a / SCALE
    dy3 = res_bem_src.transform.e / SCALE
    c0 = res_bem_src.transform.c
    f0 = res_bem_src.transform.f
    width3 = res_bem_src.width * SCALE
    height3 = res_bem_src.height * SCALE
    snapped_transform = Affine(dx3, 0.0, c0, 0.0, dy3, f0)

    # Create a temp workspace
    os.makedirs(os.path.dirname(res_output_path), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        snapped_res_vol_path = os.path.join(tmpdir, "snapped_res_volume_3ss.tif")
        snapped_nres_vol_path = os.path.join(tmpdir, "snapped_nres_volume_3ss.tif")
        snapped_res_bem_path = os.path.join(tmpdir, "snapped_res_bem_3ss.tif")
        snapped_nres_bem_path = os.path.join(tmpdir, "snapped_nres_bem_3ss.tif")

        logging.info("Snapping built volume to coarse=aligned 3\" grid (nearest).")
        volume_snap_profile = res_volume_src.profile.copy()
        volume_snap_profile.update({
            "transform": snapped_transform,
            "width": width3,
            "height": height3,
            "dtype": "float32",
            "compress": "lzw",
            "BIGTIFF": "YES",
            "count": 1,
            "tiled": True,
            "nodata": 0.0
                })
        with rasterio.open(res_volume_path) as res_volume_src, rasterio.open(nres_volume_path) as nres_volume_src:
            res_vol_vrt = WarpedVRT(
                res_volume_src,
                crs=res_bem_src.crs,           # already equal per your check
                transform=snapped_transform,
                width=width3,
                height=height3,
                resampling=Resampling.nearest,
                src_nodata=0.0, nodata=0.0,
                warp_mem_limit=512
            )
            nres_vol_vrt = WarpedVRT(
                nres_volume_src,
                crs=nres_bem_src.crs,
                transform=snapped_transform,
                width=width3,
                height=height3,
                resampling=Resampling.nearest,
                src_nodata=0.0, nodata=0.0,
                warp_mem_limit=512
            )

        logging.info("Disaggregating on snapped grid (exact 50x50 nesting)")
        res_total_out_snapped = 0.0
        nres_total_out_snapped = 0.0
        res_total_coarse = 0.0
        nres_total_coarse = 0.0
        coarse_rows_per_chunk = 400
        coarse_cols_per_chunk = 2000

        with res_vol_vrt as res_vol_snap, rasterio.open(res_bem_path) as res_bem_src, \
                nres_vol_vrt as nres_vol_snap, rasterio.open(nres_bem_path) as nres_bem_src:

            # Prepare output (snapped 3" grid)
            out_snap_profile = res_vol_snap.profile.copy()
            out_snap_profile.update(dtype="float32", count=1, compress="lzw", BIGTIFF="YES", tiled=True, nodata=0.0)

            block_x = out_snap_profile.get("blockxsize", 256)
            block_y = out_snap_profile.get("blockysize", 256)

            # use large, block-aligned chunks to minimize partial tiles
            coarse_rows_per_chunk = max(block_y // SCALE, 256)    # coarse cells per chunk row
            coarse_cols_per_chunk = max((block_x * 8) // SCALE, 512)

            with rasterio.open(snapped_res_bem_path, "w", **out_snap_profile) as res_bem_snap_dst, \
                rasterio.open(snapped_nres_bem_path, "w", **out_snap_profile) as nres_bem_snap_dst:

                for row in range(0, res_bem_src.height, coarse_rows_per_chunk):
                    h_c = min(coarse_rows_per_chunk, res_bem_src.height - row)
                    r0_f = row * SCALE
                    h_f = h_c * SCALE

                    for col in range(0, res_bem_src.width, coarse_cols_per_chunk):
                        w_c = min(coarse_cols_per_chunk, res_bem_src.width - col)
                        c0_f = col * SCALE
                        w_f = w_c * SCALE

                        # Defin the windows
                        win_c = Window(col, row, w_c, h_c)
                        win_f = Window(c0_f, r0_f, w_f, h_f)    

                        # Read snapped fine volume window (masked are zero)
                        res_vol_fine = read_f32(res_vol_snap, window=win_f)
                        nres_vol_fine = read_f32(nres_vol_snap, window=win_f)
                        res_vol_fine_data = np.asarray(res_vol_fine.filled(0.0), dtype=np.float64)
                        nres_vol_fine_data = np.asarray(nres_vol_fine.filled(0.0), dtype=np.float64)

                        # Ensure block exactness
                        res_vol_fine_data = clip_array_to_factor(res_vol_fine_data, SCALE)
                        nres_vol_fine_data = clip_array_to_factor(nres_vol_fine_data, SCALE)
                        if res_vol_fine_data.size == 0 and nres_vol_fine_data.size == 0:
                            continue

                        h_f_eff, w_f_eff = res_vol_fine_data.shape
                        h_c_eff = h_f_eff // SCALE
                        w_c_eff = w_f_eff // SCALE
                        win_f_eff = Window(c0_f, r0_f, w_f_eff, h_f_eff)
                        win_c_eff = Window(col, row, w_c_eff, h_c_eff)

                        # Aggregate fine volume to coarse
                        res_vol_coarse = resample_sum_block(res_vol_fine_data, SCALE)
                        nres_vol_coarse = resample_sum_block(nres_vol_fine_data, SCALE)

                        # Read matching coarse BEM exposure
                        res_bem_coarse_ma = read_f32(res_bem_src, window=win_c_eff)
                        nres_bem_coarse_ma = read_f32(nres_bem_src, window=win_c_eff)
                        res_bem_coarse = np.asarray(res_bem_coarse_ma.filled(0.0), dtype=np.float64)
                        nres_bem_coarse = np.asarray(nres_bem_coarse_ma.filled(0.0), dtype=np.float64)

                        if res_bem_coarse.shape != res_vol_coarse.shape or nres_bem_coarse.shape != nres_vol_coarse.shape:
                            raise ValueError(f"Mismatch in coarse BEM and volume shapes: {res_bem_coarse.shape} vs {res_vol_coarse.shape} "
                                             f"and {nres_bem_coarse.shape} vs {nres_vol_coarse.shape} at row={row}, col={col}.")
                        
                        # Coarse exposure-per-volume (safe divide)
                        nres_exp_per_vol_coarse = np.divide(
                            nres_bem_coarse, nres_vol_coarse, 
                            out=np.zeros_like(nres_bem_coarse, dtype=np.float64), 
                            where=nres_vol_coarse != 0.0
                        )
                        res_exp_per_vol_coarse = np.divide(
                            res_bem_coarse, res_vol_coarse, 
                            out=np.zeros_like(res_bem_coarse, dtype=np.float64), 
                            where=res_vol_coarse != 0.0
                        )

                        # Nearest-repeat the coarse exposure-per-volume to fine grid
                        nres_exp_per_vol_fine = repeat_2d(nres_exp_per_vol_coarse, SCALE)
                        nres_exp = (nres_exp_per_vol_fine * nres_vol_fine_data).astype(np.float32)
                        res_exp_per_vol_fine = repeat_2d(res_exp_per_vol_coarse, SCALE)
                        res_exp = (res_exp_per_vol_fine * res_vol_fine_data).astype(np.float32)

                        # Write snapped 3" window
                        res_bem_snap_dst.write(res_exp, 1, window=win_f_eff)
                        nres_bem_snap_dst.write(nres_exp, 1, window=win_f_eff)

                        # Running sums (conservation on snapped grid)
                        res_total_out_snapped += float(np.sum(res_exp, dtype="float64"))
                        nres_total_out_snapped += float(np.sum(nres_exp, dtype="float64"))
                        res_total_coarse += float(np.sum(res_bem_coarse), dtype="float64")
                        nres_total_coarse += float(np.sum(nres_bem_coarse), dtype="float64")

        logging.info("Snapped-grid conservation check.")
        res_diff_snap = res_total_out_snapped - res_total_coarse
        nres_diff_snap = nres_total_out_snapped - nres_total_coarse
        res_rel_snap = 0.0 if res_total_coarse == 0 else abs(res_diff_snap) / res_total_coarse
        nres_rel_snap = 0.0 if nres_total_coarse == 0 else abs(nres_diff_snap) / nres_total_coarse
        logging.info(f"Check Residential (snapped): sum(out) = {res_total_out_snapped:,.0f} vs sum(coarse) = {res_total_coarse:,.0f} "
                    f"(diff={res_diff_snap:,.2f}, rel={res_rel_snap:.2e}).")
        logging.info(f"Check Non-Residential (snapped): sum(out) = {nres_total_out_snapped:,.0f} vs sum(coarse) = {nres_total_coarse:,.0f} "
                    f"(diff={nres_diff_snap:,.2f}, rel={nres_rel_snap:.2e}).")
        

        logging.info("Remapping disaggregated (snapped) result back to ORIGINAL 3\" grid (sum-conserving).")
        with rasterio.open(snapped_res_bem_path) as res_bem_snap, rasterio.open(res_volume_path) as res_volume, \
            rasterio.open(snapped_nres_bem_path) as nres_bem_snap, rasterio.open(nres_volume_path) as nres_volume:
            final_profile = res_volume.profile.copy()
            final_profile.update(
                dtype="float32",
                count=1, 
                compress="lzw",
                BIGTIFF="YES",
                tiled=True,
                nodata = 0.0
            )
            with rasterio.open(res_output_path, "w", **final_profile) as res_dst, rasterio.open(nres_output_path, "w", **final_profile) as nres_dst:
                reproject(
                    source=rasterio.band(res_bem_snap, 1),
                    destination=rasterio.band(res_dst, 1),
                    src_transform=res_bem_snap.transform,
                    src_crs=res_bem_snap.crs,
                    dst_transform=res_volume.transform,
                    dst_crs=res_volume.crs,
                    resampling=Resampling.sum, # mass-conserving back to original geolocation
                    num_threads=2
                )
                reproject(
                    source=rasterio.band(nres_bem_snap, 1),
                    destination=rasterio.band(nres_dst, 1),
                    src_transform=nres_bem_snap.transform,
                    src_crs=nres_bem_snap.crs,
                    dst_transform=nres_volume.transform,
                    dst_crs=nres_volume.crs,
                    resampling=Resampling.sum, # mass-conserving back to original geolocation
                    num_threads=2
                )

logging.info("Done.")