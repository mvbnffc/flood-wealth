"""
Disaggregate GIRI BEM capital stock data at the ADM2 across an building volume layer from GHSL
This is used to estimate the RES / NRES capital stock at risk from flooding.
"""

import argparse
import logging
import math
import os
import tempfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio
from rasterio import features
import numpy as np



if __name__ == "__main__":

    try:
        adm_path: str = snakemake.input["adm2_file"]
        bem_path: str = snakemake.input["giri_bem"]
        res_volume_path: str = snakemake.input["res_volume"]
        nres_volume_path: str = snakemake.input["nres_volume"]
        res_output_path: str = snakemake.output["res_capstock"]
        nres_output_path: str = snakemake.output["nres_capstock"]
        country: str = snakemake.wildcards["ISO3"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Disaggregating GIRI Building Exposure Model (BEM) across GHS Built Volume Layers for {country}.")

logging.info("Loading ADM2 layer.")
gdf = gpd.read_file(adm_path, layer="ADM2")

logging.info(f"Load GIRI BEM CSV and filter for {country}")
df = pd.read_csv(bem_path)
df = df[df["shapeGroup"] == country].copy()

logging.info("Joining relevant CSV columns with admin layer.")
common_column = "shapeID"
# What columns do we need
df_use = df[[common_column, "res_sum", "nres_sum"]].copy()
gdf_use = gdf[[common_column, "geometry"]].copy()
gdfm = gdf_use.merge(df_use, on=common_column, how='left')
# Fill missing totals with 0
gdfm[["res_sum", "nres_sum"]] = gdfm[["res_sum", "nres_sum"]].fillna(0.0)

print(gdfm.head())

logging.info("Building look-up arrays for raster disaggregation")
# Build an intiger index per admin
gdfm = gdfm.reset_index(drop=True)
gdfm["_idx"] = np.arange(1, len(gdfm) + 1, dtype=np.int32)
# Lookup arrays
idx_to_res = np.zeros(len(gdfm) + 1, dtype=np.float64)
idx_to_nres = np.zeros(len(gdfm) + 1, dtype=np.float64)
idx_to_res[gdfm["_idx"].values] = gdfm["res_sum"].values
idx_to_nres[gdfm["_idx"].values] = gdfm["nres_sum"].values

# ---------------- HELPER FUNCTION --------------------------
def allocate_to_volume(volume_path: str, totals_by_idx: np.ndarray, label: str, out_path: str):
    with rasterio.open(volume_path) as src:
        # reproject ADM2 to raster CRS if needed
        g_adm = gdfm.to_crs(src.crs) if gdfm.crs != src.crs else gdfm

        # rasterize admin indices onto this grid
        logging.info(f"Rasterizing ADM2 indices for {label} grid.")
        shapes = zip(g_adm.geometry, g_adm["_idx"].astype(int))
        admin_idx = features.rasterize(
            shapes=shapes,
            out_shape=(src.height, src.width),
            transform=src.transform,
            fill=0,
            dtype="int32",
            all_touched=False,  # = True for inclusive borders
        )

        # read volume, sanitize
        vol = src.read(1, masked=True).filled(0.0).astype(np.float64)
        vol[vol < 0] = 0.0

        # per-admin total built volume (vectorized)
        vol_sums = np.bincount(admin_idx.ravel(), weights=vol.ravel(), minlength=len(totals_by_idx))

        # scale per admin: total_cap / total_volume
        scale = np.zeros_like(totals_by_idx, dtype=np.float64)
        nz = vol_sums > 0
        scale[nz] = totals_by_idx[nz] / vol_sums[nz]

        # allocate per pixel
        alloc = (vol * scale[admin_idx]).astype("float32")

        # write output aligned to this volume raster
        meta = src.meta.copy()
        meta.update(dtype="float32", count=1, compress="lzw", BIGTIFF="YES", nodata=0.0)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with rasterio.open(out_path, "w", **meta) as dst:
            dst.write(alloc, 1)

        # diagnostics
        allocated_sum = float(np.sum(alloc, dtype=np.float64))
        requested_sum = float(np.sum(totals_by_idx, dtype=np.float64))
        zero_vol_admins = int(np.sum((vol_sums[1:] == 0) & (totals_by_idx[1:] > 0)))
        logging.info(f"{label}: sum(output)={allocated_sum:,.0f}, sum(adm totals)={requested_sum:,.0f}, "
                     f"zero-volume admins={zero_vol_admins}")
        
logging.info("Disaggregating Residential Capital Stock")
allocate_to_volume(res_volume_path, idx_to_res, "RES", res_output_path)

logging.info("Disaggregating Non-Residential Capital Stock")
allocate_to_volume(nres_volume_path, idx_to_nres, "NRES", nres_output_path)

logging.info("Done.")