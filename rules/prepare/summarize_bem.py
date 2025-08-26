"""
Use rasterio zonalstatistics to summarize BEM data within each admin unit.
"""

import logging
import sys
import glob
import os
import gc
from tqdm import tqdm
import rasterio
import exactextract as ee
import pandas as pd
import geopandas as gpd
import warnings
import shapely
from shapely.geometry import MultiPolygon
from shapely.ops import polygonize

# Suppress GDAL warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='osgeo')
warnings.filterwarnings('ignore', message='.*gdal.UseExceptions.*')

if __name__ == "__main__":
    try:
        adm2_path: str = snakemake.input["adm2"]
        adm1_path: str = snakemake.input["adm1"]
        adm0_path: str = snakemake.input["adm0"]
        res_path: str = snakemake.input["res_raster"]
        nres_path: str = snakemake.input["nres_raster"]
        adm2_output_path: str = snakemake.output["adm2"]
        adm1_output_path: str = snakemake.output["adm1"]
        adm0_output_path: str = snakemake.output["adm0"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# ----------------------------- Common functions -----------------------------

def summarize_bem(adm_path: str, res_path: str, nres_path: str, output_path: str, ADM_level: str):
    logging.info(f"Summarizing BEM data for {ADM_level}")
    gdf = gpd.read_file(adm_path)
    logging.info(f"Loaded {len(gdf)} features")

    with rasterio.open(res_path) as res_src, rasterio.open(nres_path) as nres_src:
        if gdf.crs != res_src.crs:
            logging.info("Reprojecting admin to match raster CRS")
            gdf = gdf.to_crs(res_src.crs)


        # Single-shot exactextract calls (no chunking)
        res_stats  = ee.exact_extract(res_src,  gdf, ['sum'], include_geom=False)
        nres_stats = ee.exact_extract(nres_src, gdf, ['sum'], include_geom=False)

        res_df  = pd.DataFrame([f.get('properties', {}) for f in res_stats]).add_prefix('res_')
        nres_df = pd.DataFrame([f.get('properties', {}) for f in nres_stats]).add_prefix('nres_')

        out = pd.concat([gdf.reset_index(drop=True)[["shapeName"]], res_df, nres_df], axis=1)

    # Collapse to one row per shapeName with sums
    out = out.groupby("shapeName", as_index=False)[["res_sum", "nres_sum"]].sum(min_count=1)

    logging.info(f"Writing {len(out)} rows → {output_path}")
    out.to_csv(output_path, index=False)



# -----------------------------------------------------------------------

# Run the analysis for each admin level
summarize_bem(adm2_path, res_path, nres_path, adm2_output_path, 'ADM2')
summarize_bem(adm1_path, res_path, nres_path, adm1_output_path, 'ADM1')
summarize_bem(adm0_path, res_path, nres_path, adm0_output_path, 'ADM0')

logging.info("Done.")
