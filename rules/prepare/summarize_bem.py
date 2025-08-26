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
def summarize_bem(adm_path: str, bem_res_raster_path: str, bem_nres_raster_path: str, output_path: str, ADM_level: str):
    logging.info(f"Summarizing BEM data for ADM level {ADM_level}")
    gdf = gpd.read_file(adm_path)

    # Open once (ExactExtract's Raster) and reuse
    res_r = ee.Raster(bem_res_raster_path)
    nres_r = ee.Raster(bem_nres_raster_path)
    try:
        # Use rasterio only to fetch CRS, but don't keep it open
        with rasterio.open(bem_res_raster_path) as res_src:
            if gdf.crs != res_src.crs:
                gdf = gdf.to_crs(res_src.crs)

        chunk_size = 200  # can be larger now; not tied to file handles
        all_results = []

        for i in tqdm(range(0, len(gdf), chunk_size), desc=f"Processing {ADM_level}"):
            chunk = gdf.iloc[i:i + chunk_size]

            # Pass the pre-opened Raster objects here
            res_stats = ee.exact_extract(res_r,  chunk, ['sum'], include_geom=False)
            nres_stats = ee.exact_extract(nres_r, chunk, ['sum'], include_geom=False)

            res_df = pd.DataFrame([f['properties'] for f in res_stats]).add_prefix('res_')
            nres_df = pd.DataFrame([f['properties'] for f in nres_stats]).add_prefix('nres_')

            chunk_results = pd.concat(
                [chunk.reset_index(drop=True), res_df.reset_index(drop=True), nres_df.reset_index(drop=True)],
                axis=1
            )
            all_results.append(chunk_results)

            # Free Python objects promptly
            del res_stats, nres_stats, res_df, nres_df, chunk_results, chunk
            gc.collect()

        if not all_results:
            logging.error("No chunks processed successfully!")
            return

        results_gdf = pd.concat(all_results, ignore_index=True)
        results_gdf.to_csv(
            output_path,
            index=False,
            columns=[c for c in ['shapeName', 'res_sum', 'nres_sum'] if c in results_gdf.columns]
        )
    finally:
        # Make sure GDAL handles are closed
        try: res_r.close()
        except: pass
        try: nres_r.close()
        except: pass


# -----------------------------------------------------------------------

# Run the analysis for each admin level
summarize_bem(adm2_path, res_path, nres_path, adm2_output_path, 'ADM2')
summarize_bem(adm1_path, res_path, nres_path, adm1_output_path, 'ADM1')
summarize_bem(adm0_path, res_path, nres_path, adm0_output_path, 'ADM0')

logging.info("Done.")
