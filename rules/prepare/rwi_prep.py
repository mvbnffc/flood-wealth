"""
Script for preparing the releative wealth index (RWI)

This script clips, reprojects, and resamples the RWI using rasterio
The population dataset is used as reference (so that all datasets 
align during subsuqent analysis)
"""

import logging
import sys
import os

import rasterio
import fiona
from rasterio.warp import reproject, Resampling
import rasterio.mask


if __name__ == "__main__":
    try:
        pop_path: str = snakemake.input["pop_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        boundary_path: str = snakemake.input['boundary_file']
        output_path: str = snakemake.output["rwi_resampled"]
        country: str = snakemake.wildcards["ISO3"]
    except:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Clipping, reprojecting, and resampling relative wealth index for {country}.")

logging.info("Opening population file and extracting its grid properties")
with rasterio.open(pop_path) as pop_src:
        pop_transform = pop_src.transform
        pop_width = pop_src.width
        pop_height = pop_src.height
        pop_crs = pop_src.crs

logging.info("Loading boundary file for clipping")
# Read the boundary file and extract the geometries
with fiona.open(boundary_path, "r") as shapefile:
    boundary_shapes = [feature["geometry"] for feature in shapefile]

logging.info("Open the RWI file and reproject onto the population grid")
with rasterio.open(rwi_path) as rwi_src:
    # Copy metadata and update with population grid
    kwargs = rwi_src.meta.copy()
    kwargs.update({
        'crs': pop_crs,
        'transform': pop_transform,
        'width': pop_width,
        'height': pop_height,
        'compress': 'lzw'
    })

    # Create an intermediate output path (that will be deleted)
    temp_path = output_path.replace(".tif", "_temp.tif")

    logging.info("Saving reprojected RWI file (temporary).")
    with rasterio.open(temp_path, 'w', **kwargs) as dst:
        reproject(
            source=rasterio.band(rwi_src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=rwi_src.transform,
            src_crs=rwi_src.crs,
            dst_transform=pop_transform,
            dst_crs=pop_crs,
            resampling=Resampling.nearest,
            bigtiff='YES'
        )


logging.info("Clipping the reprojected RWI file using the boundary file")
with rasterio.open(temp_path) as src:
    # mask() returns both the clipped image and the new transform.
    out_image, out_transform = rasterio.mask.mask(src, boundary_shapes, crop=True)
    out_meta = src.meta.copy()
    out_meta.update({
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform,
        "bigtiff": "YES",
        'compress': 'lzw'
    })

    
    logging.info("Saving clipped RWI file.")
    with rasterio.open(output_path, "w", **out_meta) as dst:
        dst.write(out_image)
    
logging.info("Deleting temporary files.")
os.remove(temp_path)

logging.info("Done.")
