"""
Script for preparing the gridded gross domestic product (GDP) data

This script clips, reprojects, and resamples the GDP using rasterio
The population dataset is used as reference (so that all datasets 
align during subsuqent analysis)
"""

import logging
import sys
import os

import rasterio
import fiona
import numpy as np
from rasterio.warp import reproject, Resampling
import rasterio.mask


if __name__ == "__main__":
    try:
        pop_path: str = snakemake.input["pop_file"]
        gdp_path: str = snakemake.input["gdp_file"]
        boundary_path: str = snakemake.input['boundary_file']
        output_path: str = snakemake.output["gdp_resampled"]
        country: str = snakemake.wildcards["ISO3"]
    except:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Clipping, reprojecting, and resampling GDP data for {country}.")

logging.info("Calculating national GDP from original projection (EPSG:54009)")
# Reproject boundary for this calculation
with fiona.open(boundary_path, "r") as shapefile:
    boundary_crs = shapefile.crs

    with rasterio.open(gdp_path) as gdp_src:
        gdp_crs = gdp_src.crs

        # Reproject boundary geometries to match GDP CRS
        if boundary_crs != gdp_crs:
            boundary_shapes_54009 = []
            for feature in shapefile:
                geom_54009 = rasterio.warp.transform_geom(boundary_crs, gdp_crs, feature["geometry"])
                boundary_shapes_54009.append(geom_54009)
        else:
            boundary_shapes_54009 = [feature["geometry"] for feature in shapefile]

# Calculate national GDP using original projection.
with rasterio.open(gdp_path) as gdp_src:
    orig_clip, _ = rasterio.mask.mask(gdp_src, boundary_shapes_54009, crop=True)
    orig_band = orig_clip[0] if len(orig_clip.shape) == 3 else orig_clip
    nodata = gdp_src.nodata

    # build a mask for valid GDP values
    if nodata is not None:
        valid_mask = (orig_band != nodata) & np.isfinite(orig_band)
    else:
        valid_mask = np.isfinite(orig_band)

    # Calculate the national GDP by summing valid GDP values
    national_gdp = float(orig_band[valid_mask].sum()) if valid_mask.sum() > 0 else 0.0

logging.info(f"National GDP for {country} (original projection): {national_gdp}")

logging.info("Opening population file and extracting its grid properties")
with rasterio.open(pop_path) as pop_src:
        pop_transform = pop_src.transform
        pop_width = pop_src.width
        pop_height = pop_src.height
        pop_crs = pop_src.crs

logging.info("Open the GDP file and reproject onto the population grid")
temp_gdp_reprojected = output_path.replace(".tif", "_temp_reprojected.tif")

with rasterio.open(gdp_path) as gdp_src:
    # Copy metadata and update with population grid
    kwargs = gdp_src.meta.copy()
    kwargs.update({
        'crs': pop_crs,
        'transform': pop_transform,
        'width': pop_width,
        'height': pop_height,
        'compress': 'lzw',
        'BIGTIFF': 'YES'
    })

    logging.info("Saving reprojected GDP file (temporary).")
    with rasterio.open(temp_gdp_reprojected, 'w', **kwargs) as dst:
        reproject(
            source=rasterio.band(gdp_src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=gdp_src.transform,
            src_crs=gdp_src.crs,
            dst_transform=pop_transform,
            dst_crs=pop_crs,
            resampling=Resampling.nearest,
        )

logging.info("Loading boundary file for clipping")
# Read the boundary file and extract the geometries
with fiona.open(boundary_path, "r") as shapefile:
    boundary_shapes = [feature["geometry"] for feature in shapefile]

logging.info("Clipping the reprojected GDP file using the boundary file")
with rasterio.open(temp_gdp_reprojected) as src:
    # mask() returns both the clipped image and the new transform.
    out_image, out_transform = rasterio.mask.mask(src, boundary_shapes, crop=True)

    out_meta = src.meta.copy()
    out_meta.update({
        "dtype": 'float32',
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform,
        'compress': 'lzw'
    })

logging.info("Adjust the GDP to match national GDP totals.")

# Calculate the new GDP sum
clipped_band = out_image[0]
dst_nodata = src.nodata
valid_proj = (clipped_band != dst_nodata) if dst_nodata is not None else ~np.isnan(clipped_band)
reproj_GDP_sum = float(clipped_band[valid_proj].sum())

# Compute and apply the scaling factor
if reproj_GDP_sum > 0:
    scaling_factor = national_gdp / reproj_GDP_sum
    out_image = out_image * scaling_factor

    # Verify the scaling worked
    final_sum = float(out_image[0][valid_proj].sum())
    if not np.isclose(final_sum, national_gdp, atol=10):
        logging.warning(f"Final GDP sum {final_sum} does not match national GDP {national_gdp} after scaling.")
else:
    logging.warning("Projected GDP sum is zero or negative, skipping scaling.")

out_image = out_image.astype(np.float32)  # Ensure the output is in float32 format
out_image = np.where(out_image < 0, 0, out_image) # Remove any negative values. 

logging.info("Saving the final clipped and resampled GDP file.")
with rasterio.open(output_path, "w", **out_meta) as dst:
    dst.write(out_image)

logging.info("Deleting temporary files.")
os.remove(temp_gdp_reprojected)

logging.info("Done.")
