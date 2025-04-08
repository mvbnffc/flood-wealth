"""
Script for creating a merged raster of the Google historical inundation datasets

Raw data comes in km^2 GeoJSON tiles classfied as high, medium, and low risk:
    High: wet at least 5% of the time
    Medium: wet at least 1% of the time
    Low: wet at least 0.5% of the time

This script merges the GeoJSONs into one global (km^2 resolution) GeoTiff with cell values
corresponding to the probability of flooding in each grid cell (high=0.05, medium=0.01, low=0.005)
"""

import logging
import glob
import os

import numpy as np
import rasterio
from tqdm import tqdm

if __name__ == "__main__":
    try:
        input_path: str = snakemake.input["merge_gfd_folder"]
        output_path: str = snakemake.output["merge_gfd_file"]
    except:
        raise ValueError("Must be run via snakemake.")
    
# Set some analysis parameters
raster_resolution = 0.002245788210298803843 # degrees

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Merging Global Flood Database files into global raster.")

#### Define functions for the analysis #####################
def get_global_extent(files):
    """Calculate the combined extent of all rasters."""
    min_left, min_bottom, max_right, max_top = float('inf'), float('inf'), float('-inf'), float('-inf')
    for file in files:
        with rasterio.open(file) as src:
            bounds = src.bounds
            min_left = min(min_left, bounds.left)
            min_bottom = min(min_bottom, bounds.bottom)
            max_right = max(max_right, bounds.right)
            max_top = max(max_top, bounds.top)
    return (min_left, min_bottom, max_right, max_top)

def calculate_offsets(raster_bounds, global_extent, transform):
    """Calculate the row and column offsets for the raster within the global raster."""
    row_offset = int((raster_bounds.top - global_extent[3]) / -transform[4])
    col_offset = int((raster_bounds.left - global_extent[0]) / transform[0])
    return row_offset, col_offset

def pad_and_add_raster(src, global_raster, row_offset, col_offset, global_extent, transform):
    raster_array = src.read(1, out_dtype=np.int16)

    # Calculate the necessary padding
    top_padding = -row_offset
    left_padding = col_offset
    bottom_padding = max(global_raster.shape[0] - (raster_array.shape[0] + top_padding), 0)
    right_padding = max(global_raster.shape[1] - (raster_array.shape[1] + left_padding), 0)

    # Pad the raster array
    padded_raster = np.pad(
        raster_array,
        ((top_padding, bottom_padding), (left_padding, right_padding)),
        'constant',
        constant_values=(0, 0)  # Assuming 0 is the no-data value
    )

    # Add the padded raster to the global raster
    global_raster += padded_raster

###########################################################

logging.info("Reading raster file names.")
raster_files = glob.glob(os.path.join(input_path, "*.tif"))

logging.info("Calculate global extent.")
global_extent = get_global_extent(raster_files)
global_width = int((global_extent[2] - global_extent[0]) / raster_resolution)
global_height = int((global_extent[3] - global_extent[1]) / raster_resolution)

logging.info("Initialize the global raster")
global_raster = np.zeros((global_height, global_width), dtype=np.int16)

logging.info("Processing rasters and merging into global raster.")
for file in tqdm(raster_files):
    with rasterio.open(file) as src:
        row_offset, col_offset = calculate_offsets(src.bounds, global_extent, src.transform)
        pad_and_add_raster(src, global_raster, row_offset, col_offset, global_extent, src.transform)

logging.info("Save the global raster.")
with rasterio.open(output_path,
                   'w',
                   driver='GTiff',
                   width=global_width,
                   height=global_height,
                   count=1,
                   dtype=global_raster.dtype,
                   crs=src.crs,
                   transform=rasterio.transform.from_origin(global_extent[0], global_extent[3], raster_resolution, raster_resolution)
                   ) as dst:
    dst.write(global_raster, 1)

logging.info("Done.")

