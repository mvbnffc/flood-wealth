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
import re
import os
import glob
import json
from tqdm import tqdm
from urllib.parse import urljoin

import numpy as np
import rasterio
from rasterio.transform import from_origin
import rasterio.features
from shapely.geometry import shape, box


if __name__ == "__main__":
    try:
        input_path: str = snakemake.input["raw_folder"]
        pop_path: str = snakemake.input["pop_path"]
        output_path: str = snakemake.output["merged_file"]
    except:
        raise ValueError("Must be run via snakemake.")
    
# Set GeoJSON folder
geojson_folder = os.path.join(input_path, "inundation_history", "data")
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Merging the Google historical inundation data into a Global GeoTiff")

#### Define some functions for the analysis ####
def load_geometries_from_file(file_path):
    """
    Given a GeoJSON file path, load geometries for each risk level.
    Handles both a simple dictionary format and a FeatureCollection.
    Returns a dictionary mapping risk levels to a list of shapely geometries.
    """
    geoms = {"High_risk": [], "Medium_risk": [], "Low_risk": []}
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return geoms

    if "features" in data:
        # Handle FeatureCollection format
        for feature in data["features"]:
            risk = feature.get("properties", {}).get("name")
            if risk in geoms:
                try:
                    geom = shape(feature["geometry"])
                    geoms[risk].append(geom)
                except Exception as e:
                    print(f"Error processing geometry in {file_path} for {risk}: {e}")
    else:
        # Fallback for a non-FeatureCollection structure
        for risk in geoms.keys():
            if risk in data:
                geom_data = data[risk]
                if "type" in geom_data:
                    try:
                        geom = shape(geom_data)
                        geoms[risk].append(geom)
                    except Exception as e:
                        print(f"Error processing geometry in {file_path} for {risk}: {e}")
    return geoms


# def merge_geometries(geojson_folder):
#     """
#     Scan the folder for GeoJSON files and collect geometries for each risk level.
#     """
#     high_risk_geoms = []
#     medium_risk_geoms = []
#     low_risk_geoms = []
#     files = glob.glob(os.path.join(geojson_folder, "*.geojson"))
#     print(f"Found {len(files)} GeoJSON files.")
#     for file_path in files:
#         geoms = load_geometries_from_file(file_path)
#         high_risk_geoms.extend(geoms["High_risk"])
#         medium_risk_geoms.extend(geoms["Medium_risk"])
#         low_risk_geoms.extend(geoms["Low_risk"])
#     return high_risk_geoms, medium_risk_geoms, low_risk_geoms

################################################

logging.info("Reading geometry information from population data (will be used as reference)")
with rasterio.open(pop_path) as ref:
    transform = ref.transform
    width = ref.width
    height = ref.height
    crs = ref.crs
    bounds = ref.bounds
    resolution = ref.res

# logging.info("Merging GeoJSON geometries.")
# high_risk_geoms, medium_risk_geoms, low_risk_geoms = merge_geometries(geojson_folder)

logging.info("Preparing the output raster metadata.")
meta = {
       'driver': 'GTiff',
       'height': height,
       'width': width,
       'count': 1,
       'dtype': 'float32',
       'crs': crs,
       'transform': transform,
       'bigtiff': "YES",
       'compress': 'lzw'
   }

files = glob.glob(os.path.join(geojson_folder, "*.geojson"))
print(f"Found {len(files)} GeoJSON files.")

# Create empty raster first
with rasterio.open(output_path, 'w+', **meta) as dst:
    # Get the bounds of the destination raster
    dst_bounds = dst.bounds
    # Process one tile at a time
    for file_path in tqdm(files, desc="Processing tiles"):
        # Extract bounds from filename using regex
        filename = os.path.basename(file_path)
        coords = re.findall(r'[-]?\d+\.\d+', filename)
        if len(coords) >= 4:
            try:
                min_lat = float(coords[0])
                min_lng = float(coords[1])
                max_lat = float(coords[2])
                max_lng = float(coords[3])

                # Check if the tile intersects with our destination raster
                if (min_lng > dst_bounds.right or max_lng < dst_bounds.left or
                    min_lat > dst_bounds.top or max_lat < dst_bounds.bottom):
                    print(f"Tile {filename} outside raster bounds, skipping")
                    continue
                
                # Clip tile bounds to destination raster bounds
                min_x = max(min_lng, dst_bounds.left)
                min_y = max(min_lat, dst_bounds.bottom)
                max_x = min(max_lng, dst_bounds.right)
                max_y = min(max_lat, dst_bounds.top)
                
                # Create window from bounds with clipped coordinates
                window = rasterio.windows.from_bounds(
                    min_x, min_y, max_x, max_y, dst.transform)
                
                # Ensure window has valid dimensions
                if window.width <= 0 or window.height <= 0:
                    print(f"Window for {filename} has invalid dimensions: {window}, skipping")
                    continue
                
                # Round window to integers
                window = window.round_lengths().round_offsets()
                
                # Load geometries only for this file/tile
                geoms = load_geometries_from_file(file_path)
                
                # Skip if no geometries for this tile
                if not any(len(geoms[k]) > 0 for k in geoms):
                    print(f"No geometries found in {filename}, skipping")
                    continue
                
                # Process the window
                win_transform = rasterio.windows.transform(window, dst.transform)

                # Read existing data from the window
                window_data = dst.read(1, window=window)
                
                # Rasterize risk levels
                raster_high = rasterio.features.rasterize(
                    ((geom, 0.05) for geom in geoms["High_risk"]),
                    out_shape=(window.height, window.width),
                    transform=win_transform,
                    fill=0,
                    dtype='float32'
                    )
                raster_medium = rasterio.features.rasterize(
                    ((geom, 0.01) for geom in geoms["Medium_risk"]),
                    out_shape=(window.height, window.width),
                    transform=win_transform,
                    fill=0,
                    dtype='float32'
                    )
                raster_low = rasterio.features.rasterize(
                    ((geom, 0.005) for geom in geoms["Low_risk"]),
                    out_shape=(window.height, window.width),
                    transform=win_transform,
                    fill=0,
                    dtype='float32'
                    )

                # Combine new data
                combined = np.maximum(raster_low, np.maximum(raster_medium, raster_high))

                # Merge with existing data (taking the maximum at each pixel)
                combined = np.maximum(window_data, combined)

                # Write to combined window
                dst.write(combined, 1, window=window)
                    
            except (ValueError, IndexError) as e:
                print(f"Error parsing coordinates from {filename}: {e}")


# logging.info("Rasterizing each risk level (using rasterio windowed analysis).")
# with rasterio.open(output_path, 'w', **meta) as dst:
#     for ji, window in tqdm(dst.block_windows(1)):
#         # Calculate window's transform
#         win_transform = rasterio.windows.transform(window, transform)

#         # Define bounding box of this window (in the CRS of the raster)
#         left, top = win_transform * (0,0)
#         right, bottom = win_transform * (window.width, window.height)
#         window_bbox = box(left, bottom, right, top)

#         # Filter geometries for this window
#         high_win = [geom for geom in high_risk_geoms if geom.intersects(window_bbox)]
#         med_win = [geom for geom in medium_risk_geoms if geom.intersects(window_bbox)]
#         low_win = [geom for geom in low_risk_geoms if geom.intersects(window_bbox)]

#         # Rasterize each risk level in the window
#         # High risk: 0.05
#         raster_high_win = rasterio.features.rasterize(
#             ((geom, 0.05) for geom in high_win),
#             out_shape=(window.height, window.width),
#             transform=win_transform,
#             fill=0,
#             dtype='float32'
#         )
#         # Medium risk: 0.01
#         raster_med_win = rasterio.features.rasterize(
#             ((geom, 0.01) for geom in med_win),
#             out_shape=(window.height, window.width),
#             transform=win_transform,
#             fill=0,
#             dtype='float32'
#         )
#         # Low risk: 0.005
#         raster_low_win = rasterio.features.rasterize(
#             ((geom, 0.005) for geom in low_win),
#             out_shape=(window.height, window.width),
#             transform=win_transform,
#             fill=0,
#             dtype='float32'
#         )

#         # Combine by taking maximum value per cell
#         win_combined = np.maximum(raster_low_win, np.maximum(raster_med_win, raster_high_win))

#         # Write to combined window
#         dst.write(win_combined, 1, window=window)

logging.info("Done.")

