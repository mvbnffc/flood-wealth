"""
This script creates an (adapted) river protection GeoTiff to use in the risk analysis.
This river protection GeoTiff is created by specifying a (minimum) degree of urbanization to protect
and a return period of flooding to protect agains. The cost of adaptation measures is calculated
by calculating the length of river network within each region to be protected and applying a
formula from Boulange et al (2023) to calculate the cost of adaptation measures. 

Outputs are a river protection GeoTiff and a text file with the cost of adaptation measures (low, mid, and upper bound estimates).

Reference: Boulange et al (2023) https://link.springer.com/article/
"""

import logging
import sys
import glob
import os

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.mask import mask
import geopandas as gpd
from collections import Counter
from pyproj import Geod
from shapely.geometry import LineString, MultiLineString

if __name__ == "__main__":
    try:
        flopros_path: str = snakemake.input["flopros"]
        river_network_path: str = snakemake.input["river_network"]
        urbanization_path: str = snakemake.input["urbanization"]
        output_path: str = snakemake.output["flood_protection"]
        cost_output_path: str = snakemake.output["cost_of_protection"]
        country: str = snakemake.wildcards["ISO3"]
        rp: str = snakemake.wildcards["RP"]
        urban_class: str = snakemake.wildcards["urban_class"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Creating a flood protection layer for {country} and estimating costs. Return period: {rp}, Urbanization class: {urban_class}.")

logging.info("Creating output directories (if they doesn't already exist)")
out_dir = os.path.dirname(output_path)
cost_out_dir = os.path.dirname(cost_output_path)
os.makedirs(out_dir, exist_ok=True)
os.makedirs(cost_out_dir, exist_ok=True)

logging.info("Load files.")
river_network = gpd.read_file(river_network_path)
urbanization = gpd.read_file(urbanization_path)
# Make sure both layers are in same CRS
river_network = river_network.to_crs(urbanization.crs)
# Extract the GID level for this dataset
gadm_level = urbanization.columns[0].split('_')[-1]  # Assuming it is in first column
logging.info(f"GADM level for urbanization data is level {gadm_level}.")

logging.info("Calculate the level of flood protection within each admin region.")
modes = [] # collect the modes for each region in this list
with rasterio.open(flopros_path) as src:
    for idx, row in urbanization.iterrows():
        try:
            # Extract raster values within polygon
            out_image, out_transform = mask(src, [row.geometry], crop=True, nodata=src.nodata)
            # Get non-nodata values
            valid_data = out_image[out_image != src.nodata]
            valid_data = valid_data[~np.isnan(valid_data)]

            if valid_data.size > 0:
                # Find mode 
                counter = Counter(valid_data.flatten())
                mode_value = counter.most_common(1)[0][0]
            else:
                mode_value = 0
            urbanization.at[idx, 'FLOPROS'] = mode_value

            modes.append(mode_value)
        except Exception as e:
            logging.error(f"Error processing row {idx} for GID {row[f'GID_{gadm_level}']}: {e}")
            modes.append(0)

logging.info("Clip rivers to admin regions.")
rivers_clipped = gpd.overlay(river_network, urbanization, how='intersection')

logging.info("Calculate the length of each river segment")
geod = Geod(ellps="WGS84") # Using WGS84 ellipsoid for geodesic calculations

# Function for calculating the geodetic length of a geometry
def geodetic_length(geom):
    total = 0.0
    lines = [geom] if isinstance(geom, LineString) else geom.geoms
    for line in lines:
        lons, lats = line.xy
        for i in range(len(lons) - 1):
            # geod.inv returns the forward azimuth, back azimuth, and distance
            total += geod.inv(lons[i], lats[i], lons[i + 1], lats[i + 1])[2]
    return total

# Compute the length of each clipped segment (in meters)
rivers_clipped['length_m'] = rivers_clipped.geometry.apply(geodetic_length)

logging.info("Sum river lengths per admin region.")
lengths_by_admin = (
    rivers_clipped
    .groupby(f"GID_{gadm_level}", as_index=False)["length_m"]
    .sum()
    .rename(columns={"length_m": "river_length_m"})
)

logging.info("Merge river lengths with urbanization data.")
urbanization = urbanization.merge(lengths_by_admin, on=f"GID_{gadm_level}", how="left")
urbanization['river_length_m'] = urbanization['river_length_m'].fillna(0)  # Fill NaN with 0

logging.info("Calculate the ADMIN level costs of adaptation measures.")
# Formula for cost of adaptation is C = unit cost * River Length (km) * log2(∆ flood protection) Source: Boulange et al (2023)
# 3 values from supplementary material of that paper (million USD)
min_unit_cost = 0.7996
mid_unit_cost = 2.399
max_unit_cost = 7.196

# Loop through each admin region and calculate costs
for idx, row in urbanization.iterrows():
        if row['DEGURBA_L2'] < int(urban_class):
            urbanization.at[idx, 'adaptation_cost_min'] = 0
            urbanization.at[idx, 'adaptation_cost_mid'] = 0
            urbanization.at[idx, 'adaptation_cost_max'] = 0
        else:
            # Calculate the delta in flood protection
            delta_fp = float(rp) - row['FLOPROS']
            if delta_fp <= 0:
                urbanization.at[idx, 'adaptation_cost_min'] = 0
                urbanization.at[idx, 'adaptation_cost_mid'] = 0
                urbanization.at[idx, 'adaptation_cost_max'] = 0
            else:
                # Calculate costs based on river length (in km) and unit costs
                river_length_km = row['river_length_m'] / 1000.0
                urbanization.at[idx, 'adaptation_cost_min'] = min_unit_cost * river_length_km * np.log2(delta_fp)
                urbanization.at[idx, 'adaptation_cost_mid'] = mid_unit_cost * river_length_km * np.log2(delta_fp)
                urbanization.at[idx, 'adaptation_cost_max'] = max_unit_cost * river_length_km * np.log2(delta_fp)  

logging.info("Sum costs across admin regions and save to text file.")
costs = urbanization[['GID_' + gadm_level, 'adaptation_cost_min', 'adaptation_cost_mid', 'adaptation_cost_max']].sum()
costs.index = ['GID_' + gadm_level, 'adaptation_cost_min', 'adaptation_cost_mid', 'adaptation_cost_max']
with open(cost_output_path, 'w') as f:
    f.write(f"Min cost: {costs['adaptation_cost_min']:.2f} million USD\n")
    f.write(f"Mid cost: {costs['adaptation_cost_mid']:.2f} million USD\n")
    f.write(f"Max cost: {costs['adaptation_cost_max']:.2f} million USD\n")

logging.info("Costs of adaptation:")
logging.info(f"Minimum cost: {costs['adaptation_cost_min']:.2f} million USD")
logging.info(f"Mid cost: {costs['adaptation_cost_mid']:.2f} million USD")
logging.info(f"Maximum cost: {costs['adaptation_cost_max']:.2f} million USD")

logging.info("Create a new FLOPROS layer with adaptation incorporated.")
# Create a new GeoDataFrame for the adapted FLOPROS layer
# Grab layers that fall within urbanization class, have a river inside them, and have FLOPROS below the specified return period
filtered_urbanization = urbanization[
    (urbanization['DEGURBA_L2'] >= int(urban_class)) &
    (urbanization['FLOPROS'] < float(rp)) &
    (urbanization['river_length_m'] > 0)
]
if len(filtered_urbanization) > 0:
    with rasterio.open(flopros_path) as src:
        adaptation_raster = src.read(1).astype(np.float32)  # start with a copy of original FLOPROS raster
        profile = src.profile.copy()

        # Create shapes for rasterization (geometry and value pairs)
        shapes = []
        for idx, row in filtered_urbanization.iterrows():
            shapes.append((row.geometry, float(rp)))

        # Burn the new RP values into the adaptation raster
        burned = rasterize(
            shapes,
            out_shape=adaptation_raster.shape,
            transform=src.transform,
            fill=0,
            default_value=float(rp),  # Use the specified return period as the default value
            dtype=np.float32
        )

        # Update adaptation raster (keeping original values where no new RP was burned)
        adaptation_raster = np.where(burned > 0, burned, adaptation_raster)

        # Update profile for output
        profile.update(dtype=rasterio.float32, compress='lzw')

        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(adaptation_raster, 1)
else:
    logging.info("No admin regions met the filtering criteria - write out original FLOPROS raster.")
    with rasterio.open(flopros_path) as src:
        with rasterio.open(output_path, 'w', **src.profile) as dst:
            dst.write(src.read(1), 1)

logging.info("Done.")

