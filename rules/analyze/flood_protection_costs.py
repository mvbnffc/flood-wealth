"""
This script calculates the cost of river flood protection at the admin scale, using the 
length of river to be protected in every urban area as well as the flood protection delta
(relative to baseline protection)
"""

import logging

import rasterio
from rasterio.mask import mask
from pyproj import Geod
from shapely.geometry import LineString, MultiLineString
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely
from tqdm import tqdm

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        flopros_path: str = snakemake.input["flopros"]
        river_path: str = snakemake.input["rivers"]
        urban_path: str = snakemake.input["urban"]
        gdppc_path: str = snakemake.input["gdppc"]
        output_path: str = snakemake.output["protection_cost"]
        urban_class: str = snakemake.wildcards.urban_class
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        RP: str = snakemake.wildcards.RP
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))

logging.info(f"Calculating river flood protection costs for {country} at admin level {admin_level}.")
logging.info(f"Minimum urban area to be protected is DUC{urban_class} to RP {RP}")

# -------------------HELPER FUNCTIONS -----------------------------------

def calculate_geod_length(geom):  # Changed parameter name from 'line' to 'geom'
    '''
    Function to calculate the geodetic length of a LineString or MultiLineString
    '''
    geod = Geod(ellps="WGS84") # our data is in WGS84 projection
    length = 0

    if isinstance(geom, LineString):  
        for i in range(len(geom.coords)-1):
            lon1, lat1 = geom.coords[i]  
            lon2, lat2 = geom.coords[i + 1] 
            _, _, distance = geod.inv(lon1, lat1, lon2, lat2)
            length += distance
    elif isinstance(geom, MultiLineString): 
        # Handle MultiLineString - iterate through each component LineString
        for line_segment in geom.geoms:  
            if isinstance(line_segment, LineString):
                for i in range(len(line_segment.coords)-1):
                    lon1, lat1 = line_segment.coords[i]
                    lon2, lat2 = line_segment.coords[i + 1]
                    _, _, distance = geod.inv(lon1, lat1, lon2, lat2)
                    length += distance
    else:
        # Handle unexpected geometry types
        print(f"Warning: Unexpected geometry type {type(geom)}")
        return 0

    return length / 1000 # convert to km

def zonal_mode_for_polygon(src, geom):
    '''
    Function for raster stats per polygon (will calculate mode - for FLOPROS)
    '''
    # Mask raster to polygon, crop to bbox for speed
    out, _ = mask(src, [geom], crop=True, filled=False)
    
    # Handle masked array properly
    masked_arr = out[0]
    
    # Convert masked array to regular array, getting only valid (unmasked) values
    if hasattr(masked_arr, 'compressed'):
        # This extracts only the unmasked values
        arr = masked_arr.compressed()
    else:
        # Fallback if not a masked array
        arr = masked_arr[~np.isnan(masked_arr)]
    
    # Remove explicit nodata values if they exist
    if src.nodata is not None:
        arr = arr[arr != src.nodata]
    
    if arr.size == 0:
        return 0  # Match Script 1 behavior

    # Compute mode (works for int or float)
    vals, counts = np.unique(arr, return_counts=True)
    return float(vals[np.argmax(counts)])

def calculate_protection_costs(geom, protection_level, country_cost_adjustment):
    """
    Calculate the cost of river flood protection.
    This is done using the river length in the admin region, the baseline
    level of protection and the new design protection level.
    We use the standard Boulange et al unit cost per km of river to be protected.
    We also adjust the unit costs per country based on the GDP pc relative to global average.
    We set a lower and upper limit of unit construction costs as per the Boulandge et al 
    supplementary.
    
    This function will return both the standard (globally consistent) cost of river flood protection
    and the country_adjusted one.
    """
    max_unit_cost = 7.196 # million
    standard_unit_cost = 2.399 # million
    min_unit_cost = 0.7966 # million

    river_length = geom['riv_len_km']
    baseline_protection = geom['FLOPROS']

    # Protection delta
    delta_fp = np.max((protection_level-baseline_protection), 0)

    # Guard against log2(0) or negatives
    if delta_fp <= 0:
        adaptation_cost = 0.0
        adj_adaptation_cost = 0.0
    else:
        adaptation_cost = standard_unit_cost * river_length * np.log2(delta_fp)
        adj_adaptation_cost = adjusted_unit_cost * river_length * np.log2(delta_fp)

    # Adjusted adaptation cost
    adjusted_unit_cost = country_cost_adjustment * standard_unit_cost
    if max_unit_cost < adjusted_unit_cost:
        adjusted_unit_cost = max_unit_cost
    elif min_unit_cost > adjusted_unit_cost:
        adjusted_unit_cost = min_unit_cost
    
    adj_adaptation_cost = adjusted_unit_cost * river_length * np.log2(delta_fp)

    return adaptation_cost, adj_adaptation_cost

# -----------------------------------------------------------

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = "shapeName"
admin_areas = admin_areas[[area_unique_id_col, "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

logging.info(f"Reading the river network data")
rivers = gpd.read_file(river_path)

logging.info("Reading urbanization data")
urban = gpd.read_file(urban_path)


logging.info("Calculating the flood protection per urban area")
with rasterio.open(flopros_path) as src:
    # Ensure geometries are in the raster CRS
    urban_in_raster_crs = urban.to_crs(src.crs)
    modes = [] # will be calculating most common (mode) FLOPROS value per region
    for geom in tqdm(urban_in_raster_crs.geometry, desc="calculating urban FLOPROS"):
        modes.append(zonal_mode_for_polygon(src, geom))
urban['FLOPROS'] = modes

logging.info("Calculating the river length per urban area")
urban_areas = urban[urban['DEGURBA_L2'] >= int(urban_class)] # Filter only by urban classes we are interested in
# Intersect rivers with the urbanization layer and admin layer (nested intersection)
urban_admin_rivers = gpd.overlay(rivers, gpd.overlay(urban_areas, admin_areas, how='intersection'), how='intersection')
urban_admin_rivers['riv_len_km'] = urban_admin_rivers['geometry'].apply(calculate_geod_length)

logging.info("Calculating the protection cost per urban area")
# Get the GDP adjustment
gdp_df = pd.read_csv(gdppc_path)
country_gdp = gdp_df.loc[gdp_df["ISO"].eq(country), "GDP_pc"].iloc[0]
global_gdp  = gdp_df.loc[gdp_df["ISO"].eq("WORLD"), "GDP_pc"].iloc[0]
gdp_adjustment = country_gdp / global_gdp
# Calculate the protection cost for each urban area
urban_admin_rivers[['adaptation_cost', 'adj_adaptation_cost']] = (
    urban_admin_rivers.apply(
        lambda r: pd.Series(
            calculate_protection_costs(
                r, protection_level=int(RP), country_cost_adjustment=gdp_adjustment
        )
    ), axis = 1 
    )
)

logging.info("Sum the protection costs per Admin region")
# Group by admin region and sum sugment lengths
costs_per_admin = urban_admin_rivers.groupby('shapeName')['adaptation_cost'].sum().reset_index()
adj_costs_per_admin = urban_admin_rivers.groupby('shapeName')['adj_adaptation_cost'].sum().reset_index()
# Add total costs the the admin area dataframe
admin_areas = admin_areas.merge(costs_per_admin, how='left', left_on='shapeName', right_on='shapeName')
admin_areas = admin_areas.merge(adj_costs_per_admin, how='left', left_on='shapeName', right_on='shapeName')

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(admin_areas, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

# Debug
print(results_gdf['adaptation_cost'].sum())

logging.info("Done.")