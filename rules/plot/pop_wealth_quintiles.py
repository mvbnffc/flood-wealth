"""
This script creates a population wealth quintile map for plotting / mapping
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        pop_path: str = snakemake.input["pop_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        output_path: str = snakemake.output["wealth_quintiles"]
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Mapping population wealth quintiles in {country}.")

logging.info("Load the raster data.")
with rasterio.open(rwi_path) as src_rwi:
    rwi_data = src_rwi.read(1)
    rwi_nodata = src_rwi.nodata
with rasterio.open(pop_path) as src_pop:
    pop_data = src_pop.read(1)
    profile = src_pop.meta.copy()
    profile.update(dtype=rasterio.float32, compress='lzw', nodata=0, count=1)
    nodata = src_pop.nodata
# Ensure population and RWI data align in shape
assert rwi_data.shape == pop_data.shape, "RWI and Population rasters must match in shape."

logging.info('Removing NaN and NoData values')
# Remove NaN and NoData values
valid_mask = (~np.isnan(rwi_data)) & (rwi_data != rwi_nodata) & (pop_data > 0)
weighted_rwi = rwi_data[valid_mask]
pop_weight = pop_data[valid_mask]

logging.info('Calculate quintiles')
# Only want to work on valid data
rwi_valid = rwi_data[valid_mask]
population_valid = pop_data[valid_mask]
# Sort RWI values by population weight
sorted_indices = np.argsort(rwi_valid)
rwi_sorted = rwi_valid[sorted_indices]
population_sorted = population_valid[sorted_indices]
# Calculate cumulative population
cumulative_population = np.cumsum(population_sorted)
total_population = cumulative_population[-1]
# Determine quintile thresholds
quintile_thresholds = []
for q in [0.2, 0.4, 0.6, 0.8]:
    idx = np.searchsorted(cumulative_population, q * total_population)
    quintile_thresholds.append(rwi_sorted[idx])

logging.info("Applying population-weighted quintile thresholds to population map")
first_quintile = np.where(((rwi_data<=quintile_thresholds[0]) & (pop_data>0)), 1, 0)
second_quintile = np.where(((rwi_data>quintile_thresholds[0]) & (rwi_data<=quintile_thresholds[1]) & (pop_data>0)), 2, 0)
third_quintile = np.where(((rwi_data>quintile_thresholds[1]) & (rwi_data<=quintile_thresholds[2]) & (pop_data>0)), 3, 0)
fourth_quintile = np.where(((rwi_data>quintile_thresholds[2]) & (rwi_data<=quintile_thresholds[3]) & (pop_data>0)), 4, 0)
fifth_quintile = np.where(((rwi_data>quintile_thresholds[3]) & (pop_data>0)), 5, 0)
# Combine these into one map
quintile_map = first_quintile + second_quintile + third_quintile + fourth_quintile + fifth_quintile

logging.info("Writing wealth quintile raster.")
with rasterio.open(output_path, "w", **profile) as dst:
    dst.write(quintile_map.astype(np.int16), 1)

logging.info("Done.")