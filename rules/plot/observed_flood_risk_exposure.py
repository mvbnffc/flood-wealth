"""
This script creates a observed flood risk GeoTiff for plotting / mapping.
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        pop_path: str = snakemake.input["pop_file"]
        risk_path: str = snakemake.input["risk_file"]
        output_path: str = snakemake.output["risk_exposure"]
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Creating observed flood risk exposure GeoTiff in {country}.")

logging.info("Load the raster data.")
pop = rasterio.open(pop_path)
risk = rasterio.open(risk_path)
# Get raster info
profile = pop.meta.copy()
profile.update(dtype=rasterio.float32, compress='lzw', nodata=0, count=1)
nodata = pop.nodata

logging.info("Calculating exposure (windowed analysis) and writing to file")
with rasterio.open(output_path, 'w', **profile) as dst:
    i = 0
    for ji, window in pop.block_windows(1):
        i += 1

        affine = rasterio.windows.transform(window, pop.transform)
        height, width = rasterio.windows.shape(window)
        bbox = rasterio.windows.bounds(window, pop.transform)

        profile.update({
            'height': height,
            'width': width,
            'affine': affine
        })

        pop_array = pop.read(1, window=window)
        risk_array = risk.read(1, window=window)
        # Prepare arrays (make sure all non-positive cells are zero)
        pop_array = np.where(pop_array>0, pop_array, 0)
        risk_array = np.where(risk_array>0, risk_array, 0)
        risk_exposure = pop_array * risk_array # multiplying the two datasets

        # Write to file
        dst.write(risk_exposure.astype(rasterio.float32), window=window, indexes=1)

logging.info("Done.")