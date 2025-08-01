"""
Rasterize the L2 DUC classifications from the urbanization geopackage
Going to use the population dataset as a reference grid 
"""

import logging

import rasterio
from rasterio.features import rasterize
from shapely.geometry import shape, mapping
from shapely.ops import transform as shp_transform
import pyproj
import fiona



if __name__ == "__main__":

    try:
        urbanization_path: str = snakemake.input["urbanization"]
        pop_path: str = snakemake.input["pop_file"]
        output_path: str = snakemake.output["rasterized"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Rasterizing the GHS urbanization dataset.")

logging.info("Reading reference Population data.")
with rasterio.open(pop_path) as ref:
    meta = ref.meta.copy()
    meta.update({
        "count": 1,
        "dtype": "float32",
        "compress": "lzw",
        "BIGTIFF": "YES",
    })
    dst_crs = ref.crs

logging.info("Reading urbanization geopackage. Building a list of shapes.")
with fiona.open(urbanization_path) as gpkg:
    src_crs = gpkg.crs
    project = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True).transform
    shapes = [
        (mapping(shp_transform(project, shape(feat["geometry"]))),
         feat["properties"]["DEGURBA_L2"])
        for feat in gpkg
        if feat["geometry"] is not None
            and  feat["properties"]["DEGURBA_L2"] is not None
    ]

logging.info("Writing output raster in blocks.")

with rasterio.open(output_path, "w", **meta) as dst:
    for _, window in dst.block_windows(1):
        win_transform = dst.window_transform(window)
        tile = rasterize(
            shapes,
            out_shape=(window.height, window.width),
            transform=win_transform,
            fill=0,
            dtype="float32"
        )
        dst.write(tile, 1, window=window)

logging.info("Done.")