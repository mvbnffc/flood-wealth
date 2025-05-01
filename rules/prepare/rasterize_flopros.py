"""
Rasterize the merged protection layer from the FLOPROS shapefile
Going to use the RWI dataset as a reference grid 
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
        flopros_path: str = snakemake.input["flopros"]
        pop_path: str = snakemake.input["pop_file"]
        output_path: str = snakemake.output["rasterized_flopros"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Rasterizing the FLOPROS dataset.")

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

logging.info("Reading FLOPROS shapefile. Building a list of shapes.")
with fiona.open(flopros_path) as shp:
    src_crs = shp.crs
    project = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True).transform
    shapes = [
        (mapping(shp_transform(project, shape(feat["geometry"]))),
         feat["properties"]["MerL_Riv"])
        for feat in shp
        if feat["geometry"] is not None
            and  feat["properties"]["MerL_Riv"] is not None
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

# # Function for tiled rasterization
# def stream_shapes(shp_path):
#     # this is a generator, so we never build the full list in memory
#     with fiona.open(shp_path) as shp:
#         for feat in shp:
#             yield feat["geometry"], feat["properties"]["MerL_Riv"]

# logging.info("Rasterizing the FLOPROS dataset.")

# logging.info("Reading reference Population data.")
# with rasterio.open(pop_path) as ref:
#     meta = ref.meta.copy()
#     meta.update({
#         "count": 1,
#         "dtype": "float32",
#         "compress": "lzw",
#         "BIGTIFF": "YES",
#     })

# with 

# logging.info("Writing output raster in blocks.")
# with rasterio.open(output_path, "w", **meta) as dst:
#     shapes = stream_shapes(flopros_path)
#     # iterate over each “block” (tile) that the GeoTIFF is going to write
#     for _, window in dst.block_windows(1):
#         # window_transform maps pixel‐coords within this block to CRS
#         win_transform = dst.window_transform(window)
#         # rasterize just this block
#         tile = rasterize(
#             shapes,
#             out_shape=(window.height, window.width),
#             transform=win_transform,
#             fill=0,
#             dtype='float32'
#         )
#         # write it back
#         dst.write(tile, window=window, indexes=1)

# logging.info("Done.")