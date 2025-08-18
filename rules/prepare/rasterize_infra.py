"""
Rasterize the OSM road data for a given country.
Will output an infrastructure density layer to approximate infra losses. 
Going to use the population dataset as a reference grid 
"""

import logging
import os

import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
import rasterio.enums as rio_enums
from shapely.geometry import LineString, MultiLineString
from shapely.ops import transform as shp_transform
import geopandas as gpd
import numpy as np
import pyproj
import fiona



if __name__ == "__main__":

    try:
        osm_path: str = snakemake.input["osm_folder"]
        pop_path: str = snakemake.input["pop_file"]
        output_path: str = snakemake.output["infrastructure_raster"]
        country: str = snakemake.wildcards["ISO3"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Rasterizing the OSM infrastructure layer for {country}.")

logging.info("Reading reference Population data.")
with rasterio.open(pop_path) as ref:
    meta = ref.meta.copy()
    meta.update({
        "count": 1,
        "dtype": "float32",
        "compress": "lzw",
        "BIGTIFF": "YES",
    })
    height, width = ref.height, ref.width
    transform: Affine = ref.transform
    res_x  = float(abs(transform.a))
    res_y = float(abs(transform.e))
    dst_crs = ref.crs

logging.info("Reading road network data.")
road_file = os.path.join(osm_path, "gis_osm_roads_free_1.shp")
roads = gpd.read_file(road_file)
roads = roads[['osm_id', 'fclass', 'name', 'geometry']]
roads = roads[roads.geometry.notnull()]  # Remove rows with null geometries
roads = roads[roads.geometry.type.isin(['LineString', 'MultiLineString'])]  # Keep only LineString and MultiLineString
if roads.crs is None:
    roads = roads.set_crs(dst_crs)
if roads.crs != dst_crs:
    roads = roads.To_crs(dst_crs)
# Filter road classes
keep = {'motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'residential', 'unclassified', 'service'}
roads = roads[roads['fclass'].isin(keep)]

# Densify road lines to ~1/2 pixel so we can rasterize them properly
step_deg = 0.5 * min(res_x, res_y)  # Half pixel step in degrees
set_deg = max(step_deg, 1e-6)  # Ensure step is not zero
logging.info("Densifying road geometries to half pixel")

def densify_line(line, step):
    L = line.length
    if L == 0:
        return []
    n = max(int(np.ceil(L / step)), 1)  # Ensure at least one point
    dists = np.linspace(0, L, n+1)
    segs = []
    for i in range(n):
        p0 = line.interpolate(dists[i])
        p1 = line.interpolate(dists[i+1])
        if not p0.equals(p1):
            segs.append(LineString([p0, p1]))
    return segs

def explode_and_densify(geom, step):
    if isinstance(geom, LineString):
        return densify_line(geom, step)
    elif isinstance(geom, MultiLineString):
        out = []
        for ln in geom.geoms:
            out.extend(densify_line(ln, step))
        return out
    return []

# Rasterization parameters. Will be countring segments per pixel (as a proxy for road length). The normalize by pixel area (to account for issues around projections)
acc = np.zeros((height, width), dtype=np.float32)
CHUNK_SEGS = 200_000  # Number of segments to process in each chunk
buf = []  # Buffer to hold segments

logging.info("Rasterizing road geometries")
for geom in roads.geometry:
    for seg in explode_and_densify(geom, step_deg):
        # Burn a constant (1) per short segment
        buf.append((seg, 1.0))
        if len(buf) >= CHUNK_SEGS:
            # Rasterize the current buffer
            tile = rasterize(
                shapes=buf,
                out_shape=(height, width),
                transform=transform,
                fill=0.0,
                dtype="float32",
                all_touched=True,
                merge_alg=rio_enums.MergeAlg.add  # Use add to accumulate values
            )
            acc += tile
            buf.clear()

if buf:
    tile = rasterize(
        shapes=buf,
        out_shape=(height, width),
        transform=transform,
        fill=0.0,
        all_touched=True,
        dtype="float32",
        merge_alg=rio_enums.MergeAlg.add  # Use add to accumulate values    
    )
    acc += tile
    buf.clear()

approx_deg_length = acc * step_deg # Convert to approximate length in degrees

logging.info("Area normalization") # reduce latitude bias if in WGS84
row_idx = np.arange(height, dtype='float64')
lat_center = transform.f + (row_idx + 0.5) * transform.e
coslat = np.cos(np.deg2rad(lat_center))
coslat = np.clip(coslat, 1e-6, None)  # Avoid division by zero
# Broadcast to full array
approx_density = approx_deg_length / coslat[:, None]

logging.info("Rescale to 0-1 range")
m = approx_density.max()
if m > 0:
    approx_density = approx_density / m  # Normalize to 0-1 range


logging.info("Writing output raster.")

meta.update(nodata=0.0, dtype="float32")
os.makedirs(os.path.dirname(output_path), exist_ok=True)
with rasterio.open(output_path, "w", **meta) as dst:
    dst.write(approx_density.astype(np.float32), 1)

logging.info("Done.")