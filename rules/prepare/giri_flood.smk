"""
Prepare (clip) GIRI Flood Data
"""

rule align_giri_flood:
    """
    Fixes alignment issue (GIRI is half a cell misaligned with JRC). Note: temporary fix (takes forever)
    """
    input:
        raw_flood_file="data/inputs/flood/GIRI/global_pc_h{RP}glob.tif",
    output:
        aligned_flood_file="data/inputs/flood/GIRI/aligned/global_pc_h{RP}glob.tif",
    wildcard_constraints:
        RP="5|10|25|50|100|200|500|1000"
    shell:
        """
        gdal_translate -a_ullr -179.9995833333333333 90.0004166666666667 180.0004166666666667 -89.9995833333333333 \
               -co "COMPRESS=LZW" \
               -co "TILED=YES" \
               -co "BIGTIFF=YES" \
               {input.raw_flood_file} {output.aligned_flood_file}
        """

rule clip_giri_flood:
    """
    Clip GIRI flood raster to country boundary.
    """
    input:
        raw_flood_file="data/inputs/flood/GIRI/aligned/global_pc_h{RP}glob.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_flood_file="data/inputs/analysis/{ISO3}/{ISO3}_giri-flood_RP{RP}.tif",
    wildcard_constraints:
        RP="5|10|25|50|100|200|500|1000"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flood_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co compress=lzw \
            {input.raw_flood_file} \
            {output.trimmed_flood_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_giri-flood_RP10.tif
"""

