"""
Prepare (clip) JRC Flood Data
"""

rule clip_jrc_flood:
    """
    Clip JRC flood raster to country boundary. 
    """
    input:
        raw_flood_file="data/inputs/flood/JRC/merged/jrc_global_flood_RP{RP}.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_flood_file="data/inputs/analysis/{ISO3}/{ISO3}_jrc_flood_RP{RP}.tif",
    wildcard_constraints:
        RP="10|20|50|75|100|200|500"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flood_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -dstalpha \
            -of GTiff \
            -co compress=lzw \
            {input.raw_flood_file} \
            {output.trimmed_flood_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_jrc_flood_RP10.tif
"""