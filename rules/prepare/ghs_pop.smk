"""
Prepare (clip) GHSL Population Data
"""

rule clip_ghs_pop:
    """
    Clip GHS-POP raster to country boundary. 
    """
    input:
        raw_pop_file="data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_{RESOLUTION}_V1_0.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_pop_file="data/inputs/analysis/{ISO3}/ghs-pop_{RESOLUTION}.tif",
    wildcard_constraints:
        RESOLUTION="3ss|30ss"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_pop_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -dstalpha \
            -of GTiff \
            -co compress=lzw \
            {input.raw_pop_file} \
            {output.trimmed_pop_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/ghs-pop_3ss.tif
"""