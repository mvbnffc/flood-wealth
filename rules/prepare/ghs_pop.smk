"""
Prepare (clip) GHSL Population Data
"""

rule clip_ghs_pop:
    """
    Clip GHS-POP raster to country boundary. 
    """
    input:
        raw_pop_file="data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_pop_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_pop_file} \
            {output.trimmed_pop_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_ghs-pop.tif
"""