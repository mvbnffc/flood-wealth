"""
Prepare (clip) GHSL Built (residential) data
"""

rule clip_ghs_res:
    """
    Clip GHS-RES raster to country boundary. 
    """
    input:
        raw_res_file="data/inputs/ghs-built/GHS_BUILT_S_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_res_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_res_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_res_file} \
            {output.trimmed_res_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/KEN/KEN_ghs-res.tif
"""