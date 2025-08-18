"""
Prepare (clip) GHSL Built (residential) and also (non-residential). Area and volume.
"""

rule clip_ghs_res:
    """
    Clip GHS-RES and GHS-NRES raster to country boundary. 
    """
    input:
        raw_res_a_file="data/inputs/ghs-built/GHS_BUILT_S_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        raw_nres_a_file="data/inputs/ghs-built/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        raw_res_v_file="data/inputs/ghs-built/GHS_BUILT_V_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        raw_nres_v_file="data/inputs/ghs-built/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_res_a_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_a.tif",
        trimmed_nres_a_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-nres_a.tif",
        trimmed_res_v_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_v.tif",
        trimmed_nres_v_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-nres_v.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_res_a_file})
        mkdir --parents $(dirname {output.trimmed_nres_a_file})
        mkdir --parents $(dirname {output.trimmed_res_a_file})
        mkdir --parents $(dirname {output.trimmed_nres_a_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_res_a_file} \
            {output.trimmed_res_a_file}

        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_nres_a_file} \
            {output.trimmed_nres_a_file}

        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_res_v_file} \
            {output.trimmed_res_v_file}

        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -co compress=lzw \
            {input.raw_nres_v_file} \
            {output.trimmed_nres_v_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/KEN/KEN_ghs-res_a.tif
"""
