"""
Download GHSL built up area and volume data and then extract residential raster

References
---------
GHS-POP: https://human-settlement.emergency.copernicus.eu/ghs_buS2023.php
"""

rule download_ghs_built:
    output:
        a_total = "data/inputs/ghs-built/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        a_nres = "data/inputs/ghs-built/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        v_total = "data/inputs/ghs-built/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        v_nres = "data/inputs/ghs-built/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    shell:
        """
        output_dir=$(dirname {output.a_total})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir
        
        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_V_GLOBE_R2023A/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir
        
        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_V_GLOBE_R2023A/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        
        unzip -o $output_dir/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir

        unzip -o $output_dir/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        
        unzip -o $output_dir/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        """

rule extract_res_layer:
    """
    subtract total bu layer from nres to get res layer
    """
    input:
        a_total = "data/inputs/ghs-built/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        a_nres = "data/inputs/ghs-built/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
        v_total = "data/inputs/ghs-built/GHS_BUILT_V_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        v_nres = "data/inputs/ghs-built/GHS_BUILT_V_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    output:
        a_res = "data/inputs/ghs-built/GHS_BUILT_S_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
        v_res = "data/inputs/ghs-built/GHS_BUILT_V_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    resources:
        mem_mb=8000  # Adjust based on your system memory
    shell:
        """
        gdal_calc.py -A {input.a_total} -B {input.a_nres} --outfile={output.a_res} \
        --calc="A-B" --NoDataValue=-9999 \
        --co="TILED=YES" --co="BLOCKSIZE=512" \
        --co="COMPRESS=LZW" --co="BIGTIFF=YES" \
        --overwrite --quiet

        gdal_calc.py -A {input.v_total} -B {input.v_nres} --outfile={output.v_res} \
        --calc="A-B" --NoDataValue=-9999 \
        --co="TILED=YES" --co="BLOCKSIZE=512" \
        --co="COMPRESS=LZW" --co="BIGTIFF=YES" \
        --overwrite --quiet
        """