"""
Download GHSL built up data and then extract residential raster

References
---------
GHS-POP: https://human-settlement.emergency.copernicus.eu/ghs_buS2023.php
"""

rule download_ghs_built:
    output:
        total = "data/inputs/ghs-built/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        nres = "data/inputs/ghs-built/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    shell:
        """
        output_dir=$(dirname {output.total})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir
        
        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss/V1-0/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        
        unzip -o $output_dir/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        """

rule extract_res_layer:
    """
    subtract total bu layer from nres to get res layer
    """
    input:
        total = "data/inputs/ghs-built/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        nres = "data/inputs/ghs-built/GHS_BUILT_S_NRES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    output:
        res = "data/inputs/ghs-built/GHS_BUILT_S_RES_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    resources:
        mem_mb=8000  # Adjust based on your system memory
    shell:
        """
        gdal_calc.py -A {input.total} -B {input.nres} --outfile={output.res} \
        --calc="A-B" --NoDataValue=-9999 \
        --co="TILED=YES" --co="BLOCKSIZE=512" \
        --co="COMPRESS=LZW" --co="BIGTIFF=YES" \
        --overwrite --quiet
        """