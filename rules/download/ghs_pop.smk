"""
Download GHSL Population Data

References
---------
GHS-POP: https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php
"""

rule download_ghsl:
    output:
        "data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_{RESOLUTION}_V1_0.tif"
    wildcard_constraints:
        RESOLUTION="3ss|30ss"
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_4326_{wildcards.RESOLUTION}/V1-0/GHS_POP_E2020_GLOBE_R2023A_4326_{wildcards.RESOLUTION}_V1_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_POP_E2020_GLOBE_R2023A_4326_{wildcards.RESOLUTION}_V1_0.zip \
            -d $output_dir
        """