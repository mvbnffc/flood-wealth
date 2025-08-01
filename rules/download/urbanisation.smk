"""
Download GHSL Degree of Urbanisation Data
Vector (urban) and raster (ghs-mod)

References
---------
urban: https://human-settlement.emergency.copernicus.eu/download.php?ds=DUC
ghs-mod: https://human-settlement.emergency.copernicus.eu/download.php?ds=smod
"""

rule download_urban:
    """
    Download admin-level urbanization dataset (will be used in adaptation scenarios)
    """
    output:
        urban_file = "data/inputs/ghs-duc/GHS_DUC_GLOBE_R2023A_V2_0.xlsx"
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL//GHS_DUC_GLOBE_R2023A/V2-0/GHS_DUC_MT_GLOBE_R2023A_V2_0.zip --directory-prefix=$output_dir

        unzip -o "$output_dir/GHS_DUC_MT_GLOBE_R2023A_V2_0.zip" -d "$output_dir"
        """

rule download_ghs_mod:
    """
    Download gridded urbanization dataset (will be used in concentration index decomposition)
    """
    output:
        "data/inputs/ghs-mod/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.tif"
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000/V2-0/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.zip \
            -d $output_dir
        """


