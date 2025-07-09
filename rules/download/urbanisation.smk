"""
Download GHSL Degree of Urbanisation Data

References
---------
urban: https://human-settlement.emergency.copernicus.eu/download.php?ds=DUC
"""

rule download_urban:
    output:
        urban_file = "data/inputs/ghs-duc/GHS_DUC_GLOBE_R2023A_V2_0.xlsx"
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL//GHS_DUC_GLOBE_R2023A/V2-0/GHS_DUC_MT_GLOBE_R2023A_V2_0.zip --directory-prefix=$output_dir

        unzip -o "$output_dir/GHS_DUC_MT_GLOBE_R2023A_V2_0.zip" -d "$output_dir"
        """


rule download_hydrorivers_poo:
    output:
        river_file="data/inputs/rivers/hydroRIVERS/HydroRIVERS_v10.shp"
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        wget -nc https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_shp.zip \
            --directory-prefix="$output_dir"

        unzip -o "$output_dir/HydroRIVERS_v10_shp.zip" -d "$output_dir"

        mv "$output_dir/HydroRIVERS_v10_shp/"* "$output_dir/"
        rmdir "$output_dir/HydroRIVERS_v10_shp"
        """
