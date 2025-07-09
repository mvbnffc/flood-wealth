"""
Download River Network Data

References
---------
hydroRIVERS: https://www.hydrosheds.org/products/hydrorivers
"""

rule download_hydrorivers:
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
