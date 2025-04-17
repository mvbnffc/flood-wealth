"""
Download FLOPROS flood protection data

References
---------
FLOPROS: https://nhess.copernicus.org/articles/16/1049/2016/
"""

rule download_flopros:
    """
    Download and extract the FLOPROS dataset.
    This rule outputs the shapefile needed for further processing.
    """
    output:
        shp="data/inputs/flopros/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp"
    shell:
        """
        # Define output directory as the correct parent directory
        output_dir=$(realpath $(dirname {output.shp})/../..)

        # Create the FLOPROS directory (without nesting issues)
        mkdir -p $output_dir

        # Download the dataset if it does not already exist
        wget -nc http://dx.doi.org/10.5194/nhess-16-1049-2016-supplement \
            --directory-prefix=$output_dir

        # Unzip the dataset into the correct location
        unzip -o $output_dir/nhess-16-1049-2016-supplement -d $output_dir
        """