"""
Download JRC Surface Water Dataset
that we will use as a permanent water mask.

References
---------
PERM: https://global-surface-water.appspot.com/
"""

rule download_surface_water:
    output:
        directory("data/inputs/jrc-perm/occurrence/")
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        download_water_data occurrence -d $output_dir
        """