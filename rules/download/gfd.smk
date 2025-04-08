"""
Download Global Flood Database (GFD) flood maps

References
---------
GFD: https://global-flood-database.cloudtostreet.ai/
"""

rule download_gfd:
    output:
        gfd_folder=directory("data/inputs/gfd/raw")
    shell:
        """
        mkdir -p {output.gfd_folder}
        gsutil -m cp -r gs://gfd_v1_4/* {output.gfd_folder}/
        """