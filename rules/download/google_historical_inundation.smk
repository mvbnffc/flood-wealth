"""
Download Google Inundation History Data

References
---------
ggin: https://support.google.com/flood-hub/answer/15636595?hl=en&ref_topic=15637286
"""

rule download_google_inundation:
    """
    Will download all files from relevant Google Cloud bucket
    """
    output:
        raw_folder = directory("data/inputs/google_inun/raw/")
    script:
        "./download_google_inun.py"