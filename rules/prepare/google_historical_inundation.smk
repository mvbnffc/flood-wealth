"""
Merge, rasterize, and clip the Google historical inundation data
"""

rule merge_and_rasterize_google_inundation:
    input:
        raw_folder = "data/inputs/google_inun/raw/",
        pop_path = "data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    output:
        merged_file="data/inputs/google_inun/merged/google_historical_inun.tif"
    script:
        "./rasterize_google_inun.py"