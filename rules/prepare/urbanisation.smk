"""
Prepare urbanisation dataset (combine GADM with CSV)
"""

rule prepare_urbanisation:
    """
    Merge the smallest GADM level with the GHS DUC CSV values on urbanisation
    """
    input:
        duc_info="data/inputs/ghs-duc/GHS_DUC_GLOBE_R2023A_V2_0.xlsx",
        adm_file="data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
    output:
        urban_file= "data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.gpkg"
    script:
        "./prepare_urban.py"