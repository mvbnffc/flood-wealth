"""
Rule book for extracting misc stats in CSV format.
"""

rule pop_admin_stats:
    """
    This rule returns a CSV with population statistics at the chosen admin level
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        quintile_file="data/results/social_flood/countries/{ISO3}/map_layers/{ISO3}_wealth_quintiles.tif",
    output:
        pop_stats="data/results/temp/stats/countries/{ISO3}/map_layers/{ISO3}_{ADMIN_SLUG}_pop_quintiles.csv",
    script:
        "./admin_pop_stats.py"
"""
Test with
snakemake -c1 data/results/temp/stats/countries/KEN/map_layers/KEN_ADM1_pop_quintiles.csv
"""