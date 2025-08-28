"""
Calculate adaptation costs at the admin level
"""

rule flood_protection_costs:
    """
    This rule calculate the cost of river flood protection using the length of the river to be protected 
    as well as the delta in flood protection (relative to baseline protection levels)
    """
    input:
        admin_areas="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        flopros="data/inputs/analysis/countries/{ISO3}/{ISO3}_flopros.tif",
        urban="data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.gpkg",
        rivers="data/inputs/analysis/countries/{ISO3}/{ISO3}_river_network.gpkg",
        gdppc="config/gdppc_data.csv"
    output:
        protection_cost="data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban_class}_{ADMIN_SLUG}.gpkg"
    wildcard_constraints:
        urban_class="11|12|13|21|22|23|30",
        ADMIN_SLUG="ADM0|ADM1|ADM2"
    script:
        "./flood_protection_costs.py"
""" 
Test with
snakemake -c1 data/results/adaptation/costs/countries/RWA/RWA_adaptation-cost_fp_rp100_duc30_ADM2.gpkg
"""

configfile: "config/config.yaml"
ADMINS = ["ADM1", "ADM2"]
RPs = [50, 100]
DUC_protection = [21, 22, 23, 30]

rule fP_costs_for_all_countries:
    input:
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban}_{ADM}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, RP=RPs, urban=DUC_protection)