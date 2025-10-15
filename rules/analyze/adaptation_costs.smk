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

rule relocation_costs:
    """
    This rule calculates the sub-national costs of relocation adaptation scenario.
    """
    input:
        admin_areas="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        flopros_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_flopros.tif",
        flood_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_{model}-flood_RP10.tif",
        pop_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        urbanization_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-mod.tif",
        res_area="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_a.tif",
        res_capstock="data/inputs/analysis/countries/{ISO3}/{ISO3}_res_capstock.tif"
    output:
        relocation_costs="data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_rl_m-{model}_duc{urban_class}_{ADMIN_SLUG}.gpkg"
    wildcard_constraints:
        urban_class="11|12|13|21|22|23|30",
        ADMIN_SLUG="ADM0|ADM1|ADM2",
        model="giri|jrc|wri"
    script:
        "./relocation_costs.py"
""" 
Test with
snakemake -c1 data/results/adaptation/costs/countries/RWA/RWA_adaptation-cost_rl_m-jrc_duc11_ADM2.gpkg
"""

rule dry_proofing_costs:
    """
    This rule calculates the sub-national costs of dry-proofing adaptation scenario.
    """
    input:
        admin_areas="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        rp10_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_{model}-flood_RP10.tif",
        rp500_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_{model}-flood_RP500.tif",
        res_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_a.tif",
        res_unit_cost="data/inputs/analysis/countries/{ISO3}/{ISO3}_res_unit_cost.tif"
    output:
        dry_proofing_costs="data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_dp_m-{model}_{ADMIN_SLUG}.gpkg"
    wildcard_constraints:
        ADMIN_SLUG="ADM0|ADM1|ADM2",
        model="giri|jrc|wri"
    script:
        "./dry_proofing_costs.py"
""" 
Test with
snakemake -c1 data/results/adaptation/costs/countries/RWA/RWA_adaptation-cost_dp_m-jrc_ADM2.gpkg
"""

rule capstock_unit_cost:
    """
    This rule calculates the gridded unit cost of the residential capital stock.
    """
    input:
        res_capstock="data/inputs/analysis/countries/{ISO3}/{ISO3}_res_capstock.tif",
        res_area="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_a.tif"
    output:
        res_unit_cost="data/inputs/analysis/countries/{ISO3}/{ISO3}_res_unit_cost.tif"
    shell:
        """
        gdal_calc.py -A {input.res_capstock} -B {input.res_area} \
            --outfile={output.res_unit_cost} \
            --calc="A/B" \
            --NoDataValue=0 \
            --type=Float32
        """
    

configfile: "config/config.yaml"
countries = ['KEN']
ADMINS = ["ADM0", "ADM1", "ADM2"]
RPs = [10, 20, 50, 100]
fp_urban = [21, 22, 23, 30]
DUC_protection = [21, 22, 23, 30]
rl_urban = [11, 12, 13, 21, 22, 23, 30]

rule fp_costs_for_all_countries:
    input:
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban}_{ADM}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, RP=RPs, urban=DUC_protection)

rule pc_costs_for_all_countries:
    input:
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban}_{ADM}.gpkg",
                ISO3=countries, ADM=ADMINS, RP=RPs, urban=fp_urban),
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_rl_m-jrc_duc{urban}_{ADM}.gpkg",
                ISO3=countries, ADM=ADMINS, urban=rl_urban),
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_dp_m-jrc_{ADM}.gpkg",
                ISO3=countries, ADM=ADMINS)