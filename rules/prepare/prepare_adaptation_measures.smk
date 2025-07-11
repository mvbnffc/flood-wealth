"""
Prepare adaptation scenarios
"""

rule prepare_river_flood_protection:
    """
    Load national river network, national urbanization data, flood protection (FLOPROS), and specify level of protection and the degree of urbanization
    to protect. The script will then create a new flopros layer (with adaptation incorporated), it will also print a text file with the cost of 
    the adaptation measures calculated based on river length and protection delta (formula taken from Boulange et al (2023))

    Reference: Boulange et al (2023) https://link.springer.com/article/10.1007/s11069-023-06017-7
    """
    input:
        flopros="data/inputs/analysis/countries/{ISO3}/{ISO3}_flopros.tif",
        river_network="data/inputs/analysis/countries/{ISO3}/{ISO3}_river_network.gpkg",
        urbanization="data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.gpkg"
    output:
        flood_protection="data/inputs/analysis/countries/{ISO3}/{ISO3}_adaptation_fp_rp{RP}_duc{urban_class}.tif",
        cost_of_protection="data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban_class}.txt"
    wildcard_constraints:
        urban_class="11|12|13|21|22|23|30"
    script:
        "./prepare_river_flood_protection.py"