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

rule prepare_relocation:
    """
    Relocate all those people living in the 10-year flood plain with a flood depth > 1 m.
    User needs to specify the degree of urbanization classification for people that will be relocated.
    For these areas - flood protection layer will be set to 1000 yr RP (assuming that individual's position on wealth distribution remains the same - but now no risk)
    Rule also outputs the total # of people relocated in the country (can then be used to calculate costs)
    """
    input:
        flopros_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_flopros.tif",
        flood_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_{model}-flood_RP10.tif",
        pop_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        urbanization_path="data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.tif"
    output:
        flood_protection="data/inputs/analysis/countries/{ISO3}/{ISO3}_adaptation_rl_m-{model}_duc{urban_class}.tif",
        people_relocated="data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_rl_m-{model}_duc{urban_class}.txt"
    wildcard_constraints:
        urban_class="11|12|13|21|22|23|30",
        model="giri|jrc|wri"
    script:
        "./prepare_relocation.py"

