"""
Rule book for figure plotting.
"""

rule plot_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region, given a model and vulnerability curve
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        risk_file="data/results/flood_risk/{ISO3}/{ISO3}_{MODEL}-flood-risk_AAR_V-{VULN_CURVE}.tif",
    output:
        figure_directory = directory("figures/concentration_curves/{ISO3}/{ADMIN_SLUG}/{MODEL}/{VULN_CURVE}/"),
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./concentration_curves.py"
"""
Test with
snakemake -c1 figures/concentration_curves/RWA/ADM-0/jrc/JRC
"""