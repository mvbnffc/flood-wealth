"""
Rule book for figure plotting.
"""

rule plot_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region, given a model and vulnerability curve
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_AAR_V-{VULN_CURVE}.tif",
    output:
        figure_directory = directory("figures/concentration_curves/countries/{ISO3}/{ADMIN_SLUG}/{MODEL}/{VULN_CURVE}/"),
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./concentration_curves.py"
"""
Test with
snakemake -c1 figures/concentration_curves/countries/RWA/ADM-0/jrc/JRC
"""

rule plot_protected_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region, given a model and vulnerability curve
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_protected_AAR_V-{VULN_CURVE}.tif",
    output:
        figure_directory = directory("figures/concentration_curves/countries/{ISO3}/{ADMIN_SLUG}/{MODEL}/{VULN_CURVE}/protected/"),
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./concentration_curves.py"
"""
Test with
snakemake -c1 figures/concentration_curves/countries/RWA/ADM-0/jrc/JRC
"""

rule plot_observed_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region, given a historical observed flood dataset
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/analysis/countries/{ISO3}/{ISO3}_gfd-flood.tif",
    output:
        figure_directory = directory("figures/concentration_curves/countries/{ISO3}/{ADMIN_SLUG}/gfd/"),
    wildcard_constraints:
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./concentration_curves.py"
"""
Test with
snakemake -c1 figures/concentration_curves/countries/RWA/ADM-0/gfd
"""

rule plot_event_concentration_curve:
    """
    This rule plots the concentration curve for a specific flood event, and country
    """
    input:
        rwi_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        flood_file="data/inputs/analysis/events/DFO_{event_id}/{ISO3}_{event_id}.tif",
    output:
        figure_file = "figures/concentration_curves/events/DFO_{event_id}/{ISO3}_DFO_{event_id}_concentration_curve.png"
    script:
        "./event_concentration_curve.py"
"""
Test with
snakemake -c1 figures/concentration_curves/events/DFO_3136/BGD_DFO_3136_concentration_curve.png
"""