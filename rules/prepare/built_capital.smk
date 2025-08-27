"""
Disaggregate GIRI built capital on GHS volume grids and clip for countries
"""

rule summarise_giri_bem_admin:
    """
    Summarize BEM data at admin level.
    Note: ADM2 fails for SDN - but that's fine; SDN isn't needed.
    """
    input:
        adm2="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM2_fixed.gpkg",
        adm1="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM1_fixed.gpkg",
        adm0="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM0_fixed.gpkg",
        res_raster="data/inputs/giri/bem_5x5_valfis_res.tif",
        nres_raster="data/inputs/giri/bem_5x5_valfis_nres.tif"
    output:
        adm2="data/inputs/giri/bem_5x5_valfis_adm2.csv",
        adm1="data/inputs/giri/bem_5x5_valfis_adm1.csv",
        adm0="data/inputs/giri/bem_5x5_valfis_adm0.csv"
    script:
        "./summarize_bem.py"
"""
Test with
snakemake -c1 data/inputs/giri/bem_5x5_valfis_adm2.csv"
"""

rule clip_GIRI_BEM:
    """
    Clip GIRI BEM layer for a country
    """
    input:
        global_res_bem="data/inputs/giri/bem_5x5_valfis_res.tif",
        global_nres_bem="data/inputs/giri/bem_5x5_valfis_nres.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        trimmed_res_bem="data/inputs/analysis/countries/{ISO3}/{ISO3}_bem_5x5_res.tif",
        trimmed_nres_bem="data/inputs/analysis/countries/{ISO3}/{ISO3}_bem_5x5_nres.tif"
    shell:
        """
        set -ex 

        mkdir --parents $(dirname {output.trimmed_res_bem})

        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -co compress=lzw \
            {input.global_res_bem} \
            {output.trimmed_res_bem}

        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -of GTiff \
            -co BIGTIFF=YES \
            -co compress=lzw \
            {input.global_nres_bem} \
            {output.trimmed_nres_bem}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_bem_5x5_res.tif
"""
        

rule disaggregate_capital:
    """
    Disaggregates GIRI capital stock (at ADM2) across GHS volume grids for a country
    """
    input:
        adm2_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        giri_bem="data/inputs/giri/bem_5x5_valfis_adm2.csv",
        res_volume="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-res_v.tif",
        nres_volume="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-nres_v.tif",
    output:
        nres_capstock="data/inputs/analysis/countries/{ISO3}/{ISO3}_nres_capstock.tif",
        res_capstock="data/inputs/analysis/countries/{ISO3}/{ISO3}_res_capstock.tif",
    script:
        "./disaggregate_adm_bem.py"
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_nres_capstock.tif
"""