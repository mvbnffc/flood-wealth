"""
Download Admin Boundary Data

References
---------
GADM: https://gadm.org/data.html
geoBoundaries: https://github.com/wmgeolab/geoBoundaries
"""

configfile: "config/config.yaml"


# Run for all ISO3 codes
rule all_boundaries:
    input:
        expand("data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg", ISO3=config['iso_codes']),
        expand("data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg", ISO3=config['iso_codes'])


rule gadm:
    output:
        gpkg="data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
    shell:
        """
        cd data/inputs/boundaries/{wildcards.ISO3}
        wget https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_{wildcards.ISO3}.gpkg --output-document=gadm_{wildcards.ISO3}.gpkg
        """

rule geobounds_all:
    output:
        adm0="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM0.gpkg",
        adm1="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM1.gpkg",
        adm2="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM2.gpkg",
    shell:
        """
        mkdir -p data/inputs/boundaries/global
        cd data/inputs/boundaries/global

        wget --no-check-certificate "https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/CGAZ/geoBoundariesCGAZ_ADM0.gpkg"
        wget --no-check-certificate "https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/CGAZ/geoBoundariesCGAZ_ADM1.gpkg"
        wget --no-check-certificate "https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/CGAZ/geoBoundariesCGAZ_ADM2.gpkg"
        """

rule geobounds_iso:
    input:
        adm0="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM0.gpkg",
        adm1="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM1.gpkg",
        adm2="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM2.gpkg",
    output:
        gpkg="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.gpkg",
        json="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    run:
        print("Reading global files...")
        adm0=gpd.read_file(input.adm0).to_crs(4326)
        adm1=gpd.read_file(input.adm1).to_crs(4326)
        adm2=gpd.read_file(input.adm2).to_crs(4326)

        try:
            print("Filtering by ISO code...")
            adm0_iso=adm0[adm0["shapeGroup"]==wildcards.ISO3]
            adm1_iso=adm1[adm1["shapeGroup"]==wildcards.ISO3]
            adm2_iso=adm2[adm2["shapeGroup"]==wildcards.ISO3]

        except:
            print("ISO code does not exist in geoBoundaries dataset!")

        def check_and_fix(gdf, layer_name):
            if not gdf.is_valid.all():
                print(f"Warning: Some geometries in {layer_name} are invalid. Fixing them...")
                gdf["geometry"] = gdf["geometry"].buffer(0) # 0 buffer should fix issue
            else:
                print(f"All geometries in {layer_name} are valid.")
            return gdf
        
        # Check and fix geometries for each administrative level
        adm0_iso = check_and_fix(adm0_iso, "ADM0")
        adm1_iso = check_and_fix(adm1_iso, "ADM1")
        adm2_iso = check_and_fix(adm2_iso, "ADM2")

        print("Exporting files...")
        adm0_iso.to_file(output.gpkg, layer="ADM0", driver="GPKG")
        adm1_iso.to_file(output.gpkg, layer="ADM1", driver="GPKG")
        adm2_iso.to_file(output.gpkg, layer="ADM2", driver="GPKG")

        adm0_iso.to_file(output.json, driver="GeoJSON")