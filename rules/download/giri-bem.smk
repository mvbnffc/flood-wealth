"""
Rule taken and adapted from open-gira...

Download GIRI Building Exposure Model

Total Building Stock US$ (Resolution: 5km x 5km)

The Building Exposure Model provides information on the building type and the
economic value of the built environment for the non-residential (employment,
health and education) sectors, for each country and territory of the world.

This dataset was produced by UNEP/GRID-Geneva in May 2023.

Source
------
https://giri.unepgrid.ch

Reference
---------
Thomas Piller, Antonio Benvenuti & Andrea De Bono (2023) The GIRI global
building exposure model (BEM)
https://giri.unepgrid.ch/sites/default/files/2023-09/GIRI_BEM_report_UNIGE.pdf
"""

rule download_giri_bem:
    output:
        res="data/inputs/giri/bem_5x5_valfis_res.tif",
        nres="data/inputs/giri/bem_5x5_valfis_nres.tif",
    shell:
        """
        output_dir=$(dirname {output.res})

        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_res.tif \
            --directory-prefix=$output_dir
        wget -nc https://hazards-data.unepgrid.ch/bem_5x5_valfis_nres.tif \
            --directory-prefix=$output_dir
        """

rule summarise_giri_bem_admin_chunked:
    """
    Summarize BEM data at admin level.
    Note: ADM2 fails for SDN - but its OK as that's not a country we have RWI for
    """ 
    input:
        adm2="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM2.gpkg",
        adm1="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM1.gpkg",
        adm0="data/inputs/boundaries/global/geoBoundariesCGAZ_ADM0.gpkg",
        res_raster="data/inputs/giri/bem_5x5_valfis_res.tif",
        nres_raster="data/inputs/giri/bem_5x5_valfis_nres.tif"
    output:
        adm2="data/inputs/giri/bem_5x5_valfis_adm2.csv",
        adm1="data/inputs/giri/bem_5x5_valfis_adm1.csv",
        adm0="data/inputs/giri/bem_5x5_valfis_adm0.csv"
    shell:
        """
        echo "=== CHUNKED PROCESSING FOR ADM2 (49k features) ==="
        
        # Step 1: Get list of unique countries
        echo "Getting list of countries..."
        ogrinfo -sql "SELECT DISTINCT shapeGroup FROM globalADM2 WHERE shapeGroup IS NOT NULL ORDER BY shapeGroup" \
            {input.adm2} | \
            grep "shapeGroup (String)" | \
            sed 's/.*= //' | \
            grep -v "^$" | \
            sort -u > temp_countries.txt
        
        total_countries=$(wc -l < temp_countries.txt)
        echo "Found $total_countries countries to process"
        
        # Step 2: Initialize output file with header
        echo "fid,shapeName,shapeGroup,sum(res),sum(nres)" > {output.adm2}
        
        # Step 3: Process each country separately
        processed_countries=0
        failed_countries=0
        
        while IFS= read -r country; do
            if [[ ! -z "$country" && "$country" != "NULL" ]]; then
                processed_countries=$((processed_countries + 1))
                echo "[$processed_countries/$total_countries] Processing: $country"
                
                # Create temporary country file
                country_file="temp_${{country//[^a-zA-Z0-9]/_}}_adm2.gpkg"
                country_csv="temp_${{country//[^a-zA-Z0-9]/_}}_adm2.csv"
                
                # Extract country data
                if ogr2ogr -sql "SELECT * FROM globalADM2 WHERE shapeGroup = '$country'" \
                    -nln country_adm2 \
                    "$country_file" \
                    {input.adm2} 2>/dev/null; then
                    
                    # Count features in this country
                    feature_count=$(ogrinfo -sql "SELECT COUNT(*) FROM country_adm2" "$country_file" 2>/dev/null | grep -oP 'Integer.*= \K\d+' || echo "unknown")
                    echo "  Features: $feature_count"
                    
                    # Run exactextract on this country chunk
                    if exactextract \
                        --strategy=feature-sequential \
                        --max-cells 2 \
                        --threads 1 \
                        -p "$country_file" \
                        -r "res:{input.res_raster}" \
                        -r "nres:{input.nres_raster}" \
                        -f fid \
                        --include-col shapeName shapeGroup \
                        -o "$country_csv" \
                        -s "sum(res)" \
                        -s "sum(nres)" 2>/dev/null; then
                        
                        # Append results (skip header)
                        if [[ -f "$country_csv" && -s "$country_csv" ]]; then
                            tail -n +2 "$country_csv" >> {output.adm2}
                            echo "  ✓ Success - $(( $(wc -l < "$country_csv") - 1 )) features processed"
                        fi
                    else
                        echo "  ✗ exactextract failed for $country"
                        failed_countries=$((failed_countries + 1))
                    fi
                else
                    echo "  ✗ Could not extract data for $country"
                    failed_countries=$((failed_countries + 1))
                fi
                
                # Cleanup temporary files
                rm -f "$country_file" "$country_csv"
                
                # Progress update every 10 countries
                if (( processed_countries % 10 == 0 )); then
                    echo "  Progress: $processed_countries/$total_countries countries processed"
                fi
            fi
        done < temp_countries.txt
        
        # Step 4: Summary
        total_features=$(( $(wc -l < {output.adm2}) - 1 ))
        echo "=== ADM2 CHUNKED PROCESSING COMPLETE ==="
        echo "Countries processed: $processed_countries"
        echo "Countries failed: $failed_countries"
        echo "Total features in output: $total_features"
        
        # Cleanup
        rm -f temp_countries.txt
        
        echo "=== PROCESSING ADM1 (normal - fewer features) ==="
        exactextract \
            --strategy=feature-sequential \
            --max-cells 10 \
            --threads 2 \
            -p {input.adm1} \
            -r "res:{input.res_raster}" \
            -r "nres:{input.nres_raster}" \
            -f fid \
            --include-col shapeName shapeGroup \
            -o {output.adm1} \
            -s "sum(res)" \
            -s "sum(nres)"
        echo "ADM1 completed: $(( $(wc -l < {output.adm1}) - 1 )) features"

        echo "=== PROCESSING ADM0 (normal - fewer features) ==="
        exactextract \
            --strategy=feature-sequential \
            --max-cells 10 \
            --threads 2 \
            -p {input.adm0} \
            -r "res:{input.res_raster}" \
            -r "nres:{input.nres_raster}" \
            -f fid \
            --include-col shapeName shapeGroup \
            -o {output.adm0} \
            -s "sum(res)" \
            -s "sum(nres)"
        echo "ADM0 completed: $(( $(wc -l < {output.adm0}) - 1 )) features"
        
        echo "=== ALL PROCESSING COMPLETE ==="
        """