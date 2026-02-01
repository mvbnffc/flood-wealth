"""
Download gridded GDP data

References
---------
GDP: https://www.nature.com/articles/s41597-022-01322-5
"""

rule download_gdp:
    output:
        "data/inputs/gdp/2019gdp.tif"
    params:
        top_zip="Real GDP.zip",
        nested_zip="updated real GDP/2019.zip",
        nested_member="2019GDP.tif",
        url="https://figshare.com/ndownloader/files/31456837"
    shell:
        r"""
        set -euo pipefail

        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        # Download robustly (Figshare may return 202/HTML initially)
        rm -f "$output_dir/{params.top_zip}"
        for i in $(seq 1 20); do
          curl -L -o "$output_dir/{params.top_zip}" "{params.url}" || true
          if head -c 2 "$output_dir/{params.top_zip}" | grep -q "PK"; then
            break
          fi
          echo "Download not ready (attempt $i). Retrying..." >&2
          sleep 2
        done
        head -c 2 "$output_dir/{params.top_zip}" | grep -q "PK"

        tmp_nested="$output_dir/tmp_2019.zip"
        unzip -p "$output_dir/{params.top_zip}" "{params.nested_zip}" > "$tmp_nested"
        unzip -p "$tmp_nested" "{params.nested_member}" > "{output}"
        rm -f "$tmp_nested"
        """
