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
        """
        output_dir=$(dirname {output})
        mkdir -p $output_dir

        wget -nc "{params.url}" -O "$output_dir/{params.top_zip}"

        tmp_nested="$output_dir/tmp_2019.zip"
        unzip -p "$output_dir/{params.top_zip}" "{params.nested_zip}" \
            > "$tmp_nested"

        unzip -p "$tmp_nested" "{params.nested_member}" \
            > "{output}"

        rm "$tmp_nested"
        """