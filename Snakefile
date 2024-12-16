sizes = [103]

# output dir is specidied in _quarto.yml
rule all:
    input:
        expand("intermediate_data/filtered_small_seurat_{size}_obj.rds", size=sizes),
        expand("output/reports/filtered_small_seurat_{size}_obj.html", size=sizes)
        

rule downsample:
    input:
        "../../raw_data/wendisch_et_al_final/oli_github/Monocytes.Rds"
    output: 
        "intermediate_data/small_seurat_{size}_obj.rds"
    shell:
        "quarto render 1_first.qmd -P downsample_size={wildcards.size}"

rule create_first:
    input: 
        user="user_input.xlsx",
        precomputed="intermediate_data/small_seurat_{size}_obj.rds"
    output:
        rds="intermediate_data/filtered_small_seurat_{size}_obj.rds",
        html_name="output/reports/filtered_small_seurat_{size}_obj.html",
        #temp_dir=temp(directory("foo_{size}"))

    shell: 
        """
        quarto render 2_second.qmd -P downsample_size={wildcards.size} \
        -o filtered_small_seurat_{wildcards.size}_obj.html
        
        """

