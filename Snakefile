sizes = [101, 102, 103]

# output dir is specidied in _quarto.yml
rule all:
    input:
        expand("intermediate_data/filtered_small_seurat_{size}_obj.rds", size=sizes),
        expand("output/reports/filtered_small_seurat_{size}_obj.html", size=sizes)
        

rule downsample:
    input:
        "../../raw_data/mock_data/Monocytes_{size}.Rds"
    output: 
        "intermediate_data/small_seurat_{size}_obj.rds"
    shell:
      """
      cp 1_first.qmd {wildcards.size}_first_temp.qmd && \
      quarto render {wildcards.size}_first_temp.qmd -P downsample_size={wildcards.size} \
      -o small_seurat_{wildcards.size}_obj.html && \
      rm {wildcards.size}_first_temp.qmd 

      """

rule create_first:
    input: 
        user="user_input.xlsx",
        precomputed="intermediate_data/small_seurat_{size}_obj.rds"
    output:
        rds="intermediate_data/filtered_small_seurat_{size}_obj.rds",
        html_name="output/reports/filtered_small_seurat_{size}_obj.html"
        #temp_dir=temp(directory("foo_{size}"))

    shell: 
        """
        cp 2_second.qmd {wildcards.size}_second_temp.qmd && \
        quarto render {wildcards.size}_second_temp.qmd  -P downsample_size={wildcards.size} \
        -o filtered_small_seurat_{wildcards.size}_obj.html && \
        rm {wildcards.size}_second_temp.qmd 
        """


