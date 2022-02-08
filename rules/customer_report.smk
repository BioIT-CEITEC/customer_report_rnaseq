rule final_report:
    input:  tsv = expand("results/DE_{analysis_type}/{comparison}/{biotype}/DESeq2.tsv", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list),
    output: html = "customer_report/customer_report.html"
    params: config = "customer_report.json",
            author = config["author"],
            email = config["email"]
    conda:  "../wrappers/final_report/env.yaml"
    log:    "logs/customer_report.log"
    script: "../wrappers/final_report/script.py"

