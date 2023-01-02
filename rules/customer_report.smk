rule final_report:
    input:  tsv = expand("DE_{analysis_type}/{comparison}/DESeq2.tsv", analysis_type=analysis, comparison=comparison_dir_list, biotype=biotype_dir_list),
    output: html = "customer_report/customer_report.html"
    params: config = "customer_report.json",
            author = config["author"],
            email = config["email"],
            mqc_json = "qc_reports/all_samples/multiqc_data/multiqc_data.json"
    conda:  "../wrappers/final_report/env.yaml"
    log:    "logs/customer_report.log"
    script: "../wrappers/final_report/script.py"

