import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()
config = BR.load_organism()

##### Comparison processing #####

def get_comparison_dir_list(condition_list):
    comparison_dir_list = list()
    for condition1 in condition_list:
        if ':' in condition1:
            conditions = condition1.split(":")
            comparison_dir_list.append(conditions[0] + "_vs_" + conditions[1])
        else:
            for condition2 in condition_list[condition_list.index(condition1):]:
                if ':' not in condition2 and condition2 != condition1:
                    comparison_dir_list.append(condition1 + "_vs_" + condition2)
    return comparison_dir_list


## create list of conditions
if config['conditions_to_compare'] == "all":
    condition_list = sorted(sample_tab.condition.unique())
    condition_list_first = [condition for condition in condition_list if
                            not re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list_second = [condition for condition in condition_list if
                             re.search("ctrl|control|wildtype|wt|normal",condition,flags=re.IGNORECASE)]
    condition_list = condition_list_first + condition_list_second
    comparison_dir_list = get_comparison_dir_list(condition_list)
else:
    comparison_dir_list = get_comparison_dir_list(config['conditions_to_compare'].split(","))
    condition_list = set(config['conditions_to_compare'].replace(':',',').split(","))
    if not config['keep_not_compared_samples_for_normalization']:
        sample_tab = sample_tab[sample_tab['condition'].isin(condition_list)]

#
analysis = []
if config["featureCount"]:
    count_over_list = config['count_over'].split(",")
    if ("exon" in count_over_list):
        config["featureCount_exon"] = True
        analysis.append("featureCount_exon")
    if ("gene" in count_over_list):
        config["featureCount_gene"] = True
        analysis.append("featureCount_gene")
    if ("transcript" in count_over_list):
        config["featureCount_transcript"] = True
        analysis.append("featureCount_transcript")
    if ("three_prime_UTR" in count_over_list):
        config["featureCount_3pUTR"] = True
        analysis.append("featureCount_3pUTR")
    if ("five_prime_UTR" in count_over_list):
        config["featureCount_5pUTR"] = True
        analysis.append("featureCount_5pUTR")

if config["HTSeqCount"]:
    count_over_list = config['count_over'].split(",")
    if ("exon" in count_over_list):
        config["HTSeqCount_exon"] = True
        analysis.append("HTSeqCount_exon")
    if ("gene" in count_over_list):
        config["HTSeqCount_gene"] = True
        analysis.append("HTSeqCount_gene")
    if ("transcript" in count_over_list):
        config["HTSeqCount_transcript"] = True
        analysis.append("HTSeqCount_transcript")
    if ("three_prime_UTR" in count_over_list):
        config["HTSeqCount_3pUTR"] = True
        analysis.append("HTSeqCount_3pUTR")
    if ("five_prime_UTR" in count_over_list):
        config["HTSeqCount_5pUTR"] = True
        analysis.append("HTSeqCount_5pUTR")

if config["RSEM"]:
    analysis.append("RSEM")
if config["salmon_align"]:
    analysis.append("salmon_align")
if config["salmon_map"]:
    analysis.append("salmon_map")
if config["kallisto"]:
    analysis.append("kallisto")

biotype_dir_list = config['biotypes'].split(",")

config["analysis_type"] = "|".join(analysis)
config["biotype_list"] = "|".join(biotype_dir_list)
config["comparison"] = "|".join(comparison_dir_list)

if not config["is_paired"]:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"

# check what to report:
if config["qc_samtools"] and config["display_samtools"]:
    config["display_samtools"] = True
else:
    config["display_samtools"] = False

if config["qc_dupradar_RNA"] and config["display_dupradar_RNA"]:
    config["display_dupradar_RNA"] = True
else:
    config["display_dupradar_RNA"] = False

if config["qc_qualimap_RNA"] and config["display_qualimap_RNA"]:
    config["display_qualimap_RNA"] = True
else:
    config["display_qualimap_RNA"] = False

if config["qc_picard_RNA"] and config["display_picard_RNA"]:
    config["display_picard_RNA"] = True
else:
    config["display_picard_RNA"] = False

if config["qc_RSeQC_RNA"] and config["display_RSeQC_RNA"]:
    config["display_RSeQC_RNA"] = True
else:
    config["display_RSeQC_RNA"] = False

if config["qc_biotypes_RNA"] and config["display_biotypes_RNA"]:
    config["display_biotypes_RNA"] = True
else:
    config["display_biotypes_RNA"] = False

if config["qc_fastq_screen_RNA"] and config["display_fastq_screen_RNA"]:
    config["display_fastq_screen_RNA"] = True
else:
    config["display_fastq_screen_RNA"] = False

if config["biobloom"] and config["display_biobloom"]:
    config["display_biobloom"] = True
else:
    config["display_biobloom"] = False


os.makedirs("customer_report",exist_ok=True)

f=open("customer_report/customer_report.json", "w")
json.dump(config, f, indent=4)
f.close()

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name) + "|all_samples",
    read_pair_tag = "(_R.)?",
    lib_name = "[^\.\/]+",
    analysis_type = "|".join(analysis),
    condition_list = "|".join(condition_list),
    biotype = "|".join(biotype_dir_list),
    comparison = "|".join(comparison_dir_list)

##### Target rules #####

rule all:
    input:  report = "customer_report/customer_report.html"

##### Modules #####

include: "rules/customer_report.smk"
