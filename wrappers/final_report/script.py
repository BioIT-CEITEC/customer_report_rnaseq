#############################################################
# wrapper for rule: final_report
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: final_report \n##\n")
f.close()

command = " cp '"+os.path.abspath(os.path.dirname(__file__))+"/customer_report.Rmd' customer_report/"

shell(command)

command = """ for i in results/DE_*/*_vs_*/*/*ggplot2.pdf results/DE_*/*_vs_*/*/*ggpubr.pdf results/DE_*/*_vs_*/*/heatmap*.pdf ; do convert -units PixelsPerInch -density 200 $i $i.png ; done"""

f = open(log_filename, 'a+')
f.write("## Converting pdf image in png.\n")
f.write("## COMMAND: "+command+"\n")
shell(command)

command = """ Rscript -e "rmarkdown::render('customer_report/customer_report.Rmd', params=list(config = '""" + snakemake.params.config + """', set_author =  '""" + snakemake.params.author + """', set_email = '""" + snakemake.params.email + """'))" """ +\
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " ls " + snakemake.output.html

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " rm customer_report/customer_report.Rmd"

shell(command)