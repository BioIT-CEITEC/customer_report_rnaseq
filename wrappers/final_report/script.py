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

command = "sed -i 's/NaN/0/' "+ snakemake.params.mqc_json

f = open(log_filename, 'a+')
f.write("## Need to fix 'NaN' in JSON file \n")
f.write("## COMMAND: "+command+"\n")
shell(command)

command = " cp '"+os.path.abspath(os.path.dirname(__file__))+"/customer_report.Rmd' customer_report/"

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