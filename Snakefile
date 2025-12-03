# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified:

# To run this on TTU HPC:
#     spack load trf snakemake
# sample run line:
#     snakemake -np --config sample=Guy11_Final_S4
#         OR
#     snakemake --cores 12 --config sample=Guy11_Final_S4

from pathlib import Path

configfile: "config.yaml"

SAMPLES = config["samples"]
OUTPUT_DIR = config["output_dir"]

rule all:
    input:
        expand(f"{{sample}}_trf.bed", sample=SAMPLES)

rule trf1:
    input:
        "../nanopore/{sample}.fasta"
    output:
        f"{{sample}}.fasta.2.7.7.80.10.50.2000.dat"  ##### FIXME: hardcoded variables
    log:
        "logs/trf1_{sample}.log"
    run:
        # This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
        try:
            proc_output = subprocess.check_output(f"trf {input} 2 7 7 80 10 50 2000 -h", shell=True)
        # an exception is raised by check_output() for non-zero exit codes (usually returned to indicate failure)
        except subprocess.CalledProcessError as exc: 
            if exc.returncode == 7: #### FIXME: this is taken from the above hardcoded numbers
                # this exit code is OK
                pass
            else:
                # for all others, re-raise the exception
                raise

rule trf2:
    input:
        f"{{sample}}.fasta.2.7.7.80.10.50.2000.dat"     ##### FIXME: hardcoded variables
    output:
        f"{{sample}}_trf.bed"
    log:
        "logs/trf2_{sample}.log"
    shell:
        r"""
        python3 trf2bed.py \
            --dat {input} \
            --bed {output} \
            --tool repeatseq
        """
