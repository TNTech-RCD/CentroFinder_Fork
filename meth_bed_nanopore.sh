#!/bin/bash
#SBATCH --job-name=extract_methylation
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --account=phd-ssalimi42
#SBATCH --cpus-per-task=28
#SBATCH --mem=124G

# =========================================================
# MODULES / TOOLS
# =========================================================
spack load samtools/f6zgopc minimap2/m7v3j5b


# CONFIGURATION
# =========================================================
REFERENCE="/home/tntech.edu/ssalimi42/work/Thesis_Finall/sterategy_2/centromere/Guy11/bed_meth/Guy11_Final_S4.fasta"
READ_DIR="/home/tntech.edu/ssalimi42/work/Thesis_Finall/sterategy_2/centromere/Guy11/bed_meth"
WORKDIR="/home/tntech.edu/ssalimi42/work/Thesis_Finall/sterategy_2/centromere/Guy11/bed_meth/bed"
mkdir -p "$WORKDIR"
cd $WORKDIR
THREADS=${SLURM_CPUS_PER_TASK:-28}

INPUT_FASTQ="$1"



# =========================================================
# STEP 1. ALIGN READS (if no BAM yet)
# =========================================================

minimap2 -t "$THREADS" -ax map-ont -Y "$REFERENCE" "$READ_DIR/$INPUT_FASTQ.fastq" > "${INPUT_FASTQ}.sam"
samtools sort -@ "$THREADS" -o "${INPUT_FASTQ}_sorted.bam" "${INPUT_FASTQ}.sam"
samtools index "${INPUT_FASTQ}_sorted.bam"


# =========================================================
# STEP 2. BAM to BED
# =========================================================
echo "Extracting per-site methylation calls with modbam2bed..."
/home/tntech.edu/ssalimi42/work/Thesis_Finall/modbam2bed/modbam2bed "$REFERENCE" "${INPUT_FASTQ}_sorted.bam" > ${INPUT_FASTQ}.bed


