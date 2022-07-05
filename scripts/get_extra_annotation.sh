#!/bin/bash
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"


cd $script_dir

# echo "***Fetching input samples for PUM2 from encode (eclip foreground and background samples)..."
# wget https://www.encodeproject.org/files/ENCFF041KJT/@@download/ENCFF041KJT.fastq.gz -P input_files_test -q -nc
# echo "1/4 done"
# wget https://www.encodeproject.org/files/ENCFF462SCV/@@download/ENCFF462SCV.fastq.gz -P input_files_test -q -nc
# echo "2/4 done"
# wget https://www.encodeproject.org/files/ENCFF616FCF/@@download/ENCFF616FCF.fastq.gz -P input_files_test -q -nc
# echo "3/4 done"
# wget https://www.encodeproject.org/files/ENCFF495ZPY/@@download/ENCFF495ZPY.fastq.gz -P input_files_test -q -nc
# echo "4/4 done"



echo "***Fetching genome and transcript annotation from ensembl..."

wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz -P ../input_files  -q -nc
if [[ -e "../input_files/Homo_sapiens.GRCh38.98.chr.gtf" ]]; then
    echo "file already exists, 1/3 chromosome gtf"
else
    gunzip "../input_files/Homo_sapiens.GRCh38.98.chr.gtf.gz"
    echo "1/3 chromosome gtf"
fi

wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P ../input_files  -q -nc
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz -P ../input_files  -q -nc
cat ../input_files/Homo_sapiens.GRCh38.cdna.all.fa.gz ../input_files/Homo_sapiens.GRCh38.ncrna.fa.gz > ../input_files/Homo_sapiens.GRCh38.transcriptome.fa.gz
if [[ -e "../input_files/Homo_sapiens.GRCh38.transcriptome.fa" ]]; then
    echo "file already exists, 2/3 transcript annotation from cdnas and ncrnas"
else
    gunzip "../input_files/Homo_sapiens.GRCh38.transcriptome.fa.gz" 
    echo "2/3 transcript annotation from cdnas and ncrnas"
fi

wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -P ../input_files  -q -nc
if [[ -e "../input_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" ]]; then
    echo "file already exists, 3/3 chromosome gtf"
else
    gunzip "../input_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" 
    echo "3/3 chromosome gtf"
fi


echo "***Fetching attract motif database..."
wget  https://attract.cnic.es/attract/static/ATtRACT.zip -P ../input_files  -q -nc
if [[ -d "../input_files/ATtRACT" ]]; then
    echo "folder already exists, done"
else
    unzip ../input_files/ATtRACT.zip -d ../input_files/ATtRACT 
    echo "done"
fi

echo "***Fetching RNA-central database..."
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz -P ../input_files  -q -nc
if  [[ -e "../input_files/homo_sapiens.GRCh38.gff3" ]] ; then
    echo "file already exists, done"
else
    gunzip "../input_files/homo_sapiens.GRCh38.gff3.gz" 
    echo "done"
fi

