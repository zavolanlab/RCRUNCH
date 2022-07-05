from Bio import SeqIO
import re
import sys
import os


def main():
    """ Main function """
    __doc__ = "Create test subset of genome and annotation files"
    chromosome = '20'
    outpath = os.path.dirname(os.path.realpath(__file__))
    os.chdir(outpath)
    transcripts = set()
    for inp, outp in zip(["input_files_test/Homo_sapiens.GRCh38.98.chr.gtf"],
                        [f"input_files_test/Homo_sapiens.GRCh38.98.chr{chromosome}.gtf"]):
        with open(inp, 'r') as myfile:
            with open(outp, 'w') as outfile:
                for line in myfile.readlines():
                    if line.startswith('#'):
                        outfile.write(line)
                    elif line.startswith(f'{chromosome}\t'):
                        outfile.write(line)
                        if 'transcript_id' in line:
                            r1 = re.findall(r"transcript_id \"(\w+)\"", line)
                            transcripts.add(r1[0])
                    else:
                        continue
    for inp, outp in zip(["input_files_test/Homo_sapiens.GRCh38.transcriptome.fa"],
                        [f"input_files_test/Homo_sapiens.GRCh38.transcriptome.chr{chromosome}.fa"]):
        output_file = open(outp, 'w')
        for record in SeqIO.parse(inp, format="fasta"):
            if record.id.split('.')[0] in transcripts:
                output_file.write(record.format("fasta"))
        output_file.close()
    os.remove("input_files_test/Homo_sapiens.GRCh38.transcriptome.fa")
    os.remove("input_files_test/Homo_sapiens.GRCh38.98.chr.gtf")
    sys.stdout.write('done!\n')


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interruptions
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

