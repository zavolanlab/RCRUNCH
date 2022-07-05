# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import gzip
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter
from gzip import open as gzopen
import shutil

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def main():
    """ Main function """

    __doc__ = "Convert ENCODE eCLIP data into format appropriate for the UMI-tools."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--infile",
        help="input_fastq_file",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outfile",
        help="output_fastq_file",
        required=True,
        metavar="FILE")
    
    parser.add_argument(
        "--infile2",
        help="input_fastq_file",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--outfile2",
        help="output_fastq_file",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--read_format",
        help="if there are umis to be removed",
        required=False,
        default='other',
        metavar="str")

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    infiles = [options.infile]
    if options.infile2:
        infiles.append(options.infile2)

    outfiles = [options.outfile]
    if options.outfile2:
        outfiles.append(options.outfile2)

    for inp, outp in zip(infiles, outfiles):
        if options.read_format == 'encode':
            output_file = gzip.open(outp, 'wt')
            for record in SeqIO.parse(gzopen(inp, "rt"), format="fastq"):
                    seqid = record.id
                    barcode = seqid.split(':')[0]
                    record.id = f'{seqid}_{barcode}'
                    record.id = record.id.replace((barcode + ':'), '')
                    record.description = record.description.split(' ')[-1]
                    output_file.write(record.format("fastq"))
            output_file.close()
        else:
            shutil.copy(inp, outp)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
