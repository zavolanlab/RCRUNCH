# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import check_call


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Obtain sliding windows of read coverage from bam file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)
    # ----------------------------------------------------------
    parser.add_argument(
        "--fastq_in",
        dest="fastq_in",
        help="Fastq file which contains reads",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--fastq_out",
        dest="fastq_out",
        help="Fastq file which contains reads",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--fastq_in2",
        dest="fastq_in2",
        help="Fastq file which contains reads",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--fastq_out2",
        dest="fastq_out2",
        help="Fastq file which contains reads",
        required=False,
        metavar="FILE")

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

    # -------------------------------------------------------------------------
    #           Cast options
    # -------------------------------------------------------------------------
    if '.gz' not in options.fastq_out:
        sys.stderr.write("Output file should have a gz extension!\n")
        sys.exit(0)
    if options.fastq_out2:
        if '.gz' not in options.fastq_out:
            sys.stderr.write("Output file should have a gz extension!\n")
            sys.exit(0)

    tmp_fastq_out = options.fastq_out.replace('.gz', '')
    get_unique_reads(options.fastq_in, tmp_fastq_out)
    check_call(['gzip', tmp_fastq_out])

    if (options.fastq_in2 and options.fastq_out2):
        tmp_fastq_out2 = options.fastq_out2.replace('.gz', '')
        get_unique_reads(options.fastq_in2, tmp_fastq_out2)
        check_call(['gzip', tmp_fastq_out2])


def get_unique_reads(infile, outfile):
    # fastq = open(infile, "r")
    out_file = open(outfile, "w")

    seqid = ''
    seqid1 = 'no'
    with open(infile) as f:
        for line in f:
            line = line.replace('\n', '')
            if line.startswith('@'):
                seqid1 = seqid
                seqid = line
                if seqid == seqid1:
                    continue
                out_file.write(seqid + '\n')
            else:
                if seqid == seqid1:
                    continue
                out_file.write(line + '\n')
    out_file.close()
    return


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
