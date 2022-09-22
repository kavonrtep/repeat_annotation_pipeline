#!/usr/bin/env python
"""
This script is used to run dustmasker on a fasta file. It will create
a bed file from default dustmasker output.
"""
import argparse
import subprocess
import tempfile

# parse arguments from command line, and pass it to duskmasker

parser = argparse.ArgumentParser(
    description="""This script is used to run dustmasker on a fasta file. It will create
    a bed file from default dustmasker output.""",
    formatter_class=argparse.RawTextHelpFormatter, )
parser.add_argument(
    '-f', '--fasta', default=None, required=True, help="fasta file", type=str,
    action='store'
    )
parser.add_argument(
    '-o', '--output', default=None, required=True, help="output file name", type=str,
    action='store'
    )
parser.add_argument(
    '-w', '--window', default=60, required=False, help="dustmasker window size", type=int,
    action='store'
    )
parser.add_argument(
    '-l', '--level', default=20, required=False, help="dustmasker level", action='store'
    )
args = parser.parse_args()


def main(args):
    """
    run dustmasker and convert output do bed file
    """
    # temp file for dust maske output
    tmpfile = tempfile.NamedTemporaryFile().name
    # run dustmasker
    subprocess.call(
        ["dustmasker", "-in", args.fasta, "-out", tmpfile, "-window", str(args.window),
         "-level", str(args.level)]
        )
    # create bed file from dustmasker output
    # syntax of dustmasker output is:
    # >contig_name
    # start_position - end_position
    # start_position - end_position
    # ...
    # if not masked, only contig name is printed
    with open(args.output, "w") as f:
        with open(tmpfile, "r") as f2:
            for line in f2:
                if line[0] == ">":
                    contig_name = line.strip()[1:]
                    continue
                else:
                    line = line.strip()
                    line = line.split()
                    f.write(contig_name + "\t" + line[0] + "\t" + line[2] + "\n")


if __name__ == '__main__':
    main(args)
