#!/usr/bin/env python
'''
parse .aln file - output from cap3 program. Output is fasta file and
profile file
'''
import argparse
import re
import zipfile
import tempfile
import textwrap
import os

def parse_args():
    '''Argument parsin'''
    description = """
    parsing cap3 assembly aln output
    """

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-re',
                        '--re_file',
                        default=None,
                        required=True,
                        help="RepeatExlorer archive or directory",
                        type=str,
                        action='store')
    parser.add_argument('-f',
                        '--fasta',
                        default=None,
                        required=True,
                        help="fasta output file name",
                        type=str,
                        action='store')
    parser.add_argument('-m',
                        '--min_coverage',
                        default=5,
                        required=False,
                        help="minimum contig coverage",
                        type=int,
                        action="store")
    parser.add_argument('-L',
                        '--min_contig_length',
                        default=50,
                        required=False,
                        help="minimum contig length",
                        type=int,
                        action="store")
    return parser.parse_args()


def get_header(f):
    aln_header = "    .    :    .    :    .    :    .    :    .    :    .    :"
    contig_lead = "******************"
    aln_start = -1
    while True:
        line = f.readline()
        if not line:
            return None, None
        if line[0:18] == contig_lead:
            line2 = f.readline()
        else:
            continue
        if aln_header in line2:
            aln_start = line2.index(aln_header)
            break
    contig_name = line.split()[1] + line.split()[2]
    return contig_name, aln_start


def segment_start(f):
    pos = f.tell()
    line = f.readline()
    # detect next contig or end of file
    if "********" in line or line == "" or "Number of segment pairs = " in line:
        segment = False
    else:
        segment = True
    f.seek(pos)
    return segment


def get_segment(f, seq_start):
    if not segment_start(f):
        return None, None
    aln = []
    while True:
        line = f.readline()
        if ".    :    .    :" in line:
            continue
        if "__________" in line:
            consensus = f.readline().rstrip('\n')[seq_start:]
            f.readline()  # empty line
            break
        else:
            aln.append(line.rstrip('\n')[seq_start:])
    return aln, consensus


def aln2coverage(aln):
    coverage = [0] * len(aln[0])
    for a in aln:
        for i, c in enumerate(a):
            if c not in " -":
                coverage[i] += 1
    return coverage


def read_contig(f, seq_start):
    contig = ""
    coverage = []
    while True:
        aln, consensus = get_segment(f, seq_start)
        if aln:
            contig += consensus
            coverage += aln2coverage(aln)
        else:
            break
    return contig, coverage


def remove_gaps(consensus, coverage):
    if "-" not in consensus:
        return consensus, coverage
    new_coverage = [
        cov for cons, cov in zip(consensus, coverage) if cons != "-"
    ]
    new_consensus = consensus.replace("-", "")
    return new_consensus, new_coverage


def extract_contigs_from_re_archive(archive, aln_output):
    with zipfile.ZipFile(archive, 'r') as zip_object, open(aln_output,
                                                           'w') as fout:
        flist = zip_object.infolist()
        for fn in flist:
            if re.match('seqclust.+[.]aln$', fn.filename):
                with zip_object.open(fn.filename) as aln:
                    for l in aln:
                        fout.write(l.decode('utf-8'))
    return aln_output

def read_tarean_fasta(fobj):
    ids = []
    s = []
    for i in fobj:
        if isinstance(i, str):
            ii = i
        else:
            ii = i.decode('utf-8')
        if ii[0] == ">":
            ids.append(ii)
            s.append("")
        else:
            s[-1] = s[-1] + ii.strip()
    return ids, s


def extract_contigs_from_re_directory(re_dir, aln_output):
    with open(aln_output, 'w') as fout:
        for subdir, dirs, files in os.walk(re_dir):
            for fn in files:
                fn_full = subdir + os.sep + fn
                if re.match('^.+seqclust.+[.]aln$', fn_full):
                    print(fn_full)
                    with open(fn_full) as aln:
                        for l in aln:
                            fout.write(l)
    return aln_output



def extract_tarean_contigs_from_re_archive(archive):
    with zipfile.ZipFile(archive, 'r') as zip_object:
        flist = zip_object.infolist()
        seqs_all = []
        ids_all = []
        for fn in flist:
            if re.match("seqclust.+dir_CL[0-9]+[/]tarean_contigs.fasta", fn.filename):
                print(fn.filename)
                with zip_object.open(fn.filename) as fobj:
                    ids, seqs = read_tarean_fasta(fobj)
                    # wrap sequences
                    seqs = ["\n".join(textwrap.wrap(s + s, 80)) for s in seqs]
                    seqs_all += seqs
                    ids_all += ids
    return ids_all, seqs_all


def extract_tarean_contigs_from_re_directory(re_dir):
    seqs_all = []
    ids_all = []
    for subdir, dirs, files in os.walk(re_dir):
        for fn in files:
            fn_full = subdir + os.sep + fn
            if re.match("^.+seqclust.+dir_CL[0-9]+[/]tarean_contigs.fasta", fn_full):
                print (fn_full)
                with open(fn_full) as fobj:
                    ids, seqs = read_tarean_fasta(fobj)
                    # wrap sequences
                    seqs = ["\n".join(textwrap.wrap(s + s, 80)) for s in seqs]
                    seqs_all += seqs
                    ids_all += ids
    return ids_all, seqs_all

def filter_contigs(consensus, coverage, min_coverage=5):
    x = "".join([
        s if cov >= min_coverage else " "
        for s, cov in zip(consensus, coverage)
    ]).strip()
    consensus_N = "\n".join(textwrap.wrap(x.replace(" ", "N"),80))
    return consensus_N


def main():
    args = parse_args()
    if os.path.isdir(args.re_file):
        # extract from directory
        ids, seqs = extract_tarean_contigs_from_re_directory(args.re_file)
        aln_file = extract_contigs_from_re_directory(
            args.re_file,
            tempfile.NamedTemporaryFile().name)
    else:
        # extract aln from archive
        ids, seqs = extract_tarean_contigs_from_re_archive(args.re_file)
        aln_file = extract_contigs_from_re_archive(
            args.re_file,
            tempfile.NamedTemporaryFile().name)

    with open(aln_file, 'r') as f1, open(args.fasta, 'w') as ffasta:
        while True:
            contig_name, seq_start = get_header(f1)
            if contig_name:
                consensus, coverage = remove_gaps(*read_contig(f1, seq_start))
                clean_consensus = filter_contigs(consensus, coverage,
                                                 args.min_coverage)
                if len(clean_consensus) >= args.min_contig_length:
                    ffasta.write(">{}\n".format(contig_name))
                    ffasta.write("{}\n".format(clean_consensus))
            else:
                break

        # write tarean sequences:
        for i, s in zip(ids, seqs):
            ffasta.write(i)
            ffasta.write(s + "\n")



if __name__ == "__main__":

    main()
