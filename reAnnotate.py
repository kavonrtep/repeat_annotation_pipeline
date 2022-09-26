#!/usr/bin/env python
"""
parse blast output table to gff file
"""
import argparse
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

#  check version of python, must be at least 3.7
if sys.version_info < (3, 10):
    sys.exit("Python 3.7 or a more recent version is required.")


def overlap(a, b):
    """
    check if two intervals overlap
    """
    return max(a[0], b[0]) <= min(a[1], b[1])


def blast2disjoint(
        blastfile, seqid_counts=None, start_column=6, end_column=7, class_column=1,
        bitscore_column=11, pident_column=2, canonical_classification=True
        ):
    """
    find all interval beginning and ends in blast file and create bed file
    input blastfile is tab separated file with columns:
    'qaccver saccver pident length mismatch gapopen qstart qend sstart send
   evalue bitscore'  (default outfmt 6
    blast must be sorted on qseqid and qstart
    """
    # assume all in one chromosome!
    starts_ends = {}
    intervals = {}
    if canonical_classification:
        # make regular expression for canonical classification
        # to match: Name#classification
        # e.g. "Name_of_sequence#LTR/Ty1_copia/Angela"
        regex = re.compile(r"(.*)[#](.*)")
        group = 2
    else:
        # make regular expression for non-canonical classification
        # to match: Classification__Name
        # e.g. "LTR/Ty1_copia/Angela__Name_of_sequence"
        regex = re.compile(r"(.*)__(.*)")
        group = 1

    # identify continuous intervals
    with open(blastfile, "r") as f:
        for seqid in sorted(seqid_counts.keys()):
            n_lines = seqid_counts[seqid]
            starts_ends[seqid] = set()
            for i in range(n_lines):
                items = f.readline().strip().split()
                # note 1s and 2s labels are used to distinguish between start and end and
                # guarantee that with same coordinated start will be before end when
                # sorting (1s < 2e)
                starts_ends[seqid].add((int(items[start_column]), '1s'))
                starts_ends[seqid].add((int(items[end_column]), '2e'))
            intervals[seqid] = []
            for p1, p2 in itertools.pairwise(sorted(starts_ends[seqid])):
                if p1[1] == '1s':
                    sp = 0
                else:
                    sp = 1
                if p2[1] == '2e':
                    ep = 0
                else:
                    ep = 1
                intervals[seqid].append((p1[0] + sp, p2[0] - ep))
    # scan each blast hit against continuous region and record hit with best score
    with open(blastfile, "r") as f:
        disjoint_regions = []
        for seqid in sorted(seqid_counts.keys()):
            n_lines = seqid_counts[seqid]
            idx_of_overlaps = {}
            best_pident = defaultdict(lambda: 0.0)
            best_bitscore = defaultdict(lambda: 0.0)
            best_hit_name = defaultdict(lambda: "")
            i1 = 0
            for i in range(n_lines):
                items = f.readline().strip().split()
                start = int(items[start_column])
                end = int(items[end_column])
                pident = float(items[pident_column])
                bitscore = float(items[bitscore_column])
                classification = items[class_column]
                j = 0
                done = False
                while True:
                    # beginning of searched region - does it overlap?
                    c_ovl = overlap(intervals[seqid][i1], (start, end))
                    if c_ovl:
                        # if overlap is detected, add to dictionary
                        idx_of_overlaps[i] = [i1]
                        if best_bitscore[i1] < bitscore:
                            best_pident[i1] = pident
                            best_bitscore[i1] = bitscore
                            best_hit_name[i1] = classification
                        # add search also downstream
                        while True:
                            j += 1
                            if j + i1 >= len(intervals[seqid]):
                                done = True
                                break
                            c_ovl = overlap(intervals[seqid][i1 + j], (start, end))
                            if c_ovl:
                                idx_of_overlaps[i].append(i1 + j)
                                if best_bitscore[i1 + j] < bitscore:
                                    best_pident[i1 + j] = pident
                                    best_bitscore[i1 + j] = bitscore
                                    best_hit_name[i1 + j] = classification
                            else:
                                done = True
                                break

                    else:
                        # does no overlap - search next interval
                        i1 += 1
                    if done or i1 >= (len(intervals[seqid]) - 1):
                        break

            for i in sorted(best_pident.keys()):
                try:
                    classification = re.match(regex, best_hit_name[i]).group(group)
                except AttributeError:
                    classification = best_hit_name[i]
                record = (
                    seqid, intervals[seqid][i][0], intervals[seqid][i][1], best_pident[i],
                    classification)
                disjoint_regions.append(record)
    return disjoint_regions


def remove_short_interrupting_regions(regions, min_len=10, max_gap=2):
    """
    remove intervals shorter than min_len which are directly adjacent to other
    regions on both sides which are longer than min_len and has same classification
    """
    regions_to_remove = []
    for i in range(1, len(regions) - 1):
        if regions[i][2] - regions[i][1] < min_len:
            c1 = regions[i - 1][2] - regions[i - 1][1] > min_len
            c2 = regions[i + 1][2] - regions[i + 1][1] > min_len
            c3 = regions[i - 1][4] == regions[i + 1][4]  # same classification
            c4 = regions[i + 1][4] != regions[i][4]  # different classification
            c5 = regions[i][1] - regions[i - 1][2] < max_gap  # max gap between regions
            c6 = regions[i + 1][1] - regions[i][2] < max_gap  # max gap between regions
            if c1 and c2 and c3 & c4 and c5 and c6:
                regions_to_remove.append(i)
    for i in sorted(regions_to_remove, reverse=True):
        del regions[i]
    return regions


def remove_short_regions(regions, min_l_score=600):
    """
    remove intervals shorter than min_len
    min_l_score is the minimum score for a region to be considered
    l_score = length * PID
    """
    regions_to_remove = []
    for i in range(len(regions)):
        l_score = (regions[i][3] - 50) * (regions[i][2] - regions[i][1])
        if l_score < min_l_score:
            regions_to_remove.append(i)
    for i in sorted(regions_to_remove, reverse=True):
        del regions[i]
    return regions


def join_disjoint_regions_by_classification(disjoint_regions, max_gap=0):
    """
    merge neighboring intervals with same classification and calculate mean weighted score
    weight correspond to length of the interval
    """
    merged_regions = []
    for seqid, start, end, score, classification in disjoint_regions:
        score_length = (end - start + 1) * score
        if len(merged_regions) == 0:
            merged_regions.append([seqid, start, end, score_length, classification])
        else:
            cond_same_class = merged_regions[-1][4] == classification
            cond_same_seqid = merged_regions[-1][0] == seqid
            cond_neighboring = start - merged_regions[-1][2] + 1 <= max_gap
            if cond_same_class and cond_same_seqid and cond_neighboring:
                # extend region
                merged_regions[-1] = [merged_regions[-1][0], merged_regions[-1][1], end,
                                      merged_regions[-1][3] + score_length,
                                      merged_regions[-1][4]]
            else:
                merged_regions.append([seqid, start, end, score_length, classification])
    # recalculate length weighted score
    for record in merged_regions:
        record[3] = record[3] / (record[2] - record[1] + 1)
    return merged_regions


def write_merged_regions_to_gff3(merged_regions, outfile):
    """
    write merged regions to gff3 file
    """
    with open(outfile, "w") as f:
        # write header
        f.write("##gff-version 3\n")
        for seqid, start, end, score, classification in merged_regions:
            attributes = "Name={};score={}".format(classification, score)
            f.write(
                "\t".join(
                    [seqid, "blast_parsed", "repeat_region", str(start), str(end),
                     str(round(score,2)), ".", ".", attributes]
                    )
                )
            f.write("\n")


def sort_blast_table(
        blastfile, seqid_column=0, start_column=6, cpu=1
        ):
    """
    split blast table by seqid and sort by start position
    stores output in temp files
    columns are indexed from 0
    but cut uses 1-based indexing!
    """
    blast_sorted = tempfile.NamedTemporaryFile().name
    # create sorted dictionary seqid counts
    seq_id_counts = {}
    # sort blast file on disk using sort on seqid and start (numeric) position columns
    # using sort command as blast output could be very large
    cmd = "sort -k {0},{0} -k {1},{1}n --parallel {4} {2} > {3}".format(
        seqid_column + 1, start_column + 1, blastfile, blast_sorted, cpu
        )
    subprocess.check_call(cmd, shell=True)

    # count seqids using uniq command
    cmd = "cut -f {0} {1} | uniq -c > {2}".format(
        seqid_column + 1, blast_sorted, blast_sorted + ".counts"
        )
    subprocess.check_call(cmd, shell=True)
    # read counts file and create dictionary
    with open(blast_sorted + ".counts", "r") as f:
        for line in f:
            line = line.strip().split()
            seq_id_counts[line[1]] = int(line[0])
    # remove counts file
    subprocess.call(["rm", blast_sorted + ".counts"])
    # return sorted dictionary and sorted blast file
    return seq_id_counts, blast_sorted


def run_blastn(
        query, db, blastfile, evalue=1e-3, max_target_seqs=999999999, gapopen=2,
        gapextend=1, reward=1, penalty=-1, word_size=9, num_threads=1, outfmt="6"
        ):
    """
    run blastn
    """
    # create temporary blast database:
    db_formated = tempfile.NamedTemporaryFile().name
    cmd = "makeblastdb -in {0} -dbtype nucl -out {1}".format(db, db_formated)
    subprocess.check_call(cmd, shell=True)
    # if query is smaller than 1GB, run blast on single file
    size = os.path.getsize(query)
    print("query size: {} bytes".format(size))
    max_size = 5e8
    if size < max_size:
        cmd = ("blastn -task rmblastn -query {0} -db {1} -out {2} -evalue {3} "
               "-max_target_seqs {4} "
               "-gapopen {5} -gapextend {6} -word_size {7} -num_threads "
               "{8} -outfmt '{9}' -reward {10} -penalty {11} -dust no").format(
            query, db_formated, blastfile, evalue, max_target_seqs, gapopen, gapextend,
            word_size, num_threads, outfmt, reward, penalty
            )
        subprocess.check_call(cmd, shell=True)
    # if query is larger than 1GB, split query in chunks and run blast on each chunk
    else:
        N = int(size // max_size + 1)
        print("splitting query in {} chunks".format(N))
        # create temporary directory
        tmp_dir = tempfile.mkdtemp()
        f_parts = [open(os.path.join(tmp_dir, "part_{}".format(i)), "w") for i in range(N)]
        # split query to max N chuncks, split on fasta header
        with open(query, "r") as f:
            f_parts_cycle = itertools.cycle(f_parts)
            for line in f:
                if line.startswith(">"):
                    out = next(f_parts_cycle)
                    out.write(line)
                else:
                    out.write(line)

        for f in f_parts:
            f.close()
        # run blast on each chunk
        for i in range(N):
            cmd = ("blastn -task rmblastn -query {0} -db {1} -out {2} -evalue {3} "
                   "-max_target_seqs {4} "
                   "-gapopen {5} -gapextend {6} -word_size {7} -num_threads "
                   "{8} -outfmt '{9}' -reward {10} -penalty {11} -dust no").format(
                os.path.join(tmp_dir, "part_{}".format(i)),
                db_formated, blastfile + ".part_{}".format(i), evalue,
                max_target_seqs, gapopen, gapextend, word_size, num_threads, outfmt,
                reward, penalty
                )
            print("running blast on chunk {}".format(i))
            subprocess.check_call(cmd, shell=True)
        # merge blast results
        with open(blastfile, 'w') as outfile:
            for i in range(N):
                fp = blastfile + ".part_{}".format(i)
                print(fp)
                with open(fp) as infile:
                    for line in infile:
                        outfile.write(line)
        # remove temporary directory
        print(blastfile)
        shutil.rmtree(tmp_dir)
    # remove temporary blast database
    os.unlink(db_formated + ".nhr")
    os.unlink(db_formated + ".nin")
    os.unlink(db_formated + ".nsq")


def main():
    """
    main function
    """
    # get command line arguments
    parser = argparse.ArgumentParser(
        description="""This script is used to parse blast output table to gff file""",
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        '-i', '--input', default=None, required=True, help="input file", type=str,
        action='store'
        )
    parser.add_argument(
        '-d', '--db', default=None, required=False,
        help="Fasta file with repeat database", type=str, action='store'
        )
    parser.add_argument(
        '-o', '--output', default=None, required=True, help="output file name", type=str,
        action='store'
        )
    parser.add_argument(
        '-a', '--alternative_classification_coding', default=False,
        help="Use alternative classification coding", action='store_true'
        )
    parser.add_argument(
        '-f', '--fasta_input', default=False,
        help="Input is fasta file instead of blast table", action='store_true'
        )
    parser.add_argument(
        '-c', '--cpu', default=1, help="Number of cpu to use", type=int
        )

    args = parser.parse_args()

    if args.fasta_input:
        # run blast using blastn
        blastfile = tempfile.NamedTemporaryFile().name
        if args.db:
            run_blastn(args.input, args.db, blastfile, num_threads=args.cpu)
        else:
            sys.exit("No repeat database provided")
    else:
        blastfile = args.input

    # sort blast table
    seq_id_counts, blast_sorted = sort_blast_table(blastfile, cpu=args.cpu)
    disjoin_regions = blast2disjoint(
        blast_sorted, seq_id_counts,
        canonical_classification=not args.alternative_classification_coding
        )

    # remove short regions
    disjoin_regions = remove_short_interrupting_regions(disjoin_regions)

    # join neighboring regions with same classification
    merged_regions = join_disjoint_regions_by_classification(disjoin_regions)

    # remove short regions again
    merged_regions = remove_short_interrupting_regions(merged_regions)

    # merge  again neighboring regions with same classification
    merged_regions = join_disjoint_regions_by_classification(merged_regions, max_gap=10)

    # remove short weak regions
    merged_regions = remove_short_regions(merged_regions)

    # last merge
    merged_regions = join_disjoint_regions_by_classification(merged_regions, max_gap=20)
    write_merged_regions_to_gff3(merged_regions, args.output)
    # remove temporary files
    os.remove(blast_sorted)


if __name__ == "__main__":
    main()
