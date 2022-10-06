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
    sys.exit("Python 3.10 or a more recent version is required.")

def make_temp_files(number_of_files):
    """
    Make named temporary files, file will not be deleted upon exit!
    :param number_of_files:
    :return:
    filepaths
    """
    temp_files = []
    for i in range(number_of_files):
        temp_files.append(tempfile.NamedTemporaryFile(delete=False).name)
        os.remove(temp_files[-1])
    return temp_files


def split_fasta_to_chunks(fasta_file, chunk_size=100000000, overlap=100000):
    """
    Split fasta file to chunks, sequences longe than chuck size are split to overlaping
    peaces. If sequences are shorter, chunck with multiple sequences are created.
    :param fasta_file:

    :param fasta_file:
    :param chunk_size:
    :param overlap:
    :return:
    fasta_file_split
    matching_table (list of lists [header,chunk_number, start, end, new_header])
    """
    min_chunk_size = chunk_size * 2
    fasta_sizes_dict = read_fasta_sequence_size(fasta_file)
    # calculate size of items in fasta_dist dictionary
    fasta_size = sum(fasta_sizes_dict.values())

    # calculates ranges for splitting of fasta files and store them in list
    matching_table = []
    fasta_file_split = tempfile.NamedTemporaryFile(delete=False).name
    for header, size in fasta_sizes_dict.items():
        print(header, size, min_chunk_size)

        if size > min_chunk_size:
            number_of_chunks = int(size / chunk_size)
            adjusted_chunk_size = int(size / number_of_chunks)
            for i in range(number_of_chunks):
                start = i * adjusted_chunk_size
                end = ((i + 1) *
                       adjusted_chunk_size
                       + overlap) if i + 1 < number_of_chunks else size
                new_header = header + '_' + str(i)
                matching_table.append([header, i, start, end, new_header])
        else:
            new_header = header + '_0'
            matching_table.append([header, 0, 0, size, new_header])
    # read sequences from fasta files and split them to chunks according to matching table
    # open output and input files, use with statement to close files
    number_of_temp_files = len(matching_table)
    fasta_dict = read_single_fasta_to_dictionary(open(fasta_file, 'r'))
    with open(fasta_file_split, 'w') as fh_out:
        for header in fasta_dict:
            matching_table_part = [x for x in matching_table if x[0] == header]
            for header2, i, start, end, new_header in matching_table_part:
                fh_out.write('>' + new_header + '\n')
                fh_out.write(fasta_dict[header][start:end] + '\n')
    temp_files_fasta = make_temp_files(number_of_temp_files)
    file_handles = [open(temp_file, 'w') for temp_file in temp_files_fasta]
    # make dictionary seq_id_sorted as keys and values as file handles
    fasta_seq_size = read_fasta_sequence_size(fasta_file_split)
    seq_id_size_sorted = [i[0] for i in sorted(
        fasta_seq_size.items(), key=lambda x: int(x[1]), reverse=True
        )]
    seq_id_file_handle_dict = dict(zip(seq_id_size_sorted, itertools.cycle(file_handles)))

    # write sequences to temporary files
    with open(fasta_file_split, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]
                seq_id_file_handle_dict[header].write(line)
            else:
                seq_id_file_handle_dict[header].write(line)
    os.remove(fasta_file_split)
    # close file handles
    for file_handle in file_handles:
        file_handle.close()
    return temp_files_fasta, matching_table


def read_fasta_sequence_size(fasta_file):
    """Read size of sequence into dictionary"""
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]  # remove part of name after space
                fasta_dict[header] = 0
            else:
                fasta_dict[header] += len(line.strip())
    return fasta_dict


def read_single_fasta_to_dictionary(fh):
    """
    Read fasta file into dictionary
    :param fh:
    :return:
    fasta_dict
    """
    fasta_dict = {}
    for line in fh:
        if line[0] == '>':
            header = line.strip().split(' ')[0][1:]  # remove part of name after space
            fasta_dict[header] = []
        else:
            fasta_dict[header] += [line.strip()]
    fasta_dict = {k: ''.join(v) for k, v in fasta_dict.items()}
    return fasta_dict


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
    max_size = 3e6
    overlap = 100000
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
        print(f"query is larger than {max_size}, splitting query in chunks")
        query_parts, matching_table = split_fasta_to_chunks(query, max_size, overlap)
        print(query_parts)
        for i, part in enumerate(query_parts):
            print(f"running blast on chunk {i}")
            print(part)
            cmd = ("blastn -task rmblastn -query {0} -db {1} -out {2} -evalue {3} "
                   "-max_target_seqs {4} "
                   "-gapopen {5} -gapextend {6} -word_size {7} -num_threads "
                   "{8} -outfmt '{9}' -reward {10} -penalty {11} -dust no").format(
                part, db_formated, f'{blastfile}.{i}', evalue, max_target_seqs, gapopen,
                gapextend,
                word_size, num_threads, outfmt, reward, penalty
                )
            subprocess.check_call(cmd, shell=True)
            print(cmd)
            # remove part file
            # os.unlink(part)
        # merge blast results and recalculate start, end positions and header
        merge_blast_results(blastfile, matching_table, n_parts=len(query_parts))

    # remove temporary blast database
    os.unlink(db_formated + ".nhr")
    os.unlink(db_formated + ".nin")
    os.unlink(db_formated + ".nsq")

def merge_blast_results(blastfile, matching_table, n_parts):
    """
    Merge blast tables and recalculate start, end positions based on
    matching table
    """
    with open(blastfile, "w") as f:
        matching_table_dict = {i[4]: i for i in matching_table}
        print(matching_table_dict)
        for i in range(n_parts):
            with open(f'{blastfile}.{i}', "r") as f2:
                for line in f2:
                    line = line.strip().split("\t")
                    # seqid (header) is in column 1
                    seqid = line[0]
                    line[0] = matching_table_dict[seqid][0]
                    # increase coordinates by start position of chunk
                    line[6] = str(int(line[6]) + matching_table_dict[seqid][2])
                    line[7] = str(int(line[7]) + matching_table_dict[seqid][2])
                    f.write("\t".join(line) + "\n")
            # remove temporary blast file
            # os.unlink(f'{blastfile}.{i}')

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
