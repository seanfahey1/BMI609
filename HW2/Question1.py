import sys
from Bio import SeqIO


"""
For this assignment, you will write code to
    (1) take an input genome as a FASTA file,
    (2) construct a suffix array/tree from it,
    (3) take a read as input from a FASTA file, and
    (4) print all positions where that read occurs in the genome using the bisection algorithm.

Please submit your code, and output file with all indices (locations of match) for the input files provided (you'll
find this under Assignments - Assignment2_refgenome.fasta, Assignment2_read.fasta).
"""


class ReadNotFoundError(Exception):
    def __init__(self, read):
        self.read = read

    def __repr__(self):
        print(f"Read {self.read} not found in reference sequence.")


def suffix_array(seq):
    sa = [(seq[i:], i) for i in range(len(seq))]
    return sorted(sa, key=lambda x: x[0])


def bisect(sa, read):
    if len(sa) == 0:
        return None

    mid_point = len(sa) // 2
    # cur is a tuple of (seq, index)
    cur = sa[mid_point]

    if cur[0].startswith(read):
        return cur
    elif cur[0] > read:
        return bisect(sa[:mid_point], read)
    else:
        return bisect(sa[mid_point + 1 :], read)


def get_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    return record.seq


def write_positions(positions, read_len):
    with open("start_indexes.txt", "w") as output:
        for position in positions:
            output.write(f"{position}-{position + read_len} \n")


def main():
    # get the input files for the assignment
    ref_seq = get_fasta("Assignment2_refgenome.fasta")
    read_seq = get_fasta("Assignment2_read.fasta")

    sa = suffix_array(ref_seq)

    # keep repeating bisection algo
    matches = []
    while True:
        match = bisect(sa, read_seq)

        # stop once it fails to find a match
        if match is None:
            break
        matches.append(match)

        # remove the match so it isn't accessed on the next loop
        sa.pop(sa.index(match))

    if len(matches) == 0:
        raise ReadNotFoundError(read_seq)

    positions = [x[1] for x in matches]

    _ = [print(x) for x in positions]
    write_positions(positions, len(read_seq))


if __name__ == "__main__":
    sys.exit(main())
