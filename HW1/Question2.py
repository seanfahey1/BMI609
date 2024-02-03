import sys

"""
12.5 points - You will write a program in Python/R/C (again, I don't care what language you use - just use something
that you're comfortable with) to read a FASTQ file as input, along with the kmer length (user-input), then construct a
de Bruijn graph from the reads. Provide the output as two separate files - (a) a file containing a list of the nodes,
and (b) a file that contains the list of edges. You will then run your program on the provided FASTQ file, and obtain
the output from it. You'll find this FASTQ file under Files->Assignment Files->assignment1.fastq.
"""


def window(seq, window_size):
    windows = []

    for i in range(len(seq) - (window_size - 1)):
        windows.append(seq[i : i + window_size])

    return windows


def get_reads():
    with open("assignment1.fastq", "r") as fastq:
        fastq = fastq.readlines()

    reads = fastq[1::4]  # start @ line 1, grab every 4th line

    return reads


def get_nodes_and_edges(reads, kmer_size):
    edges = set()
    nodes = set()

    for read in reads:
        read = read.strip("\n")

        kmers = window(read, kmer_size)
        for kmer in kmers:
            start, end = kmer[:-1], kmer[1:]

            nodes.add(kmer)
            edges.add((start, end))

    nodes = list(nodes)
    edges = list(edges)

    with open("nodes.txt", "w") as out:
        out.write("\n".join(nodes))

    with open("edges.txt", "w") as out:
        out.write("\n".join([" -> ".join(x) for x in edges]))
    return nodes, edges


def main():
    try:
        kmer_size = int(input("kmer length: "))
    except ValueError:
        raise ValueError(
            "Input expects a whole number integer, no spaces or special characters."
        )

    reads = get_reads()
    nodes, edges = get_nodes_and_edges(reads, kmer_size)


if __name__ == "__main__":
    sys.exit(main())
