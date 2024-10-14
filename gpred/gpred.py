import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional
from textwrap import fill


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int, 
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int, 
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path, 
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    fasta_file = isfile(fasta_file)
    with open(fasta_file, 'r') as f:
        seq = [line.rstrip() for line in f.readlines()[1:]]
        seq = ''.join(seq)
        assert seq.isupper(), f"sequence isn't in all caps : {seq[0:5]}"
        return(seq)


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    matchobj = start_regex.search(sequence, start, stop)
    if matchobj==None:
        return(None)
    return(matchobj.span()[0])


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    matchobj = stop_regex.search(sequence,  start)
    if matchobj==None:
        return(None)
    return(matchobj.span()[0])


def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    startidx = start - max_shine_dalgarno_distance
    endidx = start - 6
    if start<0:
        return(False)
    if endidx>len(sequence)-1:
        endidx = len(sequence)-1
    matchobj = shine_regex.search(sequence, startidx, endidx)
    if matchobj == None:
        return(False)
    else:
        return(True)


def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern, 
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    predgene = []
    current_pos = 0
    while len(sequence) - current_pos >= min_gap:
        current_pos = find_start(start_regex, sequence, current_pos, len(sequence)-1)
        # print("start",current_pos)
        if current_pos == None:
            current_pos = len(sequence) 
        if current_pos != None:
            stop = find_stop(stop_regex, sequence, current_pos)
            # print("stop",current_pos)
            if stop != None:
                if stop - current_pos > min_gap:
                    # print("good size")
                    if has_shine_dalgarno(shine_regex, sequence, current_pos, max_shine_dalgarno_distance):
                        # print("has_shine_dalgarno")
                        predgene.append([current_pos,stop])
                        current_pos = stop + 3 + min_gap
                    else:
                        current_pos += 1
                else:
                    current_pos += 1
            else:
                current_pos += 1
    return predgene
                


def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with open(predicted_genes_file, "w") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            print(probable_genes)
            for pos in probable_genes:
                # print(pos)
                # print("test")
                predict_genes_writer.writerows(pos)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))

    # with open(predicted_genes_file, 'w', newline='') as f:
    #     writer = csv.writer(f, delimiter=",")
    #     writer.writerows(probable_genes)

def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str, 
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep, 
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG 
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    # args = get_arguments()
    # Let us do magic in 5' to 3'
    
    fastasequence = read_fasta("data/listeria.fna")
    # print(re.search(fastasequence, "\n"))
    # print(fastasequence[::-1])
    
    pred5 = predict_genes(fastasequence, start_regex, stop_regex, shine_regex, min_gene_len = 10, max_shine_dalgarno_distance = 50, min_gap = 5)
    reversesequence = reverse_complement(fastasequence)

    # print(re.search(fastasequence, "\\n"))
    pred3r = predict_genes(reversesequence, start_regex, stop_regex, shine_regex, min_gene_len = 10, max_shine_dalgarno_distance = 50, min_gap = 5)
    pred3 = []
    for pred in pred3r:
        pred3.append([pred[1],pred[0]])
    all_pred = pred5 + pred3
    all_pred = sorted(all_pred, key=lambda pos: pos[0])
    print(all_pred)
    # write_genes_pos("data/predict_gene.csv", all_pred)
    
    # Don't forget to uncomment !!!
    # Call these function in the order that you want
    # We reverse and complement
    #sequence_rc = reverse_complement(sequence)
    # Call to output functions
    #write_genes_pos(args.predicted_genes_file, probable_genes)
    #write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp)



if __name__ == '__main__':
    main()
