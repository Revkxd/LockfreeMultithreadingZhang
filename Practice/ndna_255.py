# Practice Problem 1-4
from Bio import SeqIO

def count_ACGTN(seq):
    return [seq.count("A"), seq.count("C"), seq.count("G"), seq.count("T"), seq.count("N")]

def GC_content(seq):
    percentage = (seq.count("G") + seq.count("C")) / len(seq) * 100
    return round(percentage, 4)

def transcribe_RNA(seq):
    return seq.replace("T", "U")

def reverse_complement(seq):
    seq = seq[::-1]
    seq = seq.replace("A", "t")
    seq = seq.replace("T", "a")
    seq = seq.replace("C", "g")
    seq = seq.replace("G", "c")
    return seq.upper()

if __name__ == "__main__":
    fasta_seq = SeqIO.parse("data/ndna_255.fa", "fasta")
    for fasta in fasta_seq:
        name, seq = fasta.id, str(fasta.seq)
        counts = count_ACGTN(seq)
        letters = ["A", "C", "G", "T", "N"]
        print(name)
        print(f'Length of sequence: {len(seq)}')
        for letter, count in zip(letters, counts):
            print(f'Count of {letter}: {count}')
        print(f'GC content: {GC_content(seq)}%')
        print(f'Transcribed RNA: {transcribe_RNA(seq)}')
        print(f'Reverse complement: {reverse_complement(seq)}')