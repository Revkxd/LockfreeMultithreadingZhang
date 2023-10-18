def find_subseq(seq, subseq):
    """
    Find all locations of subseq in seq.
    """
    positions = []
    for i in range(len(seq)):
        if seq[i:i+len(subseq)] == subseq:
            positions.append(i+1)
    return positions

if __name__ == '__main__':
    inputs = []
    with open('data/rosalind_subs.txt', 'r') as f:
        for line in f:
            inputs.append(line.strip())
    seq, subseq = inputs[0], inputs[1]
    print(f'Positions of subsequence: {find_subseq(seq, subseq)}')