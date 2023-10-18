def hamming_distance(seq1, seq2):
    distance = 0
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must be of equal length')
    print(f'Length of sequence: {len(seq1)}')
    for v1, v2 in zip(seq1, seq2):
        if v1 != v2:
            distance += 1

    return distance

if __name__ == '__main__':
    inputs = []
    with open('data/rosalind_hamm.txt', 'r') as f:
        for line in f:
            inputs.append(line.strip())
    seq1, seq2 = inputs[0], inputs[1]
    print(f'Hamming Distance: {hamming_distance(seq1, seq2)}')