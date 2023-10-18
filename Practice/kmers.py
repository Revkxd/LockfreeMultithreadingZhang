def k_mers(k):
    letters = ["A", "C", "G", "T"]

    def k_mer_generator(k):
        if k == 1:
            for letter in letters:
                yield letter
        else:
            for k_mer in k_mer_generator(k - 1):
                for letter in letters:
                    yield k_mer + letter

    return list(k_mer_generator(k))

def overlapping_k_mers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]

if __name__ == "__main__":
    k = 2
    print(f'{k}-mers: {k_mers(k)}')

    seq = "GATAT"
    k = 3
    print(f'Overlapping {k}-mers: {list(overlapping_k_mers(seq, k))}')