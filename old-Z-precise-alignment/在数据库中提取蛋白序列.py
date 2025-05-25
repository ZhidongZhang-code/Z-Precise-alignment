def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        current_sequence = ''
        for line in file:
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = current_sequence
                current_id = line.strip().split('_**_')[1].split()[0]
                current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_id is not None:
            sequences[current_id] = current_sequence
    return sequences

def filter_sequences(sequences, ids_to_keep):
    filtered_sequences = {}
    for id, sequence in sequences.items():
        if id in ids_to_keep:
            filtered_sequences[id] = sequence
    return filtered_sequences

def main():
    fasta_file = "blast_normailzation_excate_protein.fasta"
    ids_file = "fileA"

    sequences = read_fasta(fasta_file)

    with open(ids_file, 'r') as file:
        ids_to_keep = set(line.strip() for line in file)

    filtered_sequences = filter_sequences(sequences, ids_to_keep)

    for id, sequence in filtered_sequences.items():
        print(f'>{id}')
        print(sequence)

if __name__ == "__main__":
    main()
