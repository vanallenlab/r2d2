import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-dn', '--dna_normal', type=argparse.FileType('r'), help='Path to DNA normal MAF.')
    parser.add_argument('-dt', '--dna_tumor', type=argparse.FileType('r'), help='Path to DNA tumor MAF.')
    parser.add_argument('-rn', '--rna_normal', type=argparse.FileType('r'), help='Path to RNA normal MAF.')
    parser.add_argument('-rt', '--rna_tumor', type=argparse.FileType('r'), help='Path to RNA tumor MAF.')

    args = parser.parse_args()
