'''
Compare two fasta files on ACGT values, give two fasta.gz files.
==Test==
python3 tests/compare_fasta.py 
tests/data/test_output/ffd8e146-1fb1-4611-9f24-813447e0868f_consensus.fasta.gz 
tests/data/output_dir/bams/BaseCall/ffd8e146-1fb1-4611-9f24-813447e0868f.consensus.fasta.gz
'''

import sys
import gzip

def diff_fasta(fasta1, fasta2):
    total = 0
    exact = 0
    mismatch_other = list()
    mismatch_acgt = list()
    for ln, (l1, l2) in enumerate(zip(fasta1, fasta2)):
        if ln == 0:
            continue
        for cn, (c1, c2) in enumerate(zip(l1, l2)):
            total = total + 1
            if c1 == c2:
                exact = exact + 1
            elif c1 == 'N' or c2 == 'N' or c1 == '-' or c2 == '-':
                mismatch_other.append((ln, cn))
            else:
                mismatch_acgt.append((ln, cn))
    return total, mismatch_other, mismatch_acgt, exact

def diff_fasta_gzip(fasta1_gzip, fasta2_gzip):
    with gzip.open(fasta1_gzip, "rb") as fasta1_h:
        fasta1 = fasta1_h.read().decode().split('\n')
    with gzip.open(fasta2_gzip, "rb") as fasta2_h:
        fasta2 = fasta2_h.read().decode().split('\n')
    return diff_fasta(fasta1, fasta2)

def fasta_as_expected(fasta1_gzip,fasta2_gzip):
    total_bases, mismatch_other, mismatch_acgt, exact = diff_fasta_gzip(fasta1_gzip, fasta2_gzip)
    print("total bases: {0}".format(total_bases))
    print("mismatch_acgt: {0}".format(len(mismatch_acgt)))
    print("mismatch_other: {0}".format(len(mismatch_other)))
    print("exact match: {0}%".format(100 * (total_bases - len(mismatch_acgt) - len(mismatch_other)) / total_bases))
    print("match ignoring N/-: {0}%".format(100 * (total_bases - len(mismatch_acgt)) / total_bases))
    print("total mismatch: {0}%".format(100 * (len(mismatch_acgt) + len(mismatch_other)  / total_bases))) 
    print("Two fastas are same: {0}".format(total_bases == exact))

    return total_bases == exact

def main():
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    fasta_is_ok = fasta_as_expected(fasta1, fasta2)
    print(fasta_is_ok)
               
if __name__ == "__main__":
    main()
