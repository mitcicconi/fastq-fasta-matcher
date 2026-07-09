#!/usr/bin/env python3
"""Generates a tiny synthetic reference/read set under sample_data/ for demoing the app.

ref_A and ref_B get reads derived from them (with a few simulated sequencing
errors and reverse-complement strand mixed in). ref_C and ref_D get no reads,
to demonstrate the "no confident match" case. A few reads are pure noise,
to demonstrate reads with no confident reference match.
"""

import random

random.seed(7)

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def mutate(seq, n_subs):
    seq = list(seq)
    for _ in range(n_subs):
        i = random.randrange(len(seq))
        seq[i] = random.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


REFS = {
    "ref_A": random_seq(200),
    "ref_B": random_seq(180),
    "ref_C": random_seq(220),  # will have no matching reads
    "ref_D": random_seq(160),  # will have no matching reads
}

with open("sample_data/expected_references.fasta", "w") as fh:
    for name, seq in REFS.items():
        fh.write(f">{name}\n{seq}\n")

reads = []

for i in range(8):
    seq = mutate(REFS["ref_A"], n_subs=random.randint(0, 4))
    if random.random() < 0.4:
        seq = revcomp(seq)
    reads.append((f"read_A_{i}", seq))

for i in range(6):
    seq = mutate(REFS["ref_B"], n_subs=random.randint(0, 3))
    if random.random() < 0.4:
        seq = revcomp(seq)
    reads.append((f"read_B_{i}", seq))

for i in range(3):
    reads.append((f"read_noise_{i}", random_seq(150)))

random.shuffle(reads)

with open("sample_data/sequencing_results.fastq", "w") as fh:
    for name, seq in reads:
        qual = "I" * len(seq)
        fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

print(f"Wrote {len(REFS)} references and {len(reads)} reads to sample_data/")
