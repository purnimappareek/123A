#!/usr/bin/env python3
"""Global sequence alignment with the Needleman-Wunsch algorithm.

Rows correspond to sequence 1 and columns correspond to sequence 2.
The implementation uses a linear gap penalty and returns the full scoring
matrix, one optimal traceback, and the final alignment score.
"""

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class AlignmentResult:
    sequence1: str
    sequence2: str
    aligned1: str
    aligned2: str
    score: int
    matrix: List[List[int]]


def needleman_wunsch(
    sequence1: str,
    sequence2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> AlignmentResult:
    """Return one optimal global alignment of two sequences.

    Parameters are user-defined scores. The two input sequences are stripped
    of whitespace and converted to uppercase before alignment.
    """
    s1 = "".join(sequence1.split()).upper()
    s2 = "".join(sequence2.split()).upper()
    if not s1 or not s2:
        raise ValueError("Both sequences must contain at least one character.")

    n, m = len(s1), len(s2)
    matrix = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + gap
    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = matrix[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            matrix[i][j] = max(diagonal, up, left)

    # Deterministic tie-breaking: diagonal, then up, then left.
    i, j = n, m
    aligned1, aligned2 = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            step = match if s1[i - 1] == s2[j - 1] else mismatch
            if matrix[i][j] == matrix[i - 1][j - 1] + step:
                aligned1.append(s1[i - 1])
                aligned2.append(s2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            aligned1.append(s1[i - 1])
            aligned2.append("-")
            i -= 1
        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap:
            aligned1.append("-")
            aligned2.append(s2[j - 1])
            j -= 1
        else:
            raise RuntimeError("Traceback reached a cell without a valid predecessor.")

    a1 = "".join(reversed(aligned1))
    a2 = "".join(reversed(aligned2))
    return AlignmentResult(s1, s2, a1, a2, matrix[n][m], matrix)


def alignment_line(aligned1: str, aligned2: str) -> str:
    """Return | for matches, . for mismatches, and spaces for gaps."""
    return "".join("|" if a == b else "." if a != "-" and b != "-" else " " for a, b in zip(aligned1, aligned2))


def read_fasta(path: str) -> Tuple[str, str]:
    """Read one FASTA record and return (identifier, sequence)."""
    identifier = path
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    break
                identifier = line[1:].strip() or path
            else:
                chunks.append(line)
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"No sequence found in {path}")
    return identifier, sequence


def print_result(result: AlignmentResult, name1: str = "sequence1", name2: str = "sequence2") -> None:
    print(f"{name1}: {result.aligned1}")
    print(f"       {alignment_line(result.aligned1, result.aligned2)}")
    print(f"{name2}: {result.aligned2}")
    print(f"Score: {result.score}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Perform global Needleman-Wunsch alignment.")
    parser.add_argument("sequence1", help="Sequence text or path to a FASTA file")
    parser.add_argument("sequence2", help="Sequence text or path to a FASTA file")
    parser.add_argument("--match", type=int, default=1)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-2)
    args = parser.parse_args()

    import os
    if os.path.isfile(args.sequence1):
        name1, seq1 = read_fasta(args.sequence1)
    else:
        name1, seq1 = "sequence1", args.sequence1
    if os.path.isfile(args.sequence2):
        name2, seq2 = read_fasta(args.sequence2)
    else:
        name2, seq2 = "sequence2", args.sequence2

    result = needleman_wunsch(seq1, seq2, args.match, args.mismatch, args.gap)
    print_result(result, name1, name2)
