from pathlib import Path
from needleman_wunsch import needleman_wunsch, read_fasta, alignment_line

ROOT = Path(__file__).parent
DATA = ROOT / 'data'
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)

def summarize(a1, a2):
    matches = sum(x == y for x, y in zip(a1, a2))
    mismatches = sum(x != y and x != '-' and y != '-' for x, y in zip(a1, a2))
    gaps = sum(x == '-' or y == '-' for x, y in zip(a1, a2))
    return matches, mismatches, gaps

experiments = [
    ('dataset_1_human_HBB_vs_HBD', DATA/'human_hbb_P68871.fasta', DATA/'human_hbd_P02042.fasta'),
    ('dataset_2_human_HBB_vs_mouse_HBB1', DATA/'human_hbb_P68871.fasta', DATA/'mouse_hbb1_P02088.fasta'),
]

with open(OUT/'alignment_results.txt', 'w', encoding='utf-8') as out:
    out.write('Needleman-Wunsch experiments\n')
    out.write('Scoring: match=+1, mismatch=-1, gap=-2\n\n')
    for name, p1, p2 in experiments:
        id1, s1 = read_fasta(str(p1)); id2, s2 = read_fasta(str(p2))
        result = needleman_wunsch(s1, s2, 1, -1, -2)
        matches, mismatches, gaps = summarize(result.aligned1, result.aligned2)
        out.write(f'{name}\n')
        out.write(f'{id1}\n{result.aligned1}\n')
        out.write(f'  {alignment_line(result.aligned1, result.aligned2)}\n')
        out.write(f'{id2}\n{result.aligned2}\n')
        out.write(f'score={result.score}; matches={matches}; mismatches={mismatches}; gap_columns={gaps}; alignment_length={len(result.aligned1)}\n\n')
        with open(OUT/f'{name}.txt', 'w', encoding='utf-8') as one:
            one.write(f'{id1}\n{result.aligned1}\n  {alignment_line(result.aligned1, result.aligned2)}\n{id2}\n{result.aligned2}\n')
            one.write(f'\nScore: {result.score}\nMatches: {matches}\nMismatches: {mismatches}\nGap columns: {gaps}\nAlignment length: {len(result.aligned1)}\n')
print(OUT/'alignment_results.txt')
