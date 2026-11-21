# Global Sequence Alignment with Needleman-Wunsch

## Project goal

This project implements the Needleman-Wunsch dynamic programming algorithm from scratch in Python. The program performs a global, end-to-end alignment of two protein sequences, constructs a scoring matrix, performs traceback, and reports an alignment and score.

## Biological question

Hemoglobin proteins are related globin proteins involved in oxygen transport. We compare human hemoglobin beta with human hemoglobin delta, and then compare human hemoglobin beta with mouse hemoglobin beta-1. The first comparison tests similarity among related human globins. The second tests similarity across mammals.

## Data

The FASTA files were obtained from public UniProtKB records:

| File | Accession | Protein | Organism | Length |
|---|---|---|---|---:|
| `human_hbb_P68871.fasta` | P68871 | Hemoglobin subunit beta | Homo sapiens | 147 aa |
| `human_hbd_P02042.fasta` | P02042 | Hemoglobin subunit delta | Homo sapiens | 147 aa |
| `mouse_hbb1_P02088.fasta` | P02088 | Hemoglobin subunit beta-1 | Mus musculus | 147 aa |

Sources: [P68871](https://www.uniprot.org/uniprotkb/P68871), [P02042](https://www.uniprot.org/uniprotkb/P02042), and [P02088](https://www.uniprot.org/uniprotkb/P02088).

## Scoring scheme

- Match: `+1`
- Mismatch: `-1`
- Gap: `-2`

For every matrix cell, the program chooses the maximum of the diagonal, up, and left possibilities. Traceback uses the deterministic tie order diagonal, up, left.

## Reproduction

```bash
python run_experiments.py
python -m unittest -v
```

The alignment results are written to `results/alignment_results.txt`.

## Results

| Comparison | Score | Matches | Mismatches | Gap columns | Alignment length |
|---|---:|---:|---:|---:|---:|
| Human HBB vs. human HBD | 127 | 137 | 10 | 0 | 147 |
| Human HBB vs. mouse HBB1 | 89 | 118 | 29 | 0 | 147 |

The human HBB/HBD comparison produced the higher score and more exact matches. This is consistent with the expectation that two related human globin proteins are more similar to one another than human and mouse beta-globin proteins under the same scoring scheme. Both comparisons align without gaps because the selected canonical sequences have the same length and the mismatch penalty is less costly than introducing two gap penalties for these positions.

## Interpretation and limitations

The results support strong sequence similarity in both comparisons, with greater similarity between the two human proteins. The score depends on the chosen scoring parameters, so it is not a complete evolutionary distance or functional prediction. The program uses a linear gap penalty and returns one optimal traceback when multiple optimal alignments exist.

## Contributions

- **Arielah:** implemented and tested the Python Needleman-Wunsch algorithm, including the scoring matrix and traceback.
- **Purnima:** researched the biological purpose of global sequence alignment, selected and documented the public protein datasets, interpreted the scores and alignment differences, and prepared the README and presentation content.
