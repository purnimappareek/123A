# Purnima's 5-7 Minute Presentation Section

## Background and importance

My part of the project focuses on the biological purpose of sequence alignment and on interpreting the results. Biological sequences contain information about how proteins are related. When two proteins have many conserved amino acids, that similarity can provide evidence that they share an evolutionary history or an important biological function.

Our project focuses on global sequence alignment. A global alignment compares two sequences from beginning to end. This is useful when the sequences are expected to be related and have similar overall lengths. The method places similar characters in the same columns and can introduce gaps when a gap gives a better overall alignment.

We selected hemoglobin-related proteins because they are biologically meaningful examples. Hemoglobin globin chains are involved in oxygen transport. Comparing their sequences lets us connect the computational score to a real biological question: how similar are related globin proteins within one species and across species?

## Data

We used three public protein sequences from UniProtKB. The first sequence is human hemoglobin subunit beta, HBB, accession P68871. The second is human hemoglobin subunit delta, HBD, accession P02042. The third is mouse hemoglobin subunit beta-1, HBB1, accession P02088. Each canonical sequence contains 147 amino acids.

These records allow us to create two comparisons. Dataset 1 compares human HBB with human HBD. This tests similarity between two related globin proteins in the same species. Dataset 2 compares human HBB with mouse HBB1. This tests similarity between related beta-globin proteins across two mammals.

## Methods

The program uses the Needleman-Wunsch dynamic programming algorithm. It creates a matrix in which each cell represents the best alignment score for two sequence prefixes. The first row and first column are initialized with cumulative gap penalties. The program then calculates every interior cell by comparing three possibilities: a diagonal move for aligning two amino acids, an upward move for adding a gap to sequence 2, or a leftward move for adding a gap to sequence 1.

For this experiment, the scoring scheme was a match score of plus one, a mismatch score of negative one, and a gap penalty of negative two. The program keeps the largest of the three possible scores. After filling the matrix, it performs traceback from the bottom-right cell to the top-left cell. This reconstructs one optimal global alignment.

## Results

The human HBB versus human HBD comparison produced a score of 127. It contained 137 exact matches, 10 mismatches, and zero gap columns across an alignment length of 147 amino acids.

The human HBB versus mouse HBB1 comparison produced a score of 89. It contained 118 exact matches, 29 mismatches, and zero gap columns across the same alignment length.

The first comparison therefore produced the higher score and more exact matches. This is consistent with the biological expectation that two related human globin proteins should be more similar to one another than human and mouse beta-globin proteins under the same scoring scheme.

## Interpretation and conclusion

The results show strong sequence similarity in both comparisons, but the human HBB and HBD pair is more similar in this experiment. The differences between the human and mouse sequences are spread across the alignment as mismatches rather than gaps. Because the sequences are the same length, and because the mismatch penalty is less costly than introducing two gap penalties in these positions, the optimal alignments contain no gaps.

The score should be interpreted carefully. It measures similarity under the scoring values we selected. It is not a complete evolutionary distance and it does not prove that two proteins have exactly the same function. A different match, mismatch, or gap score could change the final numerical score and possibly the alignment.

My contribution was the biological background research, selection and documentation of the public datasets, interpretation of the results, README preparation, and this presentation section. Arielah handled the Python implementation, scoring matrix, traceback, and software testing. Together, the code and the biological interpretation show how a dynamic programming algorithm can answer a meaningful sequence-comparison question.
