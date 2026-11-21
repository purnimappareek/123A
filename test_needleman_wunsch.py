import unittest

from needleman_wunsch import alignment_line, needleman_wunsch


class NeedlemanWunschTests(unittest.TestCase):
    def test_identical_sequences(self):
        result = needleman_wunsch("ACGT", "ACGT")
        self.assertEqual(result.aligned1, "ACGT")
        self.assertEqual(result.aligned2, "ACGT")
        self.assertEqual(result.score, 4)

    def test_alignment_lengths_match(self):
        result = needleman_wunsch("GATTACA", "GCATGCU")
        self.assertEqual(len(result.aligned1), len(result.aligned2))
        self.assertEqual(result.score, result.matrix[-1][-1])

    def test_gap_can_be_selected(self):
        result = needleman_wunsch("AC", "ABC")
        self.assertEqual(len(result.aligned1), len(result.aligned2))
        self.assertIn("-", result.aligned1 + result.aligned2)

    def test_alignment_line(self):
        self.assertEqual(alignment_line("A-CG", "AT-G"), "|  |")


if __name__ == "__main__":
    unittest.main()
