import unittest
from src.hmmer_search import HMMERSearch

class TestHMMERSearch(unittest.TestCase):

    def setUp(self):
        self.hmmer_search = HMMERSearch()

    def test_search_for_crispr_cas_proteins(self):
        # Test with a sample protein sequence
        sample_sequence = ">protein1\nMKTAYIAKQRQISFVKSHFSRQDIL"
        results = self.hmmer_search.search(sample_sequence)
        self.assertIsInstance(results, list)
        self.assertGreater(len(results), 0)

    def test_empty_sequence(self):
        # Test with an empty sequence
        empty_sequence = ""
        results = self.hmmer_search.search(empty_sequence)
        self.assertEqual(results, [])

    def test_invalid_sequence(self):
        # Test with an invalid sequence
        invalid_sequence = ">protein2\nINVALIDSEQUENCE"
        with self.assertRaises(ValueError):
            self.hmmer_search.search(invalid_sequence)

if __name__ == '__main__':
    unittest.main()