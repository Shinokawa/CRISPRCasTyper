import unittest
from src.cas_typing import CASTyping

class TestCASTyping(unittest.TestCase):

    def setUp(self):
        # Initialize the CASTyping object with mock data or parameters
        self.cas_typing = CASTyping()

    def test_predict_cas_type(self):
        # Example input for testing
        test_protein_sequence = "MKTVRQERLKSIV"
        expected_output = "Type I"  # Replace with the expected output based on your model

        # Call the method to predict the CAS type
        result = self.cas_typing.predict_cas_type(test_protein_sequence)

        # Assert the expected output matches the result
        self.assertEqual(result, expected_output)

    def test_invalid_sequence(self):
        # Test with an invalid protein sequence
        invalid_sequence = "INVALID_SEQUENCE"

        with self.assertRaises(ValueError):
            self.cas_typing.predict_cas_type(invalid_sequence)

if __name__ == '__main__':
    unittest.main()