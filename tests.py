import unittest
from code import generate_string


class TestInputStringGeneration(unittest.TestCase):
    def test_input_string_generation(self):
        output = generate_string("ACTG", [3, 6, 1])
        self.assertEqual(output, "ACACTGACTACTGACTGGTGACTACTGACTGG")
        output2 = generate_string("TACG", [1, 2, 9])
        self.assertEqual(output2, "TATTATACGCTATTATACGCGAXGCGGACGCG")


if __name__ == '__main__':
    unittest.main()
