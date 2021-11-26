import unittest
from code import generate_string, read_file


class TestInputStringGeneration(unittest.TestCase):
    def test_input_string_generation(self):
        output = generate_string("ACTG", [3, 6, 1])
        self.assertEqual(output, "ACACTGACTACTGACTGGTGACTACTGACTGG")
        output2 = generate_string("TACG", [1, 2, 9])
        self.assertEqual(output2, "TATTATACGCTATTATACGCGACGCGGACGCG")

    def test_read_file(self):
        (x_base, x_indices, y_base, y_indices) = read_file(
            "./BaseTestcases_CS570FinalProject/input1.txt")
        self.assertEqual(x_base, "ACTG")
        self.assertEqual(x_indices, [3, 6, 1, 1])
        self.assertEqual(y_base, "TACG")
        self.assertEqual(y_indices, [1, 2, 9, 2])


if __name__ == '__main__':
    unittest.main()
