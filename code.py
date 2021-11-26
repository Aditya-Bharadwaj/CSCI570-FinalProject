import tracemalloc
import sys
import time
import re
from array import array

""" mismatch_penalties = \
[
 [0, 110, 48, 94],
 [110, 0, 118, 48],
 [48, 118, 0, 110],
 [94, 48, 110, 0]
]
 """

mismatch_penalties = {
    "A": {"A": 0,     "C": 110,    "G": 48,     "T": 94},
    "C": {"A": 110,   "C": 0,      "G": 118,    "T": 48},
    "G": {"A": 48,    "C": 118,    "G": 0,      "T": 110},
    "T": {"A": 94,    "C": 48,      "G": 110,    "T": 0}
}

# mismatch_penalties = {
#     "m": {    "m": 0,     "e": 3,    "a": 3,     "n": 1    },
#     "e": {    "m": 3,   "e": 0,      "a": 1,    "n": 3    },
#     "a": {    "m": 3,    "e": 1,    "a": 0,      "n": 3    },
#     "n": {    "m": 1,    "e": 3,      "a": 3,    "n": 0    }
# }
delta = 30


def generate_string(base_str, indices):
    """
    Builds string from given base string and list of indices
    """
    orig_base_str_len = len(base_str)
    for index in indices:
        base_str = base_str[:index+1] + base_str + base_str[index+1:]
    assert len(base_str) == ((2**len(indices)) *
                             orig_base_str_len), "Generated string doesn't match expected length!"
    return base_str


def read_file(filename):
    """
    Reads the file and spits out the base strings and indices
    """
    with open(filename, 'r') as input_txt:
        x_base = input_txt.readline().rstrip()
        value = None
        x_indices = list()
        for line in input_txt:
            value = line.rstrip()
            if value.isnumeric() == False:
                break
            x_indices.append(int(value))
        y_base = value.rstrip()
        y_indices = list()
        for line in input_txt:
            value = line.rstrip()
            y_indices.append(int(value))
        return (x_base, x_indices, y_base,  y_indices)


def process_input(filename):
    """
    Reads input file and generates strings
    """
    (x_base, x_indices, y_base, y_indices) = read_file(filename)
    assert x_base and y_base, "Base strings cannot be NoneType"
    valid_pattern = re.compile('[AGTCagtc]*')
    assert re.match(valid_pattern, x_base) and re.match(
        valid_pattern, y_base), "String contains invalid characters!"
    dna_str_x = generate_string(x_base, x_indices)
    dna_str_y = generate_string(y_base, y_indices)
    return (dna_str_x.upper(), dna_str_y.upper())


def verify_output(X, Y, output_file):
    with open(output_file, 'r') as of:
        values = [line.rstrip() for line in of.readlines()]
    X_o_f50 = values[0][:50]
    X_o_l50 = values[0][-50:]
    Y_o_f50 = values[1][:50]
    Y_o_l50 = values[1][-50:]


def calculate_alignment_cost(X, Y):
    x_len = len(X)
    y_len = len(Y)

    # Find if there's a better way to initialize
    alignment_cost = [[None for _ in range(y_len+1)] for _ in range(x_len+1)]
    for i in range(0, x_len+1):
        alignment_cost[i][0] = i * delta

    for j in range(y_len+1):
        alignment_cost[0][j] = j * delta

    for i in range(1, x_len+1):
        for j in range(1, y_len+1):
            alignment_cost[i][j] = \
                min(
                    alignment_cost[i-1][j-1] +
                mismatch_penalties[X[i-1]][Y[j-1]],
                    alignment_cost[i-1][j] + delta,
                    alignment_cost[i][j-1] + delta
            )
    return alignment_cost



def create_aligned_sequence(alignment_cost, X, Y):
    i = len(X)
    j = len(Y)
    aligned_X = str()
    aligned_Y = str()
    while i != 0 and j != 0:
        if alignment_cost[i][j] == (alignment_cost[i-1][j-1] + mismatch_penalties[X[i-1]][Y[j-1]]):
            aligned_X += X[i-1]
            aligned_Y += Y[j-1]
            i -= 1
            j -= 1
        elif alignment_cost[i][j] == (alignment_cost[i-1][j] + delta):
            aligned_X += X[i-1]
            aligned_Y += '_'
            i -= 1
        elif alignment_cost[i][j] == (alignment_cost[i][j-1] + delta):
            aligned_X += '_'
            aligned_Y += Y[j-1]
            j -= 1

    if i > 0:
        remaining_str = X[:i]
        aligned_X += remaining_str[::-1]
        aligned_Y += '_' * i
    elif j > 0:
        remaining_str = Y[:j]
        aligned_Y += remaining_str[::-1]
        aligned_X += '_' * j
    return aligned_X[::-1], aligned_Y[::-1]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    X_orig, Y_orig = process_input(sys.argv[1])
    # Y_orig, X_orig = process_input(sys.argv[1])
    alignment_cost_matrix = calculate_alignment_cost(X_orig, Y_orig)
    print("Optimal Alignment Cost:",
          alignment_cost_matrix[len(X_orig)][len(Y_orig)])
    X_a, Y_a = create_aligned_sequence(alignment_cost_matrix, X_orig, Y_orig)
    print(X_orig, Y_orig)
    print(X_a, Y_a)
    # verify_output(X_a, Y_a, 'BaseTestcases_CS570FinalProject/output1.txt')
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()
    print(" --- Finished in %s seconds --- " % (time.time() - start_time))
