import tracemalloc
import sys
import time
import re

MISMATCH_PENALTIES = {
    "A": {"A": 0.0,     "C": 110.0,    "G": 48.0,     "T": 94.0},
    "C": {"A": 110.0,   "C": 0.0,      "G": 118.0,    "T": 48.0},
    "G": {"A": 48.0,    "C": 118.0,    "G": 0.0,      "T": 110.0},
    "T": {"A": 94.0,    "C": 48.0,      "G": 110.0,    "T": 0.0}
}

DELTA = 30.0


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
            if value.isnumeric() is False:
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


def calculate_alignment_cost_brute(X, Y):
    """
    Computes cost of aligning sequences X and Y and subproblems within
    """
    x_len = len(X)
    y_len = len(Y)
    alignment_cost = [[None for _ in range(y_len+1)] for _ in range(x_len+1)]
    for i in range(x_len+1):
        alignment_cost[i][0] = i * DELTA

    for j in range(y_len+1):
        alignment_cost[0][j] = j * DELTA

    for i in range(1, x_len+1):
        for j in range(1, y_len+1):
            alignment_cost[i][j] = \
                min(
                    alignment_cost[i-1][j-1] +
                    MISMATCH_PENALTIES[X[i-1]][Y[j-1]],
                    alignment_cost[i-1][j] + DELTA,
                    alignment_cost[i][j-1] + DELTA
            )
    return alignment_cost


def create_aligned_sequence(alignment_cost, X, Y):
    """
    Creates aligned strings for X and Y using cost matrix for backtracking
    """
    i = len(X)
    j = len(Y)
    aligned_X = ''
    aligned_Y = ''
    while i != 0 and j != 0:
        if alignment_cost[i][j] == (alignment_cost[i-1][j-1] + MISMATCH_PENALTIES[X[i-1]][Y[j-1]]):
            aligned_X += X[i-1]
            aligned_Y += Y[j-1]
            i -= 1
            j -= 1
        elif alignment_cost[i][j] == (alignment_cost[i-1][j] + DELTA):
            aligned_X += X[i-1]
            aligned_Y += '_'
            i -= 1
        elif alignment_cost[i][j] == (alignment_cost[i][j-1] + DELTA):
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


def find_min_index(forward_alignment, backward_alignment):
    """
    Minimizes f + b and returns that index
    """
    min_index, min_value = 0, float('inf')
    for j in range(len(forward_alignment[-1])):
        if min_value >= forward_alignment[-1][j] + backward_alignment[-1][j]:
            min_index = j
            min_value = forward_alignment[-1][j] + backward_alignment[-1][j]
    return min_index


def backward_space_efficient_alignment(X, Y):
    """
    Computes cost of alignment backwards in a space efficient manner
    """
    x_len = len(X)
    y_len = len(Y)
    efficient_alignment = [[0 for _ in range(y_len+1)] for _ in range(2)]

    for j in range(y_len+1):
        efficient_alignment[0][y_len - j] = j * DELTA

    for i in range(x_len-1, -1, -1):
        efficient_alignment[1][-1] = (x_len - i + 2) * DELTA
        for j in range(y_len-1, -1, -1):
            efficient_alignment[1][j] = \
                min(
                    efficient_alignment[1][j+1] + DELTA,
                    efficient_alignment[0][j] + DELTA,
                    efficient_alignment[0][j+1] +
                MISMATCH_PENALTIES[X[i]][Y[j]]
            )
        efficient_alignment[0] = efficient_alignment[1].copy()

    return efficient_alignment


def forward_space_efficient_alignment(X, Y):
    """
    Computes cost of alignment forwards in a space efficient manner
    """
    x_len = len(X)
    y_len = len(Y)
    efficient_alignment = [[0 for _ in range(y_len+1)] for _ in range(2)]

    for j in range(y_len+1):
        efficient_alignment[0][j] = j * DELTA

    for i in range(1, x_len+1):
        efficient_alignment[1][0] = i * DELTA
        for j in range(1, y_len+1):
            efficient_alignment[1][j] = \
                min(
                efficient_alignment[0][j-1] +
                MISMATCH_PENALTIES[X[i-1]][Y[j-1]],
                efficient_alignment[0][j] + DELTA,
                efficient_alignment[1][j-1] + DELTA
            )
        efficient_alignment[0] = efficient_alignment[1].copy()
    return efficient_alignment


def divide_and_conquer_alignment(X, Y):
    """
    Returns aligned sequences for X and Y via Divide & Conquer
    """
    x_len = len(X)
    y_len = len(Y)
    if x_len < 2 or y_len < 2:
        return sequence_alignment_brute(X, Y)
    forward_alignment = forward_space_efficient_alignment(X[:x_len//2], Y)
    backward_alignment = backward_space_efficient_alignment(X[x_len//2:], Y)
    min_index = find_min_index(forward_alignment, backward_alignment)
    left_x, left_y, left_opt = divide_and_conquer_alignment(
        X[: x_len//2], Y[: min_index])
    right_x, right_y, right_opt = divide_and_conquer_alignment(
        X[x_len//2:], Y[min_index:])
    return (left_x+right_x), (left_y+right_y), (left_opt+right_opt)


def sequence_alignment_brute(X, Y):
    """
    Returns aligned sequences for X and Y
    """
    alignment_cost_matrix = calculate_alignment_cost_brute(X, Y)
    aligned_X, aligned_Y = create_aligned_sequence(alignment_cost_matrix, X, Y)
    return aligned_X, aligned_Y, alignment_cost_matrix[-1][-1]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    X_orig, Y_orig = process_input(sys.argv[1])
    X_aligned, Y_aligned, opt_cost = divide_and_conquer_alignment(
        X_orig, Y_orig)
    print(X_aligned[: 50], X_aligned[-50:])
    print(Y_aligned[: 50], Y_aligned[-50:])
    print(opt_cost)
    print(time.time() - start_time)
    print(tracemalloc.get_traced_memory()[1]/1000)
    tracemalloc.stop()
