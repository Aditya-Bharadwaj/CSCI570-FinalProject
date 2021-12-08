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

    # Find if there's a better way to initialize?
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
    # for row in alignment_cost:
    #    print(row)
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

# def find_min_index(forward_alignment, backward_alignment):
#     min_index, min_value = 0, float('inf')
#     for i in range(len(forward_alignment)):
#         if min_value > (forward_alignment[i][-1] + backward_alignment[i][-1]):
#             min_value = forward_alignment[i][-1] + backward_alignment[i][-1]
#             min_index = i
#     print(min_value)
#     return min_index

def find_min_index(forward_alignment, backward_alignment):
    print(forward_alignment[-1])
    print(backward_alignment[-1])
    partition = [forward_alignment[-1][j] + backward_alignment[-1][j] for j in range(len(forward_alignment[-1]))]
    min_index, min_value = 0, float('inf')
    for j in range(len(forward_alignment[-1])):
        if min_value >= forward_alignment[-1][j] + backward_alignment[-1][j]:
            min_index = j
            min_value = forward_alignment[-1][j] + backward_alignment[-1][j]
    return min_index
    #return partition.index(min(partition))
    
def space_efficient_alignment(X, Y):
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
                                efficient_alignment[0][j-1] + MISMATCH_PENALTIES[X[i-1]][Y[j-1]],
                                efficient_alignment[0][j] + DELTA,
                                efficient_alignment[1][j-1] + DELTA
                            )
        # print(efficient_alignment[0])
        efficient_alignment[0] = efficient_alignment[1].copy()
    # print(efficient_alignment[0])
    return efficient_alignment

# def space_efficient_alignment(X, Y):
#     x_len = len(X)
#     y_len = len(Y)
#     efficient_alignment = [[0 for _ in range(2)] for _ in range(x_len+1)]
    
#     for i in range(x_len+1):
#         efficient_alignment[i][0] = i * DELTA
    
#     for j in range(1, y_len+1):
#         efficient_alignment[0][1] = j * DELTA
#         for i in range(1, x_len+1):
#             efficient_alignment[i][1] = \
#                             min(
#                                 efficient_alignment[i-1][0] + MISMATCH_PENALTIES[X[i-1]][Y[j-1]],
#                                 efficient_alignment[i-1][1] + DELTA,
#                                 efficient_alignment[i][0] + DELTA
#                             )
#         for i in range(x_len+1):
#             efficient_alignment[i][0] = efficient_alignment[i][1]
#     for row in efficient_alignment:
#         print(row)
#     print()
#     return efficient_alignment

p = list()
def divide_and_conquer_alignment(X, Y):
    global p
    print("---------")
    print(X,Y)
    x_len = len(X)
    y_len = len(Y)
    if x_len < 2 or y_len < 2:
        # call original alignment fn
        # maybe we need to return here but not sure what to add to global p
        print("returning %s, %s: %d" % (X, Y, calculate_alignment_cost_brute(X, Y)[-1][-1]))
        return calculate_alignment_cost_brute(X, Y)[-1][-1]

    forward_alignment = space_efficient_alignment(X[:x_len//2], Y)
    backward_alignment = space_efficient_alignment(X[x_len//2:][::-1], Y[::-1])
    # forward_alignment = space_efficient_alignment(X, Y[:int(y_len/2)])
    # backward_alignment = space_efficient_alignment(X[::-1], Y[int((y_len/2)):][::-1])
    min_index = find_min_index(forward_alignment, backward_alignment)
    # append to p
    print("min_index: ", min_index)
    p.append((min_index, x_len//2))
    # p.append((min_index, y_len/2))
    
    print("left")
    left = divide_and_conquer_alignment(X[:x_len//2], Y[:min_index])
    print("right")
    right = divide_and_conquer_alignment(X[x_len//2:], Y[min_index:])
    # left = divide_and_conquer_alignment(X[:min_index], Y[:int(y_len/2)])
    # right = divide_and_conquer_alignment(X[min_index:], Y[int((y_len/2)):])
    print("%s, %s: %d" % (X, Y, left+right))
    return left + right
    

def sequence_alignment_brute(X, Y):
    """
    Returns aligned sequences for X and Y
    """
    alignment_cost_matrix = calculate_alignment_cost_brute(X, Y)
    aligned_X, aligned_Y = create_aligned_sequence(alignment_cost_matrix, X, Y)
    return aligned_X, aligned_Y


if __name__ == '__main__':
    # output_file = 'output.txt'
    # if len(sys.argv) == 3:
    #     output_file = sys.argv[2]
    if len(sys.argv) < 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    X_orig, Y_orig = process_input(sys.argv[1])
    alignment_cost_matrix = calculate_alignment_cost_brute(X_orig, Y_orig)
    X_aligned, Y_aligned = sequence_alignment_brute(X_orig, Y_orig)
    # print(X_orig, Y_orig)
    # with open(output_file, 'w+') as of:
    #     of.write(X_aligned[:50] + " " + X_aligned[-50:])
    #     of.write(Y_aligned[:50] + " " + Y_aligned[-50:])
    #     of.write(str(alignment_cost_matrix[-1][-1]))
    #     of.write(str(tracemalloc.get_traced_memory()[1]/1024))

    print(X_aligned[:50], X_aligned[-50:])
    print(Y_aligned[:50], Y_aligned[-50:])
    print(alignment_cost_matrix[-1][-1])
    
    # space efficient approach
    X_orig, Y_orig = "AGCT", "AAACTAC"
    print("Divide and conquer", divide_and_conquer_alignment(X_orig, Y_orig))
    print("Brute force: ", calculate_alignment_cost_brute(X_orig, Y_orig)[-1][-1])

    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()
    print(" --- Finished in %s seconds --- " % (time.time() - start_time))
    print(time.time() - start_time)
    print(tracemalloc.get_traced_memory()[1]/1024)
