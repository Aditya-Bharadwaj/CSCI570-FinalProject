import tracemalloc
import sys
import time
import numpy as np
import re
""" mismatch_penalties = \
[
 [0, 110, 48, 94],
 [110, 0, 118, 48],
 [48, 118, 0, 110],
 [94, 48, 110, 0]
]
 """

mismatch_penalties = {
    "A": {    "A": 0,     "C": 110,    "G": 48,     "T": 94    },
    "C": {    "A": 110,   "C": 0,      "G": 118,    "T": 48    },
    "G": {    "A": 48,    "C": 118,    "G": 0,      "T": 110    },
    "T": {    "A": 94,    "C": 48,      "G": 110,    "T": 0    }
} 


delta = 30

def generate_string(base_str, indices):
    """
    Builds string from given base string and list of indices
    """
    orig_base_str_len = len(base_str)
    for index in indices:
        base_str = base_str[:index+1] + base_str + base_str[index+1:]
        # print(dna_str, len(dna_str))
    assert len(base_str) == ((2**len(indices)) * orig_base_str_len), "Generated string doesn't match expected length!"
    return base_str

def process_input(filename):
    """
    Reads input file and generates strings
    """
    with open(filename, 'r') as input_txt:
        input_lines = [line.rstrip('\n') for line in input_txt.readlines()]
    
    x_base = None
    x_indices = list()

    y_base = None
    y_indices = list()
    
    for line in input_lines:
        if not x_base:
            if line.isalpha():
                x_base = line
        elif not y_base:
            if line.isalpha():
                y_base = line
        if line.isnumeric():
            if not y_base:
                x_indices.append(int(line))
            else:
                y_indices.append(int(line))
    assert x_base and y_base, "Base strings cannot be NoneType"
    
    valid_pattern = re.compile('[AGTCagtc]*')
    assert re.match(valid_pattern, x_base) and re.match(valid_pattern, y_base), "String contains invalid characters!" 
    dna_str_x = generate_string(x_base, x_indices)
    dna_str_y = generate_string(y_base, y_indices)
    return (dna_str_x.upper(), dna_str_y.upper())

def calculate_alignment_cost(X,Y):
    x_len = len(X)
    y_len = len(Y)
    alignment_cost = np.empty((x_len, y_len))
    
    for i in range(0, x_len):
        alignment_cost[i][0] = i * delta

    for j in range(y_len):
        alignment_cost[0][j] = j * delta

    for i in range(1,x_len):
        for j in range(1, y_len):
            # print(X[:i], Y[:j], mismatch_penalties[X[i]][Y[j]])
            alignment_cost[i][j] = \
                min(
                    alignment_cost[i-1][j-1] + mismatch_penalties[X[i]][Y[j]],
                    alignment_cost[i-1][j] + delta,
                    alignment_cost[i][j-1] + delta 
                )
    return alignment_cost


def create_aligned_sequence(alignment_cost, X, Y):
    i = len(X)-1
    j = len(Y)-1
    aligned_X = str()
    aligned_Y = str()
    
    while i >= 0 and j >= 0:
        # print(X[i], Y[j],i , j)
        if(X[i] == Y[j]):
            aligned_X += X[i]
            aligned_Y += Y[j]
            i -= 1
            j -= 1
        elif(alignment_cost[i][j] == alignment_cost[i-1][j-1] + mismatch_penalties[X[i]][Y[j]]):
            aligned_X += X[i]
            aligned_Y += Y[j]
            i -= 1
            j -= 1
        elif(alignment_cost[i][j] == alignment_cost[i-1][j] + delta):
            aligned_X += X[i]
            aligned_Y += '_'
            i -=  1
        elif(alignment_cost[i][j] == alignment_cost[i][j-1] + delta):
            aligned_X += '_'
            aligned_Y += Y[j]
            j -=  1
    
    if j > 0:
        aligned_X += '_' * (j+1)
    if i > 0:
        aligned_Y += '_' * (i+1)

    return aligned_X[::-1], aligned_Y[::-1]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 <filename.py> <input.txt>")
        sys.exit()
    start_time = time.time()
    tracemalloc.start()
    # X,Y = process_input(sys.argv[1])
    X = 'CTGA'
    Y = 'AGCT'
    print(X,Y)
    alc = calculate_alignment_cost(X,Y)
    print(alc, alc[len(X)-1][len(Y)-1])
    print(create_aligned_sequence(alc, X, Y))
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()
    print(" --- Finished in %s seconds --- " % (time.time() - start_time))